# Adapted from NeoFox python implementation of Hex pathogenicity scoring algorithm
# This was heavily vibe-coded

from math import ceil, floor
import argparse
from collections import defaultdict
from multiprocessing import Pool, cpu_count

import numpy as np
from Bio import SeqIO
from Bio.Align import substitution_matrices
from Bio.Data.IUPACData import protein_letters


# ---------------------------------------------------------------------------
# BLOSUM62 lookup table as a NumPy matrix (20x20, indexed by amino acid order)
# ---------------------------------------------------------------------------

_BLOSUM62_BIO = substitution_matrices.load("BLOSUM62")
AA_ORDER = protein_letters  # 20-character string of standard amino acids
AA_INDEX = {aa: i for i, aa in enumerate(AA_ORDER)}

# Build a (20, 20) integer matrix so scoring becomes array indexing
BLOSUM62_MATRIX = np.zeros((20, 20), dtype=np.int32)
for i, aa1 in enumerate(AA_ORDER):
    for j, aa2 in enumerate(AA_ORDER):
        BLOSUM62_MATRIX[i, j] = _BLOSUM62_BIO[aa1, aa2]


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def read_iedb(fasta_file):
    sequences = []
    with open(fasta_file, "r") as f:
        for line in f:
            line = line.strip()
            if line:
                if not any(aa not in protein_letters for aa in line):
                    sequences.append(line)
        # for record in SeqIO.parse(handle, "fasta"):
        #     if not any(aa not in protein_letters for aa in record.seq):
        #         sequences.append(str(record.seq))
    return sequences


def read_peptides(peptide_file):
    with open(peptide_file, "r") as f:
        return [line.strip().upper() for line in f if line.strip()]


# ---------------------------------------------------------------------------
# Weights
# ---------------------------------------------------------------------------
# Weights serve to give more importance to TCR interacting residues with regards to similarity
# Potential to introduce iedb's masking technique for different HLA allele anchors

def get_sequence_weights(length, magic_number=4):
    """Return a positional weight array for a peptide of the given length."""
    mid_score = ceil(length / 2) * magic_number
    weights = list(range(1, mid_score, magic_number))
    weights.extend(reversed(weights[0:floor(length / 2)]))
    top_floor = floor(length / 3)
    weights[0:top_floor] = list(range(1, top_floor + 1))
    tail = length - top_floor
    weights[tail:length] = list(reversed(range(1, top_floor + 1)))
    return np.array(weights, dtype=np.int32)  # shape: (length,)


# ---------------------------------------------------------------------------
# Encoding
# ---------------------------------------------------------------------------

def encode_sequences(sequences):
    """Convert a list of equal-length peptide strings to a (N, L) int32 array."""
    n = len(sequences)
    l = len(sequences[0])
    arr = np.empty((n, l), dtype=np.int32)
    for i, seq in enumerate(sequences):
        arr[i] = [AA_INDEX[aa] for aa in seq]
    return arr


# ---------------------------------------------------------------------------
# Core vectorised scoring for one length group
# ---------------------------------------------------------------------------

def score_length_group(peptides, iedb_encoded, weights):
    """
    Score all peptides of a given length against the corresponding IEDB
    reference matrix in one vectorised pass.

    Parameters
    ----------
    peptides    : (P, L) int32 array  — encoded query peptides
    iedb_encoded: (R, L) int32 array  — encoded IEDB references
    weights     : (L,)   int32 array  — positional weights

    Returns
    -------
    (P,) float32 array of best scores, one per query peptide.
    """
    P, L = peptides.shape
    R     = iedb_encoded.shape[0]

    # BLOSUM62 scores for all (reference, peptide, position) triples
    # blosum_scores[r, p, pos] = BLOSUM62_MATRIX[iedb[r, pos], peptide[p, pos]]
    #
    # Achieved without an explicit loop:
    #   - iedb_encoded[:, None, :]  → (R, 1, L)  broadcast over peptides
    #   - peptides[None, :, :]      → (1, P, L)  broadcast over references
    # Indexing BLOSUM62_MATRIX with two (R, P, L) arrays gives (R, P, L) scores.
    blosum_scores = BLOSUM62_MATRIX[
        iedb_encoded[:, None, :],   # row indices    (R, 1, L)
        peptides[None, :, :]        # column indices (1, P, L)
    ]                               # result: (R, P, L)

    # Apply positional weights and sum across positions → (R, P)
    weighted = blosum_scores * weights[None, None, :]  # broadcast (1, 1, L)
    total_scores = weighted.sum(axis=2)                # (R, P)

    # Best IEDB match for each peptide
    best_scores = total_scores.max(axis=0)                    # (P,)

    # BLOSUM62 score of each peptide against itself at each position
    self_blosum = BLOSUM62_MATRIX[
        peptides,   # (P, L) — row index
        peptides    # (P, L) — column index, same amino acid
    ]                                           # (P, L)
    self_scores = (self_blosum * weights[None, :]).sum(axis=1)  # (P,)

    return best_scores / self_scores


# ---------------------------------------------------------------------------
# Worker function for multiprocessing
# ---------------------------------------------------------------------------

def _worker(args):
    """Process one length group; called in a subprocess."""
    length, peptide_strings, iedb_strings, magic_number = args
    weights      = get_sequence_weights(length, magic_number)
    peptides_enc = encode_sequences(peptide_strings)
    iedb_enc     = encode_sequences(iedb_strings)
    best_scores  = score_length_group(peptides_enc, iedb_enc, weights)
    return list(zip(peptide_strings, best_scores.tolist()))


# ---------------------------------------------------------------------------
# Main scoring entry point
# ---------------------------------------------------------------------------

def score_all_peptides(peptides, iedb_sequences, magic_number=4, n_workers=None):
    """
    Score all peptides against the IEDB reference, grouped by length so that
    each length group is processed as a single vectorised batch.

    Returns a dict mapping peptide string → best similarity score (or None).
    """
    if n_workers is None:
        n_workers = cpu_count()

    # Group queries and references by length
    peptides_by_length = defaultdict(list)
    for p in peptides:
        peptides_by_length[len(p)].append(p)

    iedb_by_length = defaultdict(list)
    for s in iedb_sequences:
        iedb_by_length[len(s)].append(s)

    # Build work units — only lengths that appear in both sets
    jobs = []
    no_match_lengths = set()
    for length, peps in peptides_by_length.items():
        if length not in iedb_by_length:
            no_match_lengths.add(length)
            continue
        jobs.append((length, peps, iedb_by_length[length], magic_number))

    results = {}

    # Peptides with no length-matched IEDB entries
    for length in no_match_lengths:
        for p in peptides_by_length[length]:
            results[p] = None

    # Parallelise across length groups
    with Pool(processes=min(n_workers, len(jobs)) if jobs else 1) as pool:
        for group_results in pool.imap_unordered(_worker, jobs):
            for peptide, score in group_results:
                results[peptide] = score

    return results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Score peptides by similarity to an IEDB pathogen reference database."
    )
    parser.add_argument("--peptides",    required=True, help="Input file: one peptide per line.")
    parser.add_argument("--iedb",        required=True, help="IEDB reference FASTA file.")
    parser.add_argument("--output",      required=True, help="Output CSV: peptide, similarity_score.")
    parser.add_argument("--magic-number",type=int, default=4, help="Weight scaling factor (default: 4).")
    parser.add_argument("--workers",     type=int, default=None, help="Number of CPU cores to use (default: all).")
    args = parser.parse_args()

    print(f"Loading IEDB reference: {args.iedb}")
    iedb_sequences = read_iedb(args.iedb)
    print(f"  {len(iedb_sequences)} sequences loaded.")

    print(f"Loading peptides: {args.peptides}")
    peptides = read_peptides(args.peptides)
    print(f"  {len(peptides)} peptides loaded.")

    print(f"Scoring with {args.workers or cpu_count()} worker(s)...")
    results = score_all_peptides(
        peptides,
        iedb_sequences,
        magic_number=args.magic_number,
        n_workers=args.workers,
    )

    scored = sum(1 for v in results.values() if v is not None)
    print(f"  Scored {scored}/{len(peptides)} peptides successfully.")

    with open(args.output, "w") as out:
        out.write("peptide,similarity_score\n")
        for peptide in peptides:  # preserve input order
            score = results[peptide]
            out.write(f"{peptide},{score if score is not None else 'NA'}\n")

    print(f"Results written to {args.output}")


if __name__ == "__main__":
    main()