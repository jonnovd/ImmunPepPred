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

def read_reference(fasta_file):
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

def get_sequence_weights(length, magic_number=1):
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

def score_length_group(peptides, reference_encoded, weights, chunk_size=256):
    P, L = peptides.shape

    self_blosum = BLOSUM62_MATRIX[peptides, peptides]
    self_scores = (self_blosum * weights[None, :]).sum(axis=1).astype(np.float32)

    best_scores = np.full(P, np.nan, dtype=np.float32)  # nan = "no match < 1 found"

    ref_blosum_weighted = BLOSUM62_MATRIX[reference_encoded] * weights[None, :, None]

    for start in range(0, P, chunk_size):
        chunk = peptides[start:start + chunk_size]  # (C, L)
        chunk_self = self_scores[start:start + chunk_size]  # (C,)

        scores = np.tensordot(
            ref_blosum_weighted,
            np.eye(20, dtype=np.int32)[chunk].transpose(1, 0, 2),
            axes=([1, 2], [0, 2])
        )  # (R, C) raw (unnormalized) similarity scores

        # Mask out any ref score that would normalize to >= 1 for that peptide
        # (i.e. self-matches or ties), leaving only strictly-better-than-self exclusions
        below_identity = scores < chunk_self[None, :]
        masked = np.where(below_identity, scores, -np.inf)

        chunk_best = masked.max(axis=0)  # (C,), -inf where nothing qualified
        has_match = np.isfinite(chunk_best)
        best_scores[start:start + chunk_size][has_match] = chunk_best[has_match]

    return best_scores / self_scores


# ---------------------------------------------------------------------------
# Worker function for multiprocessing
# ---------------------------------------------------------------------------

def _worker(args):
    """Process one length group; called in a subprocess."""
    length, peptide_strings, ref_strings, magic_number = args
    weights      = get_sequence_weights(length, magic_number)
    peptides_enc = encode_sequences(peptide_strings)
    ref_enc     = encode_sequences(ref_strings)
    best_scores  = score_length_group(peptides_enc, ref_enc, weights)
    return list(zip(peptide_strings, best_scores.tolist()))


# ---------------------------------------------------------------------------
# Main scoring entry point
# ---------------------------------------------------------------------------

def score_all_peptides(peptides, reference_sequences, magic_number=4, n_workers=None):
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

    reference_by_length = defaultdict(list)
    for s in reference_sequences:
        reference_by_length[len(s)].append(s)

    # Build work units — only lengths that appear in both sets
    jobs = []
    no_match_lengths = set()
    for length, peps in peptides_by_length.items():
        if length not in reference_by_length:
            no_match_lengths.add(length)
            continue
        jobs.append((length, peps, reference_by_length[length], magic_number))

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
    parser.add_argument("--reference",        required=True, help="Reference FASTA file.")
    parser.add_argument("--output",      required=True, help="Output CSV: peptide, similarity_score.")
    parser.add_argument("--magic-number",type=int, default=4, help="Weight scaling factor (default: 4).")
    parser.add_argument("--workers",     type=int, default=None, help="Number of CPU cores to use (default: all).")
    args = parser.parse_args()

    print(f"Loading Reference reference: {args.reference}")
    reference_sequences = read_reference(args.reference)
    print(f"  {len(reference_sequences)} sequences loaded.")

    print(f"Loading peptides: {args.peptides}")
    peptides = read_peptides(args.peptides)
    print(f"  {len(peptides)} peptides loaded.")

    print(f"Scoring with {args.workers or cpu_count()} worker(s)...")
    results = score_all_peptides(
        peptides,
        reference_sequences,
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