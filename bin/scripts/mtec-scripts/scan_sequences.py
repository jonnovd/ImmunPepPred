"""
scan_sequences.py

Loads a pre-built Aho-Corasick automaton and counts peptide occurrences
across a sequences file.

Usage:
    python scan_sequences.py -a <automaton_file> -p <peptides_file>
                             -s <sequences_file> -o <output_file>

Input formats:
    automaton_file:  binary automaton saved by build_automaton.py
    peptides_file:   one peptide per line — same file used to build the automaton
                     (used only for writing output, not for searching)
    sequences_file:  one AA sequence per line (97 million sequences, 26 AA length)
    output_file:     TSV: peptide<tab>count

Requirements:
    pip install pyahocorasick
"""

import sys
import time
import pickle
import argparse
import ahocorasick


# ── Constants ──────────────────────────────────────────────────────────────────

REPORT_EVERY = 10_000_000


# ── Helpers ────────────────────────────────────────────────────────────────────

def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)

def get_num_peps(peptides_file: str) -> int:
    count = 0
    with open(peptides_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                count += 1

    return count


# ── Load ───────────────────────────────────────────────────────────────────────

def load_automaton(path: str) -> ahocorasick.Automaton:
    """
    Deserialise automaton from disk using pyahocorasick's native load().

    pickle.loads is passed as the value deserialiser — must match
    the serialiser used in build_automaton.py (pickle.dumps).
    """
    log(f"Loading automaton from: {path}")
    A = ahocorasick.load(path, pickle.loads)
    log(f"  Loaded — {len(A):,} peptides, {A.get_stats()['nodes_count']:,} nodes")
    return A


# ── Scan ───────────────────────────────────────────────────────────────────────

def scan_sequences(path: str, automaton: ahocorasick.Automaton, numPeps: int) -> list[int]:
    """
    Stream sequences file line by line, running each sequence through
    the automaton. Returns a counts array indexed by peptide_id (0-based).

    automaton.iter(string) yields (end_index, pep_id) for every match,
    automatically following output/dictionary links — so overlapping
    and suffix patterns are all reported correctly.
    """
    log(f"Scanning sequences from: {path}")
    counts = [0] * numPeps   # 0-based: counts[pep_id - 1]
    seq_count = 0
    match_count = 0

    with open(path, "r") as f:
        for line in f:
            seq = line.strip().upper()
            if not seq:
                continue

            for _end_idx, pep_id in automaton.iter(seq):
                counts[pep_id - 1] += 1
                match_count += 1

            seq_count += 1
            if seq_count % REPORT_EVERY == 0:
                log(f"  Scanned {seq_count // 1_000_000}M sequences "
                    f"({match_count:,} matches so far)...")

    log(f"  Scanned {seq_count:,} sequences total")
    log(f"  Total matches: {match_count:,}")
    return counts


# ── Write output ───────────────────────────────────────────────────────────────

def write_output(inpath: str, outpath: str, counts: list[int]) -> None:
    """
    Write TSV output: peptide<tab>count, one line per peptide.
    Reads peptide strings from the original peptide file to recover names,
    preserving insertion order (which matches the 1-based automaton ids).
    """
    log(f"Writing output to: {outpath}")

    with open(inpath, "r") as inFile:
        with open(outpath, "w") as f:
            id = 1
            for line in inFile:
                pep = line.strip()
                if pep:
                    f.write(f"{pep}\t{counts[id - 1]}\n")
                    id += 1

    log(f"  Wrote {id - 1:,} lines")


# ── Main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Scan sequences using a pre-built Aho-Corasick automaton"
    )
    parser.add_argument("-a", "--automaton_file",  help="Binary automaton file from build_automaton.py")
    parser.add_argument("-p", "--peptides_file",   help="Original peptide file (for output labels)")
    parser.add_argument("-s", "--sequences_file",  help="One AA sequence per line")
    parser.add_argument("-o", "--output_file",     help="TSV output: peptide<tab>count")
    args = parser.parse_args()

    numPeps = get_num_peps(args.peptides_file)

    total_start = time.time()
    log("=== Sequence Scanner ===")

    t0 = time.time()
    automaton = load_automaton(args.automaton_file)
    log(f"  Load done in {time.time() - t0:.1f}s")

    t0 = time.time()
    counts = scan_sequences(args.sequences_file, automaton, numPeps)
    log(f"  Scan done in {time.time() - t0:.1f}s")

    t0 = time.time()
    write_output(args.peptides_file, args.output_file, counts)
    log(f"  Write done in {time.time() - t0:.1f}s")

    log(f"=== Total wall time: {time.time() - total_start:.1f}s ===")


if __name__ == "__main__":
    main()
