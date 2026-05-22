"""
build_automaton.py

Builds an Aho-Corasick automaton from a peptide file and saves it to disk
for reuse across multiple scanning jobs.

Usage:
    python build_automaton.py -p <peptides_file> -o <automaton_output_file>

Output:
    A binary automaton file loadable by scan_sequences.py

Requirements:
    pip install pyahocorasick
"""

import sys
import time
import pickle
import argparse
import ahocorasick


# ── Helpers ────────────────────────────────────────────────────────────────────

def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


# ── Build ──────────────────────────────────────────────────────────────────────

def build_automaton(path: str) -> tuple[ahocorasick.Automaton, int]:
    """
    Read peptides from file and build Aho-Corasick automaton in a single pass.

    Uses STORE_INTS mode so pyahocorasick automatically stores the 1-based
    insertion index at each terminal node.

    Returns:
        A            : compiled Automaton, ready for iter()
        peptide_count: total number of peptides inserted
    """
    log(f"Loading peptides and building automaton from: {path}")

    A = ahocorasick.Automaton(ahocorasick.STORE_INTS)

    with open(path, "r") as f:
        for line in f:
            pep = line.strip()
            if not pep:
                continue
            A.add_word(pep)

    peptide_count = len(A)
    log(f"  Loaded {peptide_count:,} peptides")

    A.make_automaton()
    log(f"  Automaton built — {A.get_stats()['nodes_count']:,} nodes")

    return A, peptide_count


# ── Save ───────────────────────────────────────────────────────────────────────

def save_automaton(A: ahocorasick.Automaton, path: str) -> None:
    """
    Serialise automaton to disk using pyahocorasick's native save().

    pickle.dumps is passed as the value serialiser — required by the API
    even for STORE_INTS mode where payloads are plain ints.
    """
    log(f"Saving automaton to: {path}")
    A.save(path)
    log(f"  Saved successfully")


# ── Main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build and save Aho-Corasick automaton from peptide file"
    )
    parser.add_argument("-p", "--peptides_file",   help="One peptide per line")
    parser.add_argument("-o", "--output_file",     help="Path to save automaton binary")
    args = parser.parse_args()

    total_start = time.time()
    log("=== Automaton Builder ===")

    t0 = time.time()
    automaton, peptide_count = build_automaton(args.peptides_file)
    log(f"  Build done in {time.time() - t0:.1f}s")

    t0 = time.time()
    save_automaton(automaton, args.output_file)
    log(f"  Save done in {time.time() - t0:.1f}s")

    log(f"  Total peptides : {peptide_count:,}")
    log(f"=== Total wall time: {time.time() - total_start:.1f}s ===")


if __name__ == "__main__":
    main()
