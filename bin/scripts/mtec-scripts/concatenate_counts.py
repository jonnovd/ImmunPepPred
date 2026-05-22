"""
concatenate_counts.py

Concatenates the count columns from multiple peptide_counter TSV output files
into a single output file.

Each input file is expected to have the same peptides in the same order
(as they all derive from the same peptide input file).

Usage:
    python concatenate_counts.py -i <file1> <file2> ... -o <output_file>
    python concatenate_counts.py -i *.tsv -o combined_counts.tsv

Input format (each file):
    peptide<tab>count

Output format:
    peptide<tab>count_file1<tab>count_file2<tab>...

Requirements:
    None (stdlib only)
"""

import sys
import time
import argparse


# ── Helpers ────────────────────────────────────────────────────────────────────

def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


# ── Concatenate ────────────────────────────────────────────────────────────────

def concatenate_counts(input_paths: list[str], output_path: str) -> None:
    """
    Streams all input files simultaneously line by line, extracting the count
    column from each and writing a single merged output line per peptide.

    Uses a single open file handle per input — never loads any full file
    into memory. Peak memory usage is O(number of input files), not O(peptides).
    """
    n_files = len(input_paths)
    log(f"Opening {n_files} input files...")

    # Open all input files simultaneously
    handles = []
    try:
        for path in input_paths:
            handles.append(open(path, "r"))

        log(f"Writing merged output to: {output_path}")

        with open(output_path, "w", buffering=8 * 1024 * 1024) as out:

            # Write header: one column per input file plus a total column
            header_cols = "\t".join(f"count_{i+1}" for i in range(n_files))
            out.write(f"peptide\t{header_cols}\ttotal\n")

            line_num = 0

            # Iterate all files in lockstep — one line at a time from each
            for lines in zip(*handles):
                line_num += 1

                # Parse each file's line: peptide<tab>count
                parts = [line.rstrip("\n").split() for line in lines]

                peptide = parts[0][0]

                # Extract counts, compute total, and write merged line
                count_vals = [int(part[1]) for part in parts]
                total      = sum(count_vals)
                counts_str = "\t".join(str(c) for c in count_vals)
                out.write(f"{peptide}\t{counts_str}\t{total}\n")

                if line_num % 1_000_000 == 0:
                    log(f"  Processed {line_num // 1_000_000}M peptides...")

    finally:
        for h in handles:
            h.close()

    log(f"  Done — wrote {line_num:,} peptide rows "
        f"({n_files} count columns)")

def get_mtec_expression_file(inFilePath: str, outFilePath: str, threshold: int):
    with open(inFilePath, 'r') as inFile:
        with open(outFilePath, "w", buffering=8 * 1024 * 1024) as out:
            for line in inFile:
                line = line.strip()
                if line:
                    if line.startswith('p'):
                        out.write("peptide,mtec_expression\n")
                    else:
                        parts = line.split()
                        pep   = parts[0]
                        total = parts[-1]
                        if total > threshold:
                            out.write(f"{pep},{1}\n")
                        else:
                            out.write(f"{pep},{0}\n")
                        
# ── Main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Concatenate count columns from multiple peptide counter TSV files"
    )
    parser.add_argument("-i", "--input_files",
                        nargs="+",
                        help="Two or more TSV files (peptide<tab>count)")
    parser.add_argument("-o", "--output_file",
                        help="Output TSV: peptide<tab>count_1<tab>count_2<tab>...total")
    parser.add_argument("-c", "--output_expression_file",
                        help="Output CSV: peptide,mTEC_expression")
    parser.add_argument("-t", "--threshold",
                        help="Count Threshold to determine mTEC expression: int")
    args = parser.parse_args()

    if len(args.input_files) < 2:
        log("ERROR: at least 2 input files required")
        sys.exit(1)

    total_start = time.time()
    log("=== Peptide Count Concatenator ===")
    for i, path in enumerate(args.input_files, 1):
        log(f"  Input {i}: {path}")

    concatenate_counts(args.input_files, args.output_file)

    log(f"=== Total wall time: {time.time() - total_start:.1f}s ===")

    get_mtec_expression_file(args.output_file, args.output_expression_file, args.threshold)


if __name__ == "__main__":
    main()
