#!/usr/bin/env python3
"""
sample_peptides.py
------------------
Randomly samples N peptides from an input file, which can be either:

  1. A FASTA file of protein ORFs — in this mode, each sampled peptide is a
     random contiguous window cut from a randomly chosen source sequence.

  2. A plain text file with one peptide per line — in this mode, peptides
     are drawn (with replacement) directly from the lines that already
     match each target length. No windowing/cutting is performed.

In both modes, the LENGTH_RATIOS constant controls how many of the N
peptides are drawn at each length defined in PEPTIDE_LENGTHS.

The input type is auto-detected: if the file's first non-blank line starts
with '>', it's treated as FASTA; otherwise it's treated as a plain peptide
list (one sequence per line). You can also force the mode with --mode.

Usage:
    python sample_peptides.py <input_file> <n_peptides> [output.fasta] [--mode fasta|txt]

Arguments:
    input_file    Path to a FASTA file of ORFs, or a .txt file with one
                  peptide per line.
    n_peptides    Total number of peptides to sample.
    output.fasta  (Optional) Output file path. Defaults to 'sampled_peptides.fasta'.
    --mode        (Optional) Force 'fasta' or 'txt' input parsing instead of
                  auto-detecting from the file contents.
"""

import sys
import random
import math

# =============================================================================
# CONFIGURATION — edit these values to control peptide lengths and their ratios
# =============================================================================

# Peptide lengths to sample (in amino acids).
PEPTIDE_LENGTHS: list[int] = [8, 9, 10, 11]

# Relative sampling ratios for each length in PEPTIDE_LENGTHS (same order).
# These need not sum to 1 or 100 — they are normalised automatically.
# Example below draws equal numbers of each length.
LENGTH_RATIOS: list[float] = [350, 6100, 2400, 800]

# Random seed for reproducibility. Set to None for a truly random run.
RANDOM_SEED: int | None = 42

# =============================================================================
# END OF CONFIGURATION
# =============================================================================


def validate_config() -> None:
    if len(PEPTIDE_LENGTHS) != len(LENGTH_RATIOS):
        raise ValueError(
            f"PEPTIDE_LENGTHS has {len(PEPTIDE_LENGTHS)} entries but "
            f"LENGTH_RATIOS has {len(LENGTH_RATIOS)}. They must be the same length."
        )
    if any(r < 0 for r in LENGTH_RATIOS):
        raise ValueError("All values in LENGTH_RATIOS must be >= 0.")
    if sum(LENGTH_RATIOS) == 0:
        raise ValueError("At least one value in LENGTH_RATIOS must be > 0.")
    if any(l < 1 for l in PEPTIDE_LENGTHS):
        raise ValueError("All values in PEPTIDE_LENGTHS must be >= 1.")


def detect_input_mode(path: str) -> str:
    """Sniff the file to decide whether it's FASTA or a plain peptide list."""
    with open(path, "r") as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line:
                continue
            return "fasta" if line.startswith(">") else "txt"
    # Empty file — default to fasta parsing (will just yield 0 records).
    return "fasta"


def parse_fasta(path: str) -> list[tuple[str, str]]:
    """Parse a FASTA file, returning a list of (header, sequence) tuples."""
    records: list[tuple[str, str]] = []
    header: str | None = None
    seq_parts: list[str] = []

    with open(path, "r") as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seq = "".join(seq_parts)
                    if seq:
                        records.append((header, seq))
                    else:
                        print(f"  Warning: empty sequence for '{header}', skipping.",
                              file=sys.stderr)
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line.replace(" ", "").replace("\t", ""))

    if header is not None:
        seq = "".join(seq_parts)
        if seq:
            records.append((header, seq))

    return records


def parse_peptide_list(path: str) -> list[str]:
    """Parse a plain text file with one peptide sequence per line."""
    peptides: list[str] = []
    with open(path, "r") as fh:
        for raw_line in fh:
            line = raw_line.strip().replace(" ", "").replace("\t", "")
            if line:
                peptides.append(line)
    return peptides


def compute_sample_counts(n_total: int, ratios: list[float]) -> list[int]:
    """
    Divide n_total into per-length counts according to normalised ratios.
    Uses largest-remainder method to ensure counts sum exactly to n_total.
    """
    total_ratio = sum(ratios)
    exact = [(r / total_ratio) * n_total for r in ratios]
    floors = [math.floor(e) for e in exact]
    remainders = [(exact[i] - floors[i], i) for i in range(len(ratios))]

    shortfall = n_total - sum(floors)
    # Award remaining slots to the groups with the largest fractional parts
    for _, i in sorted(remainders, reverse=True)[:shortfall]:
        floors[i] += 1

    return floors


def sample_peptide_window(sequence: str, length: int, rng: random.Random) -> str:
    """Return a single random contiguous subsequence of the given length."""
    max_start = len(sequence) - length
    start = rng.randint(0, max_start)
    return sequence[start : start + length]


def write_fasta(records: list[tuple[str, str]], path: str, line_width: int = 60) -> None:
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for start in range(0, len(seq), line_width):
                fh.write(seq[start : start + line_width] + "\n")


def write_txt(records: list[tuple[str, str]], path: str, line_width: int = 60) -> None:
    with open(path, "w") as fh:
        for _, seq in records:
            fh.write(f'{seq}\n')


def print_length_table(sample_counts: list[int]) -> None:
    total_ratio = sum(LENGTH_RATIOS)
    norm_ratios = [r / total_ratio for r in LENGTH_RATIOS]

    print(f"\n{'Length (aa)':<14} {'Ratio':>8} {'Count':>8}")
    print("-" * 33)
    for length, ratio, count in zip(PEPTIDE_LENGTHS, norm_ratios, sample_counts):
        print(f"{length:<14} {ratio:>8.2%} {count:>8,}")
    print("-" * 33)
    print(f"{'TOTAL':<14} {'':>8} {sum(sample_counts):>8,}\n")


def sample_from_fasta(
    all_records: list[tuple[str, str]],
    sample_counts: list[int],
    rng: random.Random,
) -> list[tuple[str, str]]:
    """Sample peptides by cutting random windows out of FASTA source sequences."""
    sampled: list[tuple[str, str]] = []

    for length, count in zip(PEPTIDE_LENGTHS, sample_counts):
        # Only sequences long enough to yield a peptide of this length are eligible
        eligible = [rec for rec in all_records if len(rec[1]) >= length]

        if not eligible:
            print(
                f"  Warning: no sequences are long enough to yield peptides of "
                f"length {length} aa. Skipping {count} peptide(s).",
                file=sys.stderr,
            )
            continue

        for i in range(count):
            header, seq = rng.choice(eligible)
            peptide_seq = sample_peptide_window(seq, length, rng)
            # Label: source protein | peptide length | sample index (1-based)
            peptide_header = f"peptide_len{length}_{i+1} source={header}"
            sampled.append((peptide_header, peptide_seq))

    return sampled


def sample_from_peptide_list(
    all_peptides: list[str],
    sample_counts: list[int],
    rng: random.Random,
) -> list[tuple[str, str]]:
    """
    Sample peptides directly from an existing list, grouped by exact length.
    Sampling is WITHOUT replacement within each length group — each peptide
    from the input file can appear at most once in the output.
    """
    by_length: dict[int, list[str]] = {}
    for pep in all_peptides:
        by_length.setdefault(len(pep), []).append(pep)

    sampled: list[tuple[str, str]] = []

    for length, count in zip(PEPTIDE_LENGTHS, sample_counts):
        candidates = by_length.get(length, [])

        if not candidates:
            print(
                f"  Warning: no peptides of length {length} aa found in the input "
                f"list. Skipping {count} peptide(s).",
                file=sys.stderr,
            )
            continue

        if count > len(candidates):
            print(
                f"  Warning: requested {count} peptide(s) of length {length} aa "
                f"but only {len(candidates)} available in the input list "
                f"(sampling without replacement). Using all {len(candidates)}.",
                file=sys.stderr,
            )
            count = len(candidates)

        chosen = rng.sample(candidates, count)
        for i, peptide_seq in enumerate(chosen):
            peptide_header = f"peptide_len{length}_{i+1}"
            sampled.append((peptide_header, peptide_seq))

    return sampled


def main() -> None:
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    forced_mode = None
    if "--mode" in sys.argv:
        idx = sys.argv.index("--mode")
        if idx + 1 < len(sys.argv):
            forced_mode = sys.argv[idx + 1]
            if forced_mode not in ("fasta", "txt"):
                print(f"Error: --mode must be 'fasta' or 'txt', got '{forced_mode}'.", file=sys.stderr)
                sys.exit(1)

    if len(args) < 2:
        print(__doc__)
        sys.exit(1)

    input_path = args[0]
    try:
        n_peptides = int(args[1])
    except ValueError:
        print(f"Error: n_peptides must be an integer, got '{args[1]}'.", file=sys.stderr)
        sys.exit(1)
    if n_peptides <= 0:
        print("Error: n_peptides must be a positive integer.", file=sys.stderr)
        sys.exit(1)

    output_path = args[2] if len(args) >= 3 else "sampled_peptides.fasta"

    validate_config()

    rng = random.Random(RANDOM_SEED)
    if RANDOM_SEED is not None:
        print(f"Random seed: {RANDOM_SEED}")

    # ------------------------------------------------------------------
    # Determine input mode and parse
    # ------------------------------------------------------------------
    mode = forced_mode or detect_input_mode(input_path)
    print(f"Input mode: {mode}")
    print(f"Parsing '{input_path}' …")

    if mode == "fasta":
        all_records = parse_fasta(input_path)
        print(f"  {len(all_records)} sequences loaded.")
    else:
        all_peptides = parse_peptide_list(input_path)
        print(f"  {len(all_peptides)} peptides loaded.")

    # ------------------------------------------------------------------
    # Compute how many peptides to draw at each length
    # ------------------------------------------------------------------
    sample_counts = compute_sample_counts(n_peptides, LENGTH_RATIOS)
    print_length_table(sample_counts)

    # ------------------------------------------------------------------
    # Sample peptides
    # ------------------------------------------------------------------
    if mode == "fasta":
        sampled = sample_from_fasta(all_records, sample_counts, rng)
    else:
        sampled = sample_from_peptide_list(all_peptides, sample_counts, rng)

    # ------------------------------------------------------------------
    # Write output
    # ------------------------------------------------------------------
    write_txt(sampled, output_path)
    print(f"Wrote {len(sampled)} sampled peptides to '{output_path}'.")


if __name__ == "__main__":
    main()



# #!/usr/bin/env python3
# """
# sample_peptides.py
# ------------------
# Randomly samples N peptides (contiguous subsequences) from protein ORFs
# in a FASTA file. Each sampled peptide is a random window cut from a
# randomly chosen source sequence.

# The LENGTH_RATIOS constant controls how many of the N peptides will be
# drawn at each peptide length defined in PEPTIDE_LENGTHS.

# Usage:
#     python sample_peptides.py <input.fasta> <n_peptides> [output.fasta]

# Arguments:
#     input.fasta   Path to the input FASTA file containing protein ORFs.
#     n_peptides    Total number of peptides to sample.
#     output.fasta  (Optional) Output file path. Defaults to 'sampled_peptides.fasta'.
# """

# import sys
# import random
# import math

# # =============================================================================
# # CONFIGURATION — edit these values to control peptide lengths and their ratios
# # =============================================================================

# # Peptide lengths to sample (in amino acids).
# PEPTIDE_LENGTHS: list[int] = [8, 9, 10, 11]

# # Relative sampling ratios for each length in PEPTIDE_LENGTHS (same order).
# # These need not sum to 1 or 100 — they are normalised automatically.
# # Example below draws equal numbers of each length.
# LENGTH_RATIOS: list[float] = [350, 6100, 2400, 800]

# # Random seed for reproducibility. Set to None for a truly random run.
# RANDOM_SEED: int | None = 42

# # =============================================================================
# # END OF CONFIGURATION
# # =============================================================================


# def validate_config() -> None:
#     if len(PEPTIDE_LENGTHS) != len(LENGTH_RATIOS):
#         raise ValueError(
#             f"PEPTIDE_LENGTHS has {len(PEPTIDE_LENGTHS)} entries but "
#             f"LENGTH_RATIOS has {len(LENGTH_RATIOS)}. They must be the same length."
#         )
#     if any(r < 0 for r in LENGTH_RATIOS):
#         raise ValueError("All values in LENGTH_RATIOS must be >= 0.")
#     if sum(LENGTH_RATIOS) == 0:
#         raise ValueError("At least one value in LENGTH_RATIOS must be > 0.")
#     if any(l < 1 for l in PEPTIDE_LENGTHS):
#         raise ValueError("All values in PEPTIDE_LENGTHS must be >= 1.")


# def parse_fasta(path: str) -> list[tuple[str, str]]:
#     """Parse a FASTA file, returning a list of (header, sequence) tuples."""
#     records: list[tuple[str, str]] = []
#     header: str | None = None
#     seq_parts: list[str] = []

#     with open(path, "r") as fh:
#         for raw_line in fh:
#             line = raw_line.strip()
#             if not line:
#                 continue
#             if line.startswith(">"):
#                 if header is not None:
#                     seq = "".join(seq_parts)
#                     if seq:
#                         records.append((header, seq))
#                     else:
#                         print(f"  Warning: empty sequence for '{header}', skipping.",
#                               file=sys.stderr)
#                 header = line[1:]
#                 seq_parts = []
#             else:
#                 seq_parts.append(line.replace(" ", "").replace("\t", ""))

#     if header is not None:
#         seq = "".join(seq_parts)
#         if seq:
#             records.append((header, seq))

#     return records


# def compute_sample_counts(n_total: int, ratios: list[float]) -> list[int]:
#     """
#     Divide n_total into per-length counts according to normalised ratios.
#     Uses largest-remainder method to ensure counts sum exactly to n_total.
#     """
#     total_ratio = sum(ratios)
#     exact = [(r / total_ratio) * n_total for r in ratios]
#     floors = [math.floor(e) for e in exact]
#     remainders = [(exact[i] - floors[i], i) for i in range(len(ratios))]

#     shortfall = n_total - sum(floors)
#     # Award remaining slots to the groups with the largest fractional parts
#     for _, i in sorted(remainders, reverse=True)[:shortfall]:
#         floors[i] += 1

#     return floors


# def sample_peptide(sequence: str, length: int, rng: random.Random) -> str:
#     """Return a single random contiguous subsequence of the given length."""
#     max_start = len(sequence) - length
#     start = rng.randint(0, max_start)
#     return sequence[start : start + length]


# def write_fasta(records: list[tuple[str, str]], path: str, line_width: int = 60) -> None:
#     with open(path, "w") as fh:
#         for header, seq in records:
#             fh.write(f">{header}\n")
#             for start in range(0, len(seq), line_width):
#                 fh.write(seq[start : start + line_width] + "\n")

# def write_txt(records: list[tuple[str, str]], path: str, line_width: int = 60) -> None:
#     with open(path, "w") as fh:
#         for _, seq in records:
#             fh.write(f'{seq}\n')

# def main() -> None:
#     if len(sys.argv) < 3:
#         print(__doc__)
#         sys.exit(1)

#     input_path = sys.argv[1]
#     try:
#         n_peptides = int(sys.argv[2])
#     except ValueError:
#         print(f"Error: n_peptides must be an integer, got '{sys.argv[2]}'.", file=sys.stderr)
#         sys.exit(1)
#     if n_peptides <= 0:
#         print("Error: n_peptides must be a positive integer.", file=sys.stderr)
#         sys.exit(1)

#     output_path = sys.argv[3] if len(sys.argv) >= 4 else "sampled_peptides.fasta"

#     validate_config()

#     rng = random.Random(RANDOM_SEED)
#     if RANDOM_SEED is not None:
#         print(f"Random seed: {RANDOM_SEED}")

#     # ------------------------------------------------------------------
#     # Parse input
#     # ------------------------------------------------------------------
#     print(f"Parsing '{input_path}' …")
#     all_records = parse_fasta(input_path)
#     print(f"  {len(all_records)} sequences loaded.")

#     # ------------------------------------------------------------------
#     # Compute how many peptides to draw at each length
#     # ------------------------------------------------------------------
#     sample_counts = compute_sample_counts(n_peptides, LENGTH_RATIOS)

#     max_peptide_len = max(PEPTIDE_LENGTHS)
#     total_ratio = sum(LENGTH_RATIOS)
#     norm_ratios  = [r / total_ratio for r in LENGTH_RATIOS]

#     print(f"\n{'Length (aa)':<14} {'Ratio':>8} {'Count':>8}")
#     print("-" * 33)
#     for length, ratio, count in zip(PEPTIDE_LENGTHS, norm_ratios, sample_counts):
#         print(f"{length:<14} {ratio:>8.2%} {count:>8,}")
#     print("-" * 33)
#     print(f"{'TOTAL':<14} {'':>8} {sum(sample_counts):>8,}\n")

#     # ------------------------------------------------------------------
#     # Sample peptides for each length group
#     # ------------------------------------------------------------------
#     sampled: list[tuple[str, str]] = []

#     for length, count in zip(PEPTIDE_LENGTHS, sample_counts):
#         # Only sequences long enough to yield a peptide of this length are eligible
#         eligible = [rec for rec in all_records if len(rec[1]) >= length]

#         if not eligible:
#             print(
#                 f"  Warning: no sequences are long enough to yield peptides of "
#                 f"length {length} aa. Skipping {count} peptide(s).",
#                 file=sys.stderr,
#             )
#             continue

#         for i in range(count):
#             header, seq = rng.choice(eligible)
#             peptide_seq  = sample_peptide(seq, length, rng)
#             # Label: source protein | peptide length | sample index (1-based)
#             peptide_header = f"peptide_len{length}_{i+1} source={header}"
#             sampled.append((peptide_header, peptide_seq))

#     # Shuffle so lengths are interleaved in the output file
#     # rng.shuffle(sampled)

#     # ------------------------------------------------------------------
#     # Write output
#     # ------------------------------------------------------------------
#     write_txt(sampled, output_path)
#     print(f"Wrote {len(sampled)} sampled peptides to '{output_path}'.")


# if __name__ == "__main__":
#     main()