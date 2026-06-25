#!/usr/bin/env python3
"""
sample_peptides.py
------------------
Randomly samples N peptides from a protein/polypeptide FASTA file,
drawing from defined length groups according to user-specified ratios.

Usage:
    python sample_peptides.py <input.fasta> <n_peptides> [output.fasta]

Arguments:
    input.fasta   Path to the input FASTA file containing protein ORFs.
    n_peptides    Total number of peptides to sample.
    output.fasta  (Optional) Output file path. Defaults to 'sampled_peptides.fasta'.
"""

import sys
import random
import math
from collections import defaultdict

# =============================================================================
# CONFIGURATION — edit these values to control length groups and their ratios
# =============================================================================

# Define peptide length groups as (min_aa, max_aa) inclusive ranges.
# Sequences shorter than the first group's min or longer than the last
# group's max are silently excluded from sampling.
LENGTH_GROUPS: list[tuple[int, int]] = [
    (8,   8),   # Short peptides       (e.g. 7–11 aa)
    (9,  9),   # Medium-short         (e.g. 12–20 aa)
    (10,  10),   # Medium               (e.g. 21–50 aa)
    (11, 11),   # Medium-long          (e.g. 51–100 aa)
]

# Relative sampling ratios for each length group (must match LENGTH_GROUPS length).
# These need not sum to 1 or 100 — they are normalised automatically.
# Example below draws roughly: 10% short, 30% medium-short, 40% medium, 20% medium-long.
LENGTH_RATIOS: list[float] = [
    0.10,   # short
    0.30,   # medium-short
    0.40,   # medium
    0.20,   # medium-long
]

# Random seed for reproducibility. Set to None for a truly random run.
RANDOM_SEED: int | None = 42

# =============================================================================
# END OF CONFIGURATION
# =============================================================================


def validate_config() -> None:
    """Raise informative errors if the configuration constants are inconsistent."""
    if len(LENGTH_GROUPS) != len(LENGTH_RATIOS):
        raise ValueError(
            f"LENGTH_GROUPS has {len(LENGTH_GROUPS)} entries but "
            f"LENGTH_RATIOS has {len(LENGTH_RATIOS)}. They must be the same length."
        )
    if any(r < 0 for r in LENGTH_RATIOS):
        raise ValueError("All values in LENGTH_RATIOS must be >= 0.")
    if sum(LENGTH_RATIOS) == 0:
        raise ValueError("At least one value in LENGTH_RATIOS must be > 0.")
    for i, (lo, hi) in enumerate(LENGTH_GROUPS):
        if lo > hi:
            raise ValueError(
                f"LENGTH_GROUPS[{i}] has min={lo} > max={hi}. min must be <= max."
            )


def parse_fasta(path: str) -> list[tuple[str, str]]:
    """
    Parse a FASTA file and return a list of (header, sequence) tuples.
    Multi-line sequences are concatenated. Whitespace inside sequences
    is stripped. Empty sequences are skipped with a warning.
    """
    records: list[tuple[str, str]] = []
    header: str | None = None
    seq_parts: list[str] = []

    with open(path, "r") as fh:
        for line_no, raw_line in enumerate(fh, start=1):
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
                header = line[1:]   # strip the leading '>'
                seq_parts = []
            else:
                seq_parts.append(line.replace(" ", "").replace("\t", ""))

    # Flush the last record
    if header is not None:
        seq = "".join(seq_parts)
        if seq:
            records.append((header, seq))
        else:
            print(f"  Warning: empty sequence for '{header}', skipping.",
                  file=sys.stderr)

    return records


def assign_to_groups(
    records: list[tuple[str, str]],
    groups: list[tuple[int, int]],
) -> dict[int, list[tuple[str, str]]]:
    """
    Bucket each (header, sequence) record into the appropriate length group.

    Returns a dict keyed by group index (0-based), containing only records
    whose sequence length falls within that group's [min, max] range.
    Sequences that fall outside every group are counted and reported.
    """
    buckets: dict[int, list[tuple[str, str]]] = defaultdict(list)
    out_of_range = 0

    for header, seq in records:
        length = len(seq)
        placed = False
        for idx, (lo, hi) in enumerate(groups):
            if lo <= length <= hi:
                buckets[idx].append((header, seq))
                placed = True
                break
        if not placed:
            out_of_range += 1

    if out_of_range:
        print(
            f"  Note: {out_of_range} sequence(s) fell outside all length groups "
            "and were excluded from sampling.",
            file=sys.stderr,
        )

    return buckets


def compute_sample_counts(
    n_total: int,
    ratios: list[float],
    available: list[int],
) -> list[int]:
    """
    Allocate n_total samples across groups according to normalised ratios,
    respecting per-group availability caps.

    Uses an iterative approach: allocate by ratio, cap at availability,
    redistribute remainder to uncapped groups until stable.

    Returns a list of sample counts (same length as ratios).
    """
    n_groups = len(ratios)
    total_ratio = sum(ratios)
    norm = [r / total_ratio for r in ratios]

    counts = [0] * n_groups
    remaining = n_total
    free_groups = set(range(n_groups))   # groups not yet capped

    for _iteration in range(n_groups + 1):  # at most n_groups capping rounds
        if not free_groups or remaining == 0:
            break

        # Renormalise over free groups only
        free_ratio_sum = sum(norm[i] for i in free_groups)
        if free_ratio_sum == 0:
            break

        new_free: set[int] = set()
        allocated_this_round = 0

        for i in sorted(free_groups):
            share = (norm[i] / free_ratio_sum) * remaining
            alloc = min(math.floor(share), available[i])
            counts[i] = alloc
            allocated_this_round += alloc
            if alloc < available[i]:
                new_free.add(i)

        remaining -= allocated_this_round
        free_groups = new_free

    # Distribute any leftover (due to floor rounding) one-by-one
    # to groups that still have capacity, in ratio order (largest first)
    if remaining > 0:
        order = sorted(free_groups, key=lambda i: -norm[i])
        for i in order:
            if remaining == 0:
                break
            gap = available[i] - counts[i]
            add = min(gap, remaining)
            counts[i] += add
            remaining -= add

    if remaining > 0:
        print(
            f"  Warning: could only allocate {n_total - remaining} of {n_total} "
            "requested peptides (insufficient sequences in some length groups).",
            file=sys.stderr,
        )

    return counts


def write_fasta(records: list[tuple[str, str]], path: str, line_width: int = 60) -> None:
    """Write (header, sequence) records to a FASTA file."""
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for start in range(0, len(seq), line_width):
                fh.write(seq[start : start + line_width] + "\n")


def main() -> None:
    # ------------------------------------------------------------------
    # Argument parsing
    # ------------------------------------------------------------------
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    input_path  = sys.argv[1]
    try:
        n_peptides = int(sys.argv[2])
    except ValueError:
        print(f"Error: n_peptides must be an integer, got '{sys.argv[2]}'.", file=sys.stderr)
        sys.exit(1)
    if n_peptides <= 0:
        print("Error: n_peptides must be a positive integer.", file=sys.stderr)
        sys.exit(1)

    output_path = sys.argv[3] if len(sys.argv) >= 4 else "sampled_peptides.fasta"

    # ------------------------------------------------------------------
    # Validate configuration
    # ------------------------------------------------------------------
    validate_config()

    # ------------------------------------------------------------------
    # Seed RNG
    # ------------------------------------------------------------------
    random.seed(RANDOM_SEED)
    if RANDOM_SEED is not None:
        print(f"Random seed: {RANDOM_SEED}")

    # ------------------------------------------------------------------
    # Parse input
    # ------------------------------------------------------------------
    print(f"Parsing '{input_path}' …")
    records = parse_fasta(input_path)
    print(f"  {len(records)} sequences loaded.")

    # ------------------------------------------------------------------
    # Bucket into length groups
    # ------------------------------------------------------------------
    print("Assigning sequences to length groups …")
    buckets = assign_to_groups(records, LENGTH_GROUPS)

    total_ratio = sum(LENGTH_RATIOS)
    norm_ratios = [r / total_ratio for r in LENGTH_RATIOS]

    print(f"\n{'Group':<6} {'Range (aa)':<14} {'Ratio':>7} {'Available':>10} {'Target':>8}")
    print("-" * 50)

    available_counts = [len(buckets.get(i, [])) for i in range(len(LENGTH_GROUPS))]
    sample_counts    = compute_sample_counts(n_peptides, LENGTH_RATIOS, available_counts)

    for i, ((lo, hi), ratio, avail, target) in enumerate(
        zip(LENGTH_GROUPS, norm_ratios, available_counts, sample_counts)
    ):
        print(f"{i:<6} {f'{lo}–{hi}':<14} {ratio:>7.2%} {avail:>10,} {target:>8,}")

    total_sampled = sum(sample_counts)
    print("-" * 50)
    print(f"{'TOTAL':<6} {'':<14} {'':>7} {sum(available_counts):>10,} {total_sampled:>8,}\n")

    # ------------------------------------------------------------------
    # Sample
    # ------------------------------------------------------------------
    sampled: list[tuple[str, str]] = []
    for i, count in enumerate(sample_counts):
        pool = buckets.get(i, [])
        draw = random.sample(pool, count)
        sampled.extend(draw)

    # Shuffle the final list so groups are interleaved rather than blocked
    random.shuffle(sampled)

    # ------------------------------------------------------------------
    # Write output
    # ------------------------------------------------------------------
    write_fasta(sampled, output_path)
    print(f"Wrote {len(sampled)} sampled peptides to '{output_path}'.")


if __name__ == "__main__":
    main()