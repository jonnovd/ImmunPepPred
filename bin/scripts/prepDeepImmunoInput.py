#!/usr/bin/env python
"""
prep_deepimmuno_input.py

Builds the input CSV for DeepImmuno from a peptide file and an allele file.

Usage:
    python prep_deepimmuno_input.py -p peptide_file -a allele_file --output deepimmuno_in.csv

Input files:
    - peptide_file: one peptide per line (streamed, not loaded fully into memory)
    - allele_file:  one allele per line, e.g. "HLA-C07:04" (read into memory up front)

Output:
    - CSV with no header, one row per (peptide, allele) pair, in the form:
        pep,allele
      where allele has been reformatted from "HLA-C07:04" to "HLA-C*07:04"
"""

import argparse
import csv
import re
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Prepare DeepImmuno input CSV from peptide and allele files."
    )
    parser.add_argument(
        "-p", "--peptide_file", required=True,
        help="Path to input file containing one peptide per line."
    )
    parser.add_argument(
        "-a", "--allele_file", required=True,
        help="Path to input file containing one allele per line."
    )
    parser.add_argument(
        "--output", required=True,
        help="Path to output CSV file."
    )
    return parser.parse_args()


def read_alleles(filepath):
    """Read the allele file fully into memory and return a list of reformatted alleles."""
    alleles = []
    with open(filepath, "r") as f:
        for line in f:
            allele = line.strip()
            if allele:
                alleles.append(reformat_allele(allele))
    return alleles


def reformat_allele(allele):
    """
    Convert allele format from e.g. 'HLA-C07:04' to 'HLA-C*0704'.

    Handles gene names of varying length (e.g. A, B, C, DRB1, DQB1, etc.)
    by inserting '*' between the trailing gene letters and the leading
    digits of the allele designation, and stripping any ':' separators
    from the digit portion.
    """
    match = re.match(r"^(HLA-[A-Za-z]+)(\d.*)$", allele)
    if not match:
        raise ValueError(f"Allele '{allele}' does not match expected format 'HLA-<gene><digits>...'")
    gene_part, number_part = match.groups()
    number_part = number_part.replace(":", "")
    return f"{gene_part}*{number_part}"


def main():
    args = parse_args()

    # Read allele file fully into memory first, since it's expected to be small.
    alleles = read_alleles(args.allele_file)
    if not alleles:
        sys.exit(f"Error: no alleles found in {args.allele_file}")

    row_count = 0

    # Stream the peptide file line by line and write rows as we go,
    # so the peptide file is never fully loaded into memory.
    with open(args.peptide_file, "r") as pep_f, open(args.output, "w", newline="") as out_f:
        writer = csv.writer(out_f)
        for line in pep_f:
            pep = line.strip()
            if not pep:
                continue
            for allele in alleles:
                writer.writerow([pep, allele])
                row_count += 1

    if row_count == 0:
        sys.exit(f"Error: no peptides found in {args.peptide_file}")

    print(f"Wrote {row_count} rows to {args.output}")


if __name__ == "__main__":
    main()