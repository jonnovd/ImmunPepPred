# #!/usr/bin/env python
# """
# prep_deepimmuno_input.py

# Builds the input CSV for DeepImmuno from a peptide file and an allele file.

# Usage:
#     python prep_deepimmuno_input.py -p peptide_file -a allele_file --output deepimmuno_in.csv

# Input files:
#     - peptide_file: one peptide per line (streamed, not loaded fully into memory)
#     - allele_file:  one allele per line, e.g. "HLA-C07:04" (read into memory up front)

# Output:
#     - CSV with no header, one row per (peptide, allele) pair, in the form:
#         pep,allele
#       where allele has been reformatted from "HLA-C07:04" to "HLA-C*07:04"
# """

# import argparse
# import csv
# import re
# import sys


# def parse_args():
#     parser = argparse.ArgumentParser(
#         description="Prepare DeepImmuno input CSV from peptide and allele files."
#     )
#     parser.add_argument(
#         "-p", "--peptide_file", required=True,
#         help="Path to input file containing one peptide per line."
#     )
#     parser.add_argument(
#         "-a", "--allele_file", required=True,
#         help="Path to input file containing one allele per line."
#     )
#     parser.add_argument(
#         "--output", required=True,
#         help="Path to output CSV file."
#     )
#     return parser.parse_args()


# def read_alleles(filepath):
#     """Read the allele file fully into memory and return a list of reformatted alleles."""
#     alleles = []
#     with open(filepath, "r") as f:
#         for line in f:
#             allele = line.strip()
#             if allele:
#                 alleles.append(reformat_allele(allele))
#     return alleles


# def reformat_allele(allele):
#     """
#     Convert allele format from e.g. 'HLA-C07:04' to 'HLA-C*0704'.

#     Handles gene names of varying length (e.g. A, B, C, DRB1, DQB1, etc.)
#     by inserting '*' between the trailing gene letters and the leading
#     digits of the allele designation, and stripping any ':' separators
#     from the digit portion.
#     """
#     match = re.match(r"^(HLA-[A-Za-z]+)(\d.*)$", allele)
#     if not match:
#         raise ValueError(f"Allele '{allele}' does not match expected format 'HLA-<gene><digits>...'")
#     gene_part, number_part = match.groups()
#     number_part = number_part.replace(":", "")
#     return f"{gene_part}*{number_part}"


# def main():
#     args = parse_args()

#     # Read allele file fully into memory first, since it's expected to be small.
#     alleles = read_alleles(args.allele_file)
#     if not alleles:
#         sys.exit(f"Error: no alleles found in {args.allele_file}")

#     row_count = 0

#     # Stream the peptide file line by line and write rows as we go,
#     # so the peptide file is never fully loaded into memory.
#     with open(args.peptide_file, "r") as pep_f, open(args.output, "w", newline="") as out_f:
#         writer = csv.writer(out_f)
#         for line in pep_f:
#             pep = line.strip()
#             if not pep:
#                 continue
#             for allele in alleles:
#                 writer.writerow([pep, allele])
#                 row_count += 1

#     if row_count == 0:
#         sys.exit(f"Error: no peptides found in {args.peptide_file}")

#     print(f"Wrote {row_count} rows to {args.output}")


# if __name__ == "__main__":
#     main()

#!/usr/bin/env python
"""
prep_deepimmuno_input.py

Builds per-length input CSVs for DeepImmuno from a peptide file and an allele file.

Usage:
    python prep_deepimmuno_input.py -p peptide_file -a allele_file \\
        --minlength 9 --maxlength 11 \\
        --output9 deepimmuno_in_9mers.csv \\
        --output10 deepimmuno_in_10mers.csv \\
        --output11 deepimmuno_in_11mers.csv

Input files:
    - peptide_file: one peptide per line (streamed, not loaded fully into memory)
    - allele_file:  one allele per line, e.g. "HLA-C07:04" (read into memory up front)

Output:
    - Up to three CSVs (no header), one row per (peptide, allele) pair, in the form:
        pep,allele
      where allele has been reformatted from "HLA-C07:04" to "HLA-C*07:04".
    - Only peptides whose length is one of 9, 10, or 11 AND falls within
      [minlength, maxlength] are written, and only to the file whose
      corresponding --outputN flag was supplied.
    - A group's output file is created only if at least one peptide belongs
      to that group. If no peptides fall into a group, that group's file is
      never created (left absent, consistent with the Nextflow process's
      `optional: true` output declaration).
    - 11mer peptides are additionally routed through process_11mer_peptide()
      for special handling before being written.
"""

import argparse
import csv
import re
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Prepare per-length DeepImmuno input CSVs from peptide and allele files."
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
        "--minlength", required=True, type=int,
        help="Minimum peptide length to include."
    )
    parser.add_argument(
        "--maxlength", required=True, type=int,
        help="Maximum peptide length to include."
    )
    parser.add_argument(
        "--output9", default=None,
        help="Path to output CSV file for 9mer peptides. If omitted, 9mers are skipped."
    )
    parser.add_argument(
        "--output10", default=None,
        help="Path to output CSV file for 10mer peptides. If omitted, 10mers are skipped."
    )
    parser.add_argument(
        "--output11", default=None,
        help="Path to output CSV file for 11mer peptides. If omitted, 11mers are skipped."
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


def process_11mer_peptide(pep, allele):
    """
    Special handling for 11mer peptides before they're written out.

    TODO: fill in the actual transformation/filtering logic here. This stub
    currently just passes the peptide and allele through unchanged, so
    behavior is identical to the 9mer/10mer path until this is implemented.

    Must return either:
      - a [pep, allele] row (list of two strings) to write a single row, or
      - None to skip writing a row for this (pep, allele) pair.
    """
    return [pep, allele]


class LazyCSVWriter:
    """
    Wraps a csv writer that only opens its underlying file on the first
    call to writerow(). This ensures a group's output file is never
    created on disk unless at least one row belongs to that group.
    """

    def __init__(self, path):
        self.path = path
        self._file = None
        self._writer = None
        self.row_count = 0

    def writerow(self, row):
        if self._writer is None:
            self._file = open(self.path, "w", newline="")
            self._writer = csv.writer(self._file)
        self._writer.writerow(row)
        self.row_count += 1

    def close(self):
        if self._file is not None:
            self._file.close()


def main():
    args = parse_args()

    if args.minlength > args.maxlength:
        sys.exit(f"Error: --minlength ({args.minlength}) is greater than --maxlength ({args.maxlength})")

    # Read allele file fully into memory first, since it's expected to be small.
    alleles = read_alleles(args.allele_file)
    if not alleles:
        sys.exit(f"Error: no alleles found in {args.allele_file}")

    # Build a length -> LazyCSVWriter map, only for lengths whose --outputN
    # flag was actually supplied (and that fall within [minlength, maxlength]
    # is enforced later, per-peptide, so any mismatch just results in zero
    # rows for that group rather than an error).
    writers = {}
    if args.output9 is not None:
        writers[9] = LazyCSVWriter(args.output9)
    if args.output10 is not None:
        writers[10] = LazyCSVWriter(args.output10)
    if args.output11 is not None:
        writers[11] = LazyCSVWriter(args.output11)

    if not writers:
        sys.exit("Error: no output files specified (need at least one of --output9/--output10/--output11)")

    # Stream the peptide file line by line and write rows as we go,
    # so the peptide file is never fully loaded into memory.
    with open(args.peptide_file, "r") as pep_f:
        for line in pep_f:
            pep = line.strip()
            if not pep:
                continue

            length = len(pep)

            # Only 9/10/11mers are handled at all, and only within the
            # requested [minlength, maxlength] range.
            if length < args.minlength or length > args.maxlength:
                continue

            writer = writers.get(length)
            if writer is None:
                # Either not a 9/10/11mer, or that group's output wasn't requested.
                continue

            for allele in alleles:
                if length == 11:
                    #row = process_11mer_peptide(pep, allele)
                    writer.writerow([pep[:10], allele])
                    writer.writerow([pep[1:], allele])

                else:
                    writer.writerow([pep, allele])

    total_rows = 0
    for length, writer in writers.items():
        writer.close()
        total_rows += writer.row_count
        print(f"Wrote {writer.row_count} rows to {writer.path}" if writer.row_count else
              f"No {length}mer peptides found; {writer.path} not created")

    if total_rows == 0:
        sys.exit(f"Error: no peptides found across requested length groups in {args.peptide_file}")


if __name__ == "__main__":
    main()