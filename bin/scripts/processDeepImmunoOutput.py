#!/usr/bin/env python
"""
process_deepimmuno_output.py

Aggregates one or more DeepImmuno output CSV files into a single summary CSV,
one row per peptide.

Usage:
    python process_deepimmuno_output.py -f file1.csv file2.csv ... --output immunogenicity_deepimmuno.csv

Expected input columns (case-insensitive, order-independent):
    - peptide
    - HLA          (e.g. "HLA-A*02:01")
    - immunogenicity   (numeric score)

Output CSV (header included):
    peptide,best_immunogenicity,avg_immunogenicity,best_hla_a,best_hla_b,best_hla_c

Where:
    - best_immunogenicity: max score across all HLA alleles for that peptide
    - avg_immunogenicity:  mean score across all HLA alleles for that peptide
    - best_hla_a/b/c:      max score restricted to alleles of that HLA gene
                            (blank if the peptide has no predictions for that gene)
"""

import argparse
import csv
import re
import sys
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description="Aggregate DeepImmuno output CSV(s) into a per-peptide summary CSV."
    )
    parser.add_argument(
        "-f", "--files", required=True, nargs="+",
        help="One or more DeepImmuno output CSV files."
    )
    parser.add_argument(
        "--lengroup", required=True,
        help="Path to length group of this file"
    )
    parser.add_argument(
        "--output", required=True,
        help="Path to output summary CSV file."
    )
    return parser.parse_args()


def find_column(fieldnames, candidates):
    """Case-insensitive lookup of a column name among a set of acceptable candidates."""
    lookup = {name.strip().lower(): name for name in fieldnames}
    for candidate in candidates:
        if candidate in lookup:
            return lookup[candidate]
    return None


def parse_hla_gene(hla):
    """
    Extract the HLA gene letter(s) from an allele string.
    e.g. 'HLA-A*02:01' -> 'A', 'HLA-C*07:04' -> 'C', 'HLA-DRB1*01:01' -> 'DRB1'
    """
    match = re.match(r"^HLA-?\*?([A-Za-z]+)", hla.strip())
    if not match:
        return None
    return match.group(1).upper()

def get_paired_rows(reader, filepath):
    """
    Yields (row_a, row_b) pairs of consecutive rows from a CSV reader,
    without ever materializing the whole file in memory.
    Warns and drops a trailing unpaired row if the file has an odd row count.
    """
    for row_a in reader:
        try:
            row_b = next(reader)
        except StopIteration:
            print(f"Warning: odd number of rows in {filepath}; dropping last unpaired row",
                  file=sys.stderr)
            return
        yield row_a, row_b

def main():
    args = parse_args()

    if args.lengroup is None:
        sys.exit("Error: No length group provided")
    process11mers = int(args.lengroup) == 11

    # peptide -> list of (hla_gene, score)
    peptide_scores = defaultdict(list)

    for filepath in args.files:
        with open(filepath, "r", newline="") as f:
            reader = csv.DictReader(f, delimiter='\t')
            if reader.fieldnames is None:
                continue
		
            pep_col = "peptide"
            hla_col = "HLA"
            score_col = "immunogenicity"

            if process11mers:
                for row_a, row_b in get_paired_rows(reader, filepath):
                    pep_a, pep_b = row_a[pep_col].strip(), row_b[pep_col].strip()
                    hla_a, hla_b = row_a[hla_col].strip(), row_b[hla_col].strip()
                    score_a, score_b = float(row_a[score_col].strip()), float(row_b[score_col].strip())
                    if hla_a != hla_b:
                        print(f"Warning: mismatched alleles in pair ({hla_a} vs {hla_b}); skipping", file=sys.stderr)
                        continue  # fragments don't belong to the same peptide/allele pair — skip
                    mer11 = pep_a + pep_b[-1]
                    peptide_scores[mer11].append((parse_hla_gene(hla_a), max(score_a, score_b)))
            else:
                for row in reader:
                    pep = row[pep_col].strip()
                    hla = row[hla_col].strip()
                    score_raw = row[score_col].strip()
                    if not pep or not hla or not score_raw:
                        continue
                    try:
                        score = float(score_raw)
                    except ValueError:
                        continue

                    gene = parse_hla_gene(hla)
                    peptide_scores[pep].append((gene, score))

    if not peptide_scores:
        sys.exit("Error: no valid rows found across input files.")

    with open(args.output, "w", newline="") as out_f:
        writer = csv.writer(out_f)
        writer.writerow([
            "peptide",
            "best_immunogenicity",
            "avg_immunogenicity",
            "best_hla_a",
            "avg_hla_a",
            "best_hla_b",
            "avg_hla_b",
            "best_hla_c",
            "avg_hla_c",
        ])

        for pep, entries in peptide_scores.items():
            scores = [score for _, score in entries]
            best_overall = max(scores)
            avg_overall = sum(scores) / len(scores)

            gene_scores = {"A": [], "B": [], "C": []}
            for gene, score in entries:
                if gene in gene_scores:
                    gene_scores[gene].append(score)

            best_a = max(gene_scores["A"]) if gene_scores["A"] else ""
            avg_a = sum(gene_scores["A"]) / len(gene_scores["A"]) if gene_scores["A"] else ""
            best_b = max(gene_scores["B"]) if gene_scores["B"] else ""
            avg_b = sum(gene_scores["B"]) / len(gene_scores["B"]) if gene_scores["B"] else ""
            best_c = max(gene_scores["C"]) if gene_scores["C"] else ""
            avg_c = sum(gene_scores["C"]) / len(gene_scores["C"]) if gene_scores["C"] else ""

            writer.writerow([pep, best_overall, avg_overall, best_a, avg_a, best_b, avg_b, best_c, avg_c])

    print(f"Wrote summary for {len(peptide_scores)} peptides to {args.output}")


if __name__ == "__main__":
    main()