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


def main():
    args = parse_args()

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
            "best_hla_b",
            "best_hla_c",
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
            best_b = max(gene_scores["B"]) if gene_scores["B"] else ""
            best_c = max(gene_scores["C"]) if gene_scores["C"] else ""

            writer.writerow([pep, best_overall, avg_overall, best_a, best_b, best_c])

    print(f"Wrote summary for {len(peptide_scores)} peptides to {args.output}")


if __name__ == "__main__":
    main()