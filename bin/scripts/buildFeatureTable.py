#!/usr/bin/env python3
"""
build_feature_table.py

Builds the combined per-peptide model-input feature table for the Nextflow
BUILD_FEATURE_TABLE process, from a consolidated HLA feature CSV (output of
consolidate_hla_pred_results.py) plus a set of already-processed "default"
feature files (mtec, pathogenicity, deepimmuno, prime), merged on 'peptide'.

Usage:
    python build_feature_table.py \
        --hla hla_prediction_results.csv \
        --deepimmuno deepimmuno_file.csv \
        --prime prime_file.tsv \
        --pathogenicity pathogenicity_file.csv \
        --mtec mtec_file.csv \
        --output combined_featureTable.csv

--hla expects the consolidated output of
consolidate_hla_pred_results.py, which already carries raw per-allele rank
columns plus computed summary columns (best_rank, avg_rank, best_allele,
strong/weak_binders_count, Num_Tools_*, per-tool best_*_r/*_a,
all_strong_binders, all_weak_binders). This script selects only the subset
of those columns relevant to the model-input feature table or exploratory model output files.
"""

import argparse

import pandas as pd

# Known tool best-rank/best-allele output column pairs, as produced by
# consolidate_hla_pred_results.py. Only pairs actually present in the
# --hla file are selected.
TOOL_ALIASES = {
    'NetMHCpan': ('best_netmhc_r', 'netmhc_a'),
    'MixMHCpred': ('best_mixmhc_r', 'mixmhc_a'),
    'MHCFlurry': ('best_mhcflurry_r', 'mhcflurry_a'),
    'MHCnuggets': ('best_mhcnuggets_r', 'mhcnuggets_a'),
}

NUM_TOOLS_COLS = ['Num_Tools_0_05', 'Num_Tools_0_10', 'Num_Tools_0_50']
EXTRA_HLA_FEATS = ['all_strong_binders', 'all_weak_binders']


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _validate_peptide_column(df: pd.DataFrame, filepath: str) -> None:
    """Raise a clear error if the 'peptide' column is missing."""
    if "peptide" not in df.columns:
        raise ValueError(
            f"File '{filepath}' does not contain a 'peptide' column. "
            f"Available columns: {list(df.columns)}"
        )


# ---------------------------------------------------------------------------
# HLA processing
# ---------------------------------------------------------------------------

def process_hla_file(filepath: str) -> pd.DataFrame:
    """
    Read the consolidated HLA feature CSV and select only the columns
    needed for the model-input feature table.
    """
    df = pd.read_csv(filepath)
    _validate_peptide_column(df, filepath)

    model_cols = ['peptide', 'best_rank', 'avg_rank',
                  'strong_binders_count', 'weak_binders_count'] + NUM_TOOLS_COLS + EXTRA_HLA_FEATS
    model_cols = [c for c in model_cols if c in df.columns]

    for rank_col, allele_col in TOOL_ALIASES.values():
        if rank_col in df.columns:
            model_cols.append(rank_col)
        if allele_col in df.columns:
            model_cols.append(allele_col)

    feature_df = df[model_cols].copy()
    feature_df = feature_df.drop_duplicates(subset=['peptide'])
    return feature_df


# ---------------------------------------------------------------------------
# Default / other feature-file processing
# ---------------------------------------------------------------------------

def process_default_file(filepath: str) -> pd.DataFrame:
    """Load a file that needs no column renaming (e.g. mtec)."""
    df = pd.read_csv(filepath)
    _validate_peptide_column(df, filepath)
    df = df.drop_duplicates(subset=['peptide'])
    return df

def process_mtec_counts_file(filepath: str) -> pd.DataFrame:
    df = pd.read_csv(filepath, sep='\t', comment='#')
    df = df[['peptide', 'total']]
    df.rename(columns={'total': 'mtec_expression_count'}, inplace=True)
    _validate_peptide_column(df, filepath)
    df = df.drop_duplicates(subset=['peptide'])
    return df

def process_pathogenicity_file(filepath: str) -> pd.DataFrame:
    df = pd.read_csv(filepath)
    df.rename(columns={'similarity_score': 'pathogenicity'}, inplace=True)
    _validate_peptide_column(df, filepath)
    df = df.drop_duplicates(subset=['peptide'])
    return df


def process_deepimmuno_file(filepath: str) -> pd.DataFrame:
    df = pd.read_csv(filepath)
    df.rename(columns={
        'best_immunogenicity': 'di_best_score',
        'avg_immunogenicity': 'di_avg_score',
        'best_hla_a': 'di_best_a_score',
        'avg_hla_a': 'di_avg_a_score',
        'best_hla_b': 'di_best_b_score',
        'avg_hla_b': 'di_avg_b_score',
        'best_hla_c': 'di_best_c_score',
        'avg_hla_c': 'di_avg_c_score',
    }, inplace=True)
    _validate_peptide_column(df, filepath)
    df = df.drop_duplicates(subset=['peptide'])
    return df


def process_prime_file(filepath: str) -> pd.DataFrame:
    df = pd.read_csv(filepath, sep='\t', comment='#')
    df = df[['Peptide', '%Rank_bestAllele']]
    df.rename(columns={'Peptide': 'peptide', '%Rank_bestAllele': 'PRIME_%Rank_bestAllele'}, inplace=True)
    _validate_peptide_column(df, filepath)
    df = df.drop_duplicates(subset=['peptide'])
    return df


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build the combined per-peptide model-input feature table."
    )
    parser.add_argument('--hla', required=True, help="Consolidated HLA feature CSV (output of consolidate_hla_pred_results.py).")
    parser.add_argument('--mtec', required=True, help="MTEC feature file (no renaming needed).")
    parser.add_argument('--pathogenicity', required=True, help="Pathogenicity similarity-score file.")
    parser.add_argument('--deepimmuno', required=True, help="DeepImmuno feature file.")
    parser.add_argument('--prime', required=True, help="PRIME predictor output file.")
    parser.add_argument('--output', required=True, help="Path to write the combined feature CSV to.")

    parser.add_argument('--selfsimilarity', required=False, help="selfsimilarity output file.")
    parser.add_argument('--tap', required=False, help="tap output file.")

    return parser.parse_args(argv)


def main(argv=None) -> None:
    args = parse_args(argv)

    print(f"Processing HLA file: {args.hla}")
    hla_df = process_hla_file(args.hla)

    print(f"Processing MTEC file: {args.mtec}")
    # mtec_df = process_default_file(args.mtec)
    mtec_df = process_mtec_counts_file(args.mtec)

    print(f"Processing pathogenicity file: {args.pathogenicity}")
    pathogenicity_df = process_pathogenicity_file(args.pathogenicity)

    print(f"Processing DeepImmuno file: {args.deepimmuno}")
    deepimmuno_df = process_deepimmuno_file(args.deepimmuno)

    print(f"Processing PRIME file: {args.prime}")
    prime_df = process_prime_file(args.prime)

    selfsim_df, tap_df = None, None
    if args.selfsimilarity is not None:
        print(f"Processing selfsimilarity file: {args.selfsimilarity}")
        selfsim_df = process_default_file(args.selfsimilarity)
    
    if args.tap is not None:
        print(f"Processing TAP file: {args.tap}")
        tap_df = process_default_file(args.tap)
        tap_df.rename(columns={'pred_score' : 'tap_score'}, inplace=True)

    print("Merging feature tables on 'peptide'...")
    # 'outer' so that DeepImmuno (which may be missing peptides) contributes
    # NaN rather than dropping peptides present in the other tables.
    combined = hla_df
    for df in (mtec_df, pathogenicity_df, deepimmuno_df, prime_df, selfsim_df, tap_df):
        if df is not None:
            combined = combined.merge(df, on='peptide', how='outer')

    combined.to_csv(args.output, index=False)
    print(f"Saved {combined.shape[0]} peptides, {combined.shape[1]} columns to {args.output}")


if __name__ == '__main__':
    main()