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

--hla is expected to be the consolidated output of
consolidate_hla_pred_results.py, which already carries raw per-allele rank
columns plus computed summary columns (best_rank, avg_rank, best_allele,
strong/weak_binders_count, Num_Tools_*, per-tool best_*_r/*_a,
all_strong_binders, all_weak_binders). This script selects only the subset
of those columns relevant to the model-input feature table (dropping the
raw per-allele columns and the all_strong_binders/all_weak_binders string
columns).
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
                  'strong_binders_count', 'weak_binders_count'] + NUM_TOOLS_COLS
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

    # TODO REMOVE
    parser.add_argument('--selfsimilarity', required=True, help="selfsimilarity output file.")
    parser.add_argument('--tap', required=True, help="tap output file.")

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

    # TODO REMOVE
    print(f"Processing selfsimilarity file: {args.prime}")
    selfsim_df = process_default_file(args.selfsimilarity)

    print(f"Processing TAP file: {args.prime}")
    tap_df = process_default_file(args.tap)
    tap_df.rename(columns={'pred_score' : 'tap_score'}, inplace=True)


    # TODO DONE

    print("Merging feature tables on 'peptide'...")
    # 'outer' so that DeepImmuno (which may be missing peptides) contributes
    # NaN rather than dropping peptides present in the other tables.
    # TODO REMOVE SELFSIM AND TAP
    combined = hla_df
    for df in (mtec_df, pathogenicity_df, deepimmuno_df, prime_df, selfsim_df, tap_df):
        combined = combined.merge(df, on='peptide', how='outer')

    combined.to_csv(args.output, index=False)
    print(f"Saved {combined.shape[0]} peptides, {combined.shape[1]} columns to {args.output}")


if __name__ == '__main__':
    main()


# # #!/usr/bin/env python3
# # """
# # build_feature_table.py

# # Builds the combined per-peptide model-input feature table for the Nextflow
# # BUILD_FEATURE_TABLE process, from a consolidated HLA feature CSV (output of
# # consolidate_hla_pred_results.py) plus a set of already-processed "default"
# # feature files (mtec, pathogenicity, deepimmuno, prime), merged on 'peptide'.

# # Usage:
# #     python build_feature_table.py \
# #         --hla hla_prediction_results.csv \
# #         --deepimmuno deepimmuno_file.csv \
# #         --prime prime_file.tsv \
# #         --pathogenicity pathogenicity_file.csv \
# #         --mtec mtec_file.csv \
# #         --output combined_featureTable.csv

# # --hla is expected to be the consolidated output of
# # consolidate_hla_pred_results.py, which already carries raw per-allele rank
# # columns plus computed summary columns (best_rank, avg_rank, best_allele,
# # strong/weak_binders_count, Num_Tools_*, per-tool best_*_r/*_a,
# # all_strong_binders, all_weak_binders). This script selects only the subset
# # of those columns relevant to the model-input feature table (dropping the
# # raw per-allele columns and the all_strong_binders/all_weak_binders string
# # columns).
# # """

# # import argparse

# # import pandas as pd

# # # Known tool best-rank/best-allele output column pairs, as produced by
# # # consolidate_hla_pred_results.py. Only pairs actually present in the
# # # --hla file are selected.
# # TOOL_ALIASES = {
# #     'NetMHCpan': ('best_netmhc_r', 'netmhc_a'),
# #     'MixMHCpred': ('best_mixmhc_r', 'mixmhc_a'),
# #     'MHCFlurry': ('best_mhcflurry_r', 'mhcflurry_a'),
# #     'MHCnuggets': ('best_mhcnuggets_r', 'mhcnuggets_a'),
# # }

# # NUM_TOOLS_COLS = [
# #     'Num_Tools_Super_Strong_Binder', 'Num_Tools_Strong_Binder',
# #     'Num_Tools_Weak_Binder', 'Num_Tools_Super_Weak_Binder',
# #     'Num_Tools_0_01', 'Num_Tools_0_02', 'Num_Tools_0_03', 'Num_Tools_0_04',
# #     'Num_Tools_0_05', 'Num_Tools_0_10', 'Num_Tools_0_50',
# #     'Num_Tools_1', 'Num_Tools_2',
# # ]


# # # ---------------------------------------------------------------------------
# # # Helpers
# # # ---------------------------------------------------------------------------

# # def _validate_peptide_column(df: pd.DataFrame, filepath: str) -> None:
# #     """Raise a clear error if the 'peptide' column is missing."""
# #     if "peptide" not in df.columns:
# #         raise ValueError(
# #             f"File '{filepath}' does not contain a 'peptide' column. "
# #             f"Available columns: {list(df.columns)}"
# #         )


# # # ---------------------------------------------------------------------------
# # # HLA processing
# # # ---------------------------------------------------------------------------

# # def process_hla_file(filepath: str) -> pd.DataFrame:
# #     """
# #     Read the consolidated HLA feature CSV and select only the columns
# #     needed for the model-input feature table.
# #     """
# #     df = pd.read_csv(filepath)
# #     _validate_peptide_column(df, filepath)

# #     model_cols = ['peptide', 'best_rank', 'avg_rank',
# #                   'strong_binders_count', 'weak_binders_count'] #+ NUM_TOOLS_COLS
# #     model_cols = [c for c in model_cols if c in df.columns]

# #     # for rank_col, allele_col in TOOL_ALIASES.values():
# #     #     if rank_col in df.columns:
# #     #         model_cols.append(rank_col)
# #     #     if allele_col in df.columns:
# #     #         model_cols.append(allele_col)

# #     feature_df = df[model_cols].copy()
# #     feature_df = feature_df.drop_duplicates(subset=['peptide'])
# #     return feature_df


# # # ---------------------------------------------------------------------------
# # # Default / other feature-file processing
# # # ---------------------------------------------------------------------------

# # def process_default_file(filepath: str) -> pd.DataFrame:
# #     """Load a file that needs no column renaming (e.g. mtec)."""
# #     df = pd.read_csv(filepath)
# #     _validate_peptide_column(df, filepath)
# #     df = df.drop_duplicates(subset=['peptide'])
# #     return df


# # def process_pathogenicity_file(filepath: str) -> pd.DataFrame:
# #     df = pd.read_csv(filepath)
# #     df.rename(columns={'similarity_score': 'pathogenicity'}, inplace=True)
# #     _validate_peptide_column(df, filepath)
# #     df = df.drop_duplicates(subset=['peptide'])
# #     return df


# # def process_deepimmuno_file(filepath: str) -> pd.DataFrame:
# #     df = pd.read_csv(filepath)
# #     df.rename(columns={
# #         'best_immunogenicity': 'di_best_score',
# #         'avg_immunogenicity': 'di_avg_score',
# #         'best_hla_a': 'di_best_a_score',
# #         'avg_hla_a': 'di_avg_a_score',
# #         'best_hla_b': 'di_best_b_score',
# #         'avg_hla_b': 'di_avg_b_score',
# #         'best_hla_c': 'di_best_c_score',
# #         'avg_hla_c': 'di_avg_c_score',
# #     }, inplace=True)
# #     _validate_peptide_column(df, filepath)
# #     df = df.drop_duplicates(subset=['peptide'])
# #     return df


# # def process_prime_file(filepath: str) -> pd.DataFrame:
# #     df = pd.read_csv(filepath, sep='\t', comment='#')
# #     df = df[['Peptide', '%Rank_bestAllele']]
# #     df.rename(columns={'Peptide': 'peptide', '%Rank_bestAllele': 'PRIME_%Rank_bestAllele'}, inplace=True)
# #     _validate_peptide_column(df, filepath)
# #     df = df.drop_duplicates(subset=['peptide'])
# #     return df


# # # ---------------------------------------------------------------------------
# # # CLI entry point
# # # ---------------------------------------------------------------------------

# # def parse_args(argv=None) -> argparse.Namespace:
# #     parser = argparse.ArgumentParser(
# #         description="Build the combined per-peptide model-input feature table."
# #     )
# #     parser.add_argument('--hla', required=True, help="Consolidated HLA feature CSV (output of consolidate_hla_pred_results.py).")
# #     parser.add_argument('--mtec', required=True, help="MTEC feature file (no renaming needed).")
# #     parser.add_argument('--pathogenicity', required=True, help="Pathogenicity similarity-score file.")
# #     parser.add_argument('--deepimmuno', required=True, help="DeepImmuno feature file.")
# #     parser.add_argument('--prime', required=True, help="PRIME predictor output file.")
# #     parser.add_argument('--output', required=True, help="Path to write the combined feature CSV to.")
# #     return parser.parse_args(argv)


# # def main(argv=None) -> None:
# #     args = parse_args(argv)

# #     print(f"Processing HLA file: {args.hla}")
# #     hla_df = process_hla_file(args.hla)

# #     print(f"Processing MTEC file: {args.mtec}")
# #     mtec_df = process_default_file(args.mtec)

# #     print(f"Processing pathogenicity file: {args.pathogenicity}")
# #     pathogenicity_df = process_pathogenicity_file(args.pathogenicity)

# #     print(f"Processing DeepImmuno file: {args.deepimmuno}")
# #     deepimmuno_df = process_deepimmuno_file(args.deepimmuno)

# #     print(f"Processing PRIME file: {args.prime}")
# #     prime_df = process_prime_file(args.prime)

# #     print("Merging feature tables on 'peptide'...")
# #     # 'outer' so that DeepImmuno (which may be missing peptides) contributes
# #     # NaN rather than dropping peptides present in the other tables.
# #     combined = hla_df
# #     for df in (mtec_df, pathogenicity_df, deepimmuno_df, prime_df):
# #         combined = combined.merge(df, on='peptide', how='outer')

# #     combined.to_csv(args.output, index=False)
# #     print(f"Saved {combined.shape[0]} peptides, {combined.shape[1]} columns to {args.output}")


# # if __name__ == '__main__':
# #     main()

# #!/usr/bin/env python3
# """
# build_feature_table.py

# Builds the combined per-peptide model-input feature table for the Nextflow
# BUILD_FEATURE_TABLE process, from a single raw HLA-predictor TSV plus a set
# of already-processed "default" feature files (mtec, pathogenicity,
# deepimmuno, prime), merged on 'peptide'.

#     python build_feature_table.py \
#         --hla hla_file.tsv \
#         --deepimmuno deepimmuno_file.csv \
#         --prime prime_file.tsv \
#         --pathogenicity pathogenicity_file.csv \
#         --mtec mtec_file.csv \
#         --common_alleles common_alleles.txt \
#         --all_alleles all_alleles.txt \
#         --output combined_featureTable.csv

# HLA feature calculation
# ------------------------
# The raw --hla file is expected to contain one row per peptide with
# per-(allele, tool) '%Rank' columns, e.g. 'HLA-A01:09_MHCFlurry_%Rank',
# covering ALL alleles that were run through the predictor (single file,
# no batching).

# Two different allele subsets drive different parts of the output, per
# allele lists supplied as one-allele-per-line text files:

#   --all_alleles     Used for: best_rank, avg_rank, best_allele, and the
#                      four per-tool best-rank/best-allele columns
#                      (best_netmhc_r/netmhc_a, best_mixmhc_r/mixmhc_a,
#                      best_mhcflurry_r/mhcflurry_a, best_mhcnuggets_r/
#                      mhcnuggets_a).

#   --common_alleles  Used for: strong_binders_count, weak_binders_count,
#                      and all Num_Tools_* threshold-agreement columns.
# """

# import argparse
# import re

# import pandas as pd

# DEFAULT_TOOLS = ['NetMHCpan', 'MixMHCpred', 'MHCFlurry', 'MHCnuggets']

# # Friendly short column-name aliases for known tools' best-rank/allele
# # output columns.
# TOOL_ALIASES = {
#     'NetMHCpan': ('best_netmhc_r', 'netmhc_a'),
#     'MixMHCpred': ('best_mixmhc_r', 'mixmhc_a'),
#     'MHCFlurry': ('best_mhcflurry_r', 'mhcflurry_a'),
#     'MHCnuggets': ('best_mhcnuggets_r', 'mhcnuggets_a'),
# }


# def _tool_column_names(tool: str) -> tuple[str, str]:
#     """Return (rank_col, allele_col) output names for a tool, using the
#     friendly alias if known, otherwise deriving one from the tool name."""
#     if tool in TOOL_ALIASES:
#         return TOOL_ALIASES[tool]
#     slug = tool.lower()
#     return f'best_{slug}_r', f'{slug}_a'


# # ---------------------------------------------------------------------------
# # Helpers
# # ---------------------------------------------------------------------------

# def _validate_peptide_column(df: pd.DataFrame, filepath: str) -> None:
#     """Raise a clear error if the 'peptide' column is missing."""
#     if "peptide" not in df.columns:
#         raise ValueError(
#             f"File '{filepath}' does not contain a 'peptide' column. "
#             f"Available columns: {list(df.columns)}"
#         )


# def read_allele_list(filepath: str) -> set[str]:
#     """Read a one-allele-per-line text file into a set of allele names."""
#     with open(filepath) as f:
#         alleles = {line.strip() for line in f if line.strip()}
#     return alleles


# def _rank_cols_for_alleles(df: pd.DataFrame, tools: list[str], alleles: set[str]) -> list[str]:
#     """
#     Return the '%Rank' columns in df whose allele prefix (e.g. 'HLA-A01:09'
#     extracted from 'HLA-A01:09_MHCFlurry_%Rank') is in `alleles`, and whose
#     tool matches one of `tools`.
#     """
#     cols = []
#     for col in df.columns:
#         m = re.match(r'^(HLA-[^_]+)_', col)
#         if not m:
#             continue
#         if m.group(1) not in alleles:
#             continue
#         if any(f'{tool}_%Rank' in col for tool in tools):
#             cols.append(col)
#     return cols


# def _per_tool_best_rank(df: pd.DataFrame, tools: list[str], rank_cols: list[str]) -> pd.DataFrame:
#     """
#     Given a restricted set of rank_cols (already filtered to a particular
#     allele subset), compute each tool's best (minimum) rank per peptide.

#     Returns a DataFrame indexed like df, with one '{tool}_bestRank' column
#     per tool that has at least one matching column.
#     """
#     out = pd.DataFrame(index=df.index)
#     for tool in tools:
#         tool_cols = [c for c in rank_cols if f'{tool}_%Rank' in c]
#         if tool_cols:
#             out[f'{tool}_bestRank'] = df[tool_cols].min(axis=1)
#     return out


# # ---------------------------------------------------------------------------
# # HLA processing
# # ---------------------------------------------------------------------------

# def process_hla_file(
#     filepath: str,
#     tools: list[str],
#     common_alleles: set[str],
#     all_alleles: set[str],
# ) -> pd.DataFrame:
#     """
#     Build the HLA feature block from a single raw hla-predictor TSV,
#     covering ALL peptides and ALL alleles that were run.

#     best_rank / avg_rank / best_allele / per-tool best_*_r / *_a columns
#     are computed from the `all_alleles` subset of columns.

#     strong_binders_count / weak_binders_count / Num_Tools_* columns are
#     computed from the `common_alleles` subset of columns.
#     """
#     df = pd.read_csv(filepath, sep='\t')
#     _validate_peptide_column(df, filepath)

#     # ---- all_alleles-based: best_rank, avg_rank, best_allele, per-tool cols
#     all_rank_cols = _rank_cols_for_alleles(df, tools, all_alleles)
#     if not all_rank_cols:
#         raise ValueError(
#             f"No rank columns in '{filepath}' matched any allele in --all_alleles."
#         )

#     best_binding_allele = df[all_rank_cols].idxmin(axis=1).str.extract(r'^(HLA-[^_]+)')[0]

#     per_tool_best_all = _per_tool_best_rank(df, tools, all_rank_cols)
#     tool_best_rank_cols = list(per_tool_best_all.columns)

#     best_rank = per_tool_best_all[tool_best_rank_cols].min(axis=1)
#     avg_rank = per_tool_best_all[tool_best_rank_cols].mean(axis=1)

#     # Per-tool best-rank/best-allele output columns (all_alleles-based)
#     tool_output_cols = {}
#     for tool in tools:
#         tool_cols = [c for c in all_rank_cols if f'{tool}_%Rank' in c]
#         if not tool_cols:
#             continue
#         rank_col_name, allele_col_name = _tool_column_names(tool)
#         tool_output_cols[rank_col_name] = df[tool_cols].min(axis=1)
#         tool_output_cols[allele_col_name] = df[tool_cols].idxmin(axis=1).str.extract(r'^(HLA-[^_]+)')[0]

#     # ---- common_alleles-based: binder counts, Num_Tools_* columns
#     common_rank_cols = _rank_cols_for_alleles(df, tools, common_alleles)
#     if not common_rank_cols:
#         raise ValueError(
#             f"No rank columns in '{filepath}' matched any allele in --common_alleles."
#         )

#     melted = df[['peptide'] + common_rank_cols].melt(
#         id_vars='peptide', var_name='col', value_name='rank'
#     )
#     melted['rank'] = pd.to_numeric(melted['rank'], errors='coerce')

#     strong_counts = (
#         melted[melted['rank'] < 0.5]
#         .groupby('peptide')['col'].count()
#         .rename('strong_binders_count')
#     )
#     weak_counts = (
#         melted[melted['rank'] < 2]
#         .groupby('peptide')['col'].count()
#         .rename('weak_binders_count')
#     )

#     per_tool_best_common = _per_tool_best_rank(df, tools, common_rank_cols)
#     common_tool_best_cols = list(per_tool_best_common.columns)
#     common_best_ranks_numeric = per_tool_best_common[common_tool_best_cols].apply(
#         pd.to_numeric, errors='coerce'
#     )

#     num_tools = pd.DataFrame(index=df.index)
#     num_tools['Num_Tools_0_05'] = (common_best_ranks_numeric < 0.05).sum(axis=1)
#     num_tools['Num_Tools_0_10'] = (common_best_ranks_numeric < 0.1).sum(axis=1)
#     num_tools['Num_Tools_0_50'] = (common_best_ranks_numeric < 0.5).sum(axis=1)

#     # ---- assemble output ----
#     feature_df = pd.DataFrame()
#     feature_df['peptide'] = df['peptide']
#     feature_df['best_rank'] = best_rank
#     feature_df['avg_rank'] = avg_rank
#     feature_df['best_allele'] = best_binding_allele

#     feature_df = feature_df.merge(strong_counts, on='peptide', how='left')
#     feature_df = feature_df.merge(weak_counts, on='peptide', how='left')
#     feature_df['strong_binders_count'] = feature_df['strong_binders_count'].fillna(0).astype(int)
#     feature_df['weak_binders_count'] = feature_df['weak_binders_count'].fillna(0).astype(int)

#     for col in num_tools.columns:
#         feature_df[col] = num_tools[col].values

#     for col_name, values in tool_output_cols.items():
#         feature_df[col_name] = values.values

#     feature_df = feature_df.drop_duplicates(subset=['peptide'])
#     return feature_df


# # ---------------------------------------------------------------------------
# # Default / other feature-file processing
# # ---------------------------------------------------------------------------

# def process_default_file(filepath: str) -> pd.DataFrame:
#     """Load a file that needs no column renaming (e.g. mtec)."""
#     df = pd.read_csv(filepath)
#     _validate_peptide_column(df, filepath)
#     df = df.drop_duplicates(subset=['peptide'])
#     return df


# def process_pathogenicity_file(filepath: str) -> pd.DataFrame:
#     df = pd.read_csv(filepath)
#     df.rename(columns={'similarity_score': 'pathogenicity'}, inplace=True)
#     _validate_peptide_column(df, filepath)
#     df = df.drop_duplicates(subset=['peptide'])
#     return df


# def process_deepimmuno_file(filepath: str) -> pd.DataFrame:
#     df = pd.read_csv(filepath)
#     df.rename(columns={
#         'best_immunogenicity': 'di_best_score',
#         'avg_immunogenicity': 'di_avg_score',
#         'best_hla_a': 'di_best_a_score',
#         'avg_hla_a': 'di_avg_a_score',
#         'best_hla_b': 'di_best_b_score',
#         'avg_hla_b': 'di_avg_b_score',
#         'best_hla_c': 'di_best_c_score',
#         'avg_hla_c': 'di_avg_c_score',
#     }, inplace=True)
#     _validate_peptide_column(df, filepath)
#     df = df.drop_duplicates(subset=['peptide'])
#     return df


# def process_prime_file(filepath: str) -> pd.DataFrame:
#     df = pd.read_csv(filepath, sep='\t', comment='#')
#     df = df[['Peptide', '%Rank_bestAllele']]
#     df.rename(columns={'Peptide': 'peptide', '%Rank_bestAllele': 'PRIME_%Rank_bestAllele'}, inplace=True)
#     _validate_peptide_column(df, filepath)
#     df = df.drop_duplicates(subset=['peptide'])
#     return df


# # ---------------------------------------------------------------------------
# # CLI entry point
# # ---------------------------------------------------------------------------

# def parse_args(argv=None) -> argparse.Namespace:
#     parser = argparse.ArgumentParser(
#         description="Build the combined per-peptide model-input feature table."
#     )
#     parser.add_argument('--hla', required=True, help="Raw hla-predictor TSV (all peptides, all alleles).")
#     parser.add_argument('--mtec', required=True, help="MTEC feature file (no renaming needed).")
#     parser.add_argument('--pathogenicity', required=True, help="Pathogenicity similarity-score file.")
#     parser.add_argument('--deepimmuno', required=True, help="DeepImmuno feature file.")
#     parser.add_argument('--prime', required=True, help="PRIME predictor output file.")
#     parser.add_argument('--common_alleles', required=True, help="Text file, one allele per line: binder-count / Num_Tools_* allele subset.")
#     parser.add_argument('--all_alleles', required=True, help="Text file, one allele per line: best_rank/avg_rank/best_allele allele subset.")
#     parser.add_argument('--output', required=True, help="Path to write the combined feature CSV to.")
#     parser.add_argument('--tools', nargs='+', default=DEFAULT_TOOLS, help=f"Prediction tools to use (default: {DEFAULT_TOOLS}).")
#     return parser.parse_args(argv)


# def main(argv=None) -> None:
#     args = parse_args(argv)

#     common_alleles = read_allele_list(args.common_alleles)
#     all_alleles = read_allele_list(args.all_alleles)

#     print(f"Loaded {len(common_alleles)} common alleles, {len(all_alleles)} total alleles.")

#     print(f"Processing HLA file: {args.hla}")
#     hla_df = process_hla_file(args.hla, args.tools, common_alleles, all_alleles)

#     print(f"Processing DeepImmuno file: {args.deepimmuno}")
#     deepimmuno_df = process_deepimmuno_file(args.deepimmuno)
    
#     print(f"Processing PRIME file: {args.prime}")
#     prime_df = process_prime_file(args.prime)

#     print(f"Processing pathogenicity file: {args.pathogenicity}")
#     pathogenicity_df = process_pathogenicity_file(args.pathogenicity)

#     print(f"Processing MTEC file: {args.mtec}")
#     mtec_df = process_default_file(args.mtec)

#     print("Merging feature tables on 'peptide'...")
#     # 'outer' so that DeepImmuno (which may be missing peptides) contributes
#     # NaN rather than dropping peptides present in the other tables.
#     combined = hla_df
#     for df in (deepimmuno_df, prime_df, pathogenicity_df ,mtec_df):
#         combined = combined.merge(df, on='peptide', how='outer')

#     combined.to_csv(args.output, index=False)
#     print(f"Saved {combined.shape[0]} peptides, {combined.shape[1]} columns to {args.output}")


# if __name__ == '__main__':
#     main()