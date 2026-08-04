#!/usr/bin/env python3
"""
consolidate_hla_pred_results.py

Consolidates one or more raw hla-predictor TSV output files into a single
processed HLA feature CSV, for the Nextflow PROCESS_HLA_WORKFLOW_OUTPUT
process.

Input files are expected to all describe the SAME set of peptides, split
across different allele batches (i.e. the standard HLA_WORKFLOW output
layout: one file per allele batch). They are merged column-wise on
'peptide'.

Two allele subsets, supplied as one-allele-per-line text files, drive
different parts of the output:

  --all-alleles     Used for: best_rank, avg_rank, best_allele, the
                     per-tool best-rank/best-allele columns, and
                     all_strong_binders / all_weak_binders.

  --common-alleles  Used for: strong_binders_count, weak_binders_count,
                     and all Num_Tools_* threshold-agreement columns.

The output CSV carries all raw per-(allele, tool) rank columns from the
merged input, unfiltered, in addition to the computed summary columns
above. Downstream, build_feature_table.py reads this file and selects
only the columns it needs for the model-input feature table.

Usage:
    python consolidate_hla_pred_results.py file1.tsv file2.tsv ... \
        --all-alleles all_alleles.txt \
        --common-alleles common_alleles.txt \
        --output hla_prediction_results.csv
"""

import argparse
import re

import pandas as pd

DEFAULT_TOOLS = ['NetMHCpan', 'MixMHCpred', 'MHCFlurry']

# Friendly short column-name aliases for known tools' best-rank/allele
# output columns.
TOOL_ALIASES = {
    'NetMHCpan': ('best_netmhc_r', 'netmhc_a'),
    'MixMHCpred': ('best_mixmhc_r', 'mixmhc_a'),
    'MHCFlurry': ('best_mhcflurry_r', 'mhcflurry_a'),
    'MHCnuggets': ('best_mhcnuggets_r', 'mhcnuggets_a'),
}

# Matches raw per-(allele, tool) rank columns, e.g. 'HLA-A01:09_MHCFlurry_%Rank'.
RAW_RANK_COL_PATTERN = re.compile(r'^(HLA-[^_]+)_.+_%Rank$')


def _tool_column_names(tool: str) -> tuple[str, str]:
    """Return (rank_col, allele_col) output names for a tool, using the
    friendly alias if known, otherwise deriving one from the tool name."""
    if tool in TOOL_ALIASES:
        return TOOL_ALIASES[tool]
    tool = tool.lower()
    return f'best_{tool}_r', f'{tool}_a'


def read_allele_list(filepath: str) -> set[str]:
    """Read a one-allele-per-line text file into a set of allele names."""
    with open(filepath) as f:
        return {line.strip() for line in f if line.strip()}


def _rank_cols_for_alleles(df: pd.DataFrame, tools: list[str], alleles: set[str]) -> list[str]:
    """
    Return the raw '%Rank' columns in df whose allele prefix is in
    `alleles` and whose tool matches one of `tools`.
    """
    cols = []
    for col in df.columns:
        m = RAW_RANK_COL_PATTERN.match(col)
        if not m or m.group(1) not in alleles:
            continue
        if any(f'{tool}_%Rank' in col for tool in tools):
            cols.append(col)
    return cols


def _per_tool_best_rank(df: pd.DataFrame, tools: list[str], rank_cols: list[str]) -> pd.DataFrame:
    """Compute each tool's best (minimum) rank per peptide, restricted to rank_cols."""
    out = pd.DataFrame(index=df.index)
    for tool in tools:
        tool_cols = [c for c in rank_cols if f'{tool}_%Rank' in c]
        if tool_cols:
            out[f'{tool}_bestRank'] = df[tool_cols].min(axis=1)
    return out


def get_binding_alleles(df: pd.DataFrame, rank_cols: list[str]) -> tuple[pd.Series, pd.Series]:
    """
    For each peptide, find all unique alleles (restricted to `rank_cols`)
    that are strong (rank < 0.5) or weak-or-strong (rank < 2) binders
    across all tools.

    Returns two Series (indexed by peptide) of comma-separated, sorted
    allele name strings: (strong_alleles, weak_alleles).
    """
    melted = df[['peptide'] + rank_cols].melt(id_vars='peptide', var_name='col', value_name='rank')
    melted['allele'] = melted['col'].str.extract(r'^(HLA-[^_]+)')
    melted['rank'] = pd.to_numeric(melted['rank'], errors='coerce')

    strong = melted[melted['rank'] < 0.5]
    weak = melted[melted['rank'] < 2]

    strong_alleles = strong.groupby('peptide')['allele'].apply(lambda x: ','.join(sorted(x.unique())))
    weak_alleles = weak.groupby('peptide')['allele'].apply(lambda x: ','.join(sorted(x.unique())))

    return strong_alleles, weak_alleles


def merge_allele_batches(filepaths: list[str]) -> pd.DataFrame:
    """
    Read one or more allele-batched HLA prediction TSVs, each covering the
    same peptide set but a different allele batch, keep only 'peptide' plus
    each file's raw rank columns (per-file summary columns like
    Best_Binding_Allele or {tool}_bestRank are dropped, since all summary
    statistics are recomputed from the merged raw data), and merge them
    column-wise on 'peptide'.
    """
    merged = None
    seen_rank_cols: set[str] = set()

    for filepath in filepaths:
        df = pd.read_csv(filepath, sep='\t')
        if 'peptide' not in df.columns:
            raise ValueError(f"File '{filepath}' does not contain a 'peptide' column.")

        rank_cols = [c for c in df.columns if RAW_RANK_COL_PATTERN.match(c)]

        duplicate_cols = seen_rank_cols & set(rank_cols)
        if duplicate_cols:
            raise ValueError(
                f"File '{filepath}' contains rank column(s) already seen in a previous "
                f"file: {sorted(duplicate_cols)}. Input files are expected to cover "
                f"distinct allele batches for the same peptide set."
            )
        seen_rank_cols.update(rank_cols)

        df = df[['peptide'] + rank_cols]
        merged = df if merged is None else merged.merge(df, on='peptide', how='outer')

    return merged


def build_consolidated_hla_table(filepaths: list[str], tools: list[str], all_alleles: set[str], common_alleles: set[str]) -> pd.DataFrame:
    """Merge allele-batch files and compute all summary feature columns."""
    df = merge_allele_batches(filepaths)
    raw_rank_cols = [c for c in df.columns if RAW_RANK_COL_PATTERN.match(c)]

    # ---- all_alleles-based: best_rank, avg_rank, best_allele, per-tool cols, binder-allele strings
    all_rank_cols = _rank_cols_for_alleles(df, tools, all_alleles)
    if not all_rank_cols:
        raise ValueError("No rank columns matched any allele in --all-alleles.")

    best_binding_allele = df[all_rank_cols].idxmin(axis=1).str.extract(r'^(HLA-[^_]+)')[0]

    per_tool_best_all = _per_tool_best_rank(df, tools, all_rank_cols)
    tool_best_rank_cols = list(per_tool_best_all.columns)
    best_rank = per_tool_best_all[tool_best_rank_cols].min(axis=1)
    avg_rank = per_tool_best_all[tool_best_rank_cols].mean(axis=1)

    tool_output_cols = {}
    for tool in tools:
        tool_cols = [c for c in all_rank_cols if f'{tool}_%Rank' in c]
        if not tool_cols:
            continue
        rank_col_name, allele_col_name = _tool_column_names(tool)
        tool_output_cols[rank_col_name] = df[tool_cols].min(axis=1)
        tool_output_cols[allele_col_name] = df[tool_cols].idxmin(axis=1).str.extract(r'^(HLA-[^_]+)')[0]

    all_strong_binders, all_weak_binders = get_binding_alleles(df, all_rank_cols)

    # ---- common_alleles-based: binder counts, Num_Tools_* columns
    common_rank_cols = _rank_cols_for_alleles(df, tools, common_alleles)
    if not common_rank_cols:
        raise ValueError("No rank columns matched any allele in --common-alleles.")

    melted = df[['peptide'] + common_rank_cols].melt(id_vars='peptide', var_name='col', value_name='rank')
    melted['rank'] = pd.to_numeric(melted['rank'], errors='coerce')

    strong_counts = melted[melted['rank'] < 0.5].groupby('peptide')['col'].count().rename('strong_binders_count')
    weak_counts = melted[melted['rank'] < 2].groupby('peptide')['col'].count().rename('weak_binders_count')

    per_tool_best_common = _per_tool_best_rank(df, tools, common_rank_cols)
    common_tool_best_cols = list(per_tool_best_common.columns)
    common_best_ranks_numeric = per_tool_best_common[common_tool_best_cols].apply(pd.to_numeric, errors='coerce')

    num_tools = pd.DataFrame(index=df.index)
    num_tools['Num_Tools_0_01'] = (common_best_ranks_numeric < 0.01).sum(axis=1)
    num_tools['Num_Tools_0_02'] = (common_best_ranks_numeric < 0.02).sum(axis=1)
    num_tools['Num_Tools_0_03'] = (common_best_ranks_numeric < 0.03).sum(axis=1)
    num_tools['Num_Tools_0_04'] = (common_best_ranks_numeric < 0.04).sum(axis=1)
    num_tools['Num_Tools_0_05'] = (common_best_ranks_numeric < 0.05).sum(axis=1)
    num_tools['Num_Tools_0_10'] = (common_best_ranks_numeric < 0.1).sum(axis=1)
    num_tools['Num_Tools_0_50'] = (common_best_ranks_numeric < 0.5).sum(axis=1)
    num_tools['Num_Tools_1'] = (common_best_ranks_numeric < 1).sum(axis=1)
    num_tools['Num_Tools_2'] = (common_best_ranks_numeric < 2).sum(axis=1)

    # ---- assemble output: raw columns + computed summary columns
    result = df[['peptide'] + raw_rank_cols].copy()
    result['best_rank'] = best_rank
    result['avg_rank'] = avg_rank
    result['best_allele'] = best_binding_allele

    result = result.merge(strong_counts, on='peptide', how='left')
    result = result.merge(weak_counts, on='peptide', how='left')
    result['strong_binders_count'] = result['strong_binders_count'].fillna(0).astype(int)
    result['weak_binders_count'] = result['weak_binders_count'].fillna(0).astype(int)

    result = result.merge(all_strong_binders.rename('all_strong_binders'), on='peptide', how='left')
    result = result.merge(all_weak_binders.rename('all_weak_binders'), on='peptide', how='left')
    result['all_strong_binders'] = result['all_strong_binders'].fillna('')
    result['all_weak_binders'] = result['all_weak_binders'].fillna('')
    
    for col in num_tools.columns:
        result[col] = num_tools[col].values

    for col_name, values in tool_output_cols.items():
        result[col_name] = values.values

    return result


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Consolidate allele-batched HLA prediction TSV(s) into a single processed feature CSV."
    )
    parser.add_argument('tsv_files', nargs='+', help="One or more allele-batch hla-predictor TSV files (same peptide set).")
    parser.add_argument('--all-alleles', required=True, help="Text file, one allele per line.")
    parser.add_argument('--common-alleles', required=True, help="Text file, one allele per line.")
    parser.add_argument('--output', required=True, help="Path to write the consolidated CSV to.")
    parser.add_argument('--tools', nargs='+', default=DEFAULT_TOOLS, help=f"Prediction tools to use (default: {DEFAULT_TOOLS}).")
    return parser.parse_args(argv)


def main(argv=None) -> None:
    args = parse_args(argv)

    all_alleles = read_allele_list(args.all_alleles)
    common_alleles = read_allele_list(args.common_alleles)
    print(f"Loaded {len(all_alleles)} all_alleles, {len(common_alleles)} common_alleles.")

    result = build_consolidated_hla_table(args.tsv_files, args.tools, all_alleles, common_alleles)

    result.to_csv(args.output, index=False)
    print(f"Saved {result.shape[0]} peptides, {result.shape[1]} columns to {args.output}")


if __name__ == '__main__':
    main()