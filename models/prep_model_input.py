"""
prepare_data.py

Prepares and concatenates input data files for ML model input.

Usage:
    python prepare_data.py file1.csv file2.csv ... --output output.csv

File naming conventions:
    - Files containing 'PRIME' in their name → processed with process_prime_file()
    - Files containing 'hla' in their name → processed with process_hla_file()
    - All other files                       → processed with process_default_file()

All files are merged on the 'peptide' column into a single output DataFrame.
"""

import argparse
import pandas as pd
from pathlib import Path
from functools import reduce


# ---------------------------------------------------------------------------
# File-specific processing functions
# ---------------------------------------------------------------------------

def process_prime_file(filepath: str) -> pd.DataFrame:
    """
    """
    df = pd.read_csv(filepath, sep='\t', comment='#')
    df = df[['Peptide', '%Rank_bestAllele']]
    df.rename(columns={'Peptide': 'peptide', '%Rank_bestAllele': 'PRIME_%Rank_bestAllele'}, inplace=True)

    _validate_peptide_column(df, filepath)
    df = df.drop_duplicates(subset=['peptide'])
    return df


def process_hla_file(filepath: str) -> pd.DataFrame:
    """
    """
    df = pd.read_csv(filepath)
    df = df[['peptide', 'best_rank']]
    df.rename(columns={'best_rank':'hla_best_rank'}, inplace=True)

    _validate_peptide_column(df, filepath)
    df = df.drop_duplicates(subset=['peptide'])
    return df

def process_hla_file_extended_features(filepath: str) -> pd.DataFrame:
    """
    """
    df = pd.read_csv(filepath)
    df = df[['peptide', 'best_rank', 'avg_rank', 'best_netmhc_r', 'best_mixmhc_r', 'best_mhcflurry_r', 'best_mhcnuggets_r','weak_binders_count', 'strong_binders_count']]
    df.rename(columns={'best_rank':'hla_best_rank'}, inplace=True)

    _validate_peptide_column(df, filepath)
    df = df.drop_duplicates(subset=['peptide'])
    return df


def process_default_file(filepath: str) -> pd.DataFrame:
    """
    Load and process any file that is not TAP or HLA.

    Args:
        filepath: Path to the input file.

    Returns:
        Processed DataFrame with at least a 'peptide' column.
    """
    df = pd.read_csv(filepath)

    if 'pathogenicity' in Path(filepath).name.lower():
        df.rename(columns={'similarity_score':'pathogenicity'}, inplace=True)
    elif 'selfsimilarity' in Path(filepath).name.lower():
        df.rename(columns={'similarity_score':'selfsimilarity'}, inplace=True)
    elif 'tap' in Path(filepath).name.lower():
        df.rename(columns={'pred_score':'tap_score'}, inplace=True)

    _validate_peptide_column(df, filepath)
    df = df.drop_duplicates(subset=['peptide'])
    return df


# ---------------------------------------------------------------------------
# Routing logic
# ---------------------------------------------------------------------------

def load_and_process_file(filepath: str) -> pd.DataFrame:
    """
    Route a file to the correct processing function based on its name.

    Args:
        filepath: Path to the input file.

    Returns:
        Processed DataFrame.
    """
    name = Path(filepath).name.lower()

    if "hla" in name:
        print(f"  [HLA]     {filepath}")
        return process_hla_file_extended_features(filepath)
        return process_hla_file(filepath)
    elif "prime" in name:
        print(f"  [PRIME]     {filepath}")
        return process_prime_file(filepath)
    else:
        print(f"  [DEFAULT] {filepath}")
        return process_default_file(filepath)


# ---------------------------------------------------------------------------
# Merging DFs
# ---------------------------------------------------------------------------

def merge_on_peptide(dataframes: list[pd.DataFrame]) -> pd.DataFrame:
    """
    Merge a list of DataFrames on the 'peptide' column.

    Each DataFrame is joined on the 'peptide' key, producing one row
    per unique peptide with columns from all input files side by side.

    Args:
        dataframes: List of processed DataFrames, each with a 'peptide' column.

    Returns:
        Single merged DataFrame with 'peptide' as the key column.
    """
    combined = reduce(
        lambda left, right: pd.merge(left, right, on='peptide', how='outer'),
        dataframes
    )

    return combined


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
# CLI entry point
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Concatenate input files for ML model preparation."
    )
    parser.add_argument(
        "input_files",
        nargs="+",
        help="One or more input files to process. Use '+' as a separator between feature groups.",
    )
    parser.add_argument(
        "--output",
        default="prepared_data.csv",
        help="Path for the output CSV file (default: prepared_data.csv).",
    )
    return parser.parse_args()


def split_into_groups(file_list: list[str], separator: str = "+") -> list[list[str]]:
    """
    Split a flat list of file paths into groups using a separator token.

    Example:
        ["a.csv", "b.csv", "+", "c.csv", "+", "d.csv", "e.csv"]
        → [["a.csv", "b.csv"], ["c.csv"], ["d.csv", "e.csv"]]
    """
    groups, current = [], []
    for item in file_list:
        if item == separator:
            if current:
                groups.append(current)
                current = []
        else:
            current.append(item)
    if current:
        groups.append(current)
    return groups


def load_and_concat_group(filepaths: list[str]) -> pd.DataFrame:
    """
    Load all files in a feature group, concatenate them, and deduplicate.

    Args:
        filepaths: Files that all represent the same feature.

    Returns:
        A single deduplicated DataFrame for that feature.
    """
    dataframes = [load_and_process_file(fp) for fp in filepaths]
    concatenated = pd.concat(dataframes, ignore_index=True)
    before = len(concatenated)
    concatenated = concatenated.drop_duplicates(subset=["peptide"])
    after = len(concatenated)
    if before != after:
        print(f"    Dropped {before - after} duplicate peptides after concat.")
    return concatenated


def main() -> None:
    args = parse_args()
    groups = split_into_groups(args.input_files)

    print(f"Processing {len(groups)} feature group(s)…")
    feature_dataframes = []
    for i, group in enumerate(groups):
        print(f"  Group {i + 1} ({len(group)} file(s)):")
        feature_dataframes.append(load_and_concat_group(group))

    print("Merging feature groups…")
    combined = merge_on_peptide(feature_dataframes)

    combined = combined[~combined["peptide"].str.contains("pep")]
    combined.to_csv(args.output, index=False)
    print(f"Done. Output written to '{args.output}' ({len(combined)} rows).")



if __name__ == "__main__":
    main()