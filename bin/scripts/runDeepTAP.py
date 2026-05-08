#!/usr/bin/env python3
"""
Wrapper script for DeepTAP that processes large peptide files in chunks
to avoid memory issues. Runs DeepTAP serially on chunks of 200k peptides
and concatenates results sorted by descending pred_score.

Usage:
    python run_deeptap.py -t cla -f <peptides.csv> -o <output_dir> --deeptap <path/to/deeptap.py>

Mirrors the original DeepTAP interface with an additional --deeptap argument.
"""

import os
import sys
import math
import argparse
import subprocess
import pandas as pd
import tempfile

CHUNK_SIZE = 200000

def parse_args():
    parser = argparse.ArgumentParser(description="Chunked DeepTAP wrapper")
    parser.add_argument("-t", "--taskType",  required=True,  help="Task type (e.g. cla)")
    parser.add_argument("-f", "--file",      required=True,  help="Input CSV file with 'peptide' column")
    parser.add_argument("-o", "--outputDir", required=True,  help="Output directory")
    parser.add_argument("--deeptap",         required=True,  help="Path to the DeepTAP main python script")
    return parser.parse_args()


def split_input(input_file, chunk_dir, chunk_size):
    """Split input CSV into chunks of chunk_size rows, return list of chunk file paths."""
    df = pd.read_csv(input_file)

    if "peptide" not in df.columns:
        sys.exit(f"[ERROR] Input file '{input_file}' has no 'peptide' column. "
                 f"Found columns: {list(df.columns)}")

    n_chunks = math.ceil(len(df) / chunk_size)
    chunk_paths = []

    for i in range(n_chunks):
        chunk_df = df.iloc[i * chunk_size : (i + 1) * chunk_size]
        chunk_path = os.path.join(chunk_dir, f"chunk_{i:04d}.csv")
        chunk_df.to_csv(chunk_path, index=False)
        chunk_paths.append(chunk_path)

    print(f"[INFO] Split {len(df):,} peptides into {n_chunks} chunks of up to {chunk_size:,}")
    return chunk_paths


def run_deeptap_on_chunk(deeptap_script, task_type, chunk_path, chunk_out_dir):
    """Run DeepTAP on a single chunk CSV. Returns path to the ranked output CSV."""
    cmd = [
        sys.executable, deeptap_script,
        "-t", task_type,
        "-f", chunk_path,
        "-o", chunk_out_dir,
    ]

    print(f"[INFO] Running DeepTAP on {os.path.basename(chunk_path)} ...")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"[ERROR] DeepTAP failed on chunk {chunk_path}")
        print(result.stderr)
        sys.exit(1)

    # Reconstruct expected output filename using DeepTAP naming convention:
    # <input_basename>_DeepTAP_<taskType>_predresult_rank.csv
    base_name = os.path.basename(chunk_path).rsplit(".", 1)[0]
    expected_output = os.path.join(
        chunk_out_dir, f"{base_name}_DeepTAP_{task_type}_predresult_rank.csv"
    )

    if not os.path.exists(expected_output):
        sys.exit(f"[ERROR] Expected DeepTAP output not found: {expected_output}")

    return expected_output


def concatenate_and_sort(result_files, output_path):
    """
    Concatenate all chunk result CSVs, sort by pred_score descending,
    drop the last column, and write to output_path.
    """
    dfs = []
    for f in result_files:
        df = pd.read_csv(f)
        dfs.append(df)

    combined = pd.concat(dfs, ignore_index=True)

    # Sort by pred_score descending to match DeepTAP ranked output convention
    combined = combined.sort_values(by="pred_score", ascending=False).reset_index(drop=True)

    # Drop the last column (pred_label) as requested
    combined = combined.iloc[:, :-1]

    combined.to_csv(output_path, index=False)
    print(f"[INFO] Combined output written to: {output_path} ({len(combined):,} rows)")


def main():
    args = parse_args()

    # Validate inputs
    if not os.path.isfile(args.file):
        sys.exit(f"[ERROR] Input file not found: {args.file}")
    if not os.path.isfile(args.deeptap):
        sys.exit(f"[ERROR] DeepTAP script not found: {args.deeptap}")

    os.makedirs(args.outputDir, exist_ok=True)

    # Derive output filename mirroring DeepTAP convention
    input_basename = os.path.basename(args.file).rsplit(".", 1)[0]
    final_output = os.path.join(
        args.outputDir,
        f"{input_basename}_DeepTAP_{args.taskType}_predresult_rank.csv"
    )

    # Use a temporary directory for all intermediate chunk files
    with tempfile.TemporaryDirectory(prefix="deeptap_chunks_") as tmp_dir:
        chunk_input_dir = os.path.join(tmp_dir, "inputs")
        chunk_output_dir = os.path.join(tmp_dir, "outputs")
        os.makedirs(chunk_input_dir)
        os.makedirs(chunk_output_dir)

        # Split input into chunks
        chunk_paths = split_input(args.file, chunk_input_dir, CHUNK_SIZE)

        # Run DeepTAP on each chunk serially
        result_files = []
        for i, chunk_path in enumerate(chunk_paths):
            print(f"[INFO] Processing chunk {i + 1}/{len(chunk_paths)} ...")
            result_file = run_deeptap_on_chunk(
                args.deeptap, args.taskType, chunk_path, chunk_output_dir
            )
            result_files.append(result_file)

        # Concatenate, sort, and write final output
        concatenate_and_sort(result_files, final_output)

    print("[INFO] Done.")


if __name__ == "__main__":
    main()