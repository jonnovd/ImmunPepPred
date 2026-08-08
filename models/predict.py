#!/usr/bin/env python
"""
Predict peptide immunogenicity from a pre-computed feature table using a
saved scikit-learn pipeline (SVM / RF / HGB) and write a prioritised
CSV of results.
"""

import argparse
import json
import sys

import joblib
import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Predict peptide immunogenicity from a feature table using a saved model."
    )
    parser.add_argument(
        "-f", "--feature_table", required=True,
        help="Path to input feature_table.csv."
    )
    parser.add_argument(
        "-m", "--model", required=True,
        help="Path to saved model pipeline (.joblib)."
    )
    parser.add_argument(
        "-d", "--metadata", required=True,
        help="Path to model metadata JSON file (contains feature list)."
    )
    parser.add_argument(
        "-o", "--output", default="immunogenicity_prioritisation.csv",
        help="Path to output CSV (default: immunogenicity_prioritisation.csv)."
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Output is merged prediction file and feature table"
    )

    return parser.parse_args()


def load_metadata(metadata_path):
    with open(metadata_path, "r") as f:
        metadata = json.load(f)

    if "features" not in metadata:
        sys.exit(
            f"ERROR: metadata file '{metadata_path}' does not contain a 'features' key."
        )

    features = metadata["features"]
    if not isinstance(features, list) or len(features) == 0:
        sys.exit(
            f"ERROR: 'features' in metadata file '{metadata_path}' is empty or not a list."
        )

    return features


def load_feature_table(feature_table_path, features):
    df = pd.read_csv(feature_table_path)

    #if 'best_rank' in df.columns:
        #df.rename(columns={'best_rank' : 'hla_best_rank'}, inplace=True)

    if 'length' not in df.columns:
        df['length'] = df['peptide'].str.len()

    required_extra_cols = ["peptide"]#, "best_allele"]
    missing_extra = [c for c in required_extra_cols if c not in df.columns]
    if missing_extra:
        sys.exit(
            f"ERROR: feature_table '{feature_table_path}' is missing required "
            f"non-feature column(s): {missing_extra}"
        )

    missing_features = [f for f in features if f not in df.columns]
    if missing_features:
        sys.exit(
            f"ERROR: feature_table '{feature_table_path}' is missing feature "
            f"column(s) required by the model: {missing_features}"
        )

    return df


def main():
    args = parse_args()

    features = load_metadata(args.metadata)
    df = load_feature_table(args.feature_table, features)

    # Build the numpy array in the exact feature order the model expects
    X = df[features].to_numpy()

    # Load pipeline
    pipe = joblib.load(args.model)

    # Identify the classifier step's classes_ to locate the positive (1) class
    # regardless of internal class ordering.
    final_estimator = pipe.steps[-1][1]
    if not hasattr(final_estimator, "predict_proba"):
        sys.exit(
            f"ERROR: loaded model's final step ({type(final_estimator).__name__}) "
            "does not support predict_proba. If this is an SVC, re-save it with "
            "probability=True."
        )

    classes = final_estimator.classes_
    if len(classes) != 2:
        sys.exit(
            f"ERROR: this script only supports binary classification, but the "
            f"loaded model has classes: {list(classes)}"
        )
    if 1 not in classes:
        sys.exit(
            f"ERROR: expected class label '1' (immunogenic) among model classes, "
            f"but found: {list(classes)}"
        )
    positive_idx = list(classes).index(1)

    proba = pipe.predict_proba(X)
    prob_immunogenic = proba[:, positive_idx]
    prediction = pipe.predict(X)

    #Build output table
    out_df = pd.DataFrame({
        "peptide": df["peptide"],
        "prediction": prediction,
        "probability_immunogenic": prob_immunogenic,
        #"best_binding_allele": df["best_allele"],
    })
    if args.verbose:
        out_df = pd.concat([out_df, df.reset_index(drop=True)], axis=1)
    else:
        # Append feature columns
        out_df = pd.concat([out_df, df[features].reset_index(drop=True)], axis=1)

    # Sort most to least likely immunogenic
    out_df = out_df.sort_values("probability_immunogenic", ascending=False).reset_index(drop=True)

    out_df.to_csv(args.output, index=False)
    print(f"Wrote predictions for {len(out_df)} peptides to '{args.output}'")


if __name__ == "__main__":
    main()