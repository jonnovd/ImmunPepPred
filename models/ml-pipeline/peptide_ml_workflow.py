#!/usr/bin/env python3
"""
peptide_ml_workflow.py

Single entry point that consolidates the peptide-immunogenicity ML workflow
(previously spread across a notebook + predict.py + filter_peptide_predictions.py
+ visualiseClassificationThreshold.py) into one script with three subcommands:

    train     - grid-search + repeated-stratified-CV train/evaluate one or more
                sklearn models on one or more training sets, save fitted
                pipelines + metadata + metric figures.
    predict   - load a saved pipeline + metadata and score a feature table
                (equivalent to the old predict.py). No training happens here,
                so this is the piece you re-run alone on a validation/holdout
                set, or on your HPC after `train` has already produced models.
    evaluate  - take a predictions CSV (from `predict`) + ground-truth peptide
                lists and produce the annotated CSV + capture-summary figure
                (milestone bars / cumulative capture curve / per-bin
                immunogenic-count histogram) + classification-threshold
                figure. Equivalent to filter_peptide_predictions.py +
                visualiseClassificationThreshold.py combined.

A fourth convenience subcommand, `run`, chains train -> predict -> evaluate
for everything described in a config file in one go (intended for a single
HPC batch job that does the whole workflow unattended).

Every subcommand can be driven either by:
  --config some.yaml     (batch mode: many training sets / models / validation
                           sets in one file, see example_config.yaml)
or by explicit CLI flags for a single training set / single model / single
validation set (quick, interactive use).

Run `python peptide_ml_workflow.py <subcommand> --help` for the flags of each
subcommand.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import sklearn

import matplotlib
matplotlib.use("Agg")  # headless-safe for HPC / no-display nodes
import matplotlib.pyplot as plt

try:
    import seaborn as sns
except ImportError:
    sns = None

from sklearn.calibration import CalibratedClassifierCV
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import fbeta_score, make_scorer
from sklearn.model_selection import GridSearchCV, RepeatedStratifiedKFold, cross_validate
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import RobustScaler
from sklearn.svm import SVC, LinearSVC

try:
    import yaml
except ImportError:
    yaml = None


# ======================================================================
# Model registry
# ----------------------------------------------------------------------
# Every model is a Pipeline whose final step is always named "clf", so
# downstream code (save/load/predict) never needs to know which model it's
# dealing with. Two models (svm_linear, svm_rbf) need a "finalize" step
# after grid search because the estimator used *during* grid search
# (LinearSVC / plain SVC, both fast) doesn't expose predict_proba, which
# `predict` requires. finalize() rebuilds the model with the winning
# hyperparameters in a form that does expose predict_proba, and refits it
# on the full training set. Models with no finalize function use
# GridSearchCV's best_estimator_ directly (it is already refit on the full
# data because refit=<metric name> is truthy).
# ======================================================================

def _strip_prefix(params: dict, prefix: str = "clf__") -> dict:
    return {k[len(prefix):]: v for k, v in params.items() if k.startswith(prefix)}


def _search_pipeline_logreg():
    return Pipeline([
        ("imputer", SimpleImputer(strategy="median", add_indicator=False)),
        ("scaler", RobustScaler()),
        ("clf", LogisticRegression(random_state=0, max_iter=5000, class_weight="balanced")),
    ])


def _search_pipeline_svm_linear():
    return Pipeline([
        ("imputer", SimpleImputer(strategy="median", add_indicator=False)),
        ("scaler", RobustScaler()),
        ("clf", LinearSVC(random_state=0, max_iter=10000)),
    ])


def _finalize_svm_linear(best_params, X, y):
    c = best_params.get("C", 1)
    calibrated = CalibratedClassifierCV(LinearSVC(C=c, random_state=0, max_iter=10000), cv=5)
    pipe = Pipeline([
        ("imputer", SimpleImputer(strategy="median", add_indicator=False)),
        ("scaler", RobustScaler()),
        ("clf", calibrated),
    ])
    pipe.fit(X, y)
    return pipe


def _search_pipeline_svm_rbf():
    return Pipeline([
        ("imputer", SimpleImputer(strategy="median", add_indicator=False)),
        ("scaler", RobustScaler()),
        ("clf", SVC(kernel="rbf", random_state=0)),
    ])


def _finalize_svm_rbf(best_params, X, y):
    c = best_params.get("C", 1)
    gamma = best_params.get("gamma", "scale")
    pipe = Pipeline([
        ("imputer", SimpleImputer(strategy="median", add_indicator=False)),
        ("scaler", RobustScaler()),
        ("clf", SVC(kernel="rbf", C=c, gamma=gamma, probability=True, random_state=0)),
    ])
    pipe.fit(X, y)
    return pipe


def _search_pipeline_rf():
    return Pipeline([
        ("imputer", SimpleImputer(strategy="median", add_indicator=False)),
        ("clf", RandomForestClassifier(random_state=0, class_weight="balanced", n_jobs=1)),
    ])


def _search_pipeline_hgb():
    # HistGradientBoostingClassifier handles NaNs natively, no imputer needed.
    return Pipeline([("clf", HistGradientBoostingClassifier(random_state=0))])


MODEL_REGISTRY = {
    "logreg": {
        "search_pipeline": _search_pipeline_logreg,
        "param_grid": {"clf__C": [0.1, 1, 5, 10]},
        "finalize": None,
    },
    "svm_linear": {
        "search_pipeline": _search_pipeline_svm_linear,
        "param_grid": {"clf__C": [0.1, 1, 5, 10, 100]},
        "finalize": _finalize_svm_linear,
    },
    "svm_rbf": {
        "search_pipeline": _search_pipeline_svm_rbf,
        "param_grid": {"clf__C": [0.1, 1], "clf__gamma": [0.01, 0.1, 1]},
        "finalize": _finalize_svm_rbf,
    },
    "rf": {
        "search_pipeline": _search_pipeline_rf,
        "param_grid": {
            "clf__n_estimators": [200, 300],
            "clf__max_depth": [20, 30, 40],
            "clf__min_samples_leaf": [1, 5],
        },
        "finalize": None,
    },
    "hgb": {
        "search_pipeline": _search_pipeline_hgb,
        "param_grid": {
            "clf__max_iter": [100, 200],
            "clf__max_depth": [10, 20, 30],
            "clf__learning_rate": [0.05, 0.1],
        },
        "finalize": None,
    },
}

DEFAULT_MODELS = list(MODEL_REGISTRY.keys())

DEFAULT_SCORING = {
    "precision": "precision",
    "recall": "recall",
    "roc_auc": "roc_auc",
    "avg_precision": "average_precision",
    "f0.5": make_scorer(fbeta_score, beta=0.5),
}
DEFAULT_REFIT_METRIC = "precision"
DEFAULT_N_SPLITS = 5
DEFAULT_N_REPEATS = 2
DEFAULT_RANDOM_STATE = 42


# ======================================================================
# Shared data loading (feature tables, peptide label lists)
# ======================================================================

def load_peptide_set(path: str) -> set:
    with open(path) as f:
        return {line.strip() for line in f if line.strip()}


def load_feature_table(path: str, features: list[str], require_peptide_col: bool = True) -> pd.DataFrame:
    """Load a feature-table CSV and normalise it for use with this workflow.

    Applies the same light normalisation predict.py used to apply at
    inference time (best_rank -> hla_best_rank rename, auto-derive `length`
    from the peptide string) so training and prediction stay consistent,
    then checks that every column in `features` is actually present.
    """
    df = pd.read_csv(path)

    if "peptide" not in df.columns and require_peptide_col:
        sys.exit(f"ERROR: feature table '{path}' has no 'peptide' column.")

    if "best_rank" in df.columns and "hla_best_rank" not in df.columns:
        df.rename(columns={"best_rank": "hla_best_rank"}, inplace=True)

    if "length" in features and "length" not in df.columns:
        df["length"] = df["peptide"].str.len()

    missing = [f for f in features if f not in df.columns]
    if missing:
        sys.exit(f"ERROR: feature table '{path}' is missing feature column(s) required "
                  f"by the supplied feature list: {missing}")

    return df


def load_feature_list(path: str) -> list[str]:
    """One feature/column name per line, blank lines and '#' comments ignored."""
    feats = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                feats.append(line)
    if not feats:
        sys.exit(f"ERROR: feature list file '{path}' contained no feature names.")
    return feats


def build_labelled_dataset(feature_table_path: str, features: list[str],
                            immunogenic_peptides_path: str, non_immunogenic_peptides_path: str):
    """Replicates the notebook's data-prep cells: load feature table, load
    immunogenic/non-immunogenic peptide lists, restrict both to their
    intersection, and return (X_df, X, y, feature_list).
    """
    df = load_feature_table(feature_table_path, features)
    df = df.drop_duplicates(subset=["peptide"]).sort_values(by="peptide").reset_index(drop=True)

    imm_peps = load_peptide_set(immunogenic_peptides_path)
    non_imm_peps = load_peptide_set(non_immunogenic_peptides_path)

    overlap = imm_peps & non_imm_peps
    if overlap:
        print(f"WARNING: {len(overlap)} peptide(s) appear in both the immunogenic and "
              f"non-immunogenic lists; they will be excluded (ambiguous label).",
              file=sys.stderr)
        imm_peps = imm_peps - overlap
        non_imm_peps = non_imm_peps - overlap

    labels_df = pd.DataFrame(
        [{"peptide": p, "immunogenic": 1} for p in imm_peps] +
        [{"peptide": p, "immunogenic": 0} for p in non_imm_peps]
    )
    labels_df = labels_df.drop_duplicates(subset=["peptide"]).sort_values(by="peptide")

    label_peps = imm_peps | non_imm_peps
    df = df[df["peptide"].isin(label_peps)].reset_index(drop=True)
    labels_df = labels_df[labels_df["peptide"].isin(set(df["peptide"]))].reset_index(drop=True)

    # Align labels to df's peptide order exactly.
    labels_df = labels_df.set_index("peptide").loc[df["peptide"]].reset_index()

    if len(df) == 0:
        sys.exit("ERROR: no peptides in the feature table matched either ground-truth list.")

    print(f"Loaded {len(df)} labelled peptides "
          f"({int((labels_df['immunogenic'] == 1).sum())} immunogenic, "
          f"{int((labels_df['immunogenic'] == 0).sum())} non-immunogenic).")

    X = df[features].to_numpy()
    y = labels_df["immunogenic"].to_numpy().ravel()
    assert len(df) == X.shape[0] == len(y)
    return df, X, y


# ======================================================================
# Train + grid-search-CV
# ======================================================================

def build_cv_folds(X, y, lengths, n_splits, n_repeats, random_state):
    """Composite stratification key = class label + exact peptide length, so
    every model is evaluated on identical folds (fair comparison), matching
    the notebook.
    """
    strata = pd.Series(y).astype(str) + "_" + pd.Series(lengths).astype(str)
    cv = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=random_state)
    return list(cv.split(X, strata))


def run_model_grid_search(model_name, X, y, cv_folds, param_grid, scoring, refit_metric, n_jobs=-1):
    spec = MODEL_REGISTRY[model_name]
    pipe = spec["search_pipeline"]()
    grid = param_grid if param_grid is not None else spec["param_grid"]

    print("=" * 60)
    print(f"Running {model_name} (grid search, {len(cv_folds)} folds, "
          f"{len(list(_iter_grid(grid)))} param combinations)")

    clf = GridSearchCV(
        estimator=pipe, param_grid=grid, cv=cv_folds, scoring=scoring,
        refit=refit_metric, n_jobs=n_jobs, return_train_score=True,
    )
    clf.fit(X, y)
    best_index = clf.best_index_
    best_params_clf = _strip_prefix(clf.best_params_)

    rows = []
    for metric_name in scoring:
        fold_scores = [clf.cv_results_[f"split{i}_test_{metric_name}"][best_index]
                        for i in range(len(cv_folds))]
        for fold_idx, score in enumerate(fold_scores):
            rows.append({"model": model_name, "fold": fold_idx, "metric": metric_name,
                         "score": score, "params": str(clf.best_params_)})
    results_df = pd.DataFrame(rows)

    finalize = spec["finalize"]
    if finalize is not None:
        final_pipe = finalize(best_params_clf, X, y)
    else:
        final_pipe = clf.best_estimator_  # already refit on all of X,y

    print(f"Done: {model_name} (best params: {clf.best_params_})")
    return results_df, final_pipe, clf.best_params_


def _iter_grid(param_grid):
    from sklearn.model_selection import ParameterGrid
    return ParameterGrid(param_grid)


def plot_models_boxplot(results_dfs, output_path, title_suffix=""):
    if sns is None:
        print("WARNING: seaborn not installed, skipping boxplot figure.", file=sys.stderr)
        return
    combined = pd.concat(results_dfs, ignore_index=True)
    metrics = combined["metric"].unique()
    fig, axes = plt.subplots(1, len(metrics), figsize=(5 * len(metrics), 5), sharey=True)
    if len(metrics) == 1:
        axes = [axes]
    for ax, metric_name in zip(axes, metrics):
        subset = combined[combined["metric"] == metric_name]
        sns.boxplot(data=subset, x="model", y="score", hue="model", ax=ax, legend=False)
        sns.stripplot(data=subset, x="model", y="score", color="black", size=5, alpha=0.6, ax=ax)
        ax.set_title(metric_name)
        ax.set_xlabel("")
        ax.set_ylabel("score" if ax is axes[0] else "")
        ax.tick_params(axis="x", rotation=45)
    plt.suptitle(f"Model performance by metric — CV{title_suffix}")
    plt.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def plot_models_meanbar(results_dfs, output_path, error_bars="sd", capsize=0.15, title_suffix=""):
    if sns is None:
        print("WARNING: seaborn not installed, skipping mean-bar figure.", file=sys.stderr)
        return
    combined = pd.concat(results_dfs, ignore_index=True)
    metrics = combined["metric"].unique()
    y_min, y_max = combined["score"].min(), combined["score"].max()
    padding = (y_max - y_min) * 0.1 or 0.05
    fig, axes = plt.subplots(1, len(metrics), figsize=(5 * len(metrics), 5), sharey=True)
    if len(metrics) == 1:
        axes = [axes]
    for ax, metric_name in zip(axes, metrics):
        subset = combined[combined["metric"] == metric_name]
        sns.barplot(data=subset, x="model", y="score", hue="model", ax=ax, legend=False,
                    errorbar=error_bars, capsize=capsize, width=0.6)
        ax.set_title(metric_name)
        ax.set_xlabel("")
        ax.set_ylabel("mean score" if ax is axes[0] else "")
        ax.tick_params(axis="x", rotation=45)
    axes[0].set_ylim(y_min - padding, y_max + padding)
    plt.suptitle(f"Mean model performance by metric — CV{title_suffix}, error bars = {error_bars}")
    plt.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def save_model(pipe, best_params, results_df, output_dir: Path, model_name: str,
                training_set_name: str, features: list[str], n_train_samples: int):
    output_dir.mkdir(parents=True, exist_ok=True)
    model_path = output_dir / f"{model_name}_pipeline.joblib"
    meta_path = output_dir / f"{model_name}_pipeline_metadata.json"

    joblib.dump(pipe, model_path)

    metadata = {
        "sklearn_version": sklearn.__version__,
        "python_version": sys.version,
        "saved_at": datetime.now(timezone.utc).isoformat(),
        "features": list(features),
        "model_name": model_name,
        "training_set": training_set_name,
        "n_train_samples": n_train_samples,
        "best_params": best_params,
        "params": pipe.named_steps["clf"].get_params(),
        "cv_results_summary": results_df.groupby("metric")["score"].mean().to_dict(),
    }
    with open(meta_path, "w") as f:
        json.dump(metadata, f, indent=2, default=str)

    print(f"Saved {model_name}: {model_path.name}, {meta_path.name}")
    return model_path, meta_path


def train_one_training_set(ts_cfg: dict, global_cfg: dict, models_filter: list[str] | None,
                             output_root: Path):
    name = ts_cfg["name"]
    print("\n" + "#" * 70)
    print(f"# Training set: {name}")
    print("#" * 70)

    features = load_feature_list(ts_cfg["features_file"])
    df, X, y = build_labelled_dataset(
        ts_cfg["feature_table"], features,
        ts_cfg["immunogenic_peptides"], ts_cfg["non_immunogenic_peptides"],
    )

    cv_cfg = {**global_cfg.get("cv", {}), **ts_cfg.get("cv", {})}
    n_splits = cv_cfg.get("n_splits", DEFAULT_N_SPLITS)
    n_repeats = cv_cfg.get("n_repeats", DEFAULT_N_REPEATS)
    random_state = cv_cfg.get("random_state", DEFAULT_RANDOM_STATE)
    refit_metric = cv_cfg.get("refit_metric", DEFAULT_REFIT_METRIC)

    lengths = df["length"].to_numpy() if "length" in df.columns else df["peptide"].str.len().to_numpy()
    cv_folds = build_cv_folds(X, y, lengths, n_splits, n_repeats, random_state)

    models_cfg = global_cfg.get("models", {})
    ts_models_cfg = ts_cfg.get("models", {})
    model_names = models_filter or ts_cfg.get("model_order") or DEFAULT_MODELS

    out_dir = output_root / "saved_models" / name
    out_dir.mkdir(parents=True, exist_ok=True)

    all_results = []
    for model_name in model_names:
        if model_name not in MODEL_REGISTRY:
            print(f"WARNING: unknown model '{model_name}', skipping.", file=sys.stderr)
            continue
        m_cfg = {**models_cfg.get(model_name, {}), **ts_models_cfg.get(model_name, {})}
        if not m_cfg.get("enabled", True):
            print(f"Skipping {model_name} (disabled in config).")
            continue
        param_grid = m_cfg.get("param_grid")  # None -> registry default

        results_df, final_pipe, best_params = run_model_grid_search(
            model_name, X, y, cv_folds, param_grid, DEFAULT_SCORING, refit_metric,
        )
        save_model(final_pipe, best_params, results_df, out_dir, model_name, name,
                   features, len(df))
        all_results.append(results_df)

    if all_results:
        plot_models_boxplot(all_results, out_dir / "models_cv_metrics.png",
                             title_suffix=f" ({name})")
        plot_models_meanbar(all_results, out_dir / "models_mean_metrics.png",
                             title_suffix=f" ({name})")
        combined = pd.concat(all_results, ignore_index=True)
        summary = combined.groupby(["metric", "model"])["score"].mean()
        print(summary)
        summary.to_csv(out_dir / "models_cv_metrics_summary.csv")
    else:
        print(f"No models were trained for training set '{name}'.")

    return out_dir


# ======================================================================
# Predict (load saved pipeline + metadata, score a feature table)
# ======================================================================

def predict_with_saved_model(model_path: str, metadata_path: str, feature_table_path: str,
                               output_csv: str):
    with open(metadata_path) as f:
        metadata = json.load(f)
    features = metadata["features"]

    df = load_feature_table(feature_table_path, features)
    if "best_allele" not in df.columns:
        df["best_allele"] = np.nan

    X = df[features].to_numpy()
    pipe = joblib.load(model_path)

    final_estimator = pipe.steps[-1][1]
    if not hasattr(final_estimator, "predict_proba"):
        sys.exit(f"ERROR: loaded model's final step ({type(final_estimator).__name__}) does not "
                  f"support predict_proba.")
    classes = list(final_estimator.classes_)
    if len(classes) != 2 or 1 not in classes:
        sys.exit(f"ERROR: expected binary classifier with classes including 1, got: {classes}")
    positive_idx = classes.index(1)

    proba = pipe.predict_proba(X)[:, positive_idx]
    prediction = pipe.predict(X)

    out_df = pd.DataFrame({
        "peptide": df["peptide"],
        "prediction": prediction,
        "probability_immunogenic": proba,
        "best_binding_allele": df["best_allele"],
    })
    out_df = pd.concat([out_df, df[features].reset_index(drop=True)], axis=1)
    out_df = out_df.sort_values("probability_immunogenic", ascending=False).reset_index(drop=True)

    Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(output_csv, index=False)
    print(f"Wrote predictions for {len(out_df)} peptides to '{output_csv}'")
    return out_df


# ======================================================================
# Evaluate (filter/annotate predictions against ground truth + figures)
# ======================================================================

def is_predicted_positive(value):
    if pd.isna(value):
        return False
    if isinstance(value, (int, float)):
        return value == 1
    return str(value).strip().lower() in {"1", "true", "positive"}


def compute_topk_capture(df, immunogenic_col, score_col, fractions, lower_is_better=False):
    total_immunogenic = df[immunogenic_col].sum()
    ranked = df.sort_values(score_col, ascending=lower_is_better, kind="mergesort").reset_index(drop=True)
    results = {}
    n = len(ranked)
    for f in fractions:
        k = max(1, math.ceil(n * f / 100)) if n > 0 else 0
        captured = ranked.iloc[:k][immunogenic_col].sum()
        results[f] = (captured / total_immunogenic * 100) if total_immunogenic > 0 else 0.0
    return results, ranked


def compute_cumulative_capture_curve(ranked_df, immunogenic_col):
    n = len(ranked_df)
    total_immunogenic = ranked_df[immunogenic_col].sum()
    if n == 0 or total_immunogenic == 0:
        return [0], [0]
    cum = ranked_df[immunogenic_col].cumsum()
    x = [(i + 1) / n * 100 for i in range(n)]
    y = [(c / total_immunogenic * 100) for c in cum]
    return [0] + x, [0] + y


def plot_evaluation_summary(bar_data, ranked_df, immunogenic_col, score_col, prediction_pct,
                             bin_width, output_path):
    """3-panel summary figure:
      1. milestone bar chart (predicted-positive capture + top-k% capture)
      2. cumulative capture curve vs random baseline
      3. per-bin histogram of `score_col`, total-peptide bars, with the
         validated-immunogenic count (and %) for that bin labelled above
         each bar.
    """
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5))

    # --- Panel 1: milestone bars ---
    labels = ["Predicted\npositive"] + [f"Top {f}%" for f in bar_data["fractions"]]
    values = [prediction_pct] + [bar_data["topk"][f] for f in bar_data["fractions"]]
    bars = ax1.bar(labels, values, color="#4C72B0")
    ax1.set_ylabel("% of immunogenic peptides captured")
    ax1.set_ylim(0, 105)
    ax1.set_title("Capture of immunogenic peptides\nby prediction milestone")
    for bar, val in zip(bars, values):
        ax1.text(bar.get_x() + bar.get_width() / 2, val + 1.5, f"{val:.1f}%",
                  ha="center", va="bottom", fontsize=9)

    # --- Panel 2: cumulative capture curve ---
    x, y = compute_cumulative_capture_curve(ranked_df, immunogenic_col)
    ax2.plot(x, y, color="#C44E52", linewidth=2, label="Model ranking")
    ax2.plot([0, 100], [0, 100], color="gray", linestyle="--", linewidth=1, label="Random baseline")
    ax2.set_xlabel("% of peptides considered (ranked by score)")
    ax2.set_ylabel("Cumulative % of immunogenic peptides captured")
    ax2.set_title("Cumulative capture curve")
    ax2.set_xlim(0, 100)
    ax2.set_ylim(0, 105)
    ax2.legend(loc="lower right")

    # --- Panel 3: per-bin histogram, total bars with immunogenic bars overlaid in front ---
    lo = float(np.floor(ranked_df[score_col].min() / bin_width) * bin_width)
    hi = float(np.ceil(ranked_df[score_col].max() / bin_width) * bin_width)
    hi = max(hi, lo + bin_width)
    bins = np.arange(lo, hi + bin_width, bin_width)

    total_counts, edges = np.histogram(ranked_df[score_col], bins=bins)
    immuno_scores = ranked_df.loc[ranked_df[immunogenic_col] == 1, score_col]
    immuno_counts, _ = np.histogram(immuno_scores, bins=bins)

    centers = (edges[:-1] + edges[1:]) / 2
    width = bin_width * 0.9
    ax3.bar(centers, total_counts, width=width, color="#d62728", edgecolor="white", linewidth=0.5,
            label="All peptides in bin", zorder=1)
    ax3.bar(centers, immuno_counts, width=width, color="#1f77b4", edgecolor="white", linewidth=0.5,
            label="Validated immunogenic peptides in bin", zorder=2)
    ax3.set_xlabel(score_col)
    ax3.set_ylabel("Number of peptides")
    ax3.set_title("Peptides per score bin\n(overlaid = validated immunogenic peptides)")
    ax3.legend(loc="upper right")

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def plot_classification_threshold(df_labelled, score_col, label_col, threshold, bin_width,
                                   lower_is_better, output_path):
    scores_immuno = df_labelled.loc[df_labelled[label_col] == 1, score_col]
    scores_non_immuno = df_labelled.loc[df_labelled[label_col] == 0, score_col]

    if lower_is_better:
        tp = int((scores_immuno <= threshold).sum())
        fn = int((scores_immuno > threshold).sum())
        fp = int((scores_non_immuno <= threshold).sum())
        tn = int((scores_non_immuno > threshold).sum())
    else:
        tp = int((scores_immuno >= threshold).sum())
        fn = int((scores_immuno < threshold).sum())
        fp = int((scores_non_immuno >= threshold).sum())
        tn = int((scores_non_immuno < threshold).sum())

    recall = tp / (tp + fn) if (tp + fn) else float("nan")
    specificity = tn / (tn + fp) if (tn + fp) else float("nan")
    precision = tp / (tp + fp) if (tp + fp) else float("nan")
    accuracy = (tp + tn) / (tp + tn + fp + fn) if (tp + tn + fp + fn) else float("nan")

    bins = np.arange(0, 1 + bin_width, bin_width)
    scores_immuno_plot = scores_immuno.clip(lower=0, upper=1)
    scores_non_immuno_plot = scores_non_immuno.clip(lower=0, upper=1)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(scores_non_immuno_plot, bins=bins, alpha=0.6, color="#d62728",
            label=f"Actual non-immunogenic (n={len(scores_non_immuno)})", edgecolor="white", linewidth=0.5)
    ax.hist(scores_immuno_plot, bins=bins, alpha=0.6, color="#1f77b4",
            label=f"Actual immunogenic (n={len(scores_immuno)})", edgecolor="white", linewidth=0.5)
    ax.axvline(threshold, color="black", linestyle="--", linewidth=1.5, label=f"Threshold = {threshold:.2f}")
    ax.set_xlabel(score_col, fontsize=12)
    ax.set_ylabel("Number of peptides", fontsize=12)
    ax.set_title("Classification Threshold Visualisation", fontsize=14, fontweight="bold")
    ax.set_xlim(0, 1)
    ax.legend(loc="upper left")

    metrics_text = (f"At threshold = {threshold:.2f}\nTP: {tp}   FN: {fn}\nFP: {fp}   TN: {tn}\n"
                     f"Recall: {recall:.2f}   Specificity: {specificity:.2f}\n"
                     f"Precision: {precision:.2f}   Accuracy: {accuracy:.2f}")
    ax.text(0.98, 0.97, metrics_text, transform=ax.transAxes, fontsize=9, va="top", ha="right",
            bbox=dict(boxstyle="round", facecolor="white", edgecolor="gray", alpha=0.9))

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return {"tp": tp, "fn": fn, "fp": fp, "tn": tn, "recall": recall,
            "specificity": specificity, "precision": precision, "accuracy": accuracy}


def evaluate_predictions(predictions_csv: str, immunogenic_peptides_path: str,
                          non_immunogenic_peptides_path: str, output_dir: Path,
                          peptides_of_interest_path: str | None = None,
                          score_col: str = "probability_immunogenic",
                          bin_width_capture: float = 0.02, threshold: float = 0.5,
                          threshold_bin_width: float = 0.05, lower_is_better: bool = False,
                          prefix: str = "evaluation"):
    output_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(predictions_csv)
    required = {"peptide", score_col}
    missing = required - set(df.columns)
    if missing:
        sys.exit(f"ERROR: predictions CSV is missing required column(s): {missing}")

    if peptides_of_interest_path:
        interest = load_peptide_set(peptides_of_interest_path)
        df = df[df["peptide"].isin(interest)].copy()
    else:
        df = df.copy()

    imm_peps = load_peptide_set(immunogenic_peptides_path)
    non_imm_peps = load_peptide_set(non_immunogenic_peptides_path)
    overlap = imm_peps & non_imm_peps
    if overlap:
        print(f"WARNING: {len(overlap)} peptide(s) in both ground-truth lists, excluded.", file=sys.stderr)

    def label_peptide(pep):
        if pep in imm_peps and pep not in overlap:
            return 1
        if pep in non_imm_peps and pep not in overlap:
            return 0
        return np.nan

    df["validated_immunogenic"] = df["peptide"].apply(label_peptide)
    n_unlabelled = df["validated_immunogenic"].isna().sum()
    if n_unlabelled:
        print(f"NOTE: {n_unlabelled} peptide(s) have no ground-truth label and will be excluded "
              f"from the capture/threshold figures.", file=sys.stderr)

    annotated_csv = output_dir / f"{prefix}_annotated.csv"
    df.to_csv(annotated_csv, index=False)
    print(f"Wrote {len(df)} annotated rows to {annotated_csv}")

    df_labelled = df.dropna(subset=["validated_immunogenic"]).copy()
    df_labelled["validated_immunogenic"] = df_labelled["validated_immunogenic"].astype(int)

    n_immunogenic = int(df_labelled["validated_immunogenic"].sum())
    print(f"{n_immunogenic} validated immunogenic peptide(s) found among {len(df_labelled)} labelled peptides.")
    if n_immunogenic == 0:
        print("No validated immunogenic peptides in the (filtered) set; skipping figures.")
        return annotated_csv, None

    fractions = [10, 25, 50, 70]
    topk_pct, ranked = compute_topk_capture(df_labelled, "validated_immunogenic", score_col, fractions,
                                             lower_is_better=lower_is_better)

    if "prediction" in df_labelled.columns:
        is_pos = df_labelled["prediction"].apply(is_predicted_positive)
    else:
        is_pos = df_labelled[score_col] >= threshold if not lower_is_better else df_labelled[score_col] <= threshold
    n_pred_pos_and_imm = is_pos[df_labelled["validated_immunogenic"] == 1].sum()
    prediction_pct = n_pred_pos_and_imm / n_immunogenic * 100

    summary_fig_path = output_dir / f"{prefix}_capture_summary.png"
    plot_evaluation_summary(
        bar_data={"fractions": fractions, "topk": topk_pct}, ranked_df=ranked,
        immunogenic_col="validated_immunogenic", score_col=score_col, prediction_pct=prediction_pct,
        bin_width=bin_width_capture, output_path=summary_fig_path,
    )
    print(f"Wrote capture-summary figure to {summary_fig_path}")

    threshold_fig_path = output_dir / f"{prefix}_classification_threshold.png"
    metrics = plot_classification_threshold(
        df_labelled, score_col, "validated_immunogenic", threshold, threshold_bin_width,
        lower_is_better, threshold_fig_path,
    )
    print(f"Wrote classification-threshold figure to {threshold_fig_path}")
    print(metrics)

    return annotated_csv, {"summary_figure": str(summary_fig_path),
                            "threshold_figure": str(threshold_fig_path), "metrics": metrics}


# ======================================================================
# Config loading
# ======================================================================

def load_config(path: str) -> dict:
    if yaml is None:
        sys.exit("ERROR: PyYAML is required to use --config. Install with "
                  "`pip install pyyaml --break-system-packages`.")
    with open(path) as f:
        cfg = yaml.safe_load(f)
    if "training_sets" not in cfg:
        sys.exit(f"ERROR: config '{path}' has no 'training_sets' section.")
    return cfg


def resolve_output_root(cfg: dict, cli_output_root: str | None) -> Path:
    root = cli_output_root or cfg.get("output_root", "results")
    return Path(root)


# ======================================================================
# CLI subcommands
# ======================================================================

def cmd_train(args):
    if args.config:
        cfg = load_config(args.config)
        output_root = resolve_output_root(cfg, args.output_dir)
        ts_filter = set(args.training_set) if args.training_set else None
        for ts_cfg in cfg["training_sets"]:
            if ts_filter and ts_cfg["name"] not in ts_filter:
                continue
            train_one_training_set(ts_cfg, cfg, args.models, output_root)
    else:
        _require(args, ["training_set_name", "feature_table", "immunogenic", "non_immunogenic",
                         "features_file", "output_dir"], "train (without --config)")
        ts_cfg = {
            "name": args.training_set_name,
            "feature_table": args.feature_table,
            "immunogenic_peptides": args.immunogenic,
            "non_immunogenic_peptides": args.non_immunogenic,
            "features_file": args.features_file,
        }
        global_cfg = {"cv": {"n_splits": args.n_splits, "n_repeats": args.n_repeats,
                              "random_state": args.random_state, "refit_metric": args.refit_metric}}
        train_one_training_set(ts_cfg, global_cfg, args.models, Path(args.output_dir))


def cmd_predict(args):
    if args.config:
        cfg = load_config(args.config)
        output_root = resolve_output_root(cfg, args.output_dir)
        ts_filter = set(args.training_set) if args.training_set else None
        for ts_cfg in cfg["training_sets"]:
            if ts_filter and ts_cfg["name"] not in ts_filter:
                continue
            validation = ts_cfg.get("validation")
            if not validation:
                print(f"No validation set configured for training set '{ts_cfg['name']}', skipping predict.")
                continue
            model_dir = output_root / "saved_models" / ts_cfg["name"]
            model_names = args.models or ts_cfg.get("model_order") or DEFAULT_MODELS
            for model_name in model_names:
                model_path = model_dir / f"{model_name}_pipeline.joblib"
                meta_path = model_dir / f"{model_name}_pipeline_metadata.json"
                if not model_path.exists():
                    print(f"NOTE: no saved model for {ts_cfg['name']}/{model_name} at {model_path}, skipping.")
                    continue
                out_csv = output_root / "predictions" / ts_cfg["name"] / f"{model_name}_predictions.csv"
                predict_with_saved_model(str(model_path), str(meta_path), validation["feature_table"],
                                          str(out_csv))
    else:
        _require(args, ["model", "metadata", "feature_table", "output"], "predict (without --config)")
        predict_with_saved_model(args.model, args.metadata, args.feature_table, args.output)


def cmd_evaluate(args):
    if args.config:
        cfg = load_config(args.config)
        output_root = resolve_output_root(cfg, args.output_dir)
        ts_filter = set(args.training_set) if args.training_set else None
        for ts_cfg in cfg["training_sets"]:
            if ts_filter and ts_cfg["name"] not in ts_filter:
                continue
            validation = ts_cfg.get("validation")
            if not validation or not validation.get("immunogenic_peptides"):
                print(f"No validation ground-truth configured for '{ts_cfg['name']}', skipping evaluate.")
                continue
            model_names = args.models or ts_cfg.get("model_order") or DEFAULT_MODELS
            for model_name in model_names:
                pred_csv = output_root / "predictions" / ts_cfg["name"] / f"{model_name}_predictions.csv"
                if not pred_csv.exists():
                    print(f"NOTE: no predictions CSV for {ts_cfg['name']}/{model_name} at {pred_csv}, "
                          f"skipping (run `predict` first).")
                    continue
                eval_out = output_root / "evaluation" / ts_cfg["name"]
                evaluate_predictions(
                    str(pred_csv), validation["immunogenic_peptides"], validation["non_immunogenic_peptides"],
                    eval_out, peptides_of_interest_path=validation.get("peptides_of_interest"),
                    score_col=validation.get("score_col", "probability_immunogenic"),
                    bin_width_capture=validation.get("bin_width_capture", 0.02),
                    threshold=validation.get("threshold", 0.5),
                    threshold_bin_width=validation.get("threshold_bin_width", 0.05),
                    lower_is_better=validation.get("lower_is_better", False),
                    prefix=model_name,
                )
    else:
        _require(args, ["predictions_csv", "immunogenic", "non_immunogenic", "output_dir"],
                  "evaluate (without --config)")
        evaluate_predictions(
            args.predictions_csv, args.immunogenic, args.non_immunogenic, Path(args.output_dir),
            peptides_of_interest_path=args.peptides, score_col=args.score_column,
            bin_width_capture=args.bin_width_capture, threshold=args.threshold,
            threshold_bin_width=args.threshold_bin_width, lower_is_better=args.lower_is_better,
            prefix=args.prefix,
        )


def cmd_run(args):
    """Full workflow in one sitting: train -> predict -> evaluate, for
    everything in the config. Intended for a single HPC batch job.
    """
    if not args.config:
        sys.exit("ERROR: `run` requires --config (it drives the full multi-training-set workflow).")
    cmd_train(args)
    cmd_predict(args)
    cmd_evaluate(args)


def _require(args, names, context):
    missing = [n for n in names if getattr(args, n, None) is None]
    if missing:
        sys.exit(f"ERROR: --{'/--'.join(n.replace('_', '-') for n in missing)} required for `{context}` "
                  f"(or pass --config instead).")


# ======================================================================
# Argument parsing
# ======================================================================

def build_parser():
    parser = argparse.ArgumentParser(
        description="Peptide immunogenicity ML workflow: train, predict, evaluate.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # --- train ---
    p_train = sub.add_parser("train", help="Grid-search + CV train models on one or more training sets.")
    p_train.add_argument("--config", help="Batch-mode YAML config (see example_config.yaml).")
    p_train.add_argument("--training-set", action="append",
                          help="Restrict --config run to this training set name (repeatable).")
    p_train.add_argument("--models", nargs="+", choices=list(MODEL_REGISTRY),
                          help="Restrict to these models (default: all in config/registry).")
    p_train.add_argument("--output-dir", help="Root output directory.")
    # single-run (no --config) flags
    p_train.add_argument("--training-set-name")
    p_train.add_argument("--feature-table")
    p_train.add_argument("--immunogenic")
    p_train.add_argument("--non-immunogenic")
    p_train.add_argument("--features-file", help="Text file, one feature/column name per line.")
    p_train.add_argument("--n-splits", type=int, default=DEFAULT_N_SPLITS)
    p_train.add_argument("--n-repeats", type=int, default=DEFAULT_N_REPEATS)
    p_train.add_argument("--random-state", type=int, default=DEFAULT_RANDOM_STATE)
    p_train.add_argument("--refit-metric", default=DEFAULT_REFIT_METRIC)
    p_train.set_defaults(func=cmd_train)

    # --- predict ---
    p_pred = sub.add_parser("predict", help="Score a feature table with a saved pipeline.")
    p_pred.add_argument("--config", help="Batch-mode YAML config.")
    p_pred.add_argument("--training-set", action="append")
    p_pred.add_argument("--models", nargs="+", choices=list(MODEL_REGISTRY))
    p_pred.add_argument("--output-dir")
    # single-run flags
    p_pred.add_argument("--model", help="Path to saved pipeline (.joblib).")
    p_pred.add_argument("--metadata", help="Path to model metadata JSON.")
    p_pred.add_argument("--feature-table", help="Path to feature_table.csv to score.")
    p_pred.add_argument("--output", help="Output predictions CSV path.")
    p_pred.set_defaults(func=cmd_predict)

    # --- evaluate ---
    p_eval = sub.add_parser("evaluate", help="Annotate predictions with ground truth + capture/threshold figures.")
    p_eval.add_argument("--config", help="Batch-mode YAML config.")
    p_eval.add_argument("--training-set", action="append")
    p_eval.add_argument("--models", nargs="+", choices=list(MODEL_REGISTRY))
    p_eval.add_argument("--output-dir")
    # single-run flags
    p_eval.add_argument("--predictions-csv", help="CSV produced by `predict`.")
    p_eval.add_argument("--immunogenic", help="Validated immunogenic peptides, one per line.")
    p_eval.add_argument("--non-immunogenic", help="Validated non-immunogenic peptides, one per line.")
    p_eval.add_argument("--peptides", help="Optional peptides-of-interest subset, one per line "
                                            "(default: all peptides in the predictions CSV).")
    p_eval.add_argument("--score-column", default="probability_immunogenic")
    p_eval.add_argument("--bin-width-capture", type=float, default=0.02,
                         help="Bin width for the per-bin immunogenic-count histogram.")
    p_eval.add_argument("--threshold", type=float, default=0.5)
    p_eval.add_argument("--threshold-bin-width", type=float, default=0.05)
    p_eval.add_argument("--lower-is-better", action="store_true",
                         help="Set for rank-style scores where lower = more immunogenic.")
    p_eval.add_argument("--prefix", default="evaluation", help="Filename prefix for outputs.")
    p_eval.set_defaults(func=cmd_evaluate)

    # --- run (train + predict + evaluate, config-driven only) ---
    p_run = sub.add_parser("run", help="Full workflow (train -> predict -> evaluate) from a config file.")
    p_run.add_argument("--config", required=True)
    p_run.add_argument("--training-set", action="append")
    p_run.add_argument("--models", nargs="+", choices=list(MODEL_REGISTRY))
    p_run.add_argument("--output-dir")
    p_run.set_defaults(func=cmd_run)

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
