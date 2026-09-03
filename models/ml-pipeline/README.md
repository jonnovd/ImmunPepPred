# peptide_ml_workflow.py

One script that replaces the notebook + `predict.py` + `filter_peptide_predictions.py`
+ `visualiseClassificationThreshold.py` combo. Four subcommands:

| Subcommand | Function |
|---|---|
| `train` | Grid-search + repeated-stratified-CV each model, save fitted pipelines + metadata + CV metric figures |
| `predict` | Load a saved pipeline + metadata, score a feature table, write a prioritised CSV |
| `evaluate` | Annotate a predictions CSV with ground truth, write the capture-summary figure (milestone bars + cumulative capture curve + per-bin immunogenic-count histogram) and the classification-threshold figure |
| `run` | `train` -> `predict` -> `evaluate` for everything in a config file, in one go |

Every subcommand takes **either** `--config some.yaml` (batch mode, many
training sets/models/validation sets at once — see `example_config.yaml`)
**or** explicit flags for a single training set / model / validation set
(quick, interactive use, e.g. testing one new model on one training set).

Because `predict` and `evaluate` are separate subcommands from `train`, you
can re-run *only* the validation part later against any saved
`{model}_pipeline.joblib` + `{model}_pipeline_metadata.json` pair, without
retraining anything.

## Setup

```bash
pip install scikit-learn pandas numpy matplotlib seaborn joblib pyyaml
```

## Quick single-run examples

Train one model on one training set:
```bash
python peptide_ml_workflow.py train \
  --training-set-name Imm_Non-Imm \
  --feature-table data/input/formatted-input/model_input.csv \
  --immunogenic data/input/cedar_imm_peps_traintest.txt \
  --non-immunogenic data/input/cedar_non_imm_peps_traintest.txt \
  --features-file features/imm_non-imm_features.txt \
  --models logreg rf \
  --output-dir results
```
This writes `results/saved_models/Imm_Non-Imm/{logreg,rf}_pipeline.joblib` (+
`_metadata.json`) and `models_cv_metrics.png` / `models_mean_metrics.png`.

Score a validation set with a saved model (no training):
```bash
python peptide_ml_workflow.py predict \
  --model results/saved_models/Imm_Non-Imm/logreg_pipeline.joblib \
  --metadata results/saved_models/Imm_Non-Imm/logreg_pipeline_metadata.json \
  --feature-table data/input/validation_set.csv \
  --output results/predictions/logreg_predictions.csv
```

Annotate those predictions against ground truth and make the figures:
```bash
python peptide_ml_workflow.py evaluate \
  --predictions-csv results/predictions/logreg_predictions.csv \
  --immunogenic data/input/validated_immunogenic.txt \
  --non-immunogenic data/input/validated_non_immunogenic.txt \
  --output-dir results/evaluation \
  --prefix logreg
```
(`--peptides some_subset.txt` is optional — omit it to evaluate on every
peptide in the predictions CSV, or pass it to first restrict to a subset of
interest.)

## Batch / HPC use (config-driven)

Fill in `example_config.yaml` with your training sets, feature lists, and
(optionally) per-model hyperparameter grids, then:

```bash
# whole workflow, every training set, every model, unattended:
python peptide_ml_workflow.py run --config example_config.yaml

# or just retrain one training set, only two models:
python peptide_ml_workflow.py train --config example_config.yaml \
  --training-set Imm_Non-Imm --models logreg hgb

# or just re-run validation later, for everything already saved:
python peptide_ml_workflow.py predict  --config example_config.yaml
python peptide_ml_workflow.py evaluate --config example_config.yaml
```

### Minimal SLURM example
```bash
#!/bin/bash
#SBATCH --job-name=peptide_ml
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00

module load python/3.x   # adjust for your HPC
source venv/bin/activate

python peptide_ml_workflow.py run --config example_config.yaml
```
GridSearchCV uses `n_jobs=-1` internally, so it will use all cores SLURM
gives the job; no extra parallelisation setup needed for a single-node run.

## Notes on the models

- Every saved pipeline's final step is named `"clf"` and always exposes
  `predict_proba`, so `predict` never needs to know which algorithm it's
  looking at.
- `svm_linear`/`svm_rbf` are grid-searched using plain `LinearSVC`/`SVC`
  (fast, ranking-based scoring only) and then the **winning** hyperparameters
  are used to rebuild a final version that supports `predict_proba`
  (`CalibratedClassifierCV` for linear, `probability=True` for RBF), refit on
  the full training set.
- Feature lists are **fixed**, supplied by you per training set (a text file,
  one column name per line) — no automatic correlation-clustering/feature
  selection is run. Do that analysis separately

## Files

- `peptide_ml_workflow.py` — the script.
- `example_config.yaml` — annotated template config
