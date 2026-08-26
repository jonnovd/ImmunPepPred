python peptide_ml_workflow.py compare-features \
  --feature-table ../data/input/formatted-input/traintest_9kcleavedself_featureTable.csv \
  --immunogenic ../data/training_sets/cedar_imm_peps_traintest.txt \
  --non-immunogenic ../data/training_sets/negSetC_cedar-ni_selfBinders.txt \
  --champion-model hgb \
  --champion-features-file features/featureSet1.txt \
  --champion-model-path results/saved_models/augmented-neg/hgb_pipeline.joblib \
  --champion-metadata-path results/saved_models/augmented-neg/hgb_pipeline_metadata.json \
  --single-features avg_rank weak_binders_count mtec_expression_count pathogenicity di_best_score PRIME_%Rank_bestAllele \
  --output-dir results/champion_comparison_augmented-neg \
  # --palette "#A8CFE9" "#E4A6E4" "#A0DDAF" \
  # --font-family "Arial" \
  # --axis-label-fontsize 12 \
  # --tick-label-fontsize 10