python peptide_ml_workflow.py compare-features \
  --feature-table ../data/input/formatted-input/traintest_9kcleavedself_featureTable.csv \
  --immunogenic ../data/training_sets/assay-restricted/cedar_imm_peps_assay_restricted.txt \
  --non-immunogenic ../data/training_sets/assay-restricted/cedar_non_imm_peps_assay_restricted.txt \
  --champion-model hgb \
  --champion-features-file features/featureSet1.txt \
  --champion-model-path results/saved_models/assay-restricted/hgb_pipeline.joblib \
  --champion-metadata-path results/saved_models/assay-restricted/hgb_pipeline_metadata.json \
  --single-features avg_rank weak_binders_count mtec_expression_count pathogenicity di_best_score PRIME_%Rank_bestAllele \
  --output-dir results/champion_comparison_assay_restricted_SuppFig \