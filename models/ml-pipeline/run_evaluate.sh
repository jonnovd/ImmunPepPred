python peptide_ml_workflow.py evaluate \
  --predictions-csv results/predictions/augmented-neg/hgb_predictions.csv \
  --immunogenic ../data/validation_sets/iedb_imm_peps_noTrainMatch.txt \
  --non-immunogenic ../data/validation_sets/iedb_non-imm_peps_noTrainMatch.txt \
  --output-dir results/evaluation \
  --prefix augmented-neg/hgb \
  --prime-feature-table ../data/input/formatted-input/iedb_featureTable.csv \
  --prime-score-col PRIME_%Rank_bestAllele \
  --precision-at-k-fractions 5 10 15 20 30 40 50