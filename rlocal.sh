#!/usr/bin/env bash

set -euo pipefail

# ── Find all .nf files in the current directory ──────────────────────────────
mapfile -t NF_FILES < <(find . -maxdepth 1 -name "*.nf" | sort)

if [[ ${#NF_FILES[@]} -eq 0 ]]; then
  echo "No .nf files found in the current directory."
  exit 1
fi

# ── Choose a pipeline file ────────────────────────────────────────────────────
echo "Available Nextflow pipelines:"
for i in "${!NF_FILES[@]}"; do
  printf "  [%d] %s\n" "$((i + 1))" "${NF_FILES[$i]}"
done

while true; do
  read -rp "Select pipeline [1-${#NF_FILES[@]}]: " choice
  if [[ "$choice" =~ ^[0-9]+$ ]] && (( choice >= 1 && choice <= ${#NF_FILES[@]} )); then
    PIPELINE="${NF_FILES[$((choice - 1))]}"
    break
  fi
  echo "  Invalid selection, please try again."
done

# ── Choose run mode ───────────────────────────────────────────────────────────
echo ""
echo "Run mode:"
echo "  [1] Resume   (use cached results, default)"
echo "  [2] Restart  (fresh run, no -resume)"

while true; do
  read -rp "Select mode [1-2]: " mode
  case "$mode" in
    1) RESUME_FLAG="-resume"; break ;;
    2) RESUME_FLAG="";        break ;;
    *) echo "  Invalid selection, please try again." ;;
  esac
done

# ── Install Modules ───────────────────────────────────────────────────────────
echo ""
if command -v nextflow &>/dev/null; then
  echo "✔ Nextflow already loaded ($(command -v nextflow))"
else
  echo "Loading Nextflow module..."
  module load Nextflow/25.10.2
fi
 
if command -v mamba &>/dev/null; then
  echo "✔ Mamba already loaded ($(command -v mamba))"
else
  echo "Loading Mamba module..."
  module load Mamba/23.11.0-0
fi

# ── Build and run the command ─────────────────────────────────────────────────
CMD="nextflow run $PIPELINE $RESUME_FLAG -params-file local/params.json"

echo ""
echo "Running: $CMD"
echo "──────────────────────────────────────────"
eval "$CMD"