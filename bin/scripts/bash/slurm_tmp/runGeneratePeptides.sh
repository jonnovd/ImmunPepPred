#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --job-name=getMHCPeps

module load Python/3.13.5-GCCcore-14.3.0

INPUT_DIR=$1
OUTPUT_DIR=$2

mkdir -p ${OUTPUT_DIR}

for fasta in ${INPUT_DIR}/*.fa; do
    # Skip if no files match the glob
    [ -f "$fasta" ] || continue

    basename=$(basename ${fasta})
    sample_name=${basename%.*}

    python generatePeptides.py \
        -i ${fasta} \
        -f /data/rds/DBI/DUDBI/FUNCPROT/James/Databases/human_reference_sequences.fasta \
        -p ${OUTPUT_DIR}/${sample_name}_nc_peptides_8-11mers.txt \
        -l 8 9 10 11
done