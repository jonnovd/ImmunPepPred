#!/bin/bash
#SBATCH --job-name=IPP --partition=master-worker --ntasks=1 --mem=4000 --output=ipp_Log.txt --time=12:00:00

module load Nextflow

nextflow run main.nf -profile slurm -resume