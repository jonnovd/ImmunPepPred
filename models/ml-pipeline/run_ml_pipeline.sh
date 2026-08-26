#!/bin/bash
#PBS -N mlPipeline
#PBS -lselect=1:ncpus=8:mem=32gb
#PBS -lwalltime=8:00:00
#PBS -o mlPipeline.log
#PBS -j oe

cd $PBS_O_WORKDIR

source .venv/bin/activate
python peptide_ml_workflow.py run --config finalModelConfig.yaml
