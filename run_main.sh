#!/bin/bash
#PBS -N allMtecs
#PBS -lselect=1:ncpus=4:mem=32gb
#PBS -lwalltime=24:00:00
#PBS -o allMtecs.log
#PBS -j oe

cd $PBS_O_WORKDIR

module load Mamba/23.11.0-0
module load Nextflow/25.10.2

nextflow run main.nf -profile pbs -resume