#!/bin/bash
#PBS -N sreejan
#PBS -lselect=1:ncpus=4:mem=32gb
#PBS -lwalltime=12:00:00
#PBS -o sreejan.log
#PBS -j oe

cd $PBS_O_WORKDIR

module load Mamba/23.11.0-0
module load Nextflow/25.10.2

nextflow run main.nf -profile pbs -resume