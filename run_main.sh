#!/bin/bash
#PBS -N ipp
#PBS -lselect=1:ncpus=8:mem=64gb
#PBS -lwalltime=24:00:00
#PBS -o ipp.log
#PBS -j oe

cd $PBS_O_WORKDIR

module load Mamba/23.11.0-0
module load Nextflow/25.10.2

nextflow run main.nf -profile pbs -resume