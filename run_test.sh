#!/bin/bash
#PBS -N tapParallel
#PBS -o tapParallel.log
#PBS -lselect=1:ncpus=8:mem=128gb
#PBS -lwalltime=4:00:00
#PBS -j oe

cd $PBS_O_WORKDIR

module load Mamba/23.11.0-0
module load Nextflow/25.10.2

nextflow run tapWrapperTest.nf -profile pbs -resume