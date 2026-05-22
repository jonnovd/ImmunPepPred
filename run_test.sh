#!/bin/bash
#PBS -N filelogic
#PBS -o filelogic.log
#PBS -lselect=1:ncpus=4:mem=128gb
#PBS -lwalltime=1:00:00
#PBS -j oe

cd $PBS_O_WORKDIR

module load Mamba/23.11.0-0
module load Nextflow/25.10.2

nextflow run test.nf -profile pbs