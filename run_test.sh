#!/bin/bash
#PBS -N mtec3
#PBS -o mtec3.log
#PBS -lselect=1:ncpus=1:mem=128gb
#PBS -lwalltime=12:00:00
#PBS -j oe

cd $PBS_O_WORKDIR

module load Mamba/23.11.0-0
module load Nextflow/25.10.2

nextflow run test.nf -profile pbs