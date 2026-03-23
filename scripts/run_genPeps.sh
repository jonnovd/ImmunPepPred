#!/bin/bash
#PBS -lwalltime=12:00:00
#PBS -lselect=1:ncpus=8:mem=64gb
#PBS -N getMHCPeps

module load Python/3.13.5-GCCcore-14.3.0
cd ${PBS_O_WORKDIR}
# python scripts/generatePeptides.py -i tests/dbs/AS_RMATs_AN.NTC_E7_drugvsNTC_E7_vc.ExtractedAS.aa.3f.orfs.fa -f tests/dbs/human_reference_sequences.fasta -p tests/out/hpc-filter-out.txt -c tests/out/csv-hpc-filter.csv -l 9 10
python scripts/generatePeptides.py -i tests/dbs/SRR8615282_transcriptProtein.fasta -f tests/dbs/human_reference_sequences.fasta -p tests/out/udp-filter-out.txt -c tests/out/udp-csv-filter.csv -l 8 9 10 11
