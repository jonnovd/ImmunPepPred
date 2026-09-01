# ImmunPepPred

- The Nextflow pipeline is run from `main.nf` and configured in `nextflow.config`
- Confidential Nextflow workflows and processes have been hidden from the public git repo

## User-defined Input

### Pipeline Settings & Configuration

### ORF DB
- Fasta file containing ORFs of interest (eg. cancer-restricted ORFs)

### Peptides.txt
- Text file comprising 1 peptide (8-11mer) per line
- Peptides of interest could be from:
    - MS-identified peptides to prioritise for immunogenicity testing
    - Proteogenomic predicted peptides to prioritise for targeted MS
    - Candidate peptides for MS search DB to prioritise

### Both
- Peptides generated from ORF DB are combined with custom input peptides
- All unique peptides are run through the pipeline

### HLA class I alleles
- Context-specific alleles in the form `HLA-A02:01`
- Maximum 

### Input Constraints
#### Max allele input 
- 90 alleles including the common allele panel (max 56 unique context-specific alleles) 
    - Due to NetMHCpan's allele string input size of 1024 which ~90 alleles
    - TODO: Use peptide- and allele-level batching to solve this
- Some less common alleles are unsupported by different tools.
    - One solution is to disable the tool in the config file
    - WARNING: This will affect the HLA feature predictions and would require model retraining

## Retraining ML Model

## Running the Framework


## genPeps.py
