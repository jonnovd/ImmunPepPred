# ImmunPepPred

- The Nextflow pipeline is run from `main.nf` and configured in `nextflow.config`
- Confidential Nextflow workflows and processes have been hidden from the public git repo

## Running the Framework
- For PBS:   `qsub run_main.sh`
- For local: `bash rlocal.sh`
- For SLURM: `TODO`

## User-defined Input

### Pipeline Settings & Configuration
- Consult the Nextflow.config file
- User-defined parameters
    - Set input and output file paths
    - `batchSize`
    - `hla_batch_size`
    - `peptide_lengths`
    - `min_peptide_length`
    - `max_peptide_length`
    - Flags to indicate which processes to run. Can mostly be left at default values of true / false
    - ORF mode filtering (See `getPeptidesFromORFs.py` description at the end)
        - `pepsickleCleavageThreshold` to `null` to skip proteasomal cleavage filtering; or instead to the desired threshold
        - `filter`
        - DB files to set for advanced filtering:
            - `refGFF` aids functionality below
            - `GenesToRemoveFromRef` List of gene names to remove from the reference proteome, ie. self ORFs from which peptides should be retained 

### One of the following
#### ORF FASTA FILE
- Fasta file containing ORFs of interest (eg. cancer-restricted ORFs)

#### Peptides.txt
- Text file comprising 1 peptide (8-11mer) per line
- Peptides of interest could be from:
    - MS-identified peptides to prioritise for immunogenicity testing
    - Proteogenomic predicted peptides to prioritise for targeted MS
    - Candidate peptides for MS search DB to prioritise

#### Both
- Peptides generated from ORF DB are combined with custom input peptides
- All unique peptides are run through the pipeline

### HLA class I alleles
- Context-specific alleles in the form `HLA-A02:01`
- Maximum 56 unique alleles that don't overlap with the common allele panel (34 from [REF])

### Input Constraints
- Max input of 90 alleles including the common allele panel (max 56 unique context-specific alleles) 
    - Due to NetMHCpan's allele string input size of 1024 which ~90 alleles
    - TODO: Use peptide- and allele-level batching to solve this
- Some less common alleles are unsupported by different tools.
    - One solution is to disable the tool in the config file
    - WARNING: This will affect the HLA feature predictions and would require model retraining

## ML Model
- Initial models were tested and evaluated in `model_exploration.ipynb` and `model_evaluator.ipynb`
- DL model was implemented in `DL_model.ipynb`
- Ultimately developing and evaluating all models on different training configurations was restructured into an HPC-friendly ML pipeline in `models/ml-pipeline/peptide_ml_workflow.py`
    - See `models/ml-pipeline/README.md` for usage info

## Thesis Figures
- `models/feature_visualiser.ipynb`
- `models/robust_cedar_exploration.ipynb`
- `models/ml-pipeline/run_architecture_figure.sh`
- `models/ml-pipeline/run_champion_comparison_figure.sh`
- `models/ml-pipeline/run_champion_comparison_sfig6.sh`
- `models/ml-pipeline/run_evaluate.sh` -> `results/...`

## getPeptidesFromORFs.py
- Script to generate in silico peptides from a FASTA file
- Several user-defined filtering options for:
    - Peptide length
    - Cleavage probability
    - Removal of self-peptides
        - Retaining some self ORFs such as CTAs
