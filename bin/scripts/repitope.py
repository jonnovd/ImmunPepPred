# Adapted from author's original wrapper script to accept peptides as parameters

# Adapted from author's original wrapper script to accept peptides as parameters
import argparse
import os
import shutil
import glob
import pyper

# Parse command line arguments
parser = argparse.ArgumentParser(description="Run Repitope MHC-I epitope prioritization")
parser.add_argument("--input", required=True, help="Path to input peptide file (.txt)")
parser.add_argument("--home", required=True, help="Path to Repitope data directory")
parser.add_argument("--output", required=True, help="Path to output file (e.g. results.csv)")
parser.add_argument("--cpu", default="4", help="Number of CPUs to use (default: 4)")
parser.add_argument("--memory", default="60G", help="Java memory allocation (default: 60G)")
parser.add_argument("--pept_len_range", default="8:11", help="Peptide length range (default: 8:11)")
args = parser.parse_args()

# Set variables from arguments
cpu = args.cpu
memory = args.memory
pept_len_range = args.pept_len_range
home = args.home
target_peptide_file = args.input
output_file = args.output
tcr_frag_file = home + "/FragmentLibrary.fst"
mhci_feature_file = home + "/FeatureDF_MHCI_Weighted.10000.fst"

# Helper function to run R commands and surface any errors
def run_r(r, command, step_name):
    print(f"\n--- Running R step: {step_name} ---")
    result = r(command)
    if result:
        print(result)
    # Check if R has recorded any errors or warnings
    err = r("if(exists('last.warning')) print(last.warning)")
    if err:
        print(f"R warnings: {err}")

# Create a temporary working directory for Repitope's output
tmp_outdir = output_file + "_tmp"
os.makedirs(tmp_outdir, exist_ok=True)

r = pyper.R()

run_r(r, "options(java.parameters='-Xmx" + memory + "')\n \
library(tidyverse)\n \
library(data.table)\n \
library(Repitope)", "Load libraries")

run_r(r, "peptideSet_test <- Repitope::sequenceFilter(data.table::fread('" + target_peptide_file + "')$Peptide)", "Sequence filter")

run_r(r, "peptideSet_test <- peptideSet_test[nchar(peptideSet_test) %in% " + pept_len_range + "]", "Peptide length filter")

# Print how many peptides passed filtering
print(r("cat('Peptides after filtering:', length(peptideSet_test), '\n')"))

run_r(r, "fragLibDT <- fst::read_fst('" + tcr_frag_file + "', as.data.table=T)", "Load fragment library")

run_r(r, "featureDT_MHCI <- fst::read_fst('" + mhci_feature_file + "', as.data.table=T)", "Load feature matrix")

run_r(r, "res_MHCI <- EpitopePrioritization( \
  featureDF=featureDT_MHCI[Peptide%in%MHCI_Human$Peptide,], \
  metadataDF=MHCI_Human[,.(Peptide,Immunogenicity)],\
  peptideSet=peptideSet_test,\
  fragLib=fragLibDT,\
  aaIndexIDSet='all',\
  fragLenSet=3:8,\
  fragDepth=10000,\
  fragLibType='Weighted',\
  featureSet=MHCI_Human_MinimumFeatureSet,\
  seedSet=1:5,\
  coreN=" + cpu + ",\
  outDir='" + tmp_outdir + "')", "EpitopePrioritization")

# List what was actually written to the temp directory
print("\n--- Contents of output directory ---")
print(os.listdir(tmp_outdir))

# Find the output file written by Repitope and copy it to the desired output path
output_files = glob.glob(tmp_outdir + "/*.csv")
if not output_files:
    raise FileNotFoundError("Repitope did not produce any output files in " + tmp_outdir)
shutil.copy(output_files[0], output_file)

# Clean up temporary directory
shutil.rmtree(tmp_outdir)

# Fully working
# import argparse
# import os
# import pyper

# # Parse command line arguments
# parser = argparse.ArgumentParser(description="Run Repitope MHC-I epitope prioritization")
# parser.add_argument("--input", required=True, help="Path to input peptide file (.txt)")
# parser.add_argument("--home", required=True, help="Path to Repitope data directory")
# parser.add_argument("--outdir", required=True, help="Path to output directory")
# parser.add_argument("--cpu", default="4", help="Number of CPUs to use (default: 4)")
# parser.add_argument("--memory", default="60G", help="Java memory allocation (default: 60G)")
# parser.add_argument("--pept_len_range", default="8:11", help="Peptide length range (default: 8:11)")
# args = parser.parse_args()

# # Set variables from arguments
# cpu = args.cpu
# memory = args.memory
# pept_len_range = args.pept_len_range
# home = args.home
# target_peptide_file = args.input
# outdir = args.outdir
# tcr_frag_file = home + "/FragmentLibrary.fst"
# mhci_feature_file = home + "/FeatureDF_MHCI_Weighted.10000.fst"

# # Create output directory if it doesn't already exist
# os.makedirs(outdir, exist_ok=True)

# r = pyper.R()

# command = "options(java.parameters='-Xmx" + memory + "')\n \
# library(tidyverse)\n \
# library(data.table)\n \
# library(Repitope)"
# r(command)

# command = "peptideSet_test <- Repitope::sequenceFilter(data.table::fread('" + target_peptide_file + "')$Peptide)"
# r(command)

# command = "peptideSet_test <- peptideSet_test[nchar(peptideSet_test) %in% " + pept_len_range + "]"
# r(command)

# command = "fragLibDT <- fst::read_fst('" + tcr_frag_file + "', as.data.table=T)"
# r(command)

# command = "featureDT_MHCI <- fst::read_fst('" + mhci_feature_file + "', as.data.table=T)"
# r(command)

# command = "res_MHCI <- EpitopePrioritization( \
#   featureDF=featureDT_MHCI[Peptide%in%MHCI_Human$Peptide,], \
#   metadataDF=MHCI_Human[,.(Peptide,Immunogenicity)],\
#   peptideSet=peptideSet_test,\
#   fragLib=fragLibDT,\
#   aaIndexIDSet='all',\
#   fragLenSet=3:8,\
#   fragDepth=10000,\
#   fragLibType='Weighted',\
#   featureSet=MHCI_Human_MinimumFeatureSet,\
#   seedSet=1:5,\
#   coreN=" + cpu + ",\
#   outDir='" + outdir + "')"
# r(command)


# cpu = str(4)
# memory = "60G"
# pept_len_range = "8:11" ## For MHC-I prediction

# home = "C:/Repitope"
# target_peptide_file = home + "/peptide_test.txt"
# tcr_frag_file = home + "/FragmentLibrary.fst"
# mhci_feature_file = home + "/FeatureDF_MHCI_Weighted.10000.fst"

# import pyper
# r = pyper.R()

# command = "options(java.parameters='-Xmx" + memory + "')\n \
# library(tidyverse)\n \
# library(data.table)\n \
# library(Repitope)"
# r(command)

# command = "peptideSet_test <- Repitope::sequenceFilter(data.table::fread('" + target_peptide_file + "')$Peptide)"
# r(command)

# command = "peptideSet_test <- peptideSet_test[nchar(peptideSet_test) %in% " + pept_len_range + "]"
# r(command)

# command = "fragLibDT <- fst::read_fst('" + tcr_frag_file + "', as.data.table=T)"
# r(command)

# command = "featureDT_MHCI <- fst::read_fst('" + mhci_feature_file + "', as.data.table=T)"
# r(command)

# import os
# pid = str(os.getpid())
# os.mkdir(home + "/ProcessID_" + pid)

# command = "res_MHCI <- EpitopePrioritization( \
#   featureDF=featureDT_MHCI[Peptide%in%MHCI_Human$Peptide,], \
#   metadataDF=MHCI_Human[,.(Peptide,Immunogenicity)],\
#   peptideSet=peptideSet_test,\
#   fragLib=fragLibDT,\
#   aaIndexIDSet='all',\
#   fragLenSet=3:8,\
#   fragDepth=10000,\
#   fragLibType='Weighted',\
#   featureSet=MHCI_Human_MinimumFeatureSet,\
#   seedSet=1:5,\
#   coreN=" + cpu + ",\
#   outDir='" + home + "/MHCI_ProcessID_" + pid + "')"
# r(command)
