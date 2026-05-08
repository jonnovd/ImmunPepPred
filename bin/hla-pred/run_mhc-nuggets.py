# Fix for use of lr instead of learning_rate in newer keras:
from tensorflow.keras.optimizers.legacy import Adam
import tensorflow.keras.optimizers as opt_module
opt_module.Adam = Adam

import argparse
from mhcnuggets.src.predict import predict

# Set up argument parser
parser = argparse.ArgumentParser(description="Run MHC-Nuggets predictions.")
parser.add_argument('--peptides_file', required=True, help="Path to the file containing peptides.")
parser.add_argument('--hla_file', required=True, help="Path to the file containing HLA types.")
parser.add_argument('--output', default="mhcnuggets_results.txt",  help="Path to write the results.")

# Parse arguments
args = parser.parse_args()
    
def read_prediction_csv_into_dict(file_path, hla, data):
    with open(file_path, 'r') as csv_file:
        lines = csv_file.readlines()
        for line in lines[1:]:
            cols = line.strip().split(',')
            if cols[0] not in data:
                data[cols[0]] = {}

            if hla not in data[cols[0]]:
                data[cols[0]][hla] = {}

            data[cols[0]][hla] = cols[2]

# Read HLA types from file
hla_types = [h.strip() for h in open(args.hla_file)]

tmp_file = "tmp_predictions.txt"

prediction_results = {}

# Iterate over HLA types and predict for each
for hlatype in hla_types:
    hla = hlatype.replace("*", "")
    print(f"Predicting for HLA type: {hla}")
    predict(class_='I',
            peptides_path=args.peptides_file,
            mhc=hla,
            rank_output=True,
            output=tmp_file
            )
    tmptmp_file = "tmp_predictions_ranks..txt"
    read_prediction_csv_into_dict(tmptmp_file, hla, prediction_results)

# Write the results to a file
with open(args.output, 'w') as output_file:
    output_file.write("Peptide\tHLA\t%Rank\n")
    for peptide, hla_scores in prediction_results.items():
        for hla, score in hla_scores.items():
            output_file.write(f"{peptide}\t{hla}\t{score}\n")

