import csv
import argparse

def read_file(file_path):
    """Reads a file and returns a list of lines."""
    with open(file_path, 'r') as file:
        return [line.strip() for line in file if line.strip()]

def generate_pairs(peptides, hla_types):
    """Generates all pairs of peptides and HLA types."""
    return [{"peptide": peptide, "allele": hla} for peptide in peptides for hla in hla_types]

def write_csv(output_file, pairs):
    """Writes the pairs to a CSV file."""
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=["peptide", "allele"])
        writer.writeheader()
        writer.writerows(pairs)

def main():
    parser = argparse.ArgumentParser(description="Generate CSV of peptide-HLA pairs.")
    parser.add_argument("peptides_file", help="File containing list of peptides.")
    parser.add_argument("hla_file", help="File containing list of HLA types.")
    parser.add_argument("output_file", help="Output CSV file.")
    args = parser.parse_args()

    peptides = read_file(args.peptides_file)
    hla_types = read_file(args.hla_file)
    pairs = generate_pairs(peptides, hla_types)
    write_csv(args.output_file, pairs)

if __name__ == "__main__":
    main()