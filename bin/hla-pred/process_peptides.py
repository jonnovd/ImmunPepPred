
import argparse

def read_peptides(input_file):
    """Reads peptide sequences from a file."""
    with open(input_file, 'r') as file:
        peptides = [line.strip() for line in file]
    return peptides

def filter_peptides(peptides, min_length, max_length):
    """Filters peptides by length and removes those with non-standard amino acids."""
    standard_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
    filtered = [
        peptide for peptide in peptides
        if min_length <= len(peptide) <= max_length and set(peptide).issubset(standard_amino_acids)
    ]
    return filtered

def split_into_batches(peptides, batch_size):
    """Splits peptides into batches of a given size."""
    for i in range(0, len(peptides), batch_size):
        yield peptides[i:i + batch_size]

def write_batches(batches):
    """Writes peptide batches to separate files."""
    for i, batch in enumerate(batches):
        batch_file = f"peptides_batch_{i + 1}.txt"
        with open(batch_file, 'w') as file:
            file.write("\n".join(batch))

def main():
    parser = argparse.ArgumentParser(description="Process peptide sequences.")
    parser.add_argument("input_file", type=str, help="Path to the input file containing peptide sequences.")
    parser.add_argument("min_length", type=int, help="Minimum length of peptides to include.")
    parser.add_argument("max_length", type=int, help="Maximum length of peptides to include.")
    parser.add_argument("batch_size", type=int, help="Number of peptides per batch.")

    args = parser.parse_args()

    peptides = read_peptides(args.input_file)
    filtered_peptides = filter_peptides(peptides, args.min_length, args.max_length)
    batches = list(split_into_batches(filtered_peptides, args.batch_size))
    write_batches(batches)

    print(f"Processed {len(peptides)} to {len(filtered_peptides)} peptides and split into {len(batches)} batches.")

if __name__ == "__main__":
    main()