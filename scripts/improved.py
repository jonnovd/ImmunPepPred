# Feb 2026
# JVD
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--Input', help='Path/To/ Input fasta file')
parser.add_argument('-p', '--Peptides', help='Simple Peptides only Output .txt file')
parser.add_argument('-c', '--Csv', help='Output .csv file with transcript ids and peptides')
args = parser.parse_args()

def getMHCIPeptides(prot: str, pepLens = [8, 9, 10, 11]):
    peps = []

    for pepLen in pepLens:
        for i in range(len(prot) - pepLen + 1):
            peps.append(prot[i: i + pepLen])

    uniquePeptides = list(set(peps))
    return uniquePeptides

def generateMHCIPeptides(prot: str, pepLens = [8, 9, 10, 11]):
    # Generator function (returns iterable?)
    for pepLen in pepLens:
        for i in range(len(prot) - pepLen + 1):
            yield prot[i: i + pepLen]

def parseFasta(protFasta):
    header = None
    seq_tokens = []

    with open(protFasta) as f:

        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header:
                    yield header, "".join(seq_tokens)
                    seq_tokens = []
                
                header = line
            else:
                seq_tokens.append(line)        
    # Adding the last prots peptides in
    yield header, "".join(seq_tokens)

def main():

    if not args.Input:
        print('Supply input protein database fasta file. Usage: python generatePeptides.py -h')
        exit(1)

    if not args.Peptides and not args.Csv:
        print('No Output file specified. Use -p or -c')
        exit(1)

    pep_file = open(args.Peptides, 'w') if args.Peptides else None
    csv_file = open(args.Csv, 'w') if args.Csv else None

    for header, orf in parseFasta(args.Input):
        for pep in generateMHCIPeptides(orf):
            if pep_file:
                pep_file.write(f"{pep}\n")
            if csv_file:
                csv_file.write(f'{header}, {pep}')

    if pep_file:
        pep_file.close()
    if csv_file:
        csv_file.close()


if __name__ == '__main__':
    main()