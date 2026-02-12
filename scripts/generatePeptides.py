# Feb 2026
# JVD
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--File', help='Path/To/ Input fasta file')
parser.add_argument('-o', '--Output', help='Output .txt file')
args = parser.parse_args()

def getMHCIPeptides(prot: str):
    peps = []
    pepLens = [8, 9, 10, 11]

    for pepLen in pepLens:
        # Change here if you decide to start peptides at beginning of ORF
        for startPos in range(pepLen-1):
            # Creates a list of chunks of size n from the protein string; discards any trailing chunk < n
            nmers = [''.join(chunk) for chunk in zip(*[iter(prot[startPos:])]*pepLen)]
            peps.extend(nmers)

    return peps

def getAllNmers(protFasta):
    allNmers = []

    with open(protFasta) as f:
        prot = ''
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if prot:
                    allNmers.extend(getMHCIPeptides(prot))
                    prot = ''
            else:
                prot += line         
    # Adding the last prots peptides in
    allNmers.extend(getMHCIPeptides(prot))

    return allNmers

if __name__ == '__main__':

    if not args.File:
        print('Supply input protein database fasta file. Usage: python generatePeptides.py -h')
        exit(1)

    protFasta = args.File

    allPossiblePeptides = getAllNmers(protFasta)
    # DEBUG
    #   print(f'{len(allPossiblePeptides)} peptides produced')

    # This removes all duplicates
    uniquePeptides = list(set(allPossiblePeptides))
    # DEBUG
    #   print(f'{len(uniquePeptides)} unique peptides produced')

    if args.Output:
        with open(args.Output, 'w') as f:
            for pep in uniquePeptides:
                f.write(f"{pep}\n")
    else:
        print(uniquePeptides)
    
    # TODO
    # Group by Nmers or by proteins

