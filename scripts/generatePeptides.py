# Feb 2026
# JVD
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--Input', help='Path/To/ Input fasta file')
parser.add_argument('-p', '--Peptides', help='Simple Peptides only Output .txt file')
parser.add_argument('-c', '--Csv', help='Output .csv file with transcript ids and peptides')
args = parser.parse_args()

def getMHCIPeptides(prot: str):
    peps = []
    pepLens = [8, 9, 10, 11]

    for pepLen in pepLens:
        # Change here if you decide to start peptides at beginning of ORF
        # TODO
        # Do we need to do a sliding window here if we've already included transcripts for all 3 frames
        # Yes, because this has got to do with the peptide cleavage not necessarily starting at the start of the peptide, not transcription
        for startPos in range(pepLen):
            # Creates a list of chunks of size n from the protein string; discards any trailing chunk < n
            nmers = [''.join(chunk) for chunk in zip(*[iter(prot[startPos:])]*pepLen)]
            peps.extend(nmers)

    uniquePeptides = peps#list(set(peps))
    return uniquePeptides

def getAllNmers(protFasta):
    pIDs = []
    allNmers = []

    with open(protFasta) as f:
        pID = ''
        prot = ''
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if prot:
                    pIDs.append(pID)
                    # Get all possible Nmer peptides from the protein and add these to the list
                    # Data structure is a list of lists ie [[peps1], [peps2], [peps3]]
                    allNmers.append(getMHCIPeptides(prot))
                    prot = ''
                
                pID = line
            else:
                prot += line         
    # Adding the last prots peptides in
    pIDs.append(pID)
    allNmers.append(getMHCIPeptides(prot))

    return pIDs, allNmers

if __name__ == '__main__':

    if not args.Input:
        print('Supply input protein database fasta file. Usage: python generatePeptides.py -h')
        exit(1)

    protFasta = args.Input

    IDs, peptideLists = getAllNmers(protFasta)
    
    # DEBUG
    #   print(f'{len(allPossiblePeptides)} peptides produced')

    # This removes all duplicates
    # TODO - Check if there are any duplicate peptides from multiple transcripts
    # uniquePeptides = list(set(allPossiblePeptides))
    # DEBUG
    #   print(f'{len(uniquePeptides)} unique peptides produced')

    if args.Peptides:
        with open(args.Peptides, 'w') as f:
            for i in range(len(peptideLists)):
                peps = peptideLists[i]
                for pep in peps:
                    f.write(f"{pep}\n")

    if args.Csv:
        with open(args.Csv, 'w') as f:
            for i in range(len(IDs)):
                peps = peptideLists[i]
                for pep in peps:
                    f.write(f"{IDs[i]}, {pep}\n")

