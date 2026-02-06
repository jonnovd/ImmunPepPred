# Feb 2026
# JVD

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
    protFasta = 'tests/testUnique.fa'
    allPossiblePeptides = getAllNmers(protFasta)
    print(f'{len(allPossiblePeptides)} peptides produced')
    # This removes all duplicates
    uniquePeptides = list(set(allPossiblePeptides))
    print(f'{len(uniquePeptides)} unique peptides produced')
    
    # TODO
    # Write to text file or keep as list
    # Group by Nmers or by proteins
    # Use arguments in this file

