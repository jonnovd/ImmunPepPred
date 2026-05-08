# Feb 2026
# JVD
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--Input', help='Path/To/ Input fasta file')
parser.add_argument('-f', '--ReferenceFile', help='Path/To/ Reference Fasta file to filter reads from')
parser.add_argument('-p', '--Peptides', help='Simple Peptides only Output .txt file')
parser.add_argument('-c', '--Csv', help='Output .csv file with transcript ids and peptides')
args = parser.parse_args()

def getMHCIPeptidesSet(prot: str, pepLens = [9]):
    peps = set()

    for pepLen in pepLens:
        for i in range(len(prot) - pepLen + 1):
            peps.add(prot[i: i + pepLen])

    return peps

def generateMHCIPeptides(prot: str, pepLens = [8, 9, 10, 11]):
    # Generator function (returns iterable?)
    for pepLen in pepLens:
        for i in range(len(prot) - pepLen + 1):
            yield prot[i: i + pepLen]

# KMP Algorithm
def constructLps_KMP(pat, lps):
    
    # len stores the length of longest prefix which 
    # is also a suffix for the previous index
    len_ = 0
    m = len(pat)
    
    # lps[0] is always 0
    lps[0] = 0

    i = 1
    while i < m:
        
        # If characters match, increment the size of lps
        if pat[i] == pat[len_]:
            len_ += 1
            lps[i] = len_
            i += 1
        
        # If there is a mismatch
        else:
            if len_ != 0:
                
                # Update len to the previous lps value 
                # to avoid redundant comparisons
                len_ = lps[len_ - 1]
            else:
                
                # If no matching prefix found, set lps[i] to 0
                lps[i] = 0
                i += 1

# KMP algorithm
def search_KMP(pat, txt):
    n = len(txt)
    m = len(pat)

    lps = [0] * m
    res = []

    constructLps_KMP(pat, lps)

    # Pointers i and j, for traversing 
    # the text and pattern
    i = 0
    j = 0

    while i < n:
        
        # If characters match, move both pointers forward
        if txt[i] == pat[j]:
            i += 1
            j += 1

            # If the entire pattern is matched 
            # store the start index in result
            if j == m:
                #res.append(i - j)
                return True
                # Use LPS of previous index to 
                # skip unnecessary comparisons
                #j = lps[j - 1]
        
        # If there is a mismatch
        else:
            
            # Use lps value of previous index
            # to avoid redundant comparisons
            if j != 0:
                j = lps[j - 1]
            else:
                i += 1
    return False


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
    # Input validation
    if not args.Input:
        print('Supply input protein database fasta file. Usage: python generatePeptides.py -h')
        exit(1)
    if not args.Peptides and not args.Csv:
        print('No Output file specified. Use -p or -c')
        exit(1)
    FILTER = False 
    if args.ReferenceFile:
        FILTER = True

    pep_file = open(args.Peptides, 'w') if args.Peptides else None
    csv_file = open(args.Csv, 'w') if args.Csv else None
    allNMers = []

    for header, orf in parseFasta(args.Input):
        allNMers.extend(getMHCIPeptidesSet(orf))

    allNMers = set(allNMers)

    #Filtering out any peptides that are in the human proteome reference
    # TODO
    # This could be a place to add GTex data filtering
    # Only remove peptides produced from ORFs that are expressed in the healthy tissue in this site
    if FILTER:
        for header, prot in parseFasta(args.ReferenceFile):
            nmers_to_remove = set()
            for nmer in allNMers:
                if search_KMP(nmer, prot):
                    nmers_to_remove.add(nmer)
            allNMers -= nmers_to_remove

    for pep in allNMers:
        if pep_file:
            pep_file.write(f'{pep}\n')

    if pep_file:
        pep_file.close()
    if csv_file:
        csv_file.close()


if __name__ == '__main__':
    main()