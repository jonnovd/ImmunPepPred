# Feb 2026
# JVD
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--Input', help='Path/To/ Input fasta file')
parser.add_argument('-f', '--ReferenceFile', help='Path/To/ Reference Fasta file to filter reads from')
parser.add_argument('-p', '--Peptides', help='Simple Peptides only Output .txt file')
parser.add_argument('-c', '--Csv', help='Output .csv file with transcript ids and peptides')
parser.add_argument('-l', '--peptideLengths', nargs='+', type=int, default=[9],
                    help='Peptide lengths to generate (default: 9)')
args = parser.parse_args()

def getMHCIPeptidesSet(prot: str, pepLens = [8, 9, 10, 11]):
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

def buildReferenceSet(referenceFasta, pepLens=[8, 9, 10, 11]):
    ref_nmers = set()
    for _, prot in parseFasta(referenceFasta):
        ref_nmers.update(getMHCIPeptidesSet(prot, pepLens))
    return ref_nmers

def extractTranscriptId(header: str) -> str:
    
    return header.split()[0].lstrip('>')

def buildPeptideTranscriptMap(inputFasta, pepLens=[8, 9, 10, 11]):
    peptide_transcripts = defaultdict(set)
    for header, orf in parseFasta(inputFasta):
        transcript_id = extractTranscriptId(header)
        for pep in getMHCIPeptidesSet(orf, pepLens):
            peptide_transcripts[pep].add(transcript_id)
    return peptide_transcripts

def main():

    if not args.Input:
        print('Supply input protein database fasta file. Usage: python generatePeptides.py -h')
        exit(1)
    if not args.Peptides and not args.Csv:
        print('No Output file specified. Use -p or -c')
        exit(1)
    referenceFasta = None
    referencePepSet = None
    if args.ReferenceFile:
        referenceFasta = args.ReferenceFile

    pepLens = args.peptideLengths
    allNMers = set()
    for _, orf in parseFasta(args.Input):
        allNMers.update(getMHCIPeptidesSet(orf, pepLens))

    # TODO DEBUG
    test_file = open('tests/kmp-hpc/NoFilter.txt', 'w')
    for pep in allNMers:
        test_file.write(f"{pep}\n")

    if referenceFasta:
        referencePepSet = buildReferenceSet(referenceFasta, pepLens)
        allNMers -= referencePepSet

    pep_file = open(args.Peptides, 'w') if args.Peptides else None
    csv_file = open(args.Csv, 'w') if args.Csv else None

    for pep in allNMers:
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