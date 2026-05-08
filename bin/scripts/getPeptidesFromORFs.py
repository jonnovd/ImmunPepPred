# Feb 2026
# JVD
import argparse
from collections import defaultdict
import re

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='Path/To/ Input fasta file')
parser.add_argument('-r', '--reference-file', help='Optional Path/To/ Reference Fasta file to filter reads from')
parser.add_argument('-c', '--cleavage-prediction', help='Optional Path/To/ Pepsickle cleavage prediction output.txt')
parser.add_argument('-t', '--cleavage-threshold', help='Positivie cleavage probability threshold', type=float, default=0.4)
parser.add_argument('-x', '--reference-cleavage-prediction', help='Optional Path/To/ Pepsickle cleavage prediction output.txt for the reference proteome')
parser.add_argument('-o', '--output-txt', help='Simple Peptides only Output .txt file')
parser.add_argument('-O', '--output-csv', help='Output .csv file with transcript ids and peptides')
parser.add_argument('-l', '--peptide-lengths', nargs='+', type=int, default=[9],
                    help='Peptide lengths to generate (default: 9)')
# parser.add_argument('-g', '--gff', help='Path/To/ Reference GFF file', default=None)
# parser.add_argument('-e', '--exclude-genes', help='Path/To/ txt file containing 1 gene symbol per line for genes to remove from the reference set')
args = parser.parse_args()

def extractTranscriptId(header: str, type: str) -> str:
    if type == 'pepsickle':
        return header.split()[0].lstrip('>')
    elif type == 'inputFasta':
        return header.split()[0].lstrip('>')[:header.index('|')-1]
    elif type == 'referenceProteome':
        return header.split('|')[1].split('.')[0]
    else:
        return None
    
def extractGeneName(header):
    """Extract gene name from Gencode or SwissProt header."""
    if header.startswith('>ENSP'):
        return header.split('|')[6]
    elif header.startswith('>sp'):
        return header.split('=')[3].split()[0]
    else:
        print('Header not from Gencode or SwissProt')
        exit(1)

def getMHCIPeptidesSet(prot: str, pepLens = [8, 9, 10, 11]):
    peps = set()

    for pepLen in pepLens:
        for i in range(len(prot) - pepLen + 1):
            peps.add(prot[i: i + pepLen])

    return peps

def getCleavedMHCIPeptidesSet(prot: str, cleavageIndices: list, pepLens = [8, 9, 10, 11]):
    """Returns unique 8-11mers that end at predicted cleavage sites"""
    peps = set()

    for i in cleavageIndices:
        for pepLen in pepLens:
            if i - pepLen > -1:
                peps.add(prot[i - pepLen: i])

    return peps

def generateMHCIPeptides(prot: str, pepLens = [8, 9, 10, 11]):
    # Generator function
    for pepLen in pepLens:
        for i in range(len(prot) - pepLen + 1):
            yield prot[i: i + pepLen]

def getPeptides(prot, header, cleaveDict, pepLens):
    """Return peptides for a single protein, using cleavage if available."""
    if cleaveDict:
        transcript_id = extractTranscriptId(header, 'pepsickle')
        if transcript_id in cleaveDict:
            return getCleavedMHCIPeptidesSet(prot, cleaveDict[transcript_id], pepLens)
    return getMHCIPeptidesSet(prot, pepLens)

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

def getCleavageIndicesDict(cleaveFile: str, threshold) -> dict:
    indices = {}
    with open(cleaveFile) as f:
        for line in f:
            if line[0] != 'p':
                tokens = line.split()
                if float(tokens[2]) > threshold:
                    transcript = tokens[4]
                    pos = int(tokens[0]) - 1
                    if transcript not in indices:
                        indices[transcript] = []
                    indices[transcript].append(pos)
    return indices

def buildReferenceSet(referenceFasta, pepLens=[8, 9, 10, 11], cleaveFile=None, cleavageThreshold=1.0, gff=None, genesToExcludeFromReference=None, logFile=None):
    refNmers = set()
    foundGenes = set()
    cleaveDict = getCleavageIndicesDict(cleaveFile, cleavageThreshold) if cleaveFile else None

    excludeGenes = set()
    if genesToExcludeFromReference:
        with open(genesToExcludeFromReference) as txt:
            excludeGenes = {line.strip() for line in txt}

    if excludeGenes:
        for header, prot in parseFasta(referenceFasta):
            gene = extractGeneName(header)
            if gene in excludeGenes:
                foundGenes.add(gene)
            else:
                refNmers.update(getPeptides(prot, header, cleaveDict, pepLens))

        if logFile:
            with open(logFile, 'w') as log:
                if foundGenes:
                    for gene in sorted(foundGenes):
                        log.write(f"{gene}\n")
                else:
                    log.write("No genes found in ref file.")
                    log.write(f"Genes: {genesToExcludeFromReference}")
    else:
        for header, prot in parseFasta(referenceFasta):
            refNmers.update(getPeptides(prot, header, cleaveDict, pepLens))
    
    return refNmers

def buildPeptideTranscriptMap(inputFasta, pepLens=[8, 9, 10, 11], cleaveFile=None, threshold=0.2):
    peptide_transcripts = defaultdict(list)
    # Key (str): The amino acid sequence of a peptide.
    # Value (list): A list of transcript IDs (strings) containing that peptide.
    # Structure: { "Peptide1": ["transcriptId1", "TranscriptId2", ...], "Peptide2": [...], ... }

    cleaveDict = getCleavageIndicesDict(cleaveFile, threshold) if cleaveFile else None

    for header, orf in parseFasta(inputFasta):
        transcript_id = extractTranscriptId(header, 'inputFasta')
        for pep in getPeptides(orf, header, cleaveDict, pepLens):
            peptide_transcripts[pep].append(transcript_id)

    return peptide_transcripts

def main():
    # Input Validation
    if not args.input:
        print('Supply input protein database fasta file. Usage: python generatePeptides.py -h')
        exit(1)
    if not args.output_txt and not args.output_csv:
        print('No Output file specified. Usage: python generatePeptides.py -h')
        exit(1)
    # Arguments
    inputFile = args.input
    cleaveFile = args.cleavage_prediction if args.cleavage_prediction else None
    referenceFile = args.reference_file if args.reference_file else None
    refCleaveFile = args.reference_cleavage_prediction if referenceFile and args.reference_cleavage_prediction else None
    cleavageThreshold = args.cleavage_threshold
    gff = args.gff if args.gff else None
    genesToExcludeFromReference = args.exclude_genes if args.exclude_genes else None
    pepLens = args.peptide_lengths
    pepFile = open(args.output_txt, 'w') if args.output_txt else None
    csvFile = open(args.output_csv, 'w') if args.output_csv else None

    # Main
    allNmers = buildPeptideTranscriptMap(inputFile, pepLens, cleaveFile, cleavageThreshold)

    # Optional Filtering out canonical peptides
    if referenceFile:
        referencePepSet = buildReferenceSet(referenceFile, pepLens, refCleaveFile, cleavageThreshold, gff, genesToExcludeFromReference, 
                                                logFile=f"{inputFile.rsplit('.', 1)[0]}_excluded_genes.log")

        nmersToRemove = referencePepSet & set(allNmers.keys())
        for pep in nmersToRemove:
            del allNmers[pep]
    
    # Writing output
    for pep, transcripts in allNmers.items():
        if pepFile:
            pepFile.write(f"{pep}\n")
        if csvFile:
            csvFile.write(f'{pep}, {";".join(transcripts)}, {len(transcripts)}\n')
    if pepFile:
        pepFile.close()
    if csvFile:
        csvFile.close()

if __name__ == '__main__':
    main()