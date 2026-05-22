import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--cta', help='Path/To/ CTA file')
parser.add_argument('-g', '--gff', help='Path/To/ GFF reference file')
parser.add_argument('-f', '--file', help='Path/To/ input file')
parser.add_argument('-t', '--tmp', help='tmp file')
args = parser.parse_args()

def getCtaSet(ctaTxt):
    ctas = set()
    with open(ctaTxt) as f:
        for line in f:
            ctas.add(line.strip())
    return ctas

def getTranscriptIDsFromGff(genesToExcludeFromReference: set, gff: str):
    transcriptIDs = set()
    notFound = genesToExcludeFromReference

    with open(gff) as f:
        for line in f:
            if not line.startswith('#'):
                tokens = line.split('\t', 8)
                desc = tokens[8]

                if desc[0] == 'I':
                    # Get Gene name from gff description line
                    attrs = {}
                    for field in desc.strip().split(';'):
                        if '=' in field:
                            key, value = field.split('=', 1)  # maxsplit=1 handles values containing '='
                            attrs[key] = value
                    name = attrs.get('Name')
                    if name:
                        name = re.sub(r'-\d+$', '', name)  # remove trailing hyphen+digits
                        if name in genesToExcludeFromReference:
                            transcriptIDs.add(desc.split(';')[0].split(':')[1])
                            notFound.remove(name)
    return transcriptIDs, notFound

if __name__ == '__main__':
    # ARGS checking
    shouldFind = set()

    with open(args.tmp) as tmp:
        for line in tmp:
            shouldFind.add(line.strip().split('+')[1].strip())

    with open(args.file) as f:
        genes = set()
        for line in f:
            genes.add(line.strip())

        print(genes)
        print(len(genes))

        print(len(shouldFind))

        print(shouldFind-genes)