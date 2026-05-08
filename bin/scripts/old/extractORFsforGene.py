import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-g', '--genes', help='Path to txt file containing 1 gene symbol per line')
parser.add_argument('-r', '--reference-file', help='Path to Reference Fasta file to filter reads from')
parser.add_argument('-o', '--out-file', help='Output File')
args = parser.parse_args()

if __name__ == '__main__':
    gene_name = "HORMAD1"
    ref = args.reference_file if args.reference_file else None
    outFile = open(args.out_file, 'w') if args.out_file else None 

    if not ref or not outFile:
        print('reference file or out file not provided')
        exit(1)

    addLine = False
    with open(ref) as f:
        for line in f:
            if line.startswith('>'):
                if gene_name in line:
                    addLine = True
                    outFile.write(line)
                else:
                    addLine = False
            elif addLine:
                outFile.write(line)

    outFile.close()
