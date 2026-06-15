
import argparse

def get_mtec_expression_file(inFilePath: str, outFilePath: str, threshold: int):
    with open(inFilePath, 'r') as inFile:
        with open(outFilePath, "w") as out:
            for line in inFile:
                line = line.strip()
                if line:
                    if line.startswith('p'):
                        out.write("peptide,mtec_expression\n")
                    else:
                        parts = line.split('\t')
                        pep   = parts[0]
                        total = int(parts[-1])
                        if total > threshold:
                            out.write(f"{pep},{1}\n")
                        else:
                            out.write(f"{pep},{0}\n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Concatenate count columns from multiple peptide counter TSV files"
    )
    parser.add_argument("-o", "--output_file",
                        help="Output TSV: peptide<tab>count_1<tab>count_2<tab>...total")
    parser.add_argument("-c", "--output_expression_file",
                        help="Output CSV: peptide,mTEC_expression")
    parser.add_argument("-t", "--threshold",
                        help="Count Threshold to determine mTEC expression: int")
    args = parser.parse_args()

    get_mtec_expression_file(args.output_file, args.output_expression_file, int(args.threshold))