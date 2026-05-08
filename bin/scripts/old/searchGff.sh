#!/bin/bash

# Usage: ./check_ctas.sh -c ctas.csv -g reference.gff3
# CSV format: any first column, gene name in second column

usage() {
    echo "Usage: $0 -c <ctas.csv> -g <reference.gff3>"
    exit 1
}

while getopts "c:g:" opt; do
    case $opt in
        c) CSV="$OPTARG" ;;
        g) GFF="$OPTARG" ;;
        *) usage ;;
    esac
done

if [[ -z "$CSV" || -z "$GFF" ]]; then
    usage
fi

if [[ ! -f "$CSV" ]]; then
    echo "Error: CSV file not found: $CSV"
    exit 1
fi
if [[ ! -f "$GFF" ]]; then
    echo "Error: GFF file not found: $GFF"
    exit 1
fi

awk -F'\t' '
    # --- Load gene names from CSV ---
    NR == FNR {
        # Skip header
        if (NR == 1) next

        # Split on comma, extract second column
        n = split($0, cols, ",")
        gene = cols[2]

        # Strip quotes and whitespace
        gsub(/"/, "", gene)
        gsub(/^[ \t]+|[ \t]+$/, "", gene)

        if (gene != "") {
            ctas[gene] = 0      # 0 = not yet found
        }
        next
    }

    # --- Scan GFF file once ---
    /^#/ { next }               # skip comment lines

    substr($9, 1, 1) == "I" {
        # Extract Name= value from column 9
        match($9, /Name=([^;]+)/, arr)
        name = arr[1]
        # Remove trailing -digits suffix to get base gene name
        sub(/-[0-9]+$/, "", name)
        if (name in ctas) {
            ctas[name] = 1      # mark as found
        }
    }

    # --- Print summary after GFF is processed ---
    END {
        found     = 0
        not_found = 0
        for (gene in ctas) {
            if (ctas[gene] == 1) {
                found++
                found_genes[found] = gene
            } else {
                not_found++
                missing[not_found] = gene
            }
        }
        print "================================"
        print "Genes found in GFF:     " found
        print "Genes not found in GFF: " not_found
        print "================================"
        if (found > 0) {
            print "Found genes:"
            for (i = 1; i <= found; i++) {
                print "  + " found_genes[i]
            }
        }
        if (not_found > 0) {
            print "Missing genes:"
            for (i = 1; i <= not_found; i++) {
                print "  - " missing[i]
            }
        }
    }
' "$CSV" "$GFF"