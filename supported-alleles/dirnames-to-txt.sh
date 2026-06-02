#!/bin/bash

# Usage: ./list_dirs.sh [-f] [-d] [-p parent_dir] [-o output_file]

PARENT_DIR="."
OUTPUT_FILE="directories.txt"
SEARCH_FILES=false
SEARCH_DIRS=false

while getopts "fdp:o:" opt; do
    case $opt in
        f) SEARCH_FILES=true ;;
        d) SEARCH_DIRS=true ;;
        p) PARENT_DIR="$OPTARG" ;;
        o) OUTPUT_FILE="$OPTARG" ;;
        *) echo "Usage: $0 [-f] [-d] [-p parent_dir] [-o output_file]"; exit 1 ;;
    esac
done

# Default to both if neither -f nor -d is specified
if ! $SEARCH_FILES && ! $SEARCH_DIRS; then
    SEARCH_FILES=true
    SEARCH_DIRS=true
fi

# Build the find type expression
if $SEARCH_FILES && $SEARCH_DIRS; then
    TYPE_EXPR=( -type f -o -type d )
elif $SEARCH_FILES; then
    TYPE_EXPR=( -type f )
else
    TYPE_EXPR=( -type d )
fi

find "$PARENT_DIR" -maxdepth 1 -mindepth 1 \( "${TYPE_EXPR[@]}" \) -exec basename {} \; | sort > "$OUTPUT_FILE"

echo "Results written to $OUTPUT_FILE"