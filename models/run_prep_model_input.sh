#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# run_prep_model.sh  –  Interactively pick input files for prep_model_input.py
# ─────────────────────────────────────────────────────────────────────────────

#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# run_prep_model.sh  –  Interactively pick input files for prep_model_input.py
# ─────────────────────────────────────────────────────────────────────────────

#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# run_prep_model.sh  –  Interactively pick input files for prep_model_input.py
# ─────────────────────────────────────────────────────────────────────────────

#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# run_prep_model.sh  –  Interactively pick input files for prep_model_input.py
# ─────────────────────────────────────────────────────────────────────────────

#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# run_prep_model.sh  –  Interactively pick input files for prep_model_input.py
# ─────────────────────────────────────────────────────────────────────────────

set -uo pipefail

# ── Colours ──────────────────────────────────────────────────────────────────
BOLD='\033[1m'
CYAN='\033[1;36m'
GREEN='\033[1;32m'
YELLOW='\033[1;33m'
RED='\033[1;31m'
RESET='\033[0m'

# ── Data directories – each is one feature group (separated by '+') ──────────
DATA_DIRS=(
    "data/hla"
    "data/immunogenicity"
    "data/pathogenicity"
    "data/selfsimilarity"
    "data/tap"
)

# mtec: top-level files + all subdirs together form ONE group (no '+' within)
MTEC_DIR="data/mtec"

# ── Group storage ─────────────────────────────────────────────────────────────
# CMD_ARGS is built directly as we go.
# GROUP_BOUNDARY is a parallel bool array: CMD_ARGS[i]==1 means "insert + before this group".
# Simpler approach: we build CMD_ARGS in-place and track whether we need a '+' next.
CMD_ARGS=()
NEED_SEPARATOR=false   # flip to true after the first non-empty group is added

# ── Helper: list regular files in a directory (one level only) ───────────────
list_files() {
    local dir="$1"
    [[ ! -d "$dir" ]] && return
    find "$dir" -maxdepth 1 -type f | sort
}

# ── Helper: interactively pick files; result → global SELECTED_FILES ─────────
pick_files_from_dir() {
    local dir="$1"
    SELECTED_FILES=()
    local -a files=()

    while IFS= read -r f; do
        [[ -n "$f" ]] && files+=("$f")
    done < <(list_files "$dir")

    if [[ ${#files[@]} -eq 0 ]]; then
        echo -e "  ${YELLOW}⚠  No files found in ${dir} – skipping.${RESET}"
        return
    fi

    echo -e "\n${CYAN}${BOLD}── ${dir} ──────────────────────────────────${RESET}"
    echo -e "  ${BOLD}Available files:${RESET}"
    local i=1
    for f in "${files[@]}"; do
        printf "    ${GREEN}%2d)${RESET}  %s\n" "$i" "$(basename "$f")"
        i=$(( i + 1 ))
    done
    echo -e "    ${GREEN} a)${RESET}  Select ALL"
    echo -e "    ${GREEN} n)${RESET}  Select NONE (skip this directory)"

    echo -e "\n  Enter numbers separated by spaces (e.g. ${BOLD}1 3${RESET}),"
    echo -e "  or ${BOLD}a${RESET} for all, ${BOLD}n${RESET} to skip:"
    printf "  > "
    read -r choice

    if [[ "$choice" == "n" || "$choice" == "N" ]]; then
        echo -e "  ${YELLOW}Skipped.${RESET}"
        return
    fi

    if [[ "$choice" == "a" || "$choice" == "A" ]]; then
        SELECTED_FILES=("${files[@]}")
        echo -e "  ${GREEN}Selected all ${#files[@]} file(s).${RESET}"
        return
    fi

    local valid=0
    for num in $choice; do
        if [[ "$num" =~ ^[0-9]+$ ]] && (( num >= 1 && num <= ${#files[@]} )); then
            SELECTED_FILES+=("${files[$(( num - 1 ))]}")
            valid=$(( valid + 1 ))
        else
            echo -e "  ${RED}  ✗ Ignoring invalid entry: ${num}${RESET}"
        fi
    done

    if [[ $valid -eq 0 ]]; then
        echo -e "  ${YELLOW}No valid selection – skipping.${RESET}"
    else
        echo -e "  ${GREEN}Selected ${valid} file(s).${RESET}"
    fi
}

# ── Append a completed group's files into CMD_ARGS (with '+' separator) ──────
# Call this once per logical group, passing the files as arguments.
# Usage: commit_group file1 file2 ...
commit_group() {
    if [[ $# -eq 0 ]]; then
        return
    fi
    if [[ "$NEED_SEPARATOR" == true ]]; then
        CMD_ARGS+=("+")
    fi
    CMD_ARGS+=("$@")
    NEED_SEPARATOR=true
}

# ── Collect one standard feature-group directory ─────────────────────────────
collect_group() {
    local dir="$1"
    SELECTED_FILES=()
    pick_files_from_dir "$dir"
    if [[ ${#SELECTED_FILES[@]} -gt 0 ]]; then
        commit_group "${SELECTED_FILES[@]}"
    fi
}

# ── Collect mtec as ONE group (top-level + subdirs, no '+' between them) ──────
collect_mtec_group() {
    local base="$1"

    if [[ ! -d "$base" ]]; then
        echo -e "  ${YELLOW}⚠  Directory ${base} not found – skipping.${RESET}"
        return
    fi

    echo -e "\n${CYAN}${BOLD}════ ${base} (top-level + subdirectories = one group) ════${RESET}"

    local -a mtec_files=()

    # Top-level files
    SELECTED_FILES=()
    pick_files_from_dir "$base"
    if [[ ${#SELECTED_FILES[@]} -gt 0 ]]; then
        mtec_files+=("${SELECTED_FILES[@]}")
    fi

    # One level of subdirectories
    local -a subdirs=()
    while IFS= read -r sd; do
        [[ -n "$sd" ]] && subdirs+=("$sd")
    done < <(find "$base" -mindepth 1 -maxdepth 1 -type d | sort)

    if [[ ${#subdirs[@]} -eq 0 ]]; then
        echo -e "  ${YELLOW}  (No subdirectories found in ${base})${RESET}"
    else
        for sd in "${subdirs[@]}"; do
            SELECTED_FILES=()
            pick_files_from_dir "$sd"
            if [[ ${#SELECTED_FILES[@]} -gt 0 ]]; then
                mtec_files+=("${SELECTED_FILES[@]}")
            fi
        done
    fi

    if [[ ${#mtec_files[@]} -gt 0 ]]; then
        commit_group "${mtec_files[@]}"
    fi
}

# ═════════════════════════════════════════════════════════════════════════════
# MAIN
# ═════════════════════════════════════════════════════════════════════════════

clear
echo -e "${CYAN}${BOLD}"
echo "╔══════════════════════════════════════════════════════╗"
echo "║          prep_model_input.py  –  File Picker         ║"
echo "╚══════════════════════════════════════════════════════╝"
echo -e "${RESET}"

CMD_ARGS=()
NEED_SEPARATOR=false

# ── Step 1: collect files per feature group ───────────────────────────────────
for dir in "${DATA_DIRS[@]}"; do
    collect_group "$dir"
    if [[ "$dir" == "data/selfsimilarity" ]]; then
        collect_mtec_group "$MTEC_DIR"
    fi
done

# ── Guard: at least one file must have been selected ─────────────────────────
if [[ ${#CMD_ARGS[@]} -eq 0 ]]; then
    echo -e "\n${RED}${BOLD}No files selected. Aborting.${RESET}"
    exit 1
fi

# ── Step 2: confirm selection ─────────────────────────────────────────────────
echo -e "\n${CYAN}${BOLD}── Selected files by feature group ────────────────────${RESET}"
group_num=1
echo -e "  ${BOLD}Group ${group_num}:${RESET}"
for arg in "${CMD_ARGS[@]}"; do
    if [[ "$arg" == "+" ]]; then
        group_num=$(( group_num + 1 ))
        echo -e "  ${BOLD}Group ${group_num}:${RESET}"
    else
        echo -e "    ${GREEN}✔${RESET}  $arg"
    fi
done

# ── Step 3: ask for output file name ─────────────────────────────────────────
echo -e "\n${CYAN}${BOLD}── Output file ─────────────────────────────────────────${RESET}"
printf "  Output CSV filename ${BOLD}[default: prepared_data.csv]${RESET}: "
read -r output_file
output_file="${output_file:-prepared_data.csv}"
[[ "$output_file" != *.csv ]] && output_file="${output_file}.csv"
echo -e "  ${GREEN}Output:${RESET} ${output_file}"

# ── Step 4: show the command preview ─────────────────────────────────────────
echo -e "\n${CYAN}${BOLD}── Command preview ─────────────────────────────────────${RESET}"
echo -e "  ${BOLD}python prep_model_input.py \\${RESET}"
for arg in "${CMD_ARGS[@]}"; do
    if [[ "$arg" == "+" ]]; then
        echo -e "    ${YELLOW}+${RESET} \\"
    else
        echo -e "    $arg \\"
    fi
done
echo -e "    --output ${output_file}"

# ── Step 5: confirm and execute ───────────────────────────────────────────────
echo ""
printf "  ${BOLD}Run now? [Y/n]:${RESET} "
read -r confirm
confirm="${confirm:-Y}"

if [[ "$confirm" =~ ^[Yy]$ ]]; then
    echo -e "\n${GREEN}${BOLD}▶ Running…${RESET}\n"
    python prep_model_input.py "${CMD_ARGS[@]}" --output "$output_file"
    echo -e "\n${GREEN}${BOLD}✔ Done! Output written to: ${output_file}${RESET}"
else
    echo -e "\n${YELLOW}Aborted. Nothing was run.${RESET}"
    exit 0
fi

# set -euo pipefail

# # ── Colours ──────────────────────────────────────────────────────────────────
# BOLD='\033[1m'
# CYAN='\033[1;36m'
# GREEN='\033[1;32m'
# YELLOW='\033[1;33m'
# RED='\033[1;31m'
# RESET='\033[0m'

# # ── Data directories (in the order they will be presented) ───────────────────
# DATA_DIRS=(
#     "data/hla"
#     "data/immunogenicity"
#     "data/pathogenicity"
#     "data/selfsimilarity"
#     "data/tap"
# )

# # mtec is handled separately – top-level files + one level of subdirectories
# MTEC_DIR="data/mtec"

# # ── Helper: list files in a directory ────────────────────────────────────────
# list_files() {
#     local dir="$1"
#     if [[ ! -d "$dir" ]]; then
#         echo ""
#         return
#     fi
#     # List regular files, sorted
#     find "$dir" -maxdepth 1 -type f | sort
# }

# # ── Helper: prompt user to pick files from a directory ───────────────────────
# # Returns selected files (space-separated) via global SELECTED_FILES array
# pick_files_from_dir() {
#     local dir="$1"
#     local -a files=()

#     # Collect available files
#     while IFS= read -r f; do
#         [[ -n "$f" ]] && files+=("$f")
#     done < <(list_files "$dir")

#     if [[ ${#files[@]} -eq 0 ]]; then
#         echo -e "  ${YELLOW}⚠  No files found in ${dir} – skipping.${RESET}"
#         SELECTED_FILES=()
#         return
#     fi

#     echo -e "\n${CYAN}${BOLD}── ${dir} ──────────────────────────────────${RESET}"
#     echo -e "  ${BOLD}Available files:${RESET}"
#     local i=1
#     for f in "${files[@]}"; do
#         printf "    ${GREEN}%2d)${RESET}  %s\n" "$i" "$(basename "$f")"
#         (( i++ ))
#     done
#     echo -e "    ${GREEN} a)${RESET}  Select ALL"
#     echo -e "    ${GREEN} n)${RESET}  Select NONE (skip this directory)"

#     echo -e "\n  Enter numbers separated by spaces (e.g. ${BOLD}1 3${RESET}),"
#     echo -e "  or ${BOLD}a${RESET} for all, ${BOLD}n${RESET} to skip:"
#     printf "  > "

#     read -r choice
#     SELECTED_FILES=()

#     if [[ "$choice" == "n" || "$choice" == "N" ]]; then
#         echo -e "  ${YELLOW}Skipped.${RESET}"
#         return
#     fi

#     if [[ "$choice" == "a" || "$choice" == "A" ]]; then
#         SELECTED_FILES=("${files[@]}")
#         echo -e "  ${GREEN}Selected all ${#files[@]} file(s).${RESET}"
#         return
#     fi

#     # Parse individual numbers
#     local valid=0
#     for num in $choice; do
#         if [[ "$num" =~ ^[0-9]+$ ]] && (( num >= 1 && num <= ${#files[@]} )); then
#             SELECTED_FILES+=("${files[$(( num - 1 ))]}")
#             (( valid++ ))
#         else
#             echo -e "  ${RED}  ✗ Ignoring invalid entry: ${num}${RESET}"
#         fi
#     done

#     if (( valid == 0 )); then
#         echo -e "  ${YELLOW}No valid selection – skipping this directory.${RESET}"
#     else
#         echo -e "  ${GREEN}Selected ${valid} file(s).${RESET}"
#     fi
# }

# # ── Helper: handle mtec – top-level files + one level of subdirectories ──────
# pick_mtec() {
#     local base="$1"

#     if [[ ! -d "$base" ]]; then
#         echo -e "  ${YELLOW}⚠  Directory ${base} not found – skipping.${RESET}"
#         return
#     fi

#     echo -e "\n${CYAN}${BOLD}════ ${base} (top-level + subdirectories) ════${RESET}"

#     # ── top-level files first ────────────────────────────────────────────────
#     SELECTED_FILES=()
#     pick_files_from_dir "$base"
#     if [[ ${#SELECTED_FILES[@]} -gt 0 ]]; then
#         ALL_INPUT_FILES+=("${SELECTED_FILES[@]}")
#     fi

#     # ── one level of subdirectories ──────────────────────────────────────────
#     local -a subdirs=()
#     while IFS= read -r sd; do
#         [[ -n "$sd" ]] && subdirs+=("$sd")
#     done < <(find "$base" -mindepth 1 -maxdepth 1 -type d | sort)

#     if [[ ${#subdirs[@]} -eq 0 ]]; then
#         echo -e "  ${YELLOW}  (No subdirectories found in ${base})${RESET}"
#         return
#     fi

#     for sd in "${subdirs[@]}"; do
#         SELECTED_FILES=()
#         pick_files_from_dir "$sd"
#         if [[ ${#SELECTED_FILES[@]} -gt 0 ]]; then
#             ALL_INPUT_FILES+=("${SELECTED_FILES[@]}")
#         fi
#     done
# }

# # ═════════════════════════════════════════════════════════════════════════════
# # MAIN
# # ═════════════════════════════════════════════════════════════════════════════

# clear
# echo -e "${CYAN}${BOLD}"
# echo "╔══════════════════════════════════════════════════════╗"
# echo "║          prep_model_input.py  –  File Picker         ║"
# echo "╚══════════════════════════════════════════════════════╝"
# echo -e "${RESET}"

# ALL_INPUT_FILES=()

# # ── Step 1: iterate over each data directory ─────────────────────────────────
# for dir in "${DATA_DIRS[@]}"; do
#     SELECTED_FILES=()
#     pick_files_from_dir "$dir"
#     if [[ ${#SELECTED_FILES[@]} -gt 0 ]]; then
#         ALL_INPUT_FILES+=("${SELECTED_FILES[@]}")
#     fi
#     # Insert mtec (with subdirs) after selfsimilarity, before tap
#     if [[ "$dir" == "data/selfsimilarity" ]]; then
#         pick_mtec "$MTEC_DIR"
#     fi
# done

# # ── Guard: at least one file must be selected ────────────────────────────────
# if [[ ${#ALL_INPUT_FILES[@]} -eq 0 ]]; then
#     echo -e "\n${RED}${BOLD}No files selected. Aborting.${RESET}"
#     exit 1
# fi

# # ── Step 2: confirm selection ────────────────────────────────────────────────
# echo -e "\n${CYAN}${BOLD}── Selected input files ───────────────────────────────${RESET}"
# for f in "${ALL_INPUT_FILES[@]}"; do
#     echo -e "  ${GREEN}✔${RESET}  $f"
# done

# # ── Step 3: ask for output file name ─────────────────────────────────────────
# echo -e "\n${CYAN}${BOLD}── Output file ─────────────────────────────────────────${RESET}"
# printf "  Output CSV filename ${BOLD}[default: prepared_data.csv]${RESET}: "
# read -r output_file
# output_file="${output_file:-prepared_data.csv}"
# # Ensure .csv extension
# [[ "$output_file" != *.csv ]] && output_file="${output_file}.csv"
# echo -e "  ${GREEN}Output:${RESET} ${output_file}"

# # ── Step 4: show the command that will be run ────────────────────────────────
# echo -e "\n${CYAN}${BOLD}── Command preview ─────────────────────────────────────${RESET}"
# CMD="python prep_model_input.py ${ALL_INPUT_FILES[*]} --output ${output_file}"
# echo -e "  ${BOLD}${CMD}${RESET}"

# # ── Step 5: confirm and execute ───────────────────────────────────────────────
# echo ""
# printf "  ${BOLD}Run now? [Y/n]:${RESET} "
# read -r confirm
# confirm="${confirm:-Y}"

# if [[ "$confirm" =~ ^[Yy]$ ]]; then
#     echo -e "\n${GREEN}${BOLD}▶ Running…${RESET}\n"
#     python prep_model_input.py "${ALL_INPUT_FILES[@]}" --output "$output_file"
#     echo -e "\n${GREEN}${BOLD}✔ Done! Output written to: ${output_file}${RESET}"
# else
#     echo -e "\n${YELLOW}Aborted. Nothing was run.${RESET}"
#     exit 0
# fi

