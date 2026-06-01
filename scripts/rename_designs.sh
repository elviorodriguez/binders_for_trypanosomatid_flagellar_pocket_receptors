#!/usr/bin/env bash
# rename_designs.sh
# Usage: bash rename_designs.sh -i <input_dir> -o <output_dir> [-f <find_pattern>] [-r <replace_pattern>]

set -euo pipefail

# --- Defaults ---
FIND_PAT="design_"
REPLACE_PAT="HRGx_"

# --- Argument parsing ---
usage() {
    echo "Usage: $0 -i <input_dir> -o <output_dir> [-f <find_pattern>] [-r <replace_pattern>]"
    exit 1
}

while getopts "i:o:f:r:h" opt; do
    case $opt in
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        f) FIND_PAT="$OPTARG" ;;
        r) REPLACE_PAT="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Check required args
[[ -z "${INPUT_DIR:-}" || -z "${OUTPUT_DIR:-}" ]] && usage

# Check input exists
if [[ ! -d "$INPUT_DIR" ]]; then
    echo "Error: Input directory '$INPUT_DIR' does not exist."
    exit 1
fi

# Prevent overwriting input
if [[ "$(realpath "$INPUT_DIR")" == "$(realpath "$OUTPUT_DIR")" ]]; then
    echo "Error: Input and output directories must be different."
    exit 1
fi

# --- Step 1: Copy everything ---
echo "Copying '$INPUT_DIR' -> '$OUTPUT_DIR'..."
cp -r "$INPUT_DIR" "$OUTPUT_DIR"
echo "Copy done."

# --- Step 2: Rename files (deepest first to avoid broken paths) ---
echo "Renaming files containing '$FIND_PAT'..."
find "$OUTPUT_DIR" -depth -name "*${FIND_PAT}*" | while IFS= read -r path; do
    dir="$(dirname "$path")"
    base="$(basename "$path")"
    new_base="${base//$FIND_PAT/$REPLACE_PAT}"
    if [[ "$base" != "$new_base" ]]; then
        mv "$path" "$dir/$new_base"
        echo "  Renamed: $base -> $new_base"
    fi
done

# --- Step 3: Replace pattern inside file contents ---
echo "Replacing '$FIND_PAT' inside file contents..."
find "$OUTPUT_DIR" -depth -type f | while IFS= read -r file; do
    if grep -qF "$FIND_PAT" "$file" 2>/dev/null; then
        sed -i "s/${FIND_PAT}/${REPLACE_PAT}/g" "$file"
        echo "  Updated contents: $file"
    fi
done

echo ""
echo "Done! Output saved to: $OUTPUT_DIR"
