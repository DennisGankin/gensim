#!/bin/bash

# Usage: ./select_causal_snps.sh <bim_file> <output_dir> <num_snps1> [<num_snps2> ...]
# Example: ./select_causal_snps.sh genotype_dataset.bim output 80 320 1280

if [ $# -lt 3 ]; then
    echo "Usage: $0 <bim_file> <output_dir> <num_snps1> [<num_snps2> ...]"
    exit 1
fi

BIM_FILE=$1
OUTPUT_DIR=$2
shift 2

if [ ! -f "$BIM_FILE" ]; then
    echo "Error: BIM file '$BIM_FILE' does not exist."
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

for NUM_SNPS in "$@"; do
    OUTPUT_FILE="${OUTPUT_DIR}/causal_snps_${NUM_SNPS}.txt"

    # Count available causal SNPs from each category
    COMMON_CAUSAL_SNPS=$(grep "common_causal" "$BIM_FILE" | wc -l)
    RARE_CAUSAL_SNPS=$(grep "rare_causal" "$BIM_FILE" | wc -l)

    TOTAL_CAUSAL_AVAILABLE=$((COMMON_CAUSAL_SNPS + RARE_CAUSAL_SNPS))

    # Initialize output file
    > "$OUTPUT_FILE"

    if [ "$NUM_SNPS" -le "$TOTAL_CAUSAL_AVAILABLE" ]; then
        # Enough labeled causal SNPs, select proportionally
        NUM_COMMON=$(( NUM_SNPS * COMMON_CAUSAL_SNPS / TOTAL_CAUSAL_AVAILABLE ))
        NUM_RARE=$(( NUM_SNPS - NUM_COMMON ))

        grep "common_causal" "$BIM_FILE" | shuf -n $NUM_COMMON | cut -f2 >> "$OUTPUT_FILE"
        grep "rare_causal" "$BIM_FILE" | shuf -n $NUM_RARE | cut -f2 >> "$OUTPUT_FILE"
    else
        # Select all available causal SNPs first
        grep "common_causal" "$BIM_FILE" | cut -f2 >> "$OUTPUT_FILE"
        grep "rare_causal" "$BIM_FILE" | cut -f2 >> "$OUTPUT_FILE"

        REMAINING_SNPS=$(( NUM_SNPS - TOTAL_CAUSAL_AVAILABLE ))

        echo "Not enough labeled causal SNPs. Adding $REMAINING_SNPS SNPs from null categories for ${NUM_SNPS}."

        grep -E "common_null|lowfreq_null|rare_null" "$BIM_FILE" | shuf -n $REMAINING_SNPS | cut -f2 >> "$OUTPUT_FILE"
    fi

    echo "Causal SNPs saved in $OUTPUT_FILE"
done
