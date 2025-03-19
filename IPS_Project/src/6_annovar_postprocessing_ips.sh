#!/bin/bash

# Path to cellsnp_lite_output folder
#BASE_DIR="/net/noble/vol8/es1/ips_project/cellsnp_lite_output"
#BASE_DIR="/net/noble/vol8/es1/ips_project/cellsnp_lite_output_chr11"
#BASE_DIR="/net/noble/vol8/es1/ips_project/cellsnp_lite_output_chr14"
BASE_DIR="/net/noble/vol8/es1/ips_project/cellsnp_lite_output_XY"

R_SCRIPT="/net/noble/vol1/home/es1/aneuploidy/RCODE_QC.R" # Replace with the full path to the adjusted R script

# Loop through each sampleID folder
for SAMPLE_DIR in "$BASE_DIR"/*/; do
    # Extract the sample ID as the string before the first underscore
    SAMPLE_ID=$(basename "$SAMPLE_DIR") # Extract the sample ID
    SAMPLE_ID_SMALL=$(basename "$SAMPLE_DIR" | cut -d'_' -f1)
    
    RNA_SNPCELL="$SAMPLE_DIR/${SAMPLE_ID}_summary.tsv.gz"
    #RNA_ANNOVAR="$SAMPLE_DIR/${SAMPLE_ID_SMALL}.chr14.snp.call.hg38_multianno.txt"
    #RNA_ANNOVAR="$SAMPLE_DIR/${SAMPLE_ID_SMALL}.chr11.snp.call.hg38_multianno.txt"
    RNA_ANNOVAR="$SAMPLE_DIR/${SAMPLE_ID_SMALL}.chrX.snp.call.hg38_multianno.txt"
    OUTPUT_DIR="$SAMPLE_DIR" # Save results in the same directory
    
    # Check if required files exist
    if [[ -f "$RNA_SNPCELL" && -f "$RNA_ANNOVAR" ]]; then
        echo "Processing $SAMPLE_ID..."
        Rscript "$R_SCRIPT" "$RNA_SNPCELL" "$RNA_ANNOVAR" "$OUTPUT_DIR"
    else
        echo "Required files for $SAMPLE_ID are missing. Skipping..."
    fi
done

