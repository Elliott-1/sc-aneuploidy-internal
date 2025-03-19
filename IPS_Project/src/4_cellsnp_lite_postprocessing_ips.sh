#!/bin/bash

# Set the base directory where the cellSNP output folders are located
#base_dir="/net/noble/vol8/es1/ips_project/cellsnp_lite_output_chr11"
#base_dir="/net/noble/vol8/es1/ips_project/cellsnp_lite_output_chr14"
#base_dir="/net/noble/vol8/es1/ips_project/cellsnp_lite_output"
#base_dir="/net/noble/vol8/es1/ips_project/cellsnp_lite_output"
base_dir="/net/noble/vol8/es1/ips_project/cellsnp_lite_output_XY"
# Set the path to the R script
rscript_path="/net/noble/vol1/home/es1/aneuploidy/RCODE_process_cellsnp.R"

# Loop through each folder in the base directory
for folder in "$base_dir"/*; do
    if [ -d "$folder" ]; then  # Check if it's a directory
        folder_name=$(basename "$folder")
        
        # Define the file paths
        cellSNP_samples="$folder/cellSNP.samples.tsv"
        cellSNP_tag_AD="$folder/cellSNP.tag.AD.mtx"
        cellSNP_tag_DP="$folder/cellSNP.tag.DP.mtx"
        cellSNP_tag_OTH="$folder/cellSNP.tag.OTH.mtx"
        cellSNP_base_vcf="$folder/cellSNP.base.vcf"
        
        # Define the output file name
        output_file="$folder/${folder_name}_summary.tsv.gz"

        # Run the R script with the appropriate arguments
        Rscript "$rscript_path" \
            "$cellSNP_samples" \
            "$cellSNP_tag_AD" \
            "$cellSNP_tag_DP" \
            "$cellSNP_tag_OTH" \
            "$cellSNP_base_vcf" \
            "$output_file"

        echo "Processed $folder_name, output saved to $output_file"
    fi
done