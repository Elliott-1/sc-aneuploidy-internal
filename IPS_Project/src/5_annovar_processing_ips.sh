#!/bin/bash

# Base directory containing the subfolders
#base_dir="/net/noble/vol8/es1/ips_project/cellsnp_lite_output_chr11"
#base_dir="/net/noble/vol8/es1/ips_project/cellsnp_lite_output_chr14"
base_dir="/net/noble/vol8/es1/ips_project/cellsnp_lite_output_XY"
# Loop through each subfolder
for sample_dir in "$base_dir"/*; do
    if [ -d "$sample_dir" ]; then
        # Find the .tsv.gz file in the subfolder
        tsv_file=$(find "$sample_dir" -name "*_summary.tsv.gz")
        if [ -n "$tsv_file" ]; then
            # Extract sample ID from the file name
            sample_id=$(basename "$tsv_file" | cut -d'_' -f1)
            
            # Define output paths in the sample directory
            annov_input="${sample_dir}/${sample_id}.summary.annov.input"
            annovar_output_prefix="${sample_dir}/${sample_id}.chrX.snp.call"
            #annovar_output_prefix="${sample_dir}/${sample_id}.chr11.snp.call"
            #annovar_output_prefix="${sample_dir}/${sample_id}.chr14.snp.call"

            # Process the file with zcat, awk, sort, and uniq
            zcat "$tsv_file" | \
            awk -F"\t" 'NR>1 {print $2"\t"$3"\t"$3"\t"$4"\t"$5"\t"$1}' | \
            LANG=C sort | LANG=C uniq > "$annov_input"

            # Run ANNOVAR
            annovar/table_annovar.pl \
            "$annov_input" \
            annovar/humandb/ -buildver hg38 \
            -out "$annovar_output_prefix" \
            -remove -protocol refGene,ensGene,avsnp150,ALL.sites.2015_08 \
            -operation g,g,f,f -nastring NA -polish
        else
            echo "No *_summary.tsv.gz file found in $sample_dir"
        fi
    fi
done