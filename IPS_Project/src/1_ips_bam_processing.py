import os
import pysam
import pandas as pd

def process_bam_files(sample_list_file, input_dir, output_dir, region="X"):
    """
    Process BAM files based on a sample list.

    Parameters:
    - sample_list_file (str): Path to the sample list file (CSV).
    - input_dir (str): Directory containing input BAM files.
    - output_dir (str): Directory to save processed BAM files.
    - region (str): Chromosome or region to filter reads from (default is "X").
    """
    # Ensure input and output directories exist
    os.makedirs(input_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)

    # Read the sample list
    sample_df = pd.read_csv(sample_list_file, header=None)
    sample_list = sample_df.values.flatten().tolist()

    for sample in sample_list:
        print(f"Processing sample: {sample}")
        
        input_bam_path = os.path.join(input_dir, f"{sample}Aligned.sortedByCoord.out.markeddups.bam")
        temp_bam_path = os.path.join(output_dir, f"{sample}_X_cleaned.bam")
        sorted_bam_path = os.path.join(output_dir, f"{sample}_X_cleaned_sorted.bam")
        
        # Open the BAM file for reading
        bam_file = pysam.AlignmentFile(input_bam_path, "rb")
        
        # Open a new BAM file for writing with the filtered reads
        filtered_bam = pysam.AlignmentFile(temp_bam_path, 'wb', template=bam_file)
        
        # Filter and modify reads
        for read in bam_file.fetch(region):
            barcode = read.query_name
            read.set_tag("UB", barcode)
            read.set_tag("CB", sample)  # Set barcode tag
            filtered_bam.write(read)
        
        bam_file.close()
        filtered_bam.close()
        
        # Sort the BAM file
        pysam.sort("-o", sorted_bam_path, temp_bam_path)
        
        # Index the sorted BAM file
        pysam.index(sorted_bam_path)
        
        # Remove the temporary BAM file
        os.remove(temp_bam_path)

# Example usage
#process_bam_files(
    #sample_list_file="/net/noble/vol8/es1/ips_project/NPC_sample_list.txt",
    #input_dir="/net/noble/vol8/es1/ips_project/NPC/",
    #output_dir="/net/noble/vol8/es1/ips_project/NPC_processed/"
#)

#process_bam_files(
    #sample_list_file="/net/noble/vol8/es1/ips_project/Cardiomyocytes_sample_list.txt",
    #input_dir="/net/noble/vol8/es1/ips_project/Cardiomyocytes/",
    #output_dir="/net/noble/vol8/es1/ips_project/Cardiomyocytes_processed/"
#)

process_bam_files(
    sample_list_file="/net/noble/vol8/es1/ips_project/Cardiomyocytes_XY_samples.txt",
    input_dir="/net/noble/vol8/es1/ips_project/Cardiomyocytes/",
    output_dir="/net/noble/vol8/es1/ips_project/Cardiomyocytes_processed/"
)

process_bam_files(
    sample_list_file="/net/noble/vol8/es1/ips_project/NPC_XY_samples.txt",
    input_dir="/net/noble/vol8/es1/ips_project/NPC/",
    output_dir="/net/noble/vol8/es1/ips_project/NPC_processed/"
)

#process_bam_files(
#    sample_list_file="/net/noble/vol8/es1/ips_project/NPC_sample_list.txt",
#    input_dir="/net/noble/vol8/es1/ips_project/NPC/",
#    output_dir="/net/noble/vol8/es1/ips_project/NPC_processed_chr11/",
#    region="11"
#)

#process_bam_files(
#    sample_list_file="/net/noble/vol8/es1/ips_project/Cardiomyocytes_sample_list.txt",
#    input_dir="/net/noble/vol8/es1/ips_project/Cardiomyocytes/",
#    output_dir="/net/noble/vol8/es1/ips_project/Cardiomyocytes_processed_chr11/",
#    region="11"
#)

#process_bam_files(
#    sample_list_file="/net/noble/vol8/es1/ips_project/NPC_sample_list.txt",
#    input_dir="/net/noble/vol8/es1/ips_project/NPC/",
#    output_dir="/net/noble/vol8/es1/ips_project/NPC_processed_chr14/",
#    region="14"
#)

#process_bam_files(
#    sample_list_file="/net/noble/vol8/es1/ips_project/Cardiomyocytes_sample_list.txt",
#    input_dir="/net/noble/vol8/es1/ips_project/Cardiomyocytes/",
#    output_dir="/net/noble/vol8/es1/ips_project/Cardiomyocytes_processed_chr14/",
#    region="14"
#)
