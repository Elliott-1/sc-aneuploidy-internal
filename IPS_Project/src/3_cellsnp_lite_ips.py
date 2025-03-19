import os
import subprocess
import pandas as pd

def run_cellsnp_lite(sample_list_file, input_dir, cellnames_dir, output_dir, genome_vcf, region="X", threads=8, min_mapq=20):
    """
    Run cellsnp-lite for a given sample list.

    Parameters:
    - sample_list_file (str): Path to the sample list file (CSV).
    - input_dir (str): Directory containing input BAM files.
    - cellnames_dir (str): Directory containing cell names TSV files.
    - output_dir (str): Directory to save cellsnp-lite output.
    - genome_vcf (str): Path to the genome VCF file.
    - region (str): Chromosome or region to filter reads from (default is "X").
    - threads (int): Number of threads for cellsnp-lite (default is 8).
    - min_mapq (int): Minimum MAPQ value for filtering reads (default is 20).
    """
    # Read the sample list
    sample_df = pd.read_csv(sample_list_file, header=None)
    sample_list = sample_df.values.flatten().tolist()

    for sample in sample_list:
        print(f"Processing sample: {sample}")
        
        bam_rna = os.path.join(input_dir, f"{sample}_X_cleaned_sorted.bam")
        barcode_rna = os.path.join(cellnames_dir, f"{sample}.tsv")
        output_sample_dir = os.path.join(output_dir, sample)  # Create a sample-specific output folder
        os.makedirs(output_sample_dir, exist_ok=True)

        # Construct the cellsnp-lite command
        command = [
            "cellsnp-lite",
            "-s", bam_rna,
            "--minCOUNT", "5",
            "-b", barcode_rna,
            "-R", genome_vcf,
            "--UMItag", "UB",
            "-p", str(threads),
            "--minMAF", "0.05",
            "--minMAPQ", str(min_mapq),
            "--chrom", region,
            "-O", output_sample_dir
        ]

        # Run the command
        try:
            print(f"Running cellsnp-lite for sample {sample}...")
            subprocess.run(command, check=True)
            print(f"Completed cellsnp-lite for sample {sample}.")
        except subprocess.CalledProcessError as e:
            print(f"Error: cellsnp-lite failed for {sample}. Details: {e}")

# Example usage for Cardiomyocytes and NPC
output_dir = "/net/noble/vol8/es1/ips_project/cellsnp_lite_output_XY"
genome_vcf = "/net/noble/vol8/es1/for_scLinaX/Disteche/filtered_chrX_file.vcf"

#output_dir = "/net/noble/vol8/es1/ips_project/cellsnp_lite_output_chr11"
#genome_vcf = "/net/noble/vol8/es1/ips_project/filtered_chr11_file.vcf"

#output_dir = "/net/noble/vol8/es1/ips_project/cellsnp_lite_output_chr14"
#genome_vcf = "/net/noble/vol8/es1/ips_project/filtered_chr14_file.vcf"


# For Cardiomyocytes
#run_cellsnp_lite(
#    sample_list_file="/net/noble/vol8/es1/ips_project/Cardiomyocytes_sample_list.txt",
#    input_dir="/net/noble/vol8/es1/ips_project/Cardiomyocytes_processed/",
#    cellnames_dir="/net/noble/vol8/es1/ips_project/Cardiomyocytes_cellnames/",
#    output_dir=output_dir,
#    genome_vcf=genome_vcf
#)

# For NPC
#run_cellsnp_lite(
#    sample_list_file="/net/noble/vol8/es1/ips_project/NPC_sample_list.txt",
#    input_dir="/net/noble/vol8/es1/ips_project/NPC_processed/",
#    cellnames_dir="/net/noble/vol8/es1/ips_project/NPC_cellnames/",
#    output_dir=output_dir,
#    genome_vcf=genome_vcf
#)

# For Cardiomyocytes XY
run_cellsnp_lite(
    sample_list_file="/net/noble/vol8/es1/ips_project/Cardiomyocytes_XY_samples.txt",
    input_dir="/net/noble/vol8/es1/ips_project/Cardiomyocytes_processed/",
    cellnames_dir="/net/noble/vol8/es1/ips_project/Cardiomyocytes_cellnames/",
    output_dir=output_dir,
    genome_vcf=genome_vcf
)

# For NPC XY
run_cellsnp_lite(
    sample_list_file="/net/noble/vol8/es1/ips_project/NPC_XY_samples.txt",
    input_dir="/net/noble/vol8/es1/ips_project/NPC_processed/",
    cellnames_dir="/net/noble/vol8/es1/ips_project/NPC_cellnames/",
    output_dir=output_dir,
    genome_vcf=genome_vcf
)

# For Cardiomyocytes
#run_cellsnp_lite(
#    sample_list_file="/net/noble/vol8/es1/ips_project/Cardiomyocytes_sample_list.txt",
#    input_dir="/net/noble/vol8/es1/ips_project/Cardiomyocytes_processed_chr11/",
#    cellnames_dir="/net/noble/vol8/es1/ips_project/Cardiomyocytes_cellnames/",
#    output_dir=output_dir,
#    genome_vcf=genome_vcf,
#    region="11"
#)

# For NPC
#run_cellsnp_lite(
#    sample_list_file="/net/noble/vol8/es1/ips_project/NPC_sample_list.txt",
#    input_dir="/net/noble/vol8/es1/ips_project/NPC_processed_chr11/",
#    cellnames_dir="/net/noble/vol8/es1/ips_project/NPC_cellnames/",
#    output_dir=output_dir,
#    genome_vcf=genome_vcf,
#    region="11"
#)


#run_cellsnp_lite(
#    sample_list_file="/net/noble/vol8/es1/ips_project/Cardiomyocytes_sample_list.txt",
#    input_dir="/net/noble/vol8/es1/ips_project/Cardiomyocytes_processed_chr14/",
#    cellnames_dir="/net/noble/vol8/es1/ips_project/Cardiomyocytes_cellnames/",
#    output_dir=output_dir,
#    genome_vcf=genome_vcf,
#    region="14"
#)

# For NPC
#run_cellsnp_lite(
#    sample_list_file="/net/noble/vol8/es1/ips_project/NPC_sample_list.txt",
#    input_dir="/net/noble/vol8/es1/ips_project/NPC_processed_chr14/",
#    cellnames_dir="/net/noble/vol8/es1/ips_project/NPC_cellnames/",
#    output_dir=output_dir,
#    genome_vcf=genome_vcf,
#    region="14"
#)
