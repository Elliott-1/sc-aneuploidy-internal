#given a list of sample names, 
# and a directory to a folder, (Make it if it doesn't exist)
# Make .tsv files that contain the sample name as the first and only line and are the filename

import os
import pandas as pd

def create_tsv_files(sample_list_file, output_dir):
    """
    Creates .tsv files containing sample names.

    Parameters:
    - sample_list_file (str): Path to the .txt file containing sample names (one per line).
    - output_dir (str): Path to the output directory where .tsv files will be created.
    """
    # Read the sample list from the file
    sample_df = pd.read_csv(sample_list_file, header=None)
    sample_list = sample_df.values.flatten().tolist()

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Create a .tsv file for each sample
    for sample in sample_list:
        tsv_path = os.path.join(output_dir, f"{sample}.tsv")
        with open(tsv_path, 'w') as tsv_file:
            tsv_file.write(sample + "\n")
        print(f"Created {tsv_path}")

# Example usage

create_tsv_files("/net/noble/vol8/es1/ips_project/NPC_XY_samples.txt", "/net/noble/vol8/es1/ips_project/NPC_cellnames")

create_tsv_files("/net/noble/vol8/es1/ips_project/Cardiomyocytes_XY_samples.txt", "/net/noble/vol8/es1/ips_project/Cardiomyocytes_cellnames")
