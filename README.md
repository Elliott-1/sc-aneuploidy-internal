# sc-aneuploidy-internal

This repository is similar to the real repository which is unfortunately private until the work and data is publishable. As a result, some results/analysis is omitted. The files in this repo cover some preprocessing, introductory analysis, and the pipeline for an external tool (scLinaX).

It’s structured as following

## Aneuploidy project
- src: Contains all code for X chromosome project- used for data processing and analysis
  - preprocessing: The preprocessing folder is divided into several jupyter notebooks. They can be ran in the following order:
  - make_three_hvg_datasets.ipynb: This notebook takes the original datasets, divides them into three datasets split by origin, and then finds the top 2000 genes
  - run_scanorama.ipynb: This notebook integrates the three datasets used into a shared PCA space, then saves the coordinates
  - leiden_clustering.ipynb: This notebook runs leiden_clustering on the cells from the shared PCA space and top 2000 genes.
  - umi_threshold_and_hvg_to_full_conversion.ipynb: This notebook converts back to the full gene matrix 
  - post_processing.ipynb: This notebook does some brief sub clustering to extract rarer cell types
  - cluster_investigation.ipynb: Runs some validation statistics and transfers the cell type labels
- analysis:
  - differential_gene_expression_analysis.ipynb: This notebook calculates differentially expressed genes for each pairwise comparison and each cell type. This list is used in all downstream analysis
- figures: Assorted figures generated from the data

## IPS-Project: Side project done while waiting for cells - Looks at levels of X chromosome inactivation
- src: Code used to generate SNP files, run in order
  - 1_ips_bam_processing.py: Takes in BAM files and adds arbitrary UMI and cell barcodes
  - 2_make_cellname_files.py: Creates .tsv files of sample names
  - 3_cellsnp_lite_ips.py: Runs cellsnp-lite on each bam file
  - 4_cellsnp_lite_postprocessing_ips.sh: Collects the cellsnp-lite outputs in a complete file
  - 5_annovar_processing_ips.sh: Adds in ANNOVAR annotations
  - 6_annovar_postprocessing.ips.sh: Collects all SNP files into one main file
  - 7_trying_out_proportions_v1.ipynb: Combines replicate counts and adds additional metadata
- Methods PDF: Explains the method to generate these SNP tables
- Links to tables: Links to final tables for processing:
  - https://docs.google.com/spreadsheets/d/1CQfnLvq8gpEQ7XpK7WOVzOyOZ_rAD9XfBwDF_Y6vPeA/edit?usp=sharing
  - https://docs.google.com/spreadsheets/d/1PHVrxm7TMOwc83JmLr0eesgZWjxJmNN6n4aL-DhGaGQ/edit?usp=sharing

