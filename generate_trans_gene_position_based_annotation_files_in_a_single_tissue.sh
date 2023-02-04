#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-20:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=30GB                         # Memory total in MiB (for all cores)




tissue_name="$1"
expressed_genes_dir="$2"
baseld_annotation_dir="$3"
variant_annotation_dir="$4"



python3 generate_trans_gene_position_based_annotation_files_in_a_single_tissue.py $tissue_name $expressed_genes_dir $baseld_annotation_dir $variant_annotation_dir 