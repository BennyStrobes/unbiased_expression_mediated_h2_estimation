#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-20:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=10GB                         # Memory total in MiB (for all cores)





tissue_names_file="$1"
gtex_egenes_dir="$2"
expressed_genes_dir="$3"



python3 extract_list_of_expressed_genes_in_each_tissue.py $tissue_names_file $gtex_egenes_dir $expressed_genes_dir
