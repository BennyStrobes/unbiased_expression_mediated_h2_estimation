#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-20:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=10GB                         # Memory total in MiB (for all cores)









tissue_names_file="$1"
expressed_genes_dir="$2"
baseld_annotation_dir="$3"
variant_annotation_dir="$4"



# Create genomic annotation based on cis-windows and trans-windows
if false; then
sed 1d $tissue_names_file | while read tissue_name sample_size pseudotissue_name; do
	echo $tissue_name
	sbatch generate_gene_position_based_annotation_files_in_a_single_tissue.sh $tissue_name $expressed_genes_dir $baseld_annotation_dir $variant_annotation_dir
done
fi

