#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-24:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=20GB                         # Memory total in MiB (for all cores)




tissue_names_file="$1"
variant_annotation_dir="$2"
ldsc_code_dir="$3"
baseld_annotation_dir="$4"
ldsc_genotype_dir="$5"
ldsc_weights_dir="$6"
hapmap3_snplist="$7"
tissue_name="$8"


source /n/groups/price/ben/environments/sldsc/bin/activate
module load python/2.7.12

# Cis with intercept
for chrom_num in $(seq 1 22); do
	annotation_name=${tissue_name}"_cis_eqtl_gene_position_anno_1_mb_w_intercept."
	python ${ldsc_code_dir}ldsc.py --l2 --bfile ${ldsc_genotype_dir}"1000G.EUR.hg38."${chrom_num} --ld-wind-cm 1 --annot ${variant_annotation_dir}${annotation_name}${chrom_num}".annot" --out ${variant_annotation_dir}${annotation_name}${chrom_num} --print-snps $hapmap3_snplist
done

# cis w/o intercept
for chrom_num in $(seq 1 22); do
	annotation_name=${tissue_name}"_cis_eqtl_gene_position_anno_1_mb."
	python ${ldsc_code_dir}ldsc.py --l2 --bfile ${ldsc_genotype_dir}"1000G.EUR.hg38."${chrom_num} --ld-wind-cm 1 --annot ${variant_annotation_dir}${annotation_name}${chrom_num}".annot" --out ${variant_annotation_dir}${annotation_name}${chrom_num} --print-snps $hapmap3_snplist
done

# trans
for chrom_num in $(seq 1 22); do
	annotation_name=${tissue_name}"_trans_eqtl_gene_position_anno."
	python ${ldsc_code_dir}ldsc.py --l2 --bfile ${ldsc_genotype_dir}"1000G.EUR.hg38."${chrom_num} --ld-wind-cm 1 --annot ${variant_annotation_dir}${annotation_name}${chrom_num}".annot" --out ${variant_annotation_dir}${annotation_name}${chrom_num} --print-snps $hapmap3_snplist
done



