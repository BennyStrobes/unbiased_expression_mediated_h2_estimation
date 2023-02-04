#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-7:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=10GB                         # Memory total in MiB (for all cores)



chrom_num="$1"
baseld_annotation_dir="$2"
regen_baseline_annotation_dir="$3"
ldsc_genotype_dir="$4"
ldsc_weights_dir="$5"
hapmap3_snplist="$6"
ldsc_code_dir="$7"


source /n/groups/price/ben/environments/sldsc/bin/activate
module load python/2.7.12

if false; then
cp $baseld_annotation_dir"baselineLD."${chrom_num}".annot.gz" $regen_baseline_annotation_dir"baselineLD."${chrom_num}".annot.gz"

python ${ldsc_code_dir}ldsc.py --l2 --bfile ${ldsc_genotype_dir}"1000G.EUR.hg38."${chrom_num} --ld-wind-cm 1 --annot ${regen_baseline_annotation_dir}"baselineLD."${chrom_num}".annot.gz" --out ${regen_baseline_annotation_dir}"baselineLD."${chrom_num} --print-snps $hapmap3_snplist
fi



source ~/.bash_profile
python3 remove_eqtl_anno_from_baseline.py $regen_baseline_annotation_dir"baselineLD."${chrom_num}".annot.gz" $regen_baseline_annotation_dir"baselineLD_no_qtl."${chrom_num}".annot"


source /n/groups/price/ben/environments/sldsc/bin/activate
module load python/2.7.12

python ${ldsc_code_dir}ldsc.py --l2 --bfile ${ldsc_genotype_dir}"1000G.EUR.hg38."${chrom_num} --ld-wind-cm 1 --annot ${regen_baseline_annotation_dir}"baselineLD_no_qtl."${chrom_num}".annot" --out ${regen_baseline_annotation_dir}"baselineLD_no_qtl."${chrom_num} --print-snps $hapmap3_snplist
