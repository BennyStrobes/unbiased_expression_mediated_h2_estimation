


######################
# Input data
######################
# File name containing list of gtex tissue names
tissue_names_file="/n/scratch3/users/b/bes710/unbiased_expression_mediated_h2_estimation/input_data/tissue_info.txt"

# Directory containing file for each gtex tissue, where each file is a list of tested genes in that tissue (as well as top eqtl variant for that gene)
gtex_egenes_dir="/n/scratch3/users/b/bes710/unbiased_expression_mediated_h2_estimation/input_data/GTEx_Analysis_v8_eQTL/"


# Ldscore regression code
ldsc_code_dir="/n/groups/price/ldsc/ldsc/"

# Directory containing baselineLD variant annotations (in hg38)
baseld_annotation_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3_hg38/baselineLD_v2.2/"

# LDSC 1KG genotype files (hg19)
ldsc_genotype_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3_hg38/plink_files/"

# Ldsc weights
ldsc_weights_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3_hg38/weights/"

# hapmap3_snplist
hapmap3_snplist="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3_hg38/w_hm3.noMHC.snplist"

# Summary statistic directory
sumstat_dir="/n/groups/price/ldsc/sumstats_formatted_2021/"




######################
# Output data
######################
# Output root directory
output_root="/n/scratch3/users/b/bes710/unbiased_expression_mediated_h2_estimation/"

# Directory containing results of simulation
simulation_results_dir=$output_root"simulation_results/"

# Directory containing visualization of simulation run
visualize_simulation_dir=$output_root"visualize_simulation/"

# Directory containing list of expressed genes in each tissue
expressed_genes_dir=$output_root"expressed_genes/"

# Directory continaing variant annotation files
variant_annotation_dir=$output_root"variant_annotations/"

# Directory containing baseline-ld annotations
regen_baseline_annotation_dir=$output_root"baseline_annotations/"

# Directory containing ld-score results
ldsc_results_dir=$output_root"ldsc_results/"

# Directory containing ld-score results viz
ldsc_viz_dir=$output_root"ldsc_viz/"





######################
# Run code
######################

# SIMULATION
sh expression_mediated_h2_simulation_driver_key.sh $simulation_results_dir $visualize_simulation_dir



if false; then
source ~/.bash_profile
module load R/3.5.1
fi

if false; then
Rscript visualize_simulation_results.R $simulation_results_dir $visualize_simulation_dir
fi



# Extract list of expressed genes in each tissue
if false; then
sh extract_list_of_expressed_genes_in_each_tissue.sh $tissue_names_file $gtex_egenes_dir $expressed_genes_dir
fi

# Generate annotation files
if false; then
sh generate_gene_position_based_annotation_files.sh $tissue_names_file $expressed_genes_dir $baseld_annotation_dir $variant_annotation_dir
fi


# Generate ld-score weighted annotation files
if false; then
sed 1d $tissue_names_file | while read tissue_name sample_size pseudotissue_name; do
	sbatch generate_ld_score_annotation_files.sh $tissue_names_file $variant_annotation_dir $ldsc_code_dir $baseld_annotation_dir $ldsc_genotype_dir $ldsc_weights_dir $hapmap3_snplist $tissue_name
done
fi

# Re-generate ld-scores
if false; then
for chrom_num in $(seq 1 22); do
	sbatch regenerate_baseline_annotations.sh $chrom_num $baseld_annotation_dir $regen_baseline_annotation_dir $ldsc_genotype_dir $ldsc_weights_dir $hapmap3_snplist $ldsc_code_dir
done
fi

if false; then
source /n/groups/price/ben/environments/sldsc/bin/activate
module load python/2.7.12
fi
if false; then
trait_name="UKB_460K.blood_WHITE_COUNT"
trait_file=$sumstat_dir$trait_name".sumstats"
python ${ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${variant_annotation_dir}"Cells_Cultured_fibroblasts_cis_eqtl_gene_position_anno_1_mb_w_intercept." --w-ld-chr ${ldsc_weights_dir}"weights.hm3_noMHC." --overlap-annot --print-coefficients --frqfile-chr ${ldsc_genotype_dir}"1000G.EUR.hg38." --out ${ldsc_results_dir}${trait_name}"_sldsc_res_"

trait_name="UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED"
trait_file=$sumstat_dir$trait_name".sumstats"
python ${ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${variant_annotation_dir}"Cells_Cultured_fibroblasts_cis_eqtl_gene_position_anno_1_mb_w_intercept." --w-ld-chr ${ldsc_weights_dir}"weights.hm3_noMHC." --overlap-annot --print-coefficients --frqfile-chr ${ldsc_genotype_dir}"1000G.EUR.hg38." --out ${ldsc_results_dir}${trait_name}"_sldsc_res_"


trait_name="UKB_460K.disease_ASTHMA_DIAGNOSED"
trait_file=$sumstat_dir$trait_name".sumstats"
python ${ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${variant_annotation_dir}"Cells_Cultured_fibroblasts_cis_eqtl_gene_position_anno_1_mb_w_intercept." --w-ld-chr ${ldsc_weights_dir}"weights.hm3_noMHC." --overlap-annot --print-coefficients --frqfile-chr ${ldsc_genotype_dir}"1000G.EUR.hg38." --out ${ldsc_results_dir}${trait_name}"_sldsc_res_"

trait_name="UKB_460K.other_MORNINGPERSON"
trait_file=$sumstat_dir$trait_name".sumstats"
python ${ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${variant_annotation_dir}"Cells_Cultured_fibroblasts_cis_eqtl_gene_position_anno_1_mb_w_intercept." --w-ld-chr ${ldsc_weights_dir}"weights.hm3_noMHC." --overlap-annot --print-coefficients --frqfile-chr ${ldsc_genotype_dir}"1000G.EUR.hg38." --out ${ldsc_results_dir}${trait_name}"_sldsc_res_"

trait_name="UKB_460K.disease_CARDIOVASCULAR"
trait_file=$sumstat_dir$trait_name".sumstats"
python ${ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${variant_annotation_dir}"Cells_Cultured_fibroblasts_cis_eqtl_gene_position_anno_1_mb_w_intercept." --w-ld-chr ${ldsc_weights_dir}"weights.hm3_noMHC." --overlap-annot --print-coefficients --frqfile-chr ${ldsc_genotype_dir}"1000G.EUR.hg38." --out ${ldsc_results_dir}${trait_name}"_sldsc_res_"

trait_name="UKB_460K.blood_EOSINOPHIL_COUNT"
trait_file=$sumstat_dir$trait_name".sumstats"
python ${ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${variant_annotation_dir}"Cells_Cultured_fibroblasts_cis_eqtl_gene_position_anno_1_mb_w_intercept." --w-ld-chr ${ldsc_weights_dir}"weights.hm3_noMHC." --overlap-annot --print-coefficients --frqfile-chr ${ldsc_genotype_dir}"1000G.EUR.hg38." --out ${ldsc_results_dir}${trait_name}"_sldsc_res_"

trait_name="UKB_460K.bp_DIASTOLICadjMEDz"
trait_file=$sumstat_dir$trait_name".sumstats"
python ${ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${variant_annotation_dir}"Cells_Cultured_fibroblasts_cis_eqtl_gene_position_anno_1_mb_w_intercept." --w-ld-chr ${ldsc_weights_dir}"weights.hm3_noMHC." --overlap-annot --print-coefficients --frqfile-chr ${ldsc_genotype_dir}"1000G.EUR.hg38." --out ${ldsc_results_dir}${trait_name}"_sldsc_res_"

trait_name="UKB_460K.cov_EDU_COLLEGE"
trait_file=$sumstat_dir$trait_name".sumstats"
python ${ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${variant_annotation_dir}"Cells_Cultured_fibroblasts_cis_eqtl_gene_position_anno_1_mb_w_intercept." --w-ld-chr ${ldsc_weights_dir}"weights.hm3_noMHC." --overlap-annot --print-coefficients --frqfile-chr ${ldsc_genotype_dir}"1000G.EUR.hg38." --out ${ldsc_results_dir}${trait_name}"_sldsc_res_"

trait_name="UKB_460K.blood_MEAN_CORPUSCULAR_HEMOGLOBIN"
trait_file=$sumstat_dir$trait_name".sumstats"
python ${ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${variant_annotation_dir}"Cells_Cultured_fibroblasts_cis_eqtl_gene_position_anno_1_mb_w_intercept." --w-ld-chr ${ldsc_weights_dir}"weights.hm3_noMHC." --overlap-annot --print-coefficients --frqfile-chr ${ldsc_genotype_dir}"1000G.EUR.hg38." --out ${ldsc_results_dir}${trait_name}"_sldsc_res_"

trait_name="UKB_460K.body_BMIz"
trait_file=$sumstat_dir$trait_name".sumstats"
python ${ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${variant_annotation_dir}"Cells_Cultured_fibroblasts_cis_eqtl_gene_position_anno_1_mb_w_intercept." --w-ld-chr ${ldsc_weights_dir}"weights.hm3_noMHC." --overlap-annot --print-coefficients --frqfile-chr ${ldsc_genotype_dir}"1000G.EUR.hg38." --out ${ldsc_results_dir}${trait_name}"_sldsc_res_"

trait_name="UKB_460K.body_BALDING1"
trait_file=$sumstat_dir$trait_name".sumstats"
python ${ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${variant_annotation_dir}"Cells_Cultured_fibroblasts_cis_eqtl_gene_position_anno_1_mb_w_intercept." --w-ld-chr ${ldsc_weights_dir}"weights.hm3_noMHC." --overlap-annot --print-coefficients --frqfile-chr ${ldsc_genotype_dir}"1000G.EUR.hg38." --out ${ldsc_results_dir}${trait_name}"_sldsc_res_"

fi



if false; then
source ~/.bash_profile
module load R/3.5.1
fi

if false; then
Rscript visualize_ldsc_results.R $ldsc_results_dir $ldsc_viz_dir
fi



