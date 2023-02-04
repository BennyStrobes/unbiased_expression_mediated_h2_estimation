args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(hash)
library(dplyr)
library(reshape)
library(stringr)
options(warn=1)

figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}




extract_partitioned_h2_results <- function(trait_names, results_dir, num_snps, num_genes) {
	trait_arr <- c()
	h2_class_arr <- c()
	h2_arr <- c()
	h2_lower_arr <- c()
	h2_upper_arr <- c()
	for (trait_iter in 1:length(trait_names)) {
		trait_name <- trait_names[trait_iter]
		trait_file <- paste0(results_dir, "UKB_460K.", trait_name, "_sldsc_res_.results")
		df <- read.table(trait_file, header=TRUE)
		
		snp_h2 <- df$Coefficient[1]*num_snps
		snp_h2_se <- df$Coefficient_std_error[1]*num_snps
		snp_h2_lower = snp_h2 - 1.96*snp_h2_se
		snp_h2_upper = snp_h2 + 1.96*snp_h2_se

		trait_arr <- c(trait_arr, trait_name)
		h2_class_arr <- c(h2_class_arr, "h2_g")
		h2_arr <- c(h2_arr, snp_h2)
		h2_lower_arr <- c(h2_lower_arr, snp_h2_lower)
		h2_upper_arr <- c(h2_upper_arr, snp_h2_upper)

		gene_h2 <- df$Coefficient[2]*num_genes
		gene_h2_se <- df$Coefficient_std_error[2]*num_genes
		gene_h2_lower = gene_h2 - 1.96*gene_h2_se
		gene_h2_upper = gene_h2 + 1.96*gene_h2_se

		trait_arr <- c(trait_arr, trait_name)
		h2_class_arr <- c(h2_class_arr, "h2_e")
		h2_arr <- c(h2_arr, gene_h2)
		h2_lower_arr <- c(h2_lower_arr, gene_h2_lower)
		h2_upper_arr <- c(h2_upper_arr, gene_h2_upper)
	}

	df <- data.frame(trait=trait_arr, h2_class=h2_class_arr, h2=h2_arr, h2_upper=h2_upper_arr, h2_lower=h2_lower_arr)
	return(df)
}


make_standard_error_barplot_with_heritabilities <- function(df) {
	p<- ggplot(df, aes(x=trait, y=h2, fill=h2_class)) + 
  		geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  		geom_errorbar(aes(ymin=h2_lower, ymax=h2_upper), width=.2,
                 position=position_dodge(.9))  +
  		figure_theme() + 
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8))
  	return(p)
}




results_dir <- args[1]
output_dir <- args[2]



trait_names <- c("blood_WHITE_COUNT", "disease_ALLERGY_ECZEMA_DIAGNOSED", "disease_ASTHMA_DIAGNOSED", "other_MORNINGPERSON", "disease_CARDIOVASCULAR", "blood_EOSINOPHIL_COUNT", "bp_DIASTOLICadjMEDz", "cov_EDU_COLLEGE", "blood_MEAN_CORPUSCULAR_HEMOGLOBIN", "body_BMIz", "body_BALDING1")
num_snps <- 9991229.0
num_genes <- 21297.0

df <- extract_partitioned_h2_results(trait_names, results_dir, num_snps, num_genes)


p <- make_standard_error_barplot_with_heritabilities(df)
output_file <- paste0(output_dir, "standard_error_barplot.pdf")
ggsave(p, file=output_file, width=7.2, height=4.0, units="in")

