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


make_boxplot_showing_fraction_of_heritability_mediated_by_gene_expression <- function(simulation_results_dir, eqtl_sample_sizes) {
	method_arr <- c()
	fraction_mediated_arr <- c()
	sample_size_arr <- c()

	true_mediation <- 0

	for (sample_size_iter in 1:length(eqtl_sample_sizes)) {
		sample_size <- eqtl_sample_sizes[sample_size_iter]

		file_name <- paste0(simulation_results_dir, "simulation_N_expr_", sample_size, "_results.txt")
		df <- read.table(file_name, header=TRUE, sep="\t")
		fraction_med = df$est_h2_e/(df$est_h2_e + df$est_h2_g)

		method_arr <- c(method_arr, as.character(df$method_name))
		fraction_mediated_arr <- c(fraction_mediated_arr, fraction_med)
		sample_size_arr <- c(sample_size_arr, rep(sample_size, length(fraction_med)))

		true_mediation <- (df$true_h2_e/(df$true_h2_e+df$true_h2_g))[1]
	}

	df2 <- data.frame(method=method_arr, fraction_h2_mediated=fraction_mediated_arr, sample_size=factor(sample_size_arr, levels=eqtl_sample_sizes))

	p<-ggplot(df2, aes(x=sample_size, y=fraction_h2_mediated, color=method)) +
  		geom_boxplot() +
  		figure_theme() +
  		geom_hline(yintercept=true_mediation,linetype="dashed") +
  		labs(x="eQTL sample size", y="expression-mediated-h2 / h2", color="") +
  		theme(legend.position="bottom")
  	return(p)
}





simulation_results_dir <- args[1]
visualize_simulation_dir <- args[2]



eqtl_sample_sizes <- c(200, 400, 600, 800, 1000)



# Create boxplot demonstrating fraction of heritability mediated by gene expression
boxplot <- make_boxplot_showing_fraction_of_heritability_mediated_by_gene_expression(simulation_results_dir, eqtl_sample_sizes)
output_file <- paste0(visualize_simulation_dir, "fraction_mediated_boxplot.pdf")
ggsave(boxplot, file=output_file, width=7.2, height=4.0, units="in")
