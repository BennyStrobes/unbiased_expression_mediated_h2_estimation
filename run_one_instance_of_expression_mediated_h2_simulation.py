import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import numpy as np 
import os
import sys
import pdb
import sklearn.linear_model
import statsmodels.api as sm
import scipy.stats

def ld_score_regression_with_no_ld(gwas_chi_squared_statistics, num_snps, gwas_sample_size):
	h2 = np.mean((gwas_chi_squared_statistics - 1.0))*num_snps/gwas_sample_size
	return h2

def ld_score_regression_with_no_ld_and_standard_error(gwas_chi_squared_statistics, num_snps, gwas_sample_size):
	per_snp_values = (gwas_chi_squared_statistics - 1.0)*num_snps/gwas_sample_size
	h2 = np.mean(per_snp_values)
	h2_sem = scipy.stats.sem(per_snp_values)
	return h2, h2_sem



def compute_ld_scores_for_genotype_and_expression_with_no_ld_and_point_estimate_eqtl_effects(gene_sim_info_arr):
	genotype_ld_scores = []
	expr_ld_scores = []
	n_genes = len(gene_sim_info_arr)
	for gene_sim_info in gene_sim_info_arr:
		if np.isnan(gene_sim_info[0][0]):
			nvar = gene_sim_info[1]
			genotype_ld_scores.append(np.ones(nvar))
			expr_ld_scores.append(np.zeros(nvar))
		else:
			true_eqtl_effects = gene_sim_info[1]
			eqtl_geno = gene_sim_info[2]
			eqtl_expr = gene_sim_info[3]

			# Ridge regression
			clf = sklearn.linear_model.RidgeCV(alphas=[1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6]).fit(eqtl_geno, eqtl_expr)
			clf = sklearn.linear_model.LassoCV(cv=5, random_state=0).fit(eqtl_geno, eqtl_expr)
			eqtl_effects = clf.coef_

			# Standardize estimated effects
			eqtl_effects_var = np.dot(np.dot(eqtl_effects, np.eye(len(eqtl_effects))), eqtl_effects)
			if eqtl_effects_var > 0.0:
				standardized_eqtl_effects = eqtl_effects/np.sqrt(eqtl_effects_var)
			else:
				standardized_eqtl_effects = eqtl_effects

			genotype_ld_scores.append(np.ones(len(standardized_eqtl_effects)))
			expr_ld_scores.append(np.square(standardized_eqtl_effects))
	return np.hstack(genotype_ld_scores), np.hstack(expr_ld_scores)

def compute_ld_scores_for_genotype_and_expression_with_no_ld_and_true_eqtl_effects(gene_sim_info_arr):
	genotype_ld_scores = []
	expr_ld_scores = []
	n_genes = len(gene_sim_info_arr)
	n_standardized_genes = 0
	for gene_sim_info in gene_sim_info_arr:
		if np.isnan(gene_sim_info[0][0]):
			nvar = gene_sim_info[1]
			genotype_ld_scores.append(np.ones(nvar))
			expr_ld_scores.append(np.zeros(nvar))
		else:
			eqtl_effects = gene_sim_info[1]
			eqtl_effects_var = np.dot(np.dot(eqtl_effects, np.eye(len(eqtl_effects))), eqtl_effects)
			if eqtl_effects_var > 0.0:
				standardized_eqtl_effects = eqtl_effects/np.sqrt(eqtl_effects_var)
				n_standardized_genes = n_standardized_genes + 1
			else:
				standardized_eqtl_effects = eqtl_effects
			genotype_ld_scores.append(np.ones(len(standardized_eqtl_effects)))
			expr_ld_scores.append(np.square(standardized_eqtl_effects))
	return np.hstack(genotype_ld_scores), np.hstack(expr_ld_scores), n_standardized_genes


def run_simulation_in_single_gene_block(N_expr_samples, N_gwas_samples, N_var_per_gene_block, h2_expr, num_causal_eqtl_var_per_block):
	#####################################
	# Expression simulation
	#####################################
	# Simulate variant allele frequencies 
	afs = np.random.uniform(.05,.5,size=N_var_per_gene_block)
	eqtl_geno = np.zeros((N_expr_samples, N_var_per_gene_block))
	for variant_index in range(N_var_per_gene_block):
		af = afs[variant_index]
		a1 = np.random.choice([0,1], p=[1-af, af], size=N_expr_samples)
		a2 = np.random.choice([0,1], p=[1-af, af], size=N_expr_samples)
		g = a1 + a2
		standardized_g = (g - np.mean(g))/np.std(g)
		eqtl_geno[:, variant_index] = standardized_g
	
	eqtl_beta_sim = np.zeros(N_var_per_gene_block)
	causal_eqtl_variant_indices = np.random.choice(range(N_var_per_gene_block), num_causal_eqtl_var_per_block, replace=False)
	for causal_variant_index in causal_eqtl_variant_indices:
		eqtl_beta_sim[causal_variant_index] = np.random.normal(loc=0, scale = np.sqrt(h2_expr/num_causal_eqtl_var_per_block))

	eqtl_mean = np.dot(eqtl_geno, eqtl_beta_sim)
	E_sim = np.random.normal(loc=eqtl_mean, scale=np.sqrt(1.0-h2_expr))

	E_sim = (E_sim - np.mean(E_sim))/np.std(E_sim)

	return afs, eqtl_beta_sim, eqtl_geno, E_sim

def simulate_gwas_data_with_genotype_and_expression_effects_v3(N_expr_samples, N_gwas_samples, N_gene_blocks, N_null_gene_blocks, N_no_gene_blocks, N_var_per_gene_block, h2_expr, h2_e_trait, h2_g_trait, num_causal_eqtl_var_per_block, fraction_causal_var_gwas):
	## Simulate 
	mean_trait_value = np.zeros(N_gwas_samples)
	gene_sim_info_arr = []
	gwas_variant_arr = []
	for gene_iter in range(N_gene_blocks):
		# Simulate gene expression
		gene_sim_info = run_simulation_in_single_gene_block(N_expr_samples, N_gwas_samples, N_var_per_gene_block, h2_expr, num_causal_eqtl_var_per_block)
		gene_sim_info_arr.append(gene_sim_info)


		# Update mean of gwas
		afs = gene_sim_info[0]
		eqtl_beta_sim = gene_sim_info[1]
		gwas_geno = np.zeros((N_gwas_samples, N_var_per_gene_block))
		for variant_index in range(N_var_per_gene_block):
			af = afs[variant_index]
			a1 = np.random.choice([0,1], p=[1-af, af], size=N_gwas_samples)
			a2 = np.random.choice([0,1], p=[1-af, af], size=N_gwas_samples)
			g = a1 + a2
			standardized_g = (g - np.mean(g))/np.std(g)
			gwas_geno[:, variant_index] = standardized_g
	
		gwas_variant_arr.append(gwas_geno)

		gwas_beta_sim = np.zeros(N_var_per_gene_block)
		causal_gwas_variant_indices = np.random.choice(range(N_var_per_gene_block), int(np.round(fraction_causal_var_gwas*N_var_per_gene_block)), replace=False)
		for causal_variant_index in causal_gwas_variant_indices:
			gwas_beta_sim[causal_variant_index] = np.random.normal(loc=0, scale = np.sqrt(h2_g_trait/(fraction_causal_var_gwas*N_var_per_gene_block*(N_gene_blocks + N_null_gene_blocks + N_no_gene_blocks))))
		mean_trait_value = mean_trait_value + np.dot(gwas_geno, gwas_beta_sim)

		gwas_expression = np.dot(gwas_geno, eqtl_beta_sim)
		gene_trait_effect = np.random.normal(loc=0, scale=np.sqrt(h2_e_trait/N_gene_blocks))
		mean_trait_value = mean_trait_value + gene_trait_effect*(gwas_expression/np.std(gwas_expression))

		ld = np.corrcoef(np.transpose(gwas_geno))
		#gene_var = np.dot(np.dot(eqtl_beta_sim, np.eye(len(eqtl_beta_sim))), eqtl_beta_sim)
	
	for gene_iter in range(N_no_gene_blocks):
		# Simulate gene expression
		gene_sim_info = run_simulation_in_single_gene_block(N_expr_samples, N_gwas_samples, N_var_per_gene_block, h2_expr, num_causal_eqtl_var_per_block)
		gene_sim_info_arr.append(([np.nan], N_var_per_gene_block, np.nan, np.nan))


		# Update mean of gwas
		afs = gene_sim_info[0]
		eqtl_beta_sim = gene_sim_info[1]
		gwas_geno = np.zeros((N_gwas_samples, N_var_per_gene_block))
		for variant_index in range(N_var_per_gene_block):
			af = afs[variant_index]
			a1 = np.random.choice([0,1], p=[1-af, af], size=N_gwas_samples)
			a2 = np.random.choice([0,1], p=[1-af, af], size=N_gwas_samples)
			g = a1 + a2
			standardized_g = (g - np.mean(g))/np.std(g)
			gwas_geno[:, variant_index] = standardized_g
	
		gwas_variant_arr.append(gwas_geno)

		gwas_beta_sim = np.zeros(N_var_per_gene_block)
		causal_gwas_variant_indices = np.random.choice(range(N_var_per_gene_block), int(np.round(fraction_causal_var_gwas*N_var_per_gene_block)), replace=False)
		for causal_variant_index in causal_gwas_variant_indices:
			gwas_beta_sim[causal_variant_index] = np.random.normal(loc=0, scale = np.sqrt(h2_g_trait/(fraction_causal_var_gwas*N_var_per_gene_block*(N_gene_blocks + N_null_gene_blocks + N_no_gene_blocks))))
		mean_trait_value = mean_trait_value + np.dot(gwas_geno, gwas_beta_sim)


		ld = np.corrcoef(np.transpose(gwas_geno))
	
	for gene_iter in range(N_null_gene_blocks):
		# Simulate gene expression
		gene_sim_info_t = run_simulation_in_single_gene_block(N_expr_samples, N_gwas_samples, N_var_per_gene_block, h2_expr, num_causal_eqtl_var_per_block)
		#gene_sim_info = gene_sim_info_t[0], gene_sim_info_t[1]*0.0, gene_sim_info_t[2], gene_sim_info_t[3]
		gene_sim_info_arr.append(gene_sim_info_t)

		# Update mean of gwas
		afs = gene_sim_info[0]
		eqtl_beta_sim = gene_sim_info[1]
		gwas_geno = np.zeros((N_gwas_samples, N_var_per_gene_block))
		for variant_index in range(N_var_per_gene_block):
			af = afs[variant_index]
			a1 = np.random.choice([0,1], p=[1-af, af], size=N_gwas_samples)
			a2 = np.random.choice([0,1], p=[1-af, af], size=N_gwas_samples)
			g = a1 + a2
			standardized_g = (g - np.mean(g))/np.std(g)
			gwas_geno[:, variant_index] = standardized_g
	
		gwas_variant_arr.append(gwas_geno)

		gwas_beta_sim = np.zeros(N_var_per_gene_block)
		causal_gwas_variant_indices = np.random.choice(range(N_var_per_gene_block), int(np.round(fraction_causal_var_gwas*N_var_per_gene_block)), replace=False)
		for causal_variant_index in causal_gwas_variant_indices:
			gwas_beta_sim[causal_variant_index] = np.random.normal(loc=0, scale = np.sqrt(h2_g_trait/(fraction_causal_var_gwas*N_var_per_gene_block*(N_gene_blocks + N_null_gene_blocks + N_no_gene_blocks))))
		mean_trait_value = mean_trait_value + np.dot(gwas_geno, gwas_beta_sim)


		ld = np.corrcoef(np.transpose(gwas_geno))
	


	Y=np.random.normal(loc=mean_trait_value, scale=np.sqrt(1.0-h2_g_trait-h2_e_trait))
	gwas_z_scores = []
	for gwas_geno in gwas_variant_arr:
		nvar = gwas_geno.shape[1]
		for var_index in range(nvar):
			result = sm.OLS(Y, sm.add_constant(gwas_geno[:,var_index])).fit()
			z_score = result.tvalues[1]
			gwas_z_scores.append(z_score)
	gwas_z_scores = np.asarray(gwas_z_scores)


	return gene_sim_info_arr, gwas_z_scores, np.var(Y)

def h2_estimation_using_marginalized_approach(gene_sim_info_arr, gwas_z_scores, N_gwas_samples, N_gene_blocks, N_null_gene_blocks):
	geno_ld_scores, expr_ld_scores = compute_ld_scores_for_genotype_and_expression_with_no_ld_marginal_effects(gene_sim_info_arr)
	joint_ld_scores = np.hstack([geno_ld_scores.reshape(len(geno_ld_scores), 1), expr_ld_scores.reshape(len(expr_ld_scores), 1)])
	result = sm.OLS(np.square(gwas_z_scores), joint_ld_scores).fit()
	est_h2_g = (result.params[0] - 1.0)*len(gwas_z_scores)/N_gwas_samples
	est_h2_e = result.params[1]*(N_gene_blocks + N_null_gene_blocks)/N_gwas_samples
	return est_h2_g, est_h2_e

def compute_ld_scores_for_genotype_and_expression_with_no_ld_marginal_effects(gene_sim_info_arr):
	genotype_ld_scores = []
	expr_ld_scores = []
	n_genes = len(gene_sim_info_arr)
	for gene_sim_info in gene_sim_info_arr:
		if np.isnan(gene_sim_info[0][0]):
			nvar = gene_sim_info[1]
			genotype_ld_scores.append(np.ones(nvar))
			expr_ld_scores.append(np.zeros(nvar))
		else:
			nvar = len(gene_sim_info[1])
			genotype_ld_scores.append(np.ones(nvar))
			expr_ld_scores.append(np.ones(nvar)/nvar)
	return np.hstack(genotype_ld_scores), np.hstack(expr_ld_scores)

def h2_estimation_using_standard_approach_with_known_true_eqtl_effects(gene_sim_info_arr, gwas_z_scores, N_gwas_samples, N_gene_blocks, N_null_gene_blocks):
	geno_ld_scores, expr_ld_scores, n_standardized_genes = compute_ld_scores_for_genotype_and_expression_with_no_ld_and_true_eqtl_effects(gene_sim_info_arr)
	joint_ld_scores = np.hstack([geno_ld_scores.reshape(len(geno_ld_scores), 1), expr_ld_scores.reshape(len(expr_ld_scores), 1)])
	result = sm.OLS(np.square(gwas_z_scores), joint_ld_scores).fit()
	est_h2_g = (result.params[0] - 1.0)*len(gwas_z_scores)/N_gwas_samples
	est_h2_e = result.params[1]*(n_standardized_genes)/N_gwas_samples
	return est_h2_g, est_h2_e

def compute_ld_scores_for_genotype_and_expression_with_no_ld_and_point_estimate_eqtl_effects(gene_sim_info_arr):
    genotype_ld_scores = []
    expr_ld_scores = []
    n_genes = len(gene_sim_info_arr)
    n_standardized_genes = 0
    for gene_sim_info in gene_sim_info_arr:
        if np.isnan(gene_sim_info[0][0]):
            nvar = gene_sim_info[1]
            genotype_ld_scores.append(np.ones(nvar))
            expr_ld_scores.append(np.zeros(nvar))
        else:
            true_eqtl_effects = gene_sim_info[1]
            eqtl_geno = gene_sim_info[2]
            eqtl_expr = gene_sim_info[3]

            # Ridge regression
            clf = sklearn.linear_model.RidgeCV(alphas=[1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6]).fit(eqtl_geno, eqtl_expr)
            clf = sklearn.linear_model.LassoCV(cv=5, random_state=0).fit(eqtl_geno, eqtl_expr)
            eqtl_effects = clf.coef_

            eqtl_effects_var = np.dot(np.dot(eqtl_effects, np.eye(len(eqtl_effects))), eqtl_effects)
            if eqtl_effects_var > 0.0:
                standardized_eqtl_effects = eqtl_effects/np.sqrt(eqtl_effects_var)
                n_standardized_genes = n_standardized_genes + 1
            else:
                standardized_eqtl_effects = eqtl_effects

            genotype_ld_scores.append(np.ones(len(standardized_eqtl_effects)))
            expr_ld_scores.append(np.square(standardized_eqtl_effects))
    return np.hstack(genotype_ld_scores), np.hstack(expr_ld_scores), n_standardized_genes


def h2_estimation_using_standard_approach_with_point_estimated_eqtl_effects(gene_sim_info_arr, gwas_z_scores, N_gwas_samples, N_gene_blocks, N_null_gene_blocks):
	geno_ld_scores, expr_ld_scores, n_standardized_genes = compute_ld_scores_for_genotype_and_expression_with_no_ld_and_point_estimate_eqtl_effects(gene_sim_info_arr)
	joint_ld_scores = np.hstack([geno_ld_scores.reshape(len(geno_ld_scores), 1), expr_ld_scores.reshape(len(expr_ld_scores), 1)])
	result = sm.OLS(np.square(gwas_z_scores), joint_ld_scores).fit()
	est_h2_g = (result.params[0] - 1.0)*len(gwas_z_scores)/N_gwas_samples
	est_h2_e = result.params[1]*(n_standardized_genes)/N_gwas_samples
	print(n_standardized_genes)
	return est_h2_g, est_h2_e


#################################
# Command line args
#################################
simulation_results_dir = sys.argv[1]
expression_data_sample_size = int(sys.argv[2])



# Simulation parameters
N_simulations = 100
N_expr_samples = int(np.copy(expression_data_sample_size)*1.0)
N_gwas_samples=20000
N_gene_blocks = 40
N_null_gene_blocks = 20
N_no_gene_blocks = 40
N_var_per_block=200
h2_expr = .1
h2_e_trait = .25
h2_g_trait = .25
num_causal_eqtl_var_per_block = 5
fraction_causal_var_gwas = .5


# Initialize output file handle
output_file = simulation_results_dir + 'simulation_N_expr_' + str(N_expr_samples) + '_results.txt'
t = open(output_file,'w')
t.write('method_name\ttrue_h2_g\ttrue_h2_e\test_h2_g\test_h2_e\tsimulation_iter\n')



for itera in range(N_simulations):

	# Simulate genome-wide data
	gene_sim_info_arr, gwas_z_scores, variance_of_Y = simulate_gwas_data_with_genotype_and_expression_effects_v3(N_expr_samples, N_gwas_samples, N_gene_blocks, N_null_gene_blocks, N_no_gene_blocks, N_var_per_block, h2_expr, h2_e_trait, h2_g_trait, num_causal_eqtl_var_per_block, fraction_causal_var_gwas)

	# h2 estimation using true eqtl effects
	true_eqtl_est_h2_g, true_eqtl_est_h2_e = h2_estimation_using_standard_approach_with_known_true_eqtl_effects(gene_sim_info_arr, gwas_z_scores, N_gwas_samples, N_gene_blocks, N_null_gene_blocks)
	t.write('known_eqtl_effects_model\t' + str(h2_g_trait) + '\t' + str(h2_e_trait) + '\t' + str(true_eqtl_est_h2_g) + '\t' + str(true_eqtl_est_h2_e) + '\t' + str(itera) + '\n')

	# h2 estimation using point_estimated eqtl effects
	point_eqtl_est_h2_g, point_eqtl_est_h2_e = h2_estimation_using_standard_approach_with_point_estimated_eqtl_effects(gene_sim_info_arr, gwas_z_scores, N_gwas_samples, N_gene_blocks, N_null_gene_blocks)
	t.write('point_estimate_eqtl_effects_model\t' + str(h2_g_trait) + '\t' + str(h2_e_trait) + '\t' + str(point_eqtl_est_h2_g) + '\t' + str(point_eqtl_est_h2_e) + '\t' + str(itera) + '\n')

	# h2 estimation using marginalized approach
	marginal_est_h2_g, marginal_est_h2_e = h2_estimation_using_marginalized_approach(gene_sim_info_arr, gwas_z_scores, N_gwas_samples, N_gene_blocks, N_null_gene_blocks)
	t.write('marginal_eqtl_effects_model\t' + str(h2_g_trait) + '\t' + str(h2_e_trait) + '\t' + str(marginal_est_h2_g) + '\t' + str(marginal_est_h2_e) + '\t' + str(itera) + '\n')
	t.flush()

t.close()

