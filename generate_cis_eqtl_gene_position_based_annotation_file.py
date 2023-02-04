import numpy as np 
import os
import sys
import pdb
import gzip


def generate_gene_position_based_cis_eqtl_chrom_array(chrom_expressed_genes_df, distance_window):
	bp_distance_window = int(distance_window*1000000)

	chrom_arr = ['NULL']*300000000

	n_expressed_genes = chrom_expressed_genes_df.shape[0]

	for gene_iter in range(n_expressed_genes):
		ensamble_id = chrom_expressed_genes_df[gene_iter,0]
		start_pos = int(chrom_expressed_genes_df[gene_iter,3])
		end_pos = int(chrom_expressed_genes_df[gene_iter,4])

		if end_pos < start_pos:
			print('assumption eroror')
			pdb.set_trace()

		if end_pos-start_pos > bp_distance_window:
			print('big middle skip')
			continue

		true_start_pos = start_pos - bp_distance_window
		true_end_pos = end_pos + bp_distance_window
		if true_start_pos < 0:
			true_start_pos = 0


		for bp_pos in range(true_start_pos, true_end_pos):
			if chrom_arr[bp_pos] == 'NULL':
				chrom_arr[bp_pos] = ensamble_id
			else:
				chrom_arr[bp_pos] = chrom_arr[bp_pos] + ',' + ensamble_id
	return chrom_arr

def create_gene_to_nvar_mapping(chrom_arr, chrom_gene_names, chrom_baseline_ld_annotation_file, gene_names_output_file):
	# Initialize dictionary mapping
	dicti = {}
	for gene_name in chrom_gene_names:
		dicti[gene_name] = 0.0


	f = gzip.open(chrom_baseline_ld_annotation_file)
	head_count = 0
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		variant_pos = int(data[1])
		gene_string = chrom_arr[variant_pos]
		if gene_string == 'NULL':
			continue
		for gene_id in gene_string.split(','):
			dicti[gene_id] = dicti[gene_id] + 1
	f.close()

	t = open(gene_names_output_file,'w')
	t.write('gene_name\tnumber_of_cis_variants\n')
	for gene_name in chrom_gene_names:
		if dicti[gene_name] > 0.0:
			t.write(gene_name + '\t' + str(dicti[gene_name]) + '\n')
	t.close()


	return dicti

def create_new_annotation_file2(gene_to_nvar_mapping, chrom_expressed_genes_df, chrom_baseline_ld_annotation_file, chrom_annotation_output_file, tissue_name, distance_window):
	gene_starts = []
	gene_ends = []
	gene_weights = []
	
	bp_distance_window = int(distance_window*1000000)
	n_expressed_genes = chrom_expressed_genes_df.shape[0]

	for gene_iter in range(n_expressed_genes):
		ensamble_id = chrom_expressed_genes_df[gene_iter,0]
		start_pos = int(chrom_expressed_genes_df[gene_iter,3])
		end_pos = int(chrom_expressed_genes_df[gene_iter,4])

		true_start_pos = start_pos - bp_distance_window
		true_end_pos = end_pos + bp_distance_window

		if true_start_pos < 0:
			true_start_pos = 0
		if ensamble_id not in gene_to_nvar_mapping:
			continue
		if gene_to_nvar_mapping[ensamble_id] == 0.0:
			continue
		gene_starts.append(true_start_pos)
		gene_ends.append(true_end_pos)
		gene_weights.append(1.0/gene_to_nvar_mapping[ensamble_id])

	gene_starts = np.asarray(gene_starts)
	gene_ends = np.asarray(gene_ends)
	gene_weights = np.asarray(gene_weights)

	f = gzip.open(chrom_baseline_ld_annotation_file)
	t = open(chrom_annotation_output_file,'w')
	head_count = 0
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			t.write('\t'.join(np.asarray(data[:4])) + '\t' + tissue_name + '_cis_window' + '\n')
			continue

		variant_pos = int(data[1])
		gene_indices = (variant_pos >= gene_starts) & (variant_pos < gene_ends)


		t.write('\t'.join(np.asarray(data[:4])) + '\t')
		
		if np.sum(gene_indices) == 0:
			anno = 0.0
			t.write(str(anno) + '\n')
		else:
			anno = np.sum(gene_weights[gene_indices])
			t.write(str(anno) + '\n')
	f.close()
	f.close()

def create_new_annotation_file_w_intercept2(gene_to_nvar_mapping, chrom_expressed_genes_df, chrom_baseline_ld_annotation_file, chrom_annotation_output_file, tissue_name, distance_window):
	gene_starts = []
	gene_ends = []
	gene_weights = []
	
	bp_distance_window = int(distance_window*1000000)
	n_expressed_genes = chrom_expressed_genes_df.shape[0]

	for gene_iter in range(n_expressed_genes):
		ensamble_id = chrom_expressed_genes_df[gene_iter,0]
		start_pos = int(chrom_expressed_genes_df[gene_iter,3])
		end_pos = int(chrom_expressed_genes_df[gene_iter,4])

		true_start_pos = start_pos - bp_distance_window
		true_end_pos = end_pos + bp_distance_window

		if true_start_pos < 0:
			true_start_pos = 0
		if ensamble_id not in gene_to_nvar_mapping:
			continue
		if gene_to_nvar_mapping[ensamble_id] == 0.0:
			continue
		gene_starts.append(true_start_pos)
		gene_ends.append(true_end_pos)
		gene_weights.append(1.0/gene_to_nvar_mapping[ensamble_id])

	gene_starts = np.asarray(gene_starts)
	gene_ends = np.asarray(gene_ends)
	gene_weights = np.asarray(gene_weights)

	f = gzip.open(chrom_baseline_ld_annotation_file)
	t = open(chrom_annotation_output_file,'w')
	head_count = 0
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			t.write('\t'.join(np.asarray(data[:4])) + '\t' + 'base' + '\t' + tissue_name + '_cis_window' + '\n')
			continue

		variant_pos = int(data[1])
		gene_indices = (variant_pos >= gene_starts) & (variant_pos < gene_ends)


		t.write('\t'.join(np.asarray(data[:4])) + '\t' + '1.0' + '\t')
		
		if np.sum(gene_indices) == 0:
			anno = 0.0
			t.write(str(anno) + '\n')
		else:
			anno = np.sum(gene_weights[gene_indices])
			t.write(str(anno) + '\n')
	f.close()
	f.close()



def create_new_annotation_file(chrom_arr, gene_to_nvar_mapping, chrom_baseline_ld_annotation_file, chrom_annotation_output_file, tissue_name):
	f = gzip.open(chrom_baseline_ld_annotation_file)
	t = open(chrom_annotation_output_file,'w')
	head_count = 0
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			t.write('\t'.join(np.asarray(data[:4])) + '\t' + tissue_name + '_cis_window' + '\n')
			continue
		t.write('\t'.join(np.asarray(data[:4])) + '\t')
		
		variant_pos = int(data[1])
		gene_string = chrom_arr[variant_pos]
		variant_anno = 0.0
		if gene_string != 'NULL':
			for gene_id in gene_string.split(','):
				nvar_per_gene = gene_to_nvar_mapping[gene_id]
				variant_anno = variant_anno + (1.0/nvar_per_gene)
		t.write(str(variant_anno) + '\n')
	f.close()
	f.close()

def create_new_annotation_file_w_intercept(chrom_arr, gene_to_nvar_mapping, chrom_baseline_ld_annotation_file, chrom_annotation_output_file, tissue_name):
	f = gzip.open(chrom_baseline_ld_annotation_file)
	t = open(chrom_annotation_output_file,'w')
	head_count = 0
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			t.write('\t'.join(np.asarray(data[:4])) + '\t' + 'base\t' + tissue_name + '_cis_window' + '\n')
			continue
		t.write('\t'.join(np.asarray(data[:4])) + '\t')
		
		variant_pos = int(data[1])
		gene_string = chrom_arr[variant_pos]
		variant_anno = 0.0
		if gene_string != 'NULL':
			for gene_id in gene_string.split(','):
				nvar_per_gene = gene_to_nvar_mapping[gene_id]
				variant_anno = variant_anno + (1.0/nvar_per_gene)
		t.write('1.0\t' + str(variant_anno) + '\n')
	f.close()
	f.close()

def create_gene_to_nvar_mapping2(chrom_expressed_genes_df, chrom_gene_names, chrom_baseline_ld_annotation_file, gene_names_output_file, distance_window):
	# Create variant pos arr
	variant_pos_arr = []
	f = gzip.open(chrom_baseline_ld_annotation_file)
	head_count = 0
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		variant_pos = int(data[1])
		variant_pos_arr.append(variant_pos)
	f.close()
	variant_pos_arr = np.asarray(variant_pos_arr)

	
	bp_distance_window = int(distance_window*1000000)
	n_expressed_genes = chrom_expressed_genes_df.shape[0]

	gene_to_nvar_mapping = {}

	for gene_iter in range(n_expressed_genes):
		ensamble_id = chrom_expressed_genes_df[gene_iter,0]
		start_pos = int(chrom_expressed_genes_df[gene_iter,3])
		end_pos = int(chrom_expressed_genes_df[gene_iter,4])

		if end_pos-start_pos > bp_distance_window:
			print('big middle skip')
			continue

		true_start_pos = start_pos - bp_distance_window
		true_end_pos = end_pos + bp_distance_window
		if true_start_pos < 0:
			true_start_pos = 0

		cis_variant_indices = (variant_pos_arr >= true_start_pos) & (variant_pos_arr < true_end_pos)
		if ensamble_id in gene_to_nvar_mapping:
			print('asssumption eroeor')
			pdb.set_trace()
		gene_to_nvar_mapping[ensamble_id] = np.sum(cis_variant_indices)

	t = open(gene_names_output_file,'w')
	t.write('gene_name\tnumber_of_cis_variants\n')
	for gene_name in chrom_gene_names:
		if gene_name in gene_to_nvar_mapping:
			if gene_to_nvar_mapping[gene_name] > 0.0:
				t.write(gene_name + '\t' + str(gene_to_nvar_mapping[gene_name]) + '\n')
	t.close()

	return gene_to_nvar_mapping

def generate_annotation_file_for_this_chromosome(chrom_num, expressed_genes_df, chrom_baseline_ld_annotation_file, chrom_annotation_output_file, chrom_annotation_with_intercept_output_file, distance_window, gene_names_output_file, tissue_name):
	# Get expressed genes just on this chromosome
	chrom_indices = expressed_genes_df[:,2] == 'chr' + str(chrom_num)
	chrom_expressed_genes_df = expressed_genes_df[chrom_indices,:]


	# Quick error check
	if len(np.unique(chrom_expressed_genes_df[:,0])) != len(chrom_expressed_genes_df[:,0]):
		print('assumption eroror')
		pdb.set_trace()

	# Get array list of gene names on this chromosome
	chrom_gene_names = chrom_expressed_genes_df[:,0]

	# Extract number of 1KG variants mapped to each gene
	gene_to_nvar_mapping = create_gene_to_nvar_mapping2(chrom_expressed_genes_df, chrom_gene_names, chrom_baseline_ld_annotation_file, gene_names_output_file, distance_window)

	# Create new annotation file
	create_new_annotation_file2(gene_to_nvar_mapping, chrom_expressed_genes_df, chrom_baseline_ld_annotation_file, chrom_annotation_output_file, tissue_name, distance_window)
	create_new_annotation_file_w_intercept2(gene_to_nvar_mapping, chrom_expressed_genes_df, chrom_baseline_ld_annotation_file, chrom_annotation_with_intercept_output_file, tissue_name, distance_window)



tissue_name = sys.argv[1]
expressed_genes_dir = sys.argv[2]
baseld_annotation_dir = sys.argv[3]
variant_annotation_dir = sys.argv[4]
distance_window = sys.argv[5]



# Extract df containing expressed genes in this tissue
expressed_genes_file = expressed_genes_dir + tissue_name + '_expressed_genes.txt'
expressed_genes_df = np.loadtxt(expressed_genes_file, dtype=str, delimiter='\t')[1:,:]

# Generate annotation file for this chromosome
for chrom_num in range(1,23):
	chrom_annotation_output_file = variant_annotation_dir + tissue_name + '_cis_eqtl_gene_position_anno_' + distance_window + '_mb.' + str(chrom_num) + '.annot'
	chrom_annotation_with_intercept_output_file = variant_annotation_dir + tissue_name + '_cis_eqtl_gene_position_anno_' + distance_window + '_mb_w_intercept.' + str(chrom_num) + '.annot'

	gene_names_output_file = variant_annotation_dir + tissue_name + '_cis_eqtl_gene_position_anno_' + distance_window + '_mb.' + str(chrom_num) + '_gene_names.txt'
	chrom_baseline_ld_annotation_file = baseld_annotation_dir + 'baselineLD.' + str(chrom_num) + '.annot.gz'

	generate_annotation_file_for_this_chromosome(chrom_num, expressed_genes_df, chrom_baseline_ld_annotation_file, chrom_annotation_output_file, chrom_annotation_with_intercept_output_file, float(distance_window), gene_names_output_file, tissue_name)
