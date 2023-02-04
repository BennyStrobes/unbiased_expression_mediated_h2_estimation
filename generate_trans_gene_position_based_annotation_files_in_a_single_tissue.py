import numpy as np 
import os
import sys
import pdb
import gzip







def get_number_of_variants_per_chromosome(baseld_annotation_dir):
	var_per_chrom_arr =[]
	for chrom_num in range(1,23):
		file_name = baseld_annotation_dir + 'baselineLD.' + str(chrom_num) + '.annot.gz'
		f = gzip.open(file_name)
		counter = -1
		for line in f:
			counter = counter +1
		var_per_chrom_arr.append(counter)
		f.close()
	return np.asarray(var_per_chrom_arr)


def get_number_of_genes_per_chromosome(expressed_genes_file):
	arr = []
	for chrom_num in range(1,23):
		f = open(expressed_genes_file)
		head_count = 0
		chrom_gene_counter = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			if data[2] == 'chr' + str(chrom_num):
				chrom_gene_counter = chrom_gene_counter + 1
		f.close()
		arr.append(chrom_gene_counter)
	return np.asarray(arr)

def create_annotation_for_each_chromosome(number_of_genes_per_chromosome, number_of_variants_per_chromosome):
	arr = []
	for chrom_iter, chrom_num in enumerate(range(1,23)):
		anno = 0.0
		for chrom_iter2, chrom_num2 in enumerate(range(1,23)):
			if chrom_iter2 != chrom_iter:
				anno = anno + ((1.0/np.sum(np.delete(number_of_variants_per_chromosome,chrom_iter2)))*number_of_genes_per_chromosome[chrom_iter2])
		arr.append(anno)
	return np.asarray(arr)

def generate_annotation_file_for_this_chromosome(chrom_num, chrom_baseline_ld_annotation_file, chrom_annotation_output_file, tissue_name, snp_annotation):
	f = gzip.open(chrom_baseline_ld_annotation_file)
	t = open(chrom_annotation_output_file,'w')
	head_count = 0
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			t.write('\t'.join(np.asarray(data[:4])) + '\t' + tissue_name + '_trans_window' + '\n')
			continue
		t.write('\t'.join(np.asarray(data[:4])) + '\t' + str(snp_annotation) + '\n')
	f.close()
	t.close()
		

tissue_name = sys.argv[1]
expressed_genes_dir = sys.argv[2]
baseld_annotation_dir = sys.argv[3]
variant_annotation_dir = sys.argv[4]



# Extract df containing expressed genes in this tissue
expressed_genes_file = expressed_genes_dir + tissue_name + '_expressed_genes.txt'
expressed_genes_df = np.loadtxt(expressed_genes_file, dtype=str, delimiter='\t')[1:,:]

number_of_genes_per_chromosome = get_number_of_genes_per_chromosome(expressed_genes_file)

number_of_variants_per_chromosome = get_number_of_variants_per_chromosome(baseld_annotation_dir)

anno_per_chromosome = create_annotation_for_each_chromosome(number_of_genes_per_chromosome, number_of_variants_per_chromosome)


# Generate annotation file for this chromosome
for chrom_num in range(1,23):
	chrom_annotation_output_file = variant_annotation_dir + tissue_name + '_trans_eqtl_gene_position_anno.' + str(chrom_num) + '.annot'

	chrom_baseline_ld_annotation_file = baseld_annotation_dir + 'baselineLD.' + str(chrom_num) + '.annot.gz'

	generate_annotation_file_for_this_chromosome(chrom_num, chrom_baseline_ld_annotation_file, chrom_annotation_output_file, tissue_name, anno_per_chromosome[chrom_num-1])

