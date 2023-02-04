import numpy as np 
import os
import sys
import pdb
import gzip




def get_tissue_names_arr(tissue_names_file):
	f = open(tissue_names_file)
	head_count = 0
	arr = []

	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		arr.append(data[0])
	f.close()
	return np.asarray(arr)


def extract_and_print_list_of_expressed_genes_in_a_tissue(tissue_egene_file, tissue_output_file):
	# Open output file handle
	t = open(tissue_output_file, 'w')
	t.write('ensamble_id\tgene_id\tchrom_num\tgene_start\tgene_end\tstrand\n')

	# Open and stream egene file
	f = gzip.open(tissue_egene_file)
	head_count = 0
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split('\t')
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			num_fields = len(data)
			continue
		# Extract relevent fields
		ensamble_id = data[0]
		gene_id = data[1]
		chrom_num = data[2]
		start = int(data[3])
		end = int(data[4])
		strand = data[5]

		# Remove non-autosomal
		if chrom_num == 'chrX' or chrom_num == 'chrY':
			continue

		# Quick error check
		if end < start:
			print('assumption erroror')
			pdb.set_trace()
		if strand != '+' and strand != '-':
			print('assumption eroroorororor')
			pdb.set_trace()

		# Print to output
		t.write(ensamble_id + '\t' + gene_id + '\t' + chrom_num + '\t' + str(start) + '\t' + str(end) + '\t' + strand + '\n')
	f.close()
	t.close()

def get_expressed_genes_dicti(file_name):
	dicti = {}
	f = open(file_name)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count +1
			continue
		if data[0] in dicti:
			print('aassumption oeror')
			pdb.set_trace()
		dicti[data[0]] = 1
	f.close()
	return dicti

def check_overlap(t1, t2, expressed_genes_dir):
	dicti1 = get_expressed_genes_dicti(expressed_genes_dir + t1 + '_expressed_genes.txt')
	dicti2 = get_expressed_genes_dicti(expressed_genes_dir + t2 + '_expressed_genes.txt')
	overlap_dicti = {}
	for gene_id in [*dicti1]:
		if gene_id in dicti1 and gene_id in dicti2:
			overlap_dicti[gene_id] = 1

	pdb.set_trace()



tissue_names_file = sys.argv[1]
gtex_egenes_dir = sys.argv[2]
expressed_genes_dir = sys.argv[3]


# Extract list of tissues
tissue_names = get_tissue_names_arr(tissue_names_file)

# Loop through tissues
for tissue_name in tissue_names:

	# Output file name for this tissue
	tissue_output_file = expressed_genes_dir + tissue_name + '_expressed_genes.txt'

	# Tissue egene file
	tissue_egene_file = gtex_egenes_dir + tissue_name + '.v8.egenes.txt.gz'

	# Extract and print list of expressed genes for this tissue
	extract_and_print_list_of_expressed_genes_in_a_tissue(tissue_egene_file, tissue_output_file)
