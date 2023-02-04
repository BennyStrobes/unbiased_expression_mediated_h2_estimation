import numpy as np 
import os
import sys
import pdb
import gzip






input_file = sys.argv[1]
output_file = sys.argv[2]



head_count = 0
f = gzip.open(input_file)
t = open(output_file,'w')
for line in f:
	line = line.decode('utf-8').rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		valid_columns = []
		for i,val in enumerate(data):
			if val.endswith('MaxCPP') == False:
				valid_columns.append(i)
		valid_columns = np.asarray(valid_columns)
		t.write('\t'.join(np.asarray(data)[valid_columns]) + '\n')
		continue
	t.write('\t'.join(np.asarray(data)[valid_columns]) + '\n')
f.close()
t.close()