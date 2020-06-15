from utils import *
import sys,os
from scipy import sparse
import pandas as pd

data_file = sys.argv[1]
cell_file = sys.argv[2]
replicate_file = sys.argv[3]
out_dir = sys.argv[4]

data = pd.read_csv(data_file, delimiter = r"\s+", header=None)
mat_coo = sparse.coo_matrix((data[2], (data[0]-1, data[1])))
mat_coo = mat_coo.todense()

cell_list = np.genfromtxt(cell_file, dtype = 'str')
replicate_list = np.genfromtxt(cell_file, dtype = 'str')

contact_count = np.sum(mat_coo, axis=1)
target=int(np.ravel(np.median(contact_count, axis=0))[0])

outfile = '/net/noble/vol2/user/khj3017/cisTopic_dist_sparse/10/all_sparse_matrices_downsampled_10.txt'
out_cell = '/net/noble/vol5/user/khj3017/labels/all/all_downsampled.labeled'
out_replicate = '/net/noble/vol5/user/khj3017/labels/all/replicates/replicates_downsampled.labeled'

out = open(outfile, 'w')
out_cell = open(out_cell, 'w')
out_replicate = open(out_replicate, 'w')

counter = 0
for i in range(mat_coo.shape[0]):
	if np.ravel(contact_count[i])[0] >= target:
		row = mat_coo[i,:].getA()[0]
		vectMat = SubSampleVector(row, subSampleN = target)
		for j in range(len(vectMat)):
			if vectMat[j] > 0:
				print >> out, "%s\t%s\t%d" % (counter, j, vectMat[j])
		print >> out_cell, "%s\t%s" % (cell_list[i][0], cell_list[i][1])
		print >> out_replicate, "%s\t%s" % (replicate_list[i][0], replicate_list[i][1])
		counter += 1

out.close()
out_cell.close()
out_replicate.close()
