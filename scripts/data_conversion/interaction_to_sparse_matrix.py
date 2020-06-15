from utils import *
import os,sys

def main(argv):
	fileDir = sys.argv[1]
	outDir = sys.argv[2]
	libname = sys.argv[3]
	resolution = int(sys.argv[4])
	dist = int(sys.argv[5]) # number of bins
	
	if not os.path.exists(outDir):
		os.makedirs(outDir)
	
	assembly = 'hg19'
	chr_list = getChrs(assembly)[:-2]	
	valid_LPs = []

	for this_chr in chr_list:
		chrs, midPoints = generateChrMidpoints(chrs = [this_chr], assembly = assembly, resolution = resolution)
		offset = 0
		for mid_point in midPoints:
			for k in range(dist):
				if (offset+k) < len(midPoints):
					valid_LPs.append([this_chr, mid_point, this_chr, midPoints[k+offset]])
			offset += 1

	LP_list_saving = [LP[0] + ':' + str(LP[1]) + '-' + str(LP[3]) for LP in valid_LPs]
	set_lp_list = set(LP_list_saving)
	
	out_LP_name = os.path.join(outDir, libname + "_" + str(dist) + "_" + str(resolution) + "_LPnames.txt")
	np.savetxt(out_LP_name, LP_list_saving, fmt='%s', delimiter="\t", newline="\n")		


	list_files = os.popen('ls ' + fileDir + '/*int.bed').read()
	filenames = list_files.split('\n')
	filenames = filenames[:-1]
	num_cells = len(filenames)
	
	cell_number = 0
	for filename in filenames:
		mat_file = open(filename)
		cell_save = []
		print filename
		for line in mat_file:
        		temp = line.strip()
        		target = line.split()
        		if target[0] == target[2] :
            			lp = [target[0], int(target[1]), target[2], int(target[3])]
            			search_lp = target[0] + ":" + target[1] + "-" + target[3]
				if search_lp in set_lp_list:
                			cell_save.append([cell_number, LP_list_saving.index(search_lp), int(target[4])])
		mat_file.close()
    		cell_save = np.asarray(cell_save)
    		cell_save = cell_save[cell_save[:,1].argsort()]
   	
		out_mat_name = os.path.join(fileDir, os.path.split(filename)[1][:-8] + ".sparse.matrix_" + str(dist))
 		np.savetxt(out_mat_name, cell_save, fmt='%i', delimiter="\t", newline="\n")
		
    		cell_number += 1


if __name__ == "__main__":
        main(sys.argv)
