import sys,os

inputfolder=sys.argv[1]
mapbed=sys.argv[2]
resolution=int(sys.argv[3])
outputdir = sys.argv[4]
mapbedfh = open(mapbed)

bintocoord = {}
for line in mapbedfh :
	tokens = line.split()
	bintocoord[tokens[2]] = (tokens[0],str(int(tokens[1])+resolution/2))


list = os.popen('ls '+ inputfolder + '/human*matrix').read()
filenames = list.split('\n')[:-1]

for matrixfilename in filenames :
	if matrixfilename != '' :
		matrixfilefh = open(matrixfilename)
		outfilefh = open(os.path.join(outputdir, os.path.split(matrixfilename)[1][:-7]+'.int.bed'),'w')
		print matrixfilename
		for line in matrixfilefh :
			tokens = line.split()
			firstcoor = bintocoord[tokens[0]]; secondcoor = bintocoord[tokens[1]] 
			outtokens = (firstcoor[0], firstcoor[1], secondcoor[0], secondcoor[1], tokens[2], tokens[3] )
			outline = '\t'.join(outtokens)
			print >> outfilefh, outline
		outfilefh.close()


