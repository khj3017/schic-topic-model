import numpy as np
import time
import os
import sys
import random
import gzip
from decimal import *
from scipy.stats import spearmanr
import matplotlib.cm
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import transforms
from iced import normalization
from scipy import ndimage


def getChrLen ( chr = 'chr1', assembly = 'hg19' ):
	chrsizes = np.loadtxt('/net/trapnell/vol1/home/khj3017/proj/refdata/' + assembly + '.genomesize', dtype={'names': ('chr','len'),'formats':('S32','i4')} )
	for i in range(len(chrsizes)):
		if ( chrsizes[i][0] == chr ):
			return (chrsizes[i][1] )

def getChrs ( assembly  ):
	chrsizes = np.loadtxt('/net/trapnell/vol1/home/khj3017/proj/refdata/' + assembly + '.genomesize', dtype={'names': ('chr','len'),'formats':('S32','i4')} )
	chrs = []
	for i in range(len(chrsizes)):
		chrs.append( chrsizes[i][0] )
	return chrs 

#Zero based indexing for chr coordinates
def generateChrBins ( assembly = 'hg19', chrs = [], resolution = 5 * 10 ** 5, returnOffset = False, returnLabels = False, overlap = False ) :
	if chrs == [] :
		chrs = getChrs(assembly)

	if overlap :
		overlapOffset = 0
	else :
		overlapOffset = 1

	chrList = []; startList =[]; endList = []; offsets = [0]
	for chr in chrs :
		thisChr = chr
		chrLen = getChrLen( thisChr, assembly = assembly )
		if chr != chrs[len(chrs)-1] :
			offsets.append(offsets[len(offsets)-1] + np.ceil(chrLen/float(resolution)) )
		for i in range(chrLen/resolution) :
			chrList.append(thisChr)
			startList.append(i*resolution+overlapOffset)
			endList.append((i+1)*resolution)
		if (i+1) * resolution < chrLen :
			chrList.append(thisChr)
			startList.append((i+1)*resolution+overlapOffset)
			endList.append(chrLen)
	offsetDict = {}
	if returnOffset :
		for chr,offset in zip(chrs,offsets) :
			offsetDict[chr] = offset
			#print chr +" " + str(offset)
		

	if returnOffset :
		return [chrList,startList,endList,offsetDict]					
	if returnLabels :
		binStr = [ assembly+'|'+chrList[i]+':'+str(startList[i])+'-'+str(endList[i]) for i in range(len(chrList)) ]
		return binStr

	return [chrList,startList,endList]					

def generateChrLogBins ( assembly = 'hg19', chrs = [], resolution = 5 * 10 ** 5, returnOffset = False, returnLabels = False, overlap = False ) :
        if chrs == [] :
                chrs = getChrs(assembly)

        if overlap :
                overlapOffset = 0
        else :
                overlapOffset = 1

        chrList = []; startList =[]; endList = []; offsets = [0]
        bin_sizes = [0, 5000, 10000, 20000, 40000, 80000, 160000, 320000, 640000, 1280000, 2560000, 3000000]
	for chr in chrs :
                thisChr = chr
                chrLen = getChrLen( thisChr, assembly = assembly )
                if chr != chrs[len(chrs)-1] :
                        offsets.append(offsets[len(offsets)-1] + np.ceil(chrLen/float(resolution)) + 5 )
                for i in range(11) : 
			chrList.append(thisChr)
                        startList.append(bin_sizes[i]+overlapOffset)
                        endList.append(bin_sizes[i+1])
		for i in range(chrLen/resolution-6) :
                        chrList.append(thisChr)
                        if i == 0 : 
                        	startList.append(bin_sizes[-1]+overlapOffset)
                        	endList.append((i+7)*resolution)
			else :
				startList.append((i+6)*resolution+overlapOffset)
                        	endList.append((i+7)*resolution)
                if (i+1) * resolution < chrLen :
                        chrList.append(thisChr)
                        startList.append((i+7)*resolution+overlapOffset)
                        endList.append(chrLen)
        offsetDict = {}
        if returnOffset :
                for chr,offset in zip(chrs,offsets) :
                        offsetDict[chr] = offset
                        #print chr +" " + str(offset)


        if returnOffset :
                return [chrList,startList,endList,offsetDict]
        if returnLabels :
                binStr = [ assembly+'|'+chrList[i]+':'+str(startList[i])+'-'+str(endList[i]) for i in range(len(chrList)) ]
                return binStr

        return [chrList,startList,endList]

def generateChrLogBinsNew ( assembly = 'hg19', chrs = [], interval = 0.125, returnOffset = False, returnLabels = False, overlap = False ) :
        if chrs == [] :
                chrs = getChrs(assembly)

        if overlap :
                overlapOffset = 0
        else :
                overlapOffset = 1
	
	start = 10
        chrList = []; startList =[]; endList = []; offsets = [0]
        for chr in chrs :
                thisChr = chr
                chrLen = getChrLen( thisChr, assembly = assembly )
                end = np.log2(chrLen)
                binDist = 2 ** np.arange(start, end, interval)
                binDist = np.insert(binDist,0,0)
		for i in range(len(binDist)-1) :
                        chrList.append(thisChr)
                        startList.append(binDist[i]+overlapOffset)
                        endList.append(binDist[i+1])
                if binDist[i+1] < chrLen :
                        chrList.append(thisChr)
                        startList.append(binDist[i+1]+overlapOffset)
                        endList.append(chrLen)
                        #print chr +" " + str(offset)

        if returnLabels :
                binStr = [ assembly+'|'+chrList[i]+':'+str(startList[i])+'-'+str(endList[i]) for i in range(len(chrList)) ]
                return binStr

        return [chrList,startList,endList]

def binToIndex ( chr, start, resolution, binList ) :
	start = start - ( start % resolution ) + 1 
	chrInd = np.where(binList[0]==chr)[0]
	startInd = np.where(binList[1]==start)[0]
	finalInd = np.intersect1d(chrInd,startInd)
	if len (finalInd) != 1 :
		return -1
	return finalInd[0] 
	

def generateChrMidpointsLogBin ( assembly = 'hg19', chrs = [], resolution = 5 * 10 ** 5, returnLabels = False ) :
	if chrs == [] :
		chrs = getChrs(assembly)

	chrList = []; mpList =[];
	bin_sizes = [2500, 7500, 15000, 30000, 60000, 120000, 240000, 480000, 960000, 1920000, 2780000] 
	for chr in chrs :
		thisChr = chr
		chrLen = getChrLen( thisChr, assembly = assembly )
		for i in range(11) :
                        chrList.append(thisChr)
                        mpList.append(bin_sizes[i])
                for i in range(chrLen/resolution-5) :
                        chrList.append(thisChr)
			mpList.append((i+6)*resolution+resolution/2)

	if returnLabels :
		binStr = [ assembly+'|'+chrList[i]+':'+str(mpList[i]) for i in range(len(chrList)) ]
		return binStr

	return [chrList,mpList]					

def generateChrMidpointsLogBinNew ( assembly = 'hg19', chrs = [], interval = 0.125, returnLabels = False ) :
        if chrs == [] :
                chrs = getChrs(assembly)
        chrList = []; mpList =[]; start= 10
        for chr in chrs :
                thisChr = chr
                chrLen = getChrLen( thisChr, assembly = assembly )
		end = np.log2(chrLen)
                binDist = 2 ** np.arange(start, end, interval)
                binDist = np.insert(binDist,0,0)
		for i in range(len(binDist)-1) :
                        chrList.append(thisChr)
                        mpList.append((binDist[i+1] + binDist[i])/2)
		chrList.append(thisChr)
		mpList.append((binDist[i+1] + chrLen)/2)
		mpList = [float(Decimal("%.2f" % e)) for e in mpList]
        if returnLabels :
                binStr = [ assembly+'|'+chrList[i]+':'+str(mpList[i]) for i in range(len(chrList)) ]
                return binStr

        return [chrList,mpList]
		
def generateChrMidpoints ( assembly = 'hg19', chrs = [], resolution = 5 * 10 ** 5, returnLabels = False ) :
        if chrs == [] :
                chrs = getChrs(assembly)

        chrList = []; mpList =[];
        for chr in chrs :
                thisChr = chr
                chrLen = getChrLen( thisChr, assembly = assembly )
                for i in range(chrLen/resolution+1) :
                        chrList.append(thisChr)
                        mpList.append(i*resolution+resolution/2)

        if returnLabels :
                binStr = [ assembly+'|'+chrList[i]+':'+str(mpList[i]) for i in range(len(chrList)) ]
                return binStr

        return [chrList,mpList]	
		
	
def OutDekkerMatrix ( M, chr, resolution = 5 * 10**5, assembly = 'hg19', outfilename = '', printHeader = True, fieldStyle = 'Dekker' ):
	chrLen = getChrLen(chr,assembly=assembly)
	M = np.matrix(M)
	N = int( np.ceil( chrLen / float(resolution) ) )
	if ( N != len(M) ):
		print "Array not covering whole chromosome"
	if outfilename != '' :
		of = open(outfilename,'w')
	else :
		of = sys.stdout

	binCount = range( int( np.ceil( chrLen / float(resolution) ) ) )
	binStart = range(1,chrLen,resolution)
	binEnd = range(resolution,chrLen,resolution)
	if ( chrLen % resolution != 0 ): 
		binEnd.append(chrLen)

	binLabel = [ str(binCount[i])+'|'+assembly+'|'+chr+':'+str(binStart[i])+'-'+str(binEnd[i]) for i in range(len(binCount)) ]
	
	if fieldStyle == 'Dekker' :
		if printHeader is True :
			header=str(N)+'x'+str(N)+'\t'+('\t'.join(binLabel))+'\n'
			of.write(header)
		for i in range(len(M)) :
			of.write( binLabel[i] + '\t' + '\t'.join(map(str,np.matrix.tolist(M[i,])[0])) + '\n' )

	if fieldStyle == 'HiCRep' :
		for i in range(len(M)) :
			of.write( chr + '\t' + str(binStart[i]-1) + '\t' + str(binEnd[i]) + '\t' + '\t'.join(map(str,np.matrix.tolist(M[i,])[0])) + '\n' )



def ParseHiCProMatrix ( filename, assembly = 'hg19' ,resolution = 5 * 10**5 ):
	(chrs,starts, ends, offsetDict ) = generateChrBins ( assembly = assembly, resolution = resolution, returnOffset = True ) 
	uniquechrs = getChrs( assembly = assembly )

	chrDict = {}
	Matrices = []
	for chr,i in zip(uniquechrs,range(len(uniquechrs))) :
		N = int(np.ceil( float(getChrLen(chr,assembly=assembly)) / resolution ))
		Matrices.append( np.zeros((N,N)) )
		chrDict[chr] = i
		
	if filename[-3:] == '.gz' :
		file = gzip.open(filename)
	else :
		file = open(filename)
	#if header :
	#	file.readline()

	for line in file:
		line.strip()
                tokens = line.split()
		i = int(tokens[0])-1; k = int(tokens[1])-1; intValue = tokens[2]
		chr1 = chrs[i]; chr2 = chrs[k]; intValue = tokens[2]
		
		if  chr1 == chr2 :
			offset = int(offsetDict[chr1])
			Matrices[chrDict[chr1]][i-offset,k-offset] = intValue
			Matrices[chrDict[chr1]][k-offset,i-offset] = intValue
	
	#rowSums = sum(M)
	file.close()
	return(Matrices)


def ParseDekkerMatrix ( filename, chr = 'nochr', resolution = 10**6, returnLabel = False, irregularHeader= False ):
	if filename[-3:] == '.gz' :
		file = gzip.open(filename)
	else :
		file = open(filename)

	header = file.readline().strip()		
	while header[0] == '#' :
		header = file.readline().strip()

	if irregularHeader == False :
		N = len(header.split())-1
		M = np.zeros((N,N))
		i = 0
		labels = (header.split())[1:]
	else :
		N = len(header.split())
		M = np.zeros((N,N))
		i = 0
		labels = (header.split())


	for line in file:
		line.strip()
                tokens = line.split()
       	        M[i,] = map(float,tokens[1:len(tokens)])
		i = i + 1
	file.close()

	if chr != 'nochr' :
		labelchrs = [ label.split('|')[-1].split(':')[0] for label in labels]
		I = np.where(np.array(labelchrs) == chr)[0]
		labels = list(np.array(labels)[I])
		M = (M[:,I])[I,:]

	if returnLabel :
		return (M,labels)
	else :
		return(M)


def ParseDenseMatrix ( filename, resolution = 10**6, isInt = True ):
	if filename[-3:] == '.gz' :
		file = gzip.open(filename)
	else :
		file = open(filename)

	header = file.readline().strip(); header = header[2:]
	rowIDs = header.split()
	N = len(rowIDs)
	M = np.matrix( np.zeros((N,N)) )
	i = 0

	for line in file:
		line.strip()
                tokens = line.split()
		if isInt == True :
	       	        M[i,] = map(int,tokens)
		else :
	       	        M[i,] = map(float,tokens)
		i = i + 1
	
	#chrIDs = []
	#chrs = map(str,range(1,23)); chrs.append('X')
	#for chr in chrs :
	#	chrCount = rowIDs.count(chr)
	#	for i in range(chrCount) :
	#		chrIDs.append(i)
	#chrIDs = np.array(chrIDs)

	file.close()
	return(M,rowIDs)


def ParseWholeRaoMatrix ( folder , sample, resolution =10**6,assembly='hg19' ):
	#only human now
	if assembly == 'hg19' :
		#chrs = map(str,range(20,22))
		chrs = map(str,range(1,23))
		chrs.append('X')
		chrNames = [ 'chr' + chr for chr in chrs ]
	lengths = []
	fullMatrix = np.matrix([])
        for i in range(len(chrNames)):
		rowMatrix = np.matrix([])
		for k in range(len(chrNames)):
			
			chrN1 = chrNames[i]; chrN2 = chrNames[k];
			if ( i > k ) :
				n = getChrLen(chrN1,assembly=assembly) / resolution + 1
				m = getChrLen(chrN2,assembly=assembly) / resolution + 1
				thisMatrix=np.zeros((n,m))
			else :
				print chrN1 + ' ' + chrN2
				if i == k :
					# will not work for 1mb res
					resname = str(resolution/1000) + 'kb' 
					filename = folder + sample + '/' + resname + '_resolution_intrachromosomal/' + chrN1 + '/MAPQGE30/' + chrN1 + '_' + resname + '.RAWobserved'
				else :
					# will not work for 1mb res
					resname = str(resolution/1000) + 'kb' 
					filename = folder + sample + '_interchromosomal/' + resname + '_resolution_interchromosomal/' + chrN1 + '_' + chrN2 + '/MAPQGE30/' + chrN1 + '_' + chrN2[3:] + '_' + resname + '.RAWobserved'

				thisMatrix = ParseRaoMatrix(filename,chr=chrN1,chr2=chrN2,resolution=resolution,assembly=assembly)
				print ( np.shape(thisMatrix))
			if len(rowMatrix) == 1 :
				rowMatrix = thisMatrix
				lengths.append(len(thisMatrix))
			else :
				rowMatrix = np.hstack((rowMatrix,thisMatrix))

		if len(fullMatrix) == 1 :
			fullMatrix = rowMatrix
		else :
			fullMatrix = np.vstack((fullMatrix,rowMatrix))
	fullMatrix = fullMatrix + np.transpose(np.triu(fullMatrix,1))

	[chrList,startList,endList] = generateChrBins(assembly=assembly,chrs=chrs,resolution=resolution)
	binStr = [ assembly+'|'+chrList[i]+':'+str(startList[i])+'-'+str(endList[i]) for i in range(len(chrList)) ]

	return(fullMatrix,binStr)


def ParseRaoMatrix ( filename, chr, chr2 ='', resolution = 10**6, assembly = 'hg19', start = None, end = None ):
	if chr2 == '' :
		chr2 = chr

	if filename[-3:] == '.gz' :
		file = gzip.open(filename)
	else :
		file = open(filename)

	if chr == chr2 and start is not None and end is not None :
		R = ( end - start ) / resolution + 1
		C = R

		M = np.matrix( np.zeros((R,C)) )
		c = 0;
		for line in file:
			c += 1
			if c % 1000000 == 0 :
				print 'liebestraum'
			line.strip()
	                tokens = line.split()
       		        frag1 = int(tokens[0]) 
               		frag2 = int(tokens[1]) 

			if frag1 >=  start and frag1 <= end and frag2 >= start and frag2 <= end :
				i = ( frag1 - start ) / resolution; j = ( frag2 - start ) / resolution
		                c = int(float(tokens[2]))
       			        M [i,j] = c;
				if chr == chr2 :
					M [j,i] = c
		
	else :
		chrLen1 = getChrLen(chr,assembly=assembly)
		chrLen2 = getChrLen(chr2,assembly=assembly)
		R = chrLen1 / resolution + 1
		C = chrLen2 / resolution + 1

		M = np.matrix( np.zeros((R,C)) )

		for line in file:
			line.strip()
	                tokens = line.split()
       		        i = int(tokens[0]) / resolution
               		j = int(tokens[1]) / resolution
	                c = int(float(tokens[2]))
			if i < R and j < C :
	       		        M [i,j] = c;
				if chr == chr2 :
					 M [j,i] = c
	file.close()

	#rowSums = sum(M)
	#if ( rowSums[len(rowSums)-1] == 0 ):
	#	M = M[0:(len(M)-1),0:(len(M)-1)]
	return(M)

def ParseWholeInteractionBed ( filename , resolution = 5 * 10**5, assembly='hg19' ):
	chrNames = getChrs(assembly)
	lengths = []
	fullMatrix = np.matrix([])
        for i in range(len(chrNames)):
		rowMatrix = np.matrix([])
		for k in range(len(chrNames)):
			chrN1 = chrNames[i]; chrN2 = chrNames[k];
			#print chrN1 + ' ' + chrN2

			if ( i > k ) :
				n = getChrLen(chrN1,assembly=assembly) / resolution + 1
				m = getChrLen(chrN2,assembly=assembly) / resolution + 1
				thisMatrix=np.zeros((n,m))
			else :
				thisMatrix = ParseInteractionBed(filename,chr=chrN1,chr2=chrN2,resolution=resolution,assembly=assembly)
			if len(rowMatrix) == 1 :
				rowMatrix = thisMatrix
				lengths.append(len(thisMatrix))
			else :
				rowMatrix = np.hstack((rowMatrix,thisMatrix))

		if len(fullMatrix) == 1 :
			fullMatrix = rowMatrix
		else :
			fullMatrix = np.vstack((fullMatrix,rowMatrix))

	fullMatrix = fullMatrix + np.transpose(np.triu(fullMatrix,1))
	[chrList,startList,endList] = generateChrBins(assembly=assembly,chrs=chrNames,resolution=resolution)
	#print(len(chrList))
	binStr = [ assembly+'|'+chrList[i]+':'+str(startList[i])+'-'+str(endList[i]) for i in range(len(chrList)) ]

	return(fullMatrix,binStr)

def ParseWholeInteractionBedLogBin ( filename , resolution = 5 * 10**5, assembly='hg19' ):
        chrNames = getChrs(assembly)
        lengths = [];
        fullMatrix = np.matrix([])
        for i in range(len(chrNames)):
                rowMatrix = np.matrix([])
                for k in range(len(chrNames)):
                        chrN1 = chrNames[i]; chrN2 = chrNames[k];
                        #print chrN1 + ' ' + chrN2

                        if ( i > k ) :
                                n = getChrLen(chrN1,assembly=assembly) / resolution + 1 + 5
                                m = getChrLen(chrN2,assembly=assembly) / resolution + 1 + 5
                                thisMatrix=np.zeros((n,m))
                        else :
                                thisMatrix = ParseInteractionBedLogBin(filename,chr=chrN1,chr2=chrN2,resolution=resolution,assembly=assembly)
			if len(rowMatrix) == 1 :
                                rowMatrix = thisMatrix
                                lengths.append(len(thisMatrix))
                        else :
                                rowMatrix = np.hstack((rowMatrix,thisMatrix))

                if len(fullMatrix) == 1 :
                        fullMatrix = rowMatrix
                else :
                        fullMatrix = np.vstack((fullMatrix,rowMatrix))

        fullMatrix = fullMatrix + np.transpose(np.triu(fullMatrix,1))
        [chrList,startList,endList] = generateChrLogBins(assembly=assembly,chrs=chrNames,resolution=resolution)
        #print(len(chrList))
        binStr = [ assembly+'|'+chrList[i]+':'+str(startList[i])+'-'+str(endList[i]) for i in range(len(chrList)) ]
        return(fullMatrix,binStr)

def ParseWholeInteractionBedLogBinNew ( filename , interval = 0.125, assembly='hg19' ):
        chrNames = getChrs(assembly)
        lengths = []; start = 10
        fullMatrix = np.matrix([])
        for i in range(len(chrNames)):
                rowMatrix = np.matrix([])
                for k in range(len(chrNames)):
                        chrN1 = chrNames[i]; chrN2 = chrNames[k];
                        #print chrN1 + ' ' + chrN2

                        if ( i > k ) :
                                n = getChrLen(chrN1,assembly=assembly)
                                m = getChrLen(chrN2,assembly=assembly)
                		end = np.log2(n)
               	 		binDist = 2 ** np.arange(start, end, interval)
                		N1 = len(binDist) + 1
                		end = np.log2(m)
                		binDist2 = 2 ** np.arange(start, end, interval)
                		N2 = len(binDist2) + 1
                                thisMatrix=np.zeros((N1,N2))
                        else :
                                thisMatrix = ParseInteractionBedLogBinNew(filename,chr=chrN1,chr2=chrN2,interval = 0.125,assembly=assembly)
                        if len(rowMatrix) == 1 :
                                rowMatrix = thisMatrix
                                lengths.append(len(thisMatrix))
                        else :
                                rowMatrix = np.hstack((rowMatrix,thisMatrix))

                if len(fullMatrix) == 1 :
                        fullMatrix = rowMatrix
                else :
                        fullMatrix = np.vstack((fullMatrix,rowMatrix))

        fullMatrix = fullMatrix + np.transpose(np.triu(fullMatrix,1))
        [chrList,startList,endList] = generateChrLogBinsNew(assembly=assembly,chrs=chrNames,interval = 0.125)
        #print(len(chrList))
        binStr = [ assembly+'|'+chrList[i]+':'+str(startList[i])+'-'+str(endList[i]) for i in range(len(chrList)) ]
        return(fullMatrix,binStr)

def ParseInteractionBedLogBinNew ( filename, chr, chr2 = "", customN = 0, interval = 0.125, assembly = 'hg19', header = False, columnIndex = 4, isInt = True, zeros = True, chrPrefix = '' ):
	if chr2 == "" :
		chr2 = chr

	chr = chrPrefix + chr
        chr2 = chrPrefix + chr2
	
	if customN != 0 :
		N1 = customN; N2 = customN
	else :
		temp, mpList1= generateChrMidpointsLogBinNew(chrs=[chr])
		temp, mpList2= generateChrMidpointsLogBinNew(chrs=[chr2])
		N1 = len(mpList1)
		N2 = len(mpList2)
	if zeros :
		if isInt :
			M = np.matrix( np.zeros((N1,N2), dtype=np.uint32) )
		else :
			M = np.matrix( np.zeros((N1,N2), dtype=np.float32) )
	else :
		if isInt :
			M = np.matrix( np.ones((N1,N2), dtype=np.uint32) )
		else :
			M = np.matrix( np.ones((N1,N2), dtype=np.float32) )

	if filename[-3:] == '.gz' :
		file = gzip.open(filename)
	else:
		file = open(filename)
	if header :
		file.readline()


	for line in file:
		line.strip()
                tokens = line.split()
		
		if tokens[0] == chr and tokens[2] == chr2 :
	       	        i = float(Decimal(tokens[1])); j = float(Decimal(tokens[3]))
			if isInt :
		               	c = int(float(tokens[columnIndex]))
			else :
	               		c = float(tokens[columnIndex])

                        i_pos = mpList1.index(i)
                        j_pos = mpList2.index(j)

			M [ i_pos,  j_pos ] = c
        	       	
			if chr == chr2 :
				M [ j_pos, i_pos ] = c

		if tokens[2] == chr and tokens[0] == chr2 :
	       	        i = float(Decimal(tokens[3])); j = float(Decimal(tokens[1]))
			if isInt :
		               	c = int(float(tokens[columnIndex]))
			else :
	               		c = float(tokens[columnIndex])

                        i_pos = mpList1.index(i)
                        j_pos = mpList2.index(j)

                        M [ i_pos,  j_pos ] = c

                        if chr == chr2 :
                                M [ j_pos, i_pos ] = c

	#rowSums = sum(M)
	file.close()
	return(M)

def ParseInteractionBedLogBin ( filename, chr, chr2 = "", customN = 0, resolution = 5 * 10**5, assembly = 'hg19', header = False, columnIndex = 4, isInt = True, zeros = True, chrPrefix = '' ):
        if chr2 == "" :
                chr2 = chr
        if customN != 0 :
                N1 = customN; N2 = customN
        else :
                chrLen1 = getChrLen(chr,assembly=assembly)
                chrLen2 = getChrLen(chr2,assembly=assembly)
                N1 = chrLen1 / resolution + 1 + 5
                N2 = chrLen2 / resolution + 1 + 5
        if zeros :
                if isInt :
                        M = np.matrix( np.zeros((N1,N2), dtype=np.uint32) )
                else :
                        M = np.matrix( np.zeros((N1,N2), dtype=np.float32) )
        else :
                if isInt :
                        M = np.matrix( np.ones((N1,N2), dtype=np.uint32) )
                else :
                        M = np.matrix( np.ones((N1,N2), dtype=np.float32) )

        if filename[-3:] == '.gz' :
                file = gzip.open(filename)
        else :
                file = open(filename)
        if header :
                file.readline()

        chr = chrPrefix + chr
        chr2 = chrPrefix + chr2

        bin_sizes = [2500, 7500, 15000, 30000, 60000, 120000, 240000, 480000, 960000, 1920000, 2780000]

        for line in file:
                line.strip()
                tokens = line.split()

                if tokens[0] == chr and tokens[2] == chr2 :
                        i = int(tokens[1]); j = int(tokens[3])
                        if isInt :
                                c = int(float(tokens[columnIndex]))
                        else :
                                c = float(tokens[columnIndex])

                        if i < 3000000:
                                i_pos = [k for k, x in enumerate(bin_sizes) if x == i][0]
                        else :
                                i_pos = i / resolution + 5

                        if j < 3000000:
                                j_pos = [m for m, x in enumerate(bin_sizes) if x == j][0]
                        else :
                                j_pos = j / resolution + 5

                        M [ i_pos,  j_pos ] = c

                        if chr == chr2 :
                                M [ j_pos, i_pos ] = c

                if tokens[2] == chr and tokens[0] == chr2 :
                        i = int(tokens[3]); j = int(tokens[1])
                        if isInt :
                                c = int(float(tokens[columnIndex]))
                        else :
                                c = float(tokens[columnIndex])


                        if i < 3000000:
                                i_pos = [k for k, x in enumerate(bin_sizes) if x == i][0]
                        else :
                                i_pos = i / resolution + 5

                        if j < 3000000:
                                j_pos = [m for m, x in enumerate(bin_sizes) if x == j][0]
                        else :
                                j_pos = j / resolution + 5

                        M [ i_pos,  j_pos ] = c

                        if chr == chr2 :
                                M [ j_pos, i_pos ] = c

        #rowSums = sum(M)
        file.close()
        return(M)

def ParseInteractionBed ( filename, chr, chr2 = "", customN = 0, resolution = 5 * 10**5, assembly = 'hg19', header = False, columnIndex = 4, isInt = True, zeros = True, chrPrefix = '' ):
        if chr2 == "" :
                chr2 = chr
        if customN != 0 :
                N1 = customN; N2 = customN
        else :
                chrLen1 = getChrLen(chr,assembly=assembly)
                chrLen2 = getChrLen(chr2,assembly=assembly)
                N1 = chrLen1 / resolution + 1
                N2 = chrLen2 / resolution + 1
        if zeros :
                if isInt :
                        M = np.matrix( np.zeros((N1,N2), dtype=np.uint32) )
                else :
                        M = np.matrix( np.zeros((N1,N2), dtype=np.float32) )
        else :
                if isInt :
                        M = np.matrix( np.ones((N1,N2), dtype=np.uint32) )
                else :
                        M = np.matrix( np.ones((N1,N2), dtype=np.float32) )

        if filename[-3:] == '.gz' :
                file = gzip.open(filename)
        else :
                file = open(filename)
        if header :
                file.readline()

        chr = chrPrefix + chr
        chr2 = chrPrefix + chr2

        for line in file:
                line.strip()
                tokens = line.split()

                if tokens[0] == chr and tokens[2] == chr2 :
                        i = int(tokens[1]); j = int(tokens[3])
                        if isInt :
                                c = int(float(tokens[columnIndex]))
                        else :
                                c = float(tokens[columnIndex])
                        M [ i / resolution, j / resolution ] = c
                        if chr == chr2 :
                                M [ j / resolution, i / resolution ] = c

                if tokens[2] == chr and tokens[0] == chr2 :
                        i = int(tokens[3]); j = int(tokens[1])
                        if isInt :
                                c = int(float(tokens[columnIndex]))
                        else :
                                c = float(tokens[columnIndex])
                        M [ i / resolution, j / resolution ] = c
                        if chr == chr2 :
                                M [ j / resolution, i / resolution ] = c


        #rowSums = sum(M)
        file.close()
        return(M)


def ParseWholeInteractionBedByMatrix ( filename, resolution = 5 * 10**5, assembly = 'hg19', zeros = True, columnIndex = 4, isInt = True, header = False ):

	[chromosomes, starts, ends ] = generateChrBins ( assembly = assembly,resolution = resolution )
	chrs = getChrs(assembly)

	matrixDict = {}; c= 0;
	for chr, c  in zip(chrs,range(23)):
		matrixDict[chr] = c;
	
	Matrices = []
	for chr, c  in zip(chrs,range(23)):
		chrLen = getChrLen(chr,assembly=assembly)
		N = chrLen / resolution + 1
		#print N
		if zeros == True :
			Matrices.append( np.matrix( np.zeros((N,N)) ) )
		else :
			Matrices.append( np.matrix( np.ones((N,N)) ) )
		
	if filename[-3:] == '.gz' :
		file = gzip.open(filename)
	else :
		file = open(filename)
	if header :
		file.readline()

	#print 'Parsing Whole File'
	for line in file:
		line.strip()
                tokens = line.split()
		C1 = tokens[2]; C2 = tokens[0]
		if C1 == C2 and C1 != 'chrY' :
	       	        i = int(tokens[3]); j = int(tokens[1])
        	       	c = float(tokens[columnIndex])
			if isInt == True :
				c = int(c)
       	        	I = int(i / resolution); J = int(j / resolution)
	       	        Matrices[matrixDict[C1]][ I, J ] = c;
	       	        Matrices[matrixDict[C1]][ J, I ] = c;

	file.close()
	#print 'Done Parsing' 
	binStrs = []
	for chr in chrs :
		binStr = []
		for i in range(len(chromosomes)): 
			if chromosomes[i] == chr :
				binStr.append( chromosomes[i]+':'+str(starts[i])+'-'+str(ends[i]) )
		binStrs.append(binStr)

	return (Matrices,binStrs)

def ParseWholeInteractionBedByMatrixLogBin ( filename, resolution = 5 * 10**5, assembly = 'hg19', zeros = True, columnIndex = 4, isInt = True, header = False ):

        [chromosomes, starts, ends ] = generateChrLogBins ( assembly = assembly,resolution = resolution )
        chrs = getChrs(assembly)

        matrixDict = {}; c= 0;
        for chr, c  in zip(chrs,range(23)):
                matrixDict[chr] = c;

        Matrices = []
        for chr, c  in zip(chrs,range(23)):
                chrLen = getChrLen(chr,assembly=assembly)
                N = chrLen / resolution + 1 + 5
                #print N
                if zeros == True :
                        Matrices.append( np.matrix( np.zeros((N,N)) ) )
                else :
                        Matrices.append( np.matrix( np.ones((N,N)) ) )

        if filename[-3:] == '.gz' :
                file = gzip.open(filename)
        else :
                file = open(filename)
        if header :
                file.readline()

        #print 'Parsing Whole File'
        for line in file:
                line.strip()
                tokens = line.split()
                C1 = tokens[2]; C2 = tokens[0]
                if C1 == C2 and C1 != 'chrY' :
                        i = int(tokens[3]); j = int(tokens[1])
                        c = float(tokens[columnIndex])
                        if isInt == True :
                                c = int(c)
                        i_pos = i / resolution + 5
                        j_pos = j / resolution + 5

                        bin_sizes = [2500, 7500, 15000, 30000, 60000, 120000, 240000, 480000, 960000, 1920000, 2780000]
                        if i < 3000000:
                                i_pos = [k for k, x in enumerate(bin_sizes) if x == i][0]
                        if j < 3000000:
                                j_pos = [m for m, x in enumerate(bin_sizes) if x == j][0]

                        Matrices[matrixDict[C1]][ i_pos, j_pos ] = c;
                        Matrices[matrixDict[C1]][ j_pos, i_pos ] = c;

        file.close()
        #print 'Done Parsing' 
        binStrs = []
        for chr in chrs :
                binStr = []
                for i in range(len(chromosomes)):
                        if chromosomes[i] == chr :
                                binStr.append( chromosomes[i]+':'+str(starts[i])+'-'+str(ends[i]) )
                binStrs.append(binStr)

        return (Matrices,binStrs)

def ParseWholeInteractionBedByMatrixLogBinNew ( filename, interval = 0.125, assembly = 'hg19', zeros = True, columnIndex = 4, isInt = True, header = False ):
	
	[chromosomes, starts, ends ] = generateChrLogBinsNew ( assembly = assembly,interval = 0.125)

        chrs = getChrs(assembly)

        matrixDict = {}; c= 0;
        for chr, c  in zip(chrs,range(23)):
                matrixDict[chr] = c;
        Matrices = []
        for chr, c  in zip(chrs,range(23)):
                temp, mpList = generateChrMidpointsLogBinNew(chrs=[chr])
                N = len(mpList)
		#print N
                if zeros == True :
                        Matrices.append( np.matrix( np.zeros((N,N)) ) )
                else :
                        Matrices.append( np.matrix( np.ones((N,N)) ) )

        if filename[-3:] == '.gz' :
                file = gzip.open(filename)
        else :
                file = open(filename)
        if header :
                file.readline()

        #print 'Parsing Whole File'
        for line in file:
                line.strip()
                tokens = line.split()
                C1 = tokens[2]; C2 = tokens[0]
                if C1 == C2 and C1 != 'chrY' :
                        i = float(Decimal(tokens[3])); j = float(Decimal(tokens[1]));
                        c = float(tokens[columnIndex])
                        if isInt == True :
                                c = int(c)
                        
			temp, mpList = generateChrMidpointsLogBinNew(chrs=[C1])

                        i_pos = mpList.index(i)
                        j_pos = mpList.index(j)

                        Matrices[matrixDict[C1]][ i_pos, j_pos ] = c;
                        Matrices[matrixDict[C1]][ j_pos, i_pos ] = c;

        file.close()
        #print 'Done Parsing' 
        binStrs = []
        for chr in chrs :
                binStr = []
                for i in range(len(chromosomes)):
                        if chromosomes[i] == chr :
                                binStr.append( chromosomes[i]+':'+str("%.2f" % starts[i])+'-'+str("%.2f" % ends[i]) )
                binStrs.append(binStr)

        return (Matrices,binStrs)


def ParseWholeInteractionFast ( filename, resolution = 5 * 10**5, assembly = 'hg19', interchrOnly = False, columnIndex = 4, isInt = True, header = False, zeros = True ):
	[chrs, starts, ends, chrOffsets ] = generateChrBins ( assembly = assembly,resolution = resolution, returnOffset = True )
	totalBinCount = len(chrs)	
	if zeros :
		M = np.matrix( np.zeros((totalBinCount,totalBinCount)) )
	else :
		M = np.matrix( np.ones((totalBinCount,totalBinCount)) )

	if filename[-3:] == '.gz' :
		file = gzip.open(filename)
	else :
		file = open(filename)
	if header :
		file.readline()

	for line in file:
		line.strip()
                tokens = line.split()
		C1 = tokens[2]; C2 = tokens[0]
		if C1 == 'chrY' or C2 == 'chrY' :
			continue
		if interchrOnly and C1 == C2 :
			continue
       	        i = int(tokens[3]); j = int(tokens[1])
               	c = float(tokens[columnIndex])
		if isInt :
			c = int(c)
       	        I = int(chrOffsets[C1] + (i / resolution)); J = int(chrOffsets[C2] + (j / resolution))
       	        M [ I, J ] = c; M [ J, I ] = c

	file.close()

	binStr = [ assembly+'|'+chrs[i]+':'+str(starts[i])+'-'+str(ends[i]) for i in range(len(chrs)) ]

	return (M,binStr)


def ParseSegregatedInteractionBed ( filename, chr, resolution = 10**6, assembly = 'hg19', ref1 = '', ref2 ='' ):

	chrLen = getChrLen(chr,assembly=assembly)
	N = chrLen / resolution + 1
	M = np.zeros((N,N))

	if filename[-3:] == '.gz' :
		file = gzip.open(filename)
	else :
		file = open(filename)
	
	for line in file:
		line.strip()
                tokens = line.split()
		if tokens[0] == chr and tokens[2] == chr and tokens[4] == ref1 and tokens[5] == ref2 :
	       	        i = int(tokens[1])
        	       	j = int(tokens[3])
                	c = int(float(tokens[6]))
	       	        M [ i / resolution, j / resolution ] = c
        	       	M [ j / resolution, i / resolution ] = c

	rowSums = sum(M)
	file.close()
	return(M)

def OutInteractionFormat ( CM, filename, chr, resolution = 5 * 10**5, isInt = True, midPoint = True, append = False ):
	if append == False :
		out = open(filename,'w')
	else :
		out = open(filename,'a')


	for i in range(len(CM)):
		for k in range(len(CM)):
			if ( CM[i,k] > 0 ):
				if midPoint :
					firstmid = i * resolution + resolution / 2
					secondmid = k * resolution + resolution / 2
				else :
					firstmid = i * resolution
					secondmid = k * resolution

				if isInt :
					out.write ( "%s\t%d\t%s\t%d\t%d\n" % (chr,firstmid,chr,secondmid,CM[i,k]) )
				else :
					out.write ( "%s\t%d\t%s\t%d\t%f\n" % (chr,firstmid,chr,secondmid,CM[i,k]) )

	out.close()


def OutInteractionFormatLogBin ( CM, filename, chr, resolution = 5 * 10**5, isInt = True, midPoint = True, append = False ):
        if append == False :
                out = open(filename,'w')
        else :
                out = open(filename,'a')

        first_frag = [0, 5000, 10000, 20000, 40000, 80000, 160000, 320000, 640000, 1280000, 2560000]
        mid_bin  = [2500, 7500, 15000, 30000, 60000, 120000, 240000, 480000, 960000, 1920000, 2780000]

        for i in range(len(CM)):
                for k in range(len(CM)):
                        if ( CM[i,k] > 0 ):
                                if midPoint :
                                	if i < 11:
                                        	firstmid = mid_bin[i]
                                	elif i >= 11:
                                        	firstmid = (i-6)*resolution + resolution/2
                                	if k < 11:
                                        	secondmid = mid_bin[k]
                                	elif k >= 11:
                                        	secondmid = (k-6)*resolution + resolution/2
				else :
				        if i < 11:
                                                firstmid = first_frag[i]
                                        elif i >= 11 :
                                                firstmid = (i-6)*resolution
                                        if k < 11:
                                                secondmid = first_frag[k]
                                        elif k >= 11:
                                                secondmid = (k-6)*resolution
                                if isInt :
                                        out.write ( "%s\t%d\t%s\t%d\t%d\n" % (chr,firstmid,chr,secondmid,CM[i,k]) )
                                else :
                                        out.write ( "%s\t%d\t%s\t%d\t%f\n" % (chr,firstmid,chr,secondmid,CM[i,k]) )

        out.close()

def OutInteractionFormatLogBinNew ( CM, filename, chr, interval = 0.125, isInt = True, midPoint = True, append = False ):
        if append == False :
                out = open(filename,'w')
        else :
                out = open(filename,'a')


        for i in range(len(CM)):
                for k in range(len(CM)):
                        if ( CM[i,k] > 0 ):
                                if midPoint :
                                        firstmid = i * resolution + resolution / 2
                                        secondmid = k * resolution + resolution / 2
                                else :
                                        firstmid = i * resolution
                                        secondmid = k * resolution

                                if isInt :
                                        out.write ( "%s\t%d\t%s\t%d\t%d\n" % (chr,firstmid,chr,secondmid,CM[i,k]) )
                                else :
                                        out.write ( "%s\t%d\t%s\t%d\t%f\n" % (chr,firstmid,chr,secondmid,CM[i,k]) )

        out.close()

def OutInteractionFormatwithLabels ( CM, filename, labels, resolution,append = False, isInt = True ):
	if append == False :
		out = open(filename,'w')
	else :
		out = open(filename,'a')

	( I,K ) = np.shape(CM)

	for i in range(I):
		for k in range(i,K):
			if ( CM[i,k] > 0 ):
				(chr1,firstmid) = labels[i].split(':')
				(chr2,secondmid) = labels[k].split(':')
				firstmid = str(int(firstmid.split('-')[0]) + resolution/2 -1)
				secondmid = str(int(secondmid.split('-')[0]) + resolution/2 -1)

				if isInt :
					out.write ( "%s\t%s\t%s\t%s\t%d\n" % (chr1,firstmid,chr2,secondmid,CM[i,k]) )
				else :
					out.write ( "%s\t%s\t%s\t%s\t%f\n" % (chr1,firstmid,chr2,secondmid,CM[i,k]) )

	out.close()

def OutInteractionFormatwithLabelsLogBin ( CM, filename, labels, resolution,append = False, isInt = True ):
        if append == False :
                out = open(filename,'w')
        else :
                out = open(filename,'a')

        ( I,K ) = np.shape(CM)
	first_frag = [0, 5000, 10000, 20000, 40000, 80000, 160000, 320000, 640000, 1280000, 2560000]
	mid_bin  = [2500, 7500, 15000, 30000, 60000, 120000, 240000, 480000, 960000, 1920000, 2780000]
        for i in range(I):
                for k in range(i,K):
                        if ( CM[i,k] > 0 ):
                                (chr1,firstmid) = ((labels[i].split('|'))[1]).split(':')
                                (chr2,secondmid) = ((labels[k].split('|'))[1]).split(':')
				firstmid = int(firstmid.split('-')[0])-1; secondmid = int(secondmid.split('-')[0])-1
				if firstmid < 3000000:
                               		firstmid = [m for m, x in enumerate(first_frag) if x == firstmid][0]
					firstmid = str(mid_bin[firstmid])
				else :
					firstmid = str(firstmid + resolution/2)
                        	if secondmid < 3000000:
                                	secondmid = [m for m, x in enumerate(first_frag) if x == secondmid][0]
					secondmid = str(mid_bin[secondmid])
				else :
					secondmid = str(secondmid + resolution/2)

                                if isInt :
                                        out.write ( "%s\t%s\t%s\t%s\t%d\n" % (chr1,firstmid,chr2,secondmid,CM[i,k]) )
                                else :
                                        out.write ( "%s\t%s\t%s\t%s\t%f\n" % (chr1,firstmid,chr2,secondmid,CM[i,k]) )

        out.close()

def OutInteractionFormatwithLabelsLogBinNew ( CM, filename, labels, interval = 0.125, append = False, isInt = True ):
        if append == False :
                out = open(filename,'w')
        else :
                out = open(filename,'a')

        ( I,K ) = np.shape(CM)
        
	for i in range(I):
                for k in range(i,K):
                        if ( CM[i,k] > 0 ):
                                (chr1,firstmid) = ((labels[i].split('|'))[1]).split(':')
                                (chr2,secondmid) = ((labels[k].split('|'))[1]).split(':')
                                firstmid1 = float(firstmid.split('-')[0])-1; firstmid2 = float(firstmid.split('-')[1])
				secondmid1 = float(secondmid.split('-')[0])-1; secondmid2 = float(secondmid.split('-')[1])

				firstmid = (firstmid2 + firstmid1) / 2
				secondmid = (secondmid2 + secondmid1) / 2

                                if isInt :
                                        out.write ( "%s\t%.2f\t%s\t%.2f\t%d\n" % (chr1,firstmid,chr2,secondmid,CM[i,k]) )
                                else :
                                        out.write ( "%s\t%.2f\t%s\t%.2f\t%f\n" % (chr1,firstmid,chr2,secondmid,CM[i,k]) )

        out.close()

def OutHiClusterFormat ( CM, filename, resolution = 5 * 10**5,  append = False ):
        if append == False :
                out = open(filename,'w')
        else :
                out = open(filename,'a')


        for i in range(len(CM)):
                for k in range(len(CM)):
                        if ( CM[i,k] > 0 ):
                        	out.write ( "%d\t%d\t%d\n" % (i,k,CM[i,k]) )

        out.close()

def OutEpigenomeInteractionFormat ( CM, filename, chr, start = None, end = None, resolution = 10**6, maxDist = None, logTransform = False ):
        out = open(filename,'w')
        chrLen = getChrLen(chr)
        if maxDist is not None :
                distInd = int( np.ceil( maxDist / resolution ) )
        else :
                distInd = int( chrLen / resolution )

        if start is not None and end is not None :
                startInd = np.floor( start / float(resolution) )
               	endInd = np.ceil( end / float(resolution) )
              	offset = startInd * resolution
                CM = (CM[startInd:endInd,:])[:,startInd:endInd]
        else :
                offset = 0
        if logTransform :
                CM = np.log(CM)
                CM [ np.isinf(CM) ] = 0
        for i in range(len(CM)):
                for k in range(len(CM)):
                        #No Diagonal, makes for an easy life
                        if CM[i,k] > 0 and i != k and np.abs(i-k) < distInd :

                                firstStart = ( i * resolution ) + 1 + offset
                                secondStart = ( k * resolution ) + 1 + offset
                                firstEnd = ( (i+1) * resolution )  + offset
                                secondEnd = ( (k+1) * resolution ) + offset
                                if firstEnd > chrLen :
                                        firstEnd = chrLen
                                if secondEnd > chrLen :
                                        secondEnd = chrLen

                                out.write ( "%s,%d,%d %s:%d-%d %f\n" % (chr,firstStart,firstEnd,chr,secondStart,secondEnd,float(CM[i,k])) )
                                #out.write ( "%s,%d,%d %s:%d-%d %f\n" % (chr,firstStart,firstEnd,chr,secondStart,secondEnd,np.log(float(CM[i,k]))) )
        out.close()

def OutHibrowseInteractionFormat ( CM, filename, chr, start = None, end = None, resolution = 10**6 ):
	out = open(filename,'w')
	chrLen = getChrLen(chr)
	
	if start is not None and end is not None :
		startInd = np.floor( start / float(resolution) )
		endInd = np.ceil( end / float(resolution) )
		offset = startInd * resolution
		CM = (CM[startInd:endInd,:])[:,startInd:endInd]

	for i in range(len(CM)):
		for k in range(i,len(CM)):
			firstStart = ( i * resolution ) + offset
			secondStart = ( k * resolution ) + offset
			firstEnd = ( (i+1) * resolution )  + offset
			secondEnd = ( (k+1) * resolution ) + offset
			if firstEnd > chrLen :
				firstEnd = chrLen
			if secondEnd > chrLen :
				secondEnd = chrLen
			#print str(i) + ' ' + str(k)
			#print str(firstStart) + ' ' + str(firstEnd)
			#print str(secondStart) + ' ' + str(secondEnd)
				
			out.write ( "%s %d %d %s %d %d %f\n" % (chr,firstStart,firstEnd,chr,secondStart,secondEnd,float(CM[i,k])) )
	out.close()


def KillDiag(CM):

	for i in range(len(CM)):
		for k in range(len(CM)):
			if ( abs(i-k) < 2  ):
				CM[i,k] = 0
	return(CM)

def refillMatrix(CM, chr, resolution, labels=None, starts=None, assembly = 'hg19' ):
	chrLen = getChrLen(chr,assembly)
	print 'cuhcuh'
	N = int(np.ceil( float(chrLen) / resolution ))
	fullCM = np.zeros( (N, N ) )

	if labels is not None :
		starts = [ label.split('|')[1].split(':')[1].split('-')[0] for label in labels]

	for i in range(np.shape(CM)[0]) :
		for k in range(np.shape(CM)[1]) :
			fullCM[int(int(starts[i])/resolution),int(int(starts[k])/resolution)] = CM[i,k]
	
	return fullCM 

def SubSampleMatrix(CM, subSampleN = 1000000, symmetric = True ):

	if subSampleN >= np.sum(np.triu(CM)) : 
		print 'Asked for ' + str(subSampleN) + ' entries, matrix has ' + str(np.sum(np.triu(CM)))
		print 'Sampling more entries than available, returning original matrix'
		return CM

	index1 = []
	index2 = []
	subCM = np.zeros((len(CM),len(CM)))

	for i in range(len(CM)):
		for k in range(i,len(CM)):
			count = int(CM[i,k])
			v1=np.empty(count); v1.fill(i)				
			v2=np.empty(count); v2.fill(k)				
			index1.extend(v1); index2.extend(v2)
	
	index1 = np.array(index1)
	index2 = np.array(index2)
	shufIndex = range(0,len(index1));	random.shuffle(shufIndex);
	subSampleIndex = np.random.choice(shufIndex,size=subSampleN,replace=False)
	index1 = index1[subSampleIndex]
	index2 = index2[subSampleIndex]

	for i in range(len(index1)):
		a = int(index1[i]); b = int(index2[i])
		subCM[a,b] = subCM[a,b] + 1

	subCM = subCM + np.triu(subCM,1).T		
	return(np.matrix(subCM))

def SubSampleVector(CV, subSampleN = 1000000):

        if subSampleN >= sum(CV) :
                print 'Asked for ' + str(subSampleN) + ' entries, matrix has ' + str(sum(CV))
                print 'Sampling more entries than available, returning original matrix'
                return CV

        index1 = []
        subCV = np.zeros(len(CV))

        for i in range(len(CV)):
        	count = int(CV[i])
                v1=np.empty(count); v1.fill(i)
                index1.extend(v1);

        index1 = np.array(index1)
        shufIndex = range(0,len(index1));       random.shuffle(shufIndex);
        subSampleIndex = np.random.choice(shufIndex,size=subSampleN,replace=False)
        index1 = index1[subSampleIndex]

        for i in range(len(index1)):
                a = int(index1[i])
                subCV[a] = subCV[a] + 1

        return(subCV)


def countBins ( Matrices ):
	N = len(Matrices[0])
	binV = np.zeros(N)
	distV = np.zeros(N)

	countByDist = []
	for CM in Matrices:
		N = len(CM)
		for i in range(0,N):
			for k in range(i,N):
				dist = k-i
				distV[dist] = distV[dist] + CM[i,k]
				binV[dist] = binV[dist] + 1

	return( (distV,binV) )


def noiseMatrix ( CM, distanceProbVector, percentNoise = 0.1, mixProp = 0.1 ):

	CM = CM.astype(int)

	mapIndices = np.where( CM.sum(0) > 0 )[0]
	nonmapIndices = np.where( CM.sum(0) == 0 )[0]
	contactSum = sum(CM)

	#SeqDepth = np.matrix.sum(np.asmatrix(CM))
	SeqDepth = int( np.sum(CM) / 2 + np.sum(np.diagonal(CM)) / 2 )
	N = len(CM)
	countByDist = []

	for i in range(0,N):
		for k in range(i,N):
			dist = k-i
			if ( len(countByDist)-1 < dist ):
				countByDist.append( [ float(contactSum[i]) * contactSum[k] ] )
			else:
				countByDist[dist].append( float(contactSum[i]) * contactSum[k] )

	distV = distanceProbVector[0:N]

	#OMITTED NORMALIZATION
	#distV = distV / binV
	conProb = distV / sum(distV)
	uniformProb = np.ones(len(conProb)) / len(conProb)
	
	conProb = conProb * (1 - mixProp ) + uniformProb * mixProp
	sampleN = int(SeqDepth * percentNoise )
	indices = range(0,N)
	noiseMatrix = np.zeros((N,N))

	print 'Samping %d reads from noise model' % (sampleN)
	t1 = time.time()
	samplesDist = np.random.choice(indices,size=sampleN, p=conProb)

	for i in range(0,max(samplesDist)+1):
		if ( sum(countByDist[i]) != 0 ):
			noToSample = len(np.where(samplesDist==i)[0])
			samplesIndex = np.random.choice( range(0,N-i), size =noToSample, p = ( np.array( countByDist[i] ) / sum(countByDist[i]) ) )
			for k in range(len(samplesIndex)):
				a = samplesIndex[k]
				b = samplesIndex[k] + i
				noiseMatrix[a,b] = noiseMatrix[a,b] + 1
				noiseMatrix[b,a] = noiseMatrix[b,a] + 1
	
	t2 = time.time()
	print 'Time is %f' % (t2-t1)
	return ( noiseMatrix )
	

def stratifiedSample ( V, F, strataSize = 100 ):
	N = len(V)
	V = np.array(V)
	F = np.array(F)

	strataCount = int(np.ceil(float(N) / strataSize))
	sortInd = np.argsort(F)
	strata = []
	strataMax = []

	#print '%d to stratify, %d strata to be filled' % (N,strataCount)

	
	for i in range(strataCount) :
		stratum = V [ sortInd[ (strataSize*(i) ) : (strataSize*(i+1)) ] ]
		stratumF = F [ sortInd[ (strataSize*(i) ) : (strataSize*(i+1)) ] ]
		strata.append( stratum )
		strataMax.append(max(stratumF))
		#print str(strataSize*(i)) + ' ' + str(strataSize*(i+1)) + ' ' + str(len(stratum))


	sample = []
	for i in range(len(V) ):
		if ( F[i] == 0 ) :
			sample.append (0)
		else :
			stratumInd = 0
			for k in range(strataCount) :
				#if ( F[i] >= strataMax[k] ):
				if ( F[i] <= strataMax[k] ):
					stratumInd = k
					break
			if ( stratumInd == 0 ):
				stratumInd = k
			sample.append ( np.random.choice(strata[k],size=1)[0] )

	return ( sample )

def uniformMatrix ( CM, subSampleCount = 1000000, bias = False ):
	(R,C) = np.shape(CM)
	marginal = np.sum(np.array(CM),1)
	uniSampleCM = np.matrix( np.zeros((R,C)) )
	#triuSum = sum(np.arange(R)+1)
	
	indexMap = []
	indexProb = []
	for i in range(R) :
	    	for k in range(i,R) :
			if marginal[i] != 0 and marginal[k] != 0 :
	       			indexMap.append([i,k])
				if bias :
					indexProb.append(marginal[i] * marginal[k])

	if bias :
		totalProb = float(sum(indexProb))
		indexProb = [ iP / totalProb for iP in indexProb ]
		triuSample = np.random.choice(len(indexMap),subSampleCount,p=indexProb)
	else :
		triuSample = np.random.choice(len(indexMap),subSampleCount)
        	
	for s in triuSample :
	    	(i,k) = indexMap[s]
    		uniSampleCM[i,k] += 1
	uniSampleCM += np.transpose(np.triu(uniSampleCM,1))

	return (uniSampleCM)


def shuffleMatrix ( CM, stratumSize = 50 ):
	#Convert to integer
	CM = CM.astype(int)
	#Get marginals and number of rows
	contactSum = np.sum(np.array(CM),1)

	print contactSum 
	print contactSum[1]
	print contactSum[1] * contactSum[2]
	N = len(CM)

	# For matrix entry Mik, store Marginal i * Marginal k in CountByDist
	# and the Mik itself in matrixByDist
	countByDist = []
	matrixByDist = []
	for i in range(0,N):
		for k in range(i,N):
			dist = k-i
			if ( len(countByDist)-1 < dist ):
				countByDist.append( [ float(contactSum[i]) * contactSum[k] ] )
				matrixByDist.append( [ int( CM[i,k] ) ] )
			else:
				countByDist[dist].append( float(contactSum[i]) * contactSum[k] )
				matrixByDist[dist].append( int( CM[i,k] ) )
	
	noiseMatrix = np.zeros((N,N))
	t1 = time.time()

	for i in range(len(matrixByDist)):
	#for i in range(1):
		#print "dist is %d" % (i)
		thisSample = stratifiedSample(matrixByDist[i],countByDist[i],stratumSize)
		for k in range(len(thisSample)):				
			noiseMatrix[k,k+i] = thisSample[k]
		
	for i in range(0,N):
		for k in range(i,N):
			noiseMatrix[k,i] = noiseMatrix[i,k]
	
	t2 = time.time()
	print 'Time is %f' % (t2-t1)
	return ( noiseMatrix )

def plotHeatmap (CM, fname=None, title=None , perc = 98, perc2 = 0, percSim = False, cutoff = None, cutoff2 = None, axis = None, colorbar = False, show = True,
        colorbarax = None, digitize = False, labels = None, labelsloc = None, binNames = None,  colorThreshold = None, colorMapName = 'bwr' ):

        if cutoff is None:
                cutoff = np.percentile(CM,perc)
                #if cutoff < 1 :
                #       cutoff = 1
        if cutoff2 is None:
                cutoff2 = np.percentile(CM,perc2)

        if percSim and cutoff * cutoff2 < 0 :
                cutoff = np.max( ( abs(cutoff), abs(cutoff2) ) )
                cutoff2 = -np.max( ( abs(cutoff), abs(cutoff2) ) )

        CM[ cutoff <= CM ] = cutoff
        CM[ cutoff2 >= CM ] = cutoff2
        #print cutoff

        if colorThreshold is not None :
                (minT, maxT ) = colorThreshold
        else :
                #(minT, maxT ) = (np.min(CM), np.max(CM) )
                (minT, maxT ) = (np.min(CM), cutoff )

        if  digitize :
                DCM = np.digitize(CM,np.percentile(CM[CM!=0],np.arange(10,90,10)))
                DCM[CM==0] = 0
                CM = DCM
        if axis is not None :
		ref = axis.matshow(CM, origin="bottom",cmap=plt.get_cmap(colorMapName),vmax=maxT,vmin=minT)
		if title is not None :
                        axis.set_title(title)
                if labels is not None and labelsloc is not None :
                        axis.set_yticks( labelsloc )
                        axis.set_yticklabels( labels )
                        axis.set_xticks( [] )
                        #axis.set_xticks( labelsloc )
                        #axis.set_xticklabels( labels )
                        #axis.set_yticks( [] )
                elif binNames is not None :
                        xSeq = np.arange(np.shape(CM)[1], step = np.shape(CM)[1] / 4 )
                        ySeq = np.arange(np.shape(CM)[0], step = np.shape(CM)[0] / 4 )
                        print ySeq
                        #axis.set_xticks( xSeq ); axis.set_xticklabels ( binNames[xSeq] )
                        axis.set_xticks( [] ); axis.set_xticklabels ( [] )
                        axis.set_yticks( ySeq ); axis.set_yticklabels ( binNames[ySeq], fontsize = 150 )
                axis.set_xlim( 0, np.shape(CM)[1] )
                axis.set_ylim( 0, np.shape(CM)[0] )
		axis.set_xticks(np.arange(0, np.shape(CM)[0], 20))
		axis.set_yticks(np.arange(0, np.shape(CM)[1], 20))
		#labs = np.arange(0, np.shape(CM)[1], 20) / 2 # for 500k resolution
		#labs = [item + 'Mb' for item in labs.astype(str)]
		#axis.set_xticklabels(labs)
		#axis.set_yticklabels(labs)
                if colorbar and colorbarax is not None :
                        plt.colorbar(ref,cax=colorbarax,orientation='horizontal')

        #AXIS AND COLORBAR DO NOT GO TOGETHER
        else :
                #plt.figure(figsize=(12,12))
		ref = plt.matshow(CM, origin="bottom",cmap=matplotlib.cm.get_cmap(colorMapName),vmax=cutoff,vmin=minT)
                #ref.set_size_inches(c(12,12))
                if title is not None :
                        plt.title(title)
                if colorbar :
                        plt.colorbar()
                if labels is not None and labelsloc is not None :
                        plt.yticks( labelsloc, labels )
                        plt.xticks( [] )

        if  fname is not None :
                plt.savefig(fname,dpi=720)
        if show :
                plt.show()
        return(ref)
	

# for figures and powerpoint slides
# different colors and scale on color bars
def plotHeatmap2 (CM, fname=None, title=None , cutoff= 0.00002, axis = None, colorbar = False, show = True,colorbarax = None, digitize = False, alpha = 0.5, labels = None, labelsloc = None, binNames = None,  colorThreshold = None, colorMapName = 'bwr', vmax=0 ):
    	if colorThreshold is not None :
                (minT, maxT ) = colorThreshold
        else :
                #(minT, maxT ) = (np.min(CM), np.max(CM) )
                (minT, maxT ) = (np.min(CM), cutoff )
	if axis is not None :
		#cmap = matplotlib.cm.get_cmap('jet',20)
		#cmap.set_bad('white', 1.)
		#cmap = matplotlib.cm.get_cmap('Reds')
		cmap=matplotlib.cm.get_cmap(colorMapName)
		#cmap.set_bad(color='w')
		CM = np.ma.masked_where(CM <= cutoff, CM)
		ref = axis.matshow(CM, origin="bottom",cmap=matplotlib.cm.get_cmap(colorMapName),vmax=vmax*0.00001, vmin=minT, alpha = alpha)
		if title is not None :
			axis.set_title(title)
		if labels is not None and labelsloc is not None :
			axis.set_yticks( labelsloc )
			axis.set_yticklabels( labels )
			axis.set_xticks( [] )
			#axis.set_xticks( labelsloc )
			#axis.set_xticklabels( labels )
			#axis.set_yticks( [] )
		elif binNames is not None :
			xSeq = np.arange(np.shape(CM)[1], step = np.shape(CM)[1] / 4 )
			ySeq = np.arange(np.shape(CM)[0], step = np.shape(CM)[0] / 4 )
			print ySeq
			#axis.set_xticks( xSeq ); axis.set_xticklabels ( binNames[xSeq] )
			axis.set_xticks( [] ); axis.set_xticklabels ( [] )
			axis.set_yticks( ySeq ); axis.set_yticklabels ( binNames[ySeq], fontsize = 150 )
		axis.set_xlim( 0, np.shape(CM)[1] )
		axis.set_ylim( 0, np.shape(CM)[0] )
		axis.set_yticks(np.arange(0, np.shape(CM)[0], 10))
		if colorbar and colorbarax is not None :
			plt.colorbar(ref,cax=colorbarax,orientation='horizontal')

	#AXIS AND COLORBAR DO NOT GO TOGETHER
	else :
		#plt.figure(figsize=(12,12))
		ref = plt.matshow(CM, origin="bottom",cmap=matplotlib.cm.get_cmap('jet', 10).set_bad('w'),vmax=np.median(CM),vmin=minT)
		#ref = plt.matshow(CM, origin="bottom",cmap=matplotlib.cm.get_cmap(colorMapName),vmax=maxT,vmin=minT)
		#ref.set_size_inches(c(12,12))
		if title is not None :
			plt.title(title)
		if colorbar : 
			plt.colorbar()
		if labels is not None and labelsloc is not None :
			plt.yticks( labelsloc, labels )
			plt.xticks( [] )

	if  fname is not None :
		plt.savefig(fname,dpi=720)
	if show :
		plt.show()
	return(ref)

def plotHeatmapLPcomp (CM, fname=None, title=None , perc = 98, perc2 = 0, percSim = False, cutoff = None, cutoff2 = None, axis = None, colorbar = False, show = True,
        colorbarax = None, digitize = False, labels = None, labelsloc = None, binNames = None,  colorThreshold = None, colorMapName = 'bwr' ):

        if cutoff is None:
                cutoff = np.percentile(CM,perc)
                if cutoff < 1 :
                        cutoff = 1
        if cutoff2 is None:
                cutoff2 = np.percentile(CM,perc2)

        if percSim and cutoff * cutoff2 < 0 :
                cutoff = np.max( ( abs(cutoff), abs(cutoff2) ) )
                cutoff2 = -np.max( ( abs(cutoff), abs(cutoff2) ) )

        CM[ cutoff <= CM ] = cutoff
        CM[ cutoff2 >= CM ] = cutoff2

        if colorThreshold is not None :
                (minT, maxT ) = colorThreshold
        else :
                #(minT, maxT ) = (np.min(CM), np.max(CM) )
                (minT, maxT ) = (np.min(CM), cutoff )

        if  digitize :
                DCM = np.digitize(CM,np.percentile(CM[CM!=0],np.arange(10,90,10)))
                DCM[CM==0] = 0
                CM = DCM
        if axis is not None :
                cmap = matplotlib.cm.get_cmap('Reds')
                ref = axis.matshow(CM, origin="bottom",cmap=cmap,vmax=np.max(CM)*0.1,vmin=minT)
                axis.set_xlim( 0, np.shape(CM)[1] )
                axis.set_ylim( 0, np.shape(CM)[0] )
                axis.set_yticks(np.arange(0, np.shape(CM)[0], 10))
                if colorbar and colorbarax is not None :
                        plt.colorbar(ref,cax=colorbarax,orientation='horizontal')

        else :
                ref = plt.matshow(CM, origin="bottom",cmap=matplotlib.cm.get_cmap('jet', 10).set_bad('w'),vmax=np.median(CM),vmin=minT)
        
	if  fname is not None :
                plt.savefig(fname,dpi=720)
        if show :
                plt.show()
        return(ref)


def iceMatrix (CM, labels = [],  percentile = 0, outputBias = False ) :
	#THIS LINE IS FOR COMPABILITY OF NEWER VERSIONS OF NUMPY AND
	#CONVERING A MATRIX INTO NDARRAY AS ICED ASKS FOR..
	CM = np.asarray(CM)

	if ( percentile != 0 and percentile < 100 ) :
		rowSum = sum(CM)
		I = np.squeeze( np.asarray( rowSum > 1 ) )
		CM_f1 = CM[I,:]; CM_f1 = CM_f1[:,I]

		rowSum = sum(CM_f1)
		I = np.squeeze( np.asarray( rowSum > np.percentile(rowSum,percentile) ) )
		CM_f2 = CM_f1[I,:]; CM_f2 = CM_f2[:,I]
		Ind = np.where( np.squeeze( np.asarray( sum(CM) ) ) > np.percentile(sum(CM_f1),percentile) ) [0]
		CM = CM_f2
		if labels != [] :
			labels = [ labels[i] for i in Ind ]

	if outputBias :
		(NCM, biasVector) = normalization.ICE_normalization(CM, output_bias = outputBias )
	else :
		NCM = normalization.ICE_normalization(CM, output_bias = outputBias )

	if labels != [] :
		return(NCM,labels)	
	elif outputBias :
		return(NCM,biasVector)
	else :
		return(NCM)	


def filterMatrix (CM, percentile = 1, lower = True, returnIndex = False, returnBoth = False ) :
	if ( percentile >= 0 and percentile < 100 ) :
		rowSum = np.asarray(np.sum(CM,0)).reshape(-1)
		nonzeroIndex = np.where( rowSum > 0 )[0]
		#print (np.shape(CM))
		#print len(nonzeroIndex)

		#This USED TO BE EQUAL OR GREATER TO, NOW IT DOES NOT ALLOW EQUALITY
		#if lower :
		#	percIndex = np.where ( rowSum >= np.percentile(rowSum[nonzeroIndex],percentile) )[0]
		#else :
		#	percIndex = np.where ( rowSum <= np.percentile(rowSum[nonzeroIndex],percentile) )[0]

		if lower :
			percIndex = np.where ( rowSum > np.percentile(rowSum[nonzeroIndex],percentile) )[0]
		else :
			percIndex = np.where ( rowSum < np.percentile(rowSum[nonzeroIndex],percentile) )[0]

		#print np.percentile(rowSum[nonzeroIndex],percentile)
		#print np.sort(rowSum)[0:100]
		#print np.sort(rowSum[[nonzeroIndex]])[0:100]

		Ind = np.intersect1d ( nonzeroIndex, percIndex )
		
		#print len(Ind)
		#print Ind
	else :
		print ' Wrong percentile value, exiting '
		return CM

	filteredCM = (CM[Ind,:])[:,Ind]

	if returnBoth :
		return (filteredCM,Ind)
	elif returnIndex :
		return Ind
	else :
		return filteredCM 

def simpleNormalize( CM ) :
	contactSum = np.sum(CM,1)
	NCM = CM
	#FIX THIS
	SeqDepth = np.matrix.sum(np.asmatrix(CM))
	N = len(CM)
	(J, K) = np.shape(CM)
	for j in range(J) :
		for k in range(j,K) :
			normCount = CM[j,k] * SeqDepth / ( contactSum[j] * contactSum[k] )
			NCM[j,k] = normCount
			NCM[k,j] = normCount

	return(NCM)

def coverageNorm (CM, labels = [], percentile = 0, returnInd = False ) :
	rowSum = sum(CM)
        I = np.squeeze( np.asarray( rowSum > 0 ) )
        #print np.where( np.squeeze( np.asarray( rowSum == 0 ) ) == True )
        CM.f1 = CM[I,:]; CM.f1 = CM.f1[:,I]
        rowSum = sum(CM.f1)
        I = np.squeeze( np.asarray( rowSum >= np.percentile(rowSum,percentile) ) )
        CM.f2 = CM.f1[I,:]; CM.f2 = CM.f2[:,I]
        Ind = np.where( np.squeeze( np.asarray( sum(CM) ) ) > np.percentile(sum(CM.f1),percentile) ) [0]
        #print np.percentile(sum(CM),percentile)
        #print Ind
        #print len(Ind)
        #print np.shape(CM.f2)
        CM = CM.f2

	#print ( np.shape(CM))
	#print( len(labels))
	NM = np.outer(sum(CM),sum(CM))
	NM = np.sqrt(NM)
	NCM = CM / NM
	if labels != [] :
	        labels = [ labels[i] for i in Ind ]
		if returnInd :
			return(NCM,labels,Ind)	
		else :
			return(NCM,labels)	
	else :
		if returnInd :
			return(NCM,Ind)
		else :
			return(NCM)	


def plotFullMatrix ( CM, labels ) :
	chrs = map(str,range(1,23))
	#chrs = []
	chrs.append('X')
	
	#M = np.matrix ( np.zeros(23,23) )
	for chr1 in chrs :
		for chr2 in chrs :
			Ind1 = np.where( labels == str(chr1) )[0]
			Ind2 = np.where( labels == str(chr2) )[0]
			thisCM = CM[Ind1,:]
			thisCM = thisCM[:,Ind2]
			plotHeatmap(thisCM,chr1 + chr2 + '.png')
		break

def extractfromFullMatrix ( CM, labels, chr1='1',chr2='2', returnLabel= False ) :
	if chr1[0:3] == 'chr' :
		chr1 = chr1[3:len(chr1)]
	if chr2[0:3] == 'chr' :
		chr2 = chr2[3:len(chr2)]

	chrLabels = [ label.split('|')[1].split(':')[0] for label in labels]
	Ind1 = np.where( np.array(chrLabels) == 'chr'+str(chr1) )[0]
	Ind2 = np.where( np.array(chrLabels) == 'chr'+str(chr2) )[0]
	thisCM = CM[Ind1,:]
	thisCM = thisCM[:,Ind2]
	labels1 = [ labels[i] for i in Ind1 ]
	labels2 = [ labels[i] for i in Ind2 ]

	if returnLabel :
		return(thisCM,labels1,labels2)
	else :
		return thisCM


def plotChrbyChrMatrix ( CM, labels, fname, assembly = 'hg19', colorThreshold = '', returnMatrix = True, plot = True ) :
	chrs = getChrs(assembly)
	print chrs
	M = np.matrix ( np.zeros(( len(chrs),len(chrs) )) )
	binCountM = ( np.zeros((23,23)) )
	binCount = np.array([])
	intCount = np.array([])
	for i in range(len(chrs)) :
		chr1 = chrs[i]
		for k in range(i,len(chrs)) :
			if ( i != k ) :
				chr2 = chrs[k]
				thisCM = extractfromFullMatrix ( CM, labels, chr1=chr1,chr2=chr2)	
				#plotHeatmap(thisCM,'meh'+chr1+'_'+chr2+'.png')
				sum = np.matrix.sum(thisCM)
				M[i,k] = sum
				M[k,i] = M[i,k]
				binCountM[i,k] = np.shape(thisCM)[0] * np.shape(thisCM)[1]
				binCountM[k,i] = binCountM[i,k]
				binCount = np.append(binCount, [ binCountM[i,k] ] )
				intCount = np.append(intCount, [ sum ] )
				#print chr1 + ' ' + chr2 + str(M[i,k])

	#plotHeatmap(M,fname+'before.png',labels=chrs,colorbar=True,perc=100)
	
	#print (binCount)
	#print (intCount)
	
	expCountPerBin = np.sum(intCount) / float(np.sum(binCount))
	#print (expCountPerBin)

	for i in range(len(chrs)) :
		for k in range(len(chrs)) :
			if ( i != k ) :
				#print str(M[i,k]) + ' ' + str( M[i,k] / ( expCountPerBin * binCountM[i,k] ) )
				M[i,k] = M[i,k] / ( expCountPerBin * binCountM[i,k] )
			#else :
			#	M[i,k] = np.nan
				
	#print(np.sum(M,1))
	#np.savetxt('test.txt',M)
	#plotHeatmap(M,fname+'.beforeICE.png',labels=chrs,colorbar=True,perc=100,colorThreshold=colorThreshold)
	#M = iceMatrix(M)
	#print(np.sum(M,1))
	#np.savetxt('test2.txt',M)
	if plot :
		plotHeatmap(M,fname,labels=chrs,labelsloc=range(0,len(chrs)),colorbar=True,perc=100,colorThreshold=colorThreshold)
	if returnMatrix :
		return M

def inferContactProbability ( CM, labels, resolution = 10**6, counts = False, prob = True, addPseudo = True ) :
	chrs = map(str,range(1,23))
	chrs.append('X')
	chrNames = [ 'chr' + chr for chr in chrs ]	
	
	#ASSUMING CHR1 is LONGEST
	distCount = np.zeros( ( getChrLen('chr1') / resolution ) + 1 ) + 1
	binCount = np.zeros( ( getChrLen('chr1') / resolution ) + 1 )
	for i in range(len(chrs)) :
		thisCM, thisLabel, thisLabel2 = extractfromFullMatrix (CM, labels, chr1 = chrs[i], chr2 = chrs[i], returnLabel = True )
		starts = [ float( label.split('|')[1].split(':')[1].split('-')[0] ) for label in thisLabel ]
		(J, K) = np.shape(thisCM)
		for j in range(J) :
			if np.sum(thisCM[j,:]) != 0 :
				for k in range(j,K) :
					#dist = k - j		
					dist = int( np.ceil( ( starts[k] - starts[j] ) / resolution ) )
					binCount[dist] += 1
					distCount[dist] += thisCM[j,k]
	if addPseudo :
		binCount[np.where(binCount==0)] = 1
	scaledCount = distCount/binCount
	if counts :
		return(distCount)
	else :
		if prob :
			return(scaledCount / np.sum(scaledCount) )
		else :
			return(scaledCount)


#def inferContactProbabilityFromInteraction ( Mlist, labelList, resolution = 10**6, assembly = 'hg19', counts = False, prob = True, addPseudo = True ) :
def inferContactProbabilityFromInteraction ( Mlist, labelList, resolution = 10**6, assembly = 'hg19', returnProb = False ) :
	
	chrNames = getChrs(assembly)	

	#ASSUMING CHR1 is LONGEST
	distCount = np.ones( ( getChrLen('chr1') / resolution ) + 1 )
	binCount = np.zeros( ( getChrLen('chr1') / resolution ) + 1 )

	rangeIndex = len( labelList[0][0].split('|') ) - 1
	#print (rangeIndex)

	for i in range(len(Mlist)) :
		#print i
		thisCM = Mlist[i]; thisLabel = labelList[i]
		starts = [ float( label.split('|')[rangeIndex].split(':')[1].split('-')[0] ) for label in thisLabel ]
		(J, K) = np.shape(thisCM)
		for j in range(J) :
			if np.sum(thisCM[j,:]) != 0 :
				for k in range(j,K) :
					dist = int( np.ceil( ( starts[k] - starts[j] ) / resolution ) )
					binCount[dist] += 1
					distCount[dist] += thisCM[j,k]


	#if addPseudo :
	#	binCount[np.where(binCount==0)] = 1
	#scaledCount = distCount/binCount
	#if counts :
	#	return(distCount)
	#else :
	#	if prob :
	#		return(scaledCount / np.sum(scaledCount) )
	#	else :
	#		return(scaledCount)

	if returnProb :
		binCount[np.where(binCount==0)[0]] = 1
		scaledCount = distCount/binCount
		return (scaledCount / np.sum(scaledCount) )
	else :
		return (binCount,distCount)


def inferContactProbabilityFromInteractionByChr ( Mlist, labelList, resolution = 10**6, assembly = 'hg19', returnProb = False ) :
	
	chrNames = getChrs(assembly)	

	#ASSUMING CHR1 is LONGEST
	
	distCountList =  []	
	binCountList =  []	

	rangeIndex = len( labelList[0][0].split('|') ) - 1
	#print (rangeIndex)

	for i in range(len(Mlist)) :
		#print i
		distCount = np.ones( ( getChrLen('chr1') / resolution ) + 1 )
		binCount = np.zeros( ( getChrLen('chr1') / resolution ) + 1 )

		thisCM = Mlist[i]; thisLabel = labelList[i]
		starts = [ float( label.split('|')[rangeIndex].split(':')[1].split('-')[0] ) for label in thisLabel ]
		(J, K) = np.shape(thisCM)
		for j in range(J) :
			if np.sum(thisCM[j,:]) != 0 :
				for k in range(j,K) :
					dist = int( np.ceil( ( starts[k] - starts[j] ) / resolution ) )
					binCount[dist] += 1
					distCount[dist] += thisCM[j,k]
		binCountList.append(binCount)
		distCountList.append(distCount)

	if returnProb :
		probList = []
		for i in range(len(binCountList)) :
			binCount = binCountList[i];	distCount = distCountList[i]
			binCount[np.where(binCount==0)[0]] = 1
			scaledCount = distCount/binCount
			probList.append(scaledCount / np.sum(scaledCount))
		return probList
	else :
		return (binCountList,distCountList)


def inferContactProbabilityByTAD ( Mlist, labelList, tads, resolution = 40000, assembly = 'hg19', chrNames = None, byCompartment = False, TADstrata = None ) : 
	if TADstrata is None :
		TADstrata = np.array ( [1000000,1500000,2000000] ) 
	if ( chrNames ) is None :
		chrNames = getChrs(assembly)	
	if byCompartment is False :
		distCount = np.zeros ( ( len(TADstrata),max(TADstrata)/resolution ) )
		binCount = np.zeros ( ( len(TADstrata),max(TADstrata)/resolution ) )
	else :
		distCountA = np.zeros ( ( len(TADstrata),max(TADstrata)/resolution ) );
		binCountA = np.zeros ( ( len(TADstrata),max(TADstrata)/resolution ) ); 
		distCountB = np.zeros ( ( len(TADstrata),max(TADstrata)/resolution ) );
		binCountB = np.zeros ( ( len(TADstrata),max(TADstrata)/resolution ) ); 


	rangeIndex = len( labelList[0][0].split('|') ) - 1

	for i in range(len(Mlist)) :
		thisCM = Mlist[i]; thisLabel = labelList[i]
		starts = [ float( label.split('|')[rangeIndex].split(':')[1].split('-')[0] ) for label in thisLabel ]
		Ind = np.where(tads['chr']==chrNames[i])[0]
		for k in Ind :
			tadSize = tads['end'][k] - tads['start'][k]
			if byCompartment :
				compValue = tads['compave'][k]
			if tadSize > np.max(TADstrata) :
				continue
			tadSizeIndex = np.min( ( np.sum( tadSize > TADstrata ), len(TADstrata)-1 ) )
			startI = tads['start'][k] /resolution; endI = tads['end'][k] /resolution; 
			for j in range(startI,endI) :
				if np.sum(thisCM[j,:]) != 0 :
					for k in range(j,endI) :
						dist = int( np.ceil( ( starts[k] - starts[j] ) / resolution ) )
						if byCompartment is False :
							binCount[tadSizeIndex,dist] += 1
							distCount[tadSizeIndex,dist] += thisCM[j,k]
						else :
							if compValue > 0 :
								binCountA[tadSizeIndex,dist] += 1
								distCountA[tadSizeIndex,dist] += thisCM[j,k]
							if compValue < 0 :
								binCountB[tadSizeIndex,dist] += 1
								distCountB[tadSizeIndex,dist] += thisCM[j,k]
	if byCompartment is True :
		return(binCountA,distCountA,binCountB,distCountB)
	else :
		return(binCount,distCount)


def inferContactProbabilityByCompartmentFromInteraction ( Mlist, labelList, compartment, resolution = 10**6, assembly = 'hg19', returnProb = False  ) :
	
	chrNames = getChrs(assembly)	

	#ASSUMING CHR1 is LONGEST
	distCount_A = np.ones( ( getChrLen('chr1') / resolution ) + 1 ); binCount_A = np.zeros( ( getChrLen('chr1') / resolution ) + 1 )
	distCount_B = np.ones( ( getChrLen('chr1') / resolution ) + 1 ); binCount_B = np.zeros( ( getChrLen('chr1') / resolution ) + 1 )
	distCount_AB = np.ones( ( getChrLen('chr1') / resolution ) + 1 ); binCount_AB = np.zeros( ( getChrLen('chr1') / resolution ) + 1 )

	rangeIndex = len( labelList[0][0].split('|') ) - 1

	for i in range(len(Mlist)) :
		thisCM = Mlist[i]; thisLabel = labelList[i]
		starts = [ float( label.split('|')[rangeIndex].split(':')[1].split('-')[0] ) for label in thisLabel ]
		compV = compartment['comp'] [ np.where(compartment['chr'] == chrNames[i])[0] ]
		(J, K) = np.shape(thisCM)
		for j in range(J) :
			if np.sum(thisCM[j,:]) != 0 :
				for k in range(j,K) :
					dist = int( np.ceil( ( starts[k] - starts[j] ) / resolution ) )
					comp1 = compV[int(starts[k] / resolution)];	comp2 = compV[int(starts[j] / resolution)]
					if ( comp1 == 0 or comp2 == 0 ) :
						continue;
					elif ( comp1 > 0 and comp2 > 0 ) :
						binCount_A[dist] += 1; distCount_A[dist] += thisCM[j,k]
					elif ( comp1 < 0 and comp2 < 0 ) :
						binCount_B[dist] += 1; distCount_B[dist] += thisCM[j,k]
					else :
						binCount_AB[dist] += 1; distCount_AB[dist] += thisCM[j,k]

	
	if returnProb :
		binCount_A[np.where(binCount_A==0)] = 1; scaledCount_A = distCount_A/binCount_A
		binCount_B[np.where(binCount_B==0)] = 1; scaledCount_B = distCount_B/binCount_B
		binCount_AB[np.where(binCount_AB==0)] = 1; scaledCount_AB = distCount_AB/binCount_AB
		return (scaledCount_A / np.sum(scaledCount_A),
			scaledCount_B / np.sum(scaledCount_B),
			scaledCount_AB / np.sum(scaledCount_AB) )
	else :
		return (binCount_A,binCount_B,binCount_AB,distCount_A,distCount_B,distCount_AB)



def plotTAD (m, chr, resolution, mStarts, TADs, subStart = None, subEnd = None, plotName = None ) :
        m = np.zeros(np.shape(m)) + np.triu(m)
        rM = refillMatrix(m,chr,resolution,starts=mStarts)
        (chrs,starts,ends) = generateChrBins(chrs=[chr],resolution=resolution,overlap=True)
        if subStart is None :
                subStart = 0
        if subEnd is None :
                subEnd = np.shape(rM)[0]

        binNames= np.array( [ s + 'kb' for s in map( str, np.array(mStarts)[subStart:subEnd]/1000 ) ] )
        fig,ax = plt.subplots()
        plotHeatmap(rM[subStart:subEnd,subStart:subEnd],show=False,axis=ax,perc=95,binNames=binNames)

        colors = ['green','orange','cyan','magenta']
        for i in range(len(TADs['start'])-1) :
                tadStart = np.where(starts == TADs['start'][i])[0][0]
                tadEnd = np.where(ends == TADs['end'][i])[0][0]
                plt.plot([tadStart+.5, tadEnd-.5], [tadStart+.5+1, tadEnd-.5+1], color=colors[i%4], linestyle='-', linewidth=2)

        if plotName is not None	:
                plt.savefig(plotName)




def removeDistanceEffect ( CM, labels, prob, resolution = 1000000 ) :
	distanceRemovedCM = CM
	CMsum = np.sum(np.triu(CM))
	rangeIndex = len( labels[0].split('|') ) - 1
	starts = [ float( label.split('|')[rangeIndex].split(':')[1].split('-')[0] ) for label in labels ]

	(J, K) = np.shape(CM)
	for j in range(J) :
		if np.sum(CM[j,:]) != 0 :
			for k in range(K) :
				dist = int( np.ceil( ( np.absolute( starts[k] - starts[j] )  ) / resolution ) )
				#print str(k) + ' ' + str(j) + ' ' + str(dist)
				#print CM[j,k]; print prob[dist]; print CMsum
				distanceRemovedCM[j,k] = CM[j,k] / ( prob[dist] * CMsum )

	return distanceRemovedCM

def generateExpected ( CM, labels, chr, probability, returnLabel = False ) :
	(thisCM,label1,label2) = extractfromFullMatrix (CM, labels, chr1 = chr, chr2 = chr, returnLabel = True )	
	expCM = thisCM
	sumCounts = np.sum( np.triu(thisCM) )
	
	(J, K) = np.shape(thisCM)
	for j in range(J) :
		for k in range(j,K) :
			dist = k - j
			binCount = J - dist		
			expCount = probability[dist] * sumCounts
			#expCount = probability[dist] * sumCounts / binCount
			expCM[j,k] = expCount; expCM[k,j] = expCount; 
			# TEMPORARY MEASURE FOR 0 DIAGONAL
			#if dist == 0 :
			#	expCM[j,k] = 1
			
			

	if returnLabel :
		return (expCM,label1)
	else :
		return (expCM)

def decomposeContactMatrix ( CM, labels, contactProb, chr, corrMethod = 'spearman', resolution = 1000000, validationV = None, assembly = 'hg19', returnDist = False, EigI = 0 ) :
	distanceRemovedCM = removeDistanceEffect(CM,labels,contactProb,resolution=resolution)
	rangeIndex = len( labels[0].split('|') ) - 1
	starts = [ float( label.split('|')[rangeIndex].split(':')[1].split('-')[0] ) for label in labels ]
	print getChrLen(chr,assembly)
	print chr
	N = int( np.ceil ( getChrLen(chr,assembly) / float(resolution) ) )
	print 'N ' + str(N)

	#FILTER ZERO SUM ROWS
	rowSum = np.asarray(np.sum(distanceRemovedCM,0)).reshape(-1)
        nonzeroIndex = np.where( rowSum > 0 )[0]
	distanceRemovedCM = (distanceRemovedCM [nonzeroIndex,: ])[:,nonzeroIndex]
	starts = (np.array(starts)[nonzeroIndex]).tolist()

	subInd = np.array( [ int( start / resolution ) for start in starts ] )
	print 'len of subind ' + str(len(subInd))

	if corrMethod == 'pearson' :
		corrM = np.corrcoef(distanceRemovedCM)	
	elif corrMethod == 'spearman' :
		corrM, pval = spearmanr(distanceRemovedCM)
	else :
		print 'Invalid Corr Method'
		exit(0)
	(w,v) = np.linalg.eig(corrM)
	fullEig = np.zeros(N)
	
	#print fullEig[subInd]
	#print v[:,EigI]
	fullEig[subInd] = v[:,EigI] 



	if validationV  is None :
		Flip = 1;
	else :
		if np.corrcoef(fullEig,validationV)[0,1] < 0 :
			Flip = -1
		else :
			Flip = 1
		print 'Corrcoef ' + str(np.corrcoef(fullEig,validationV)[0,1])

	fullEig = fullEig * Flip
		
	fullCorrM = np.matrix(np.zeros( (N,N) ) )
	for i in range(np.shape(corrM)[0]) :
		for k in range(np.shape(corrM)[1]) :
			fullCorrM[subInd[i],subInd[k]] = corrM[i,k]

	#np.save('matrix.npy',fullCorrM)
	#np.save('vector.npy',fullEig)
	if returnDist :
		return (fullCorrM, fullEig, distanceRemovedCM )
	else :
		return (fullCorrM, fullEig )

	
def plotEigenVector ( filename, chr = '1' ) :

	(m,labels,lenlabels)=ParseDenseMatrix(filename)
	labels = np.array(labels)
	(nm,Ind)=iceMatrix(m,percentile=5)
	ice_labels = labels[Ind]
	ice_lenlabels = lenlabels[Ind]

	prob = inferContactProbability ( m, labels )
	em = generateExpected( m, labels, chr, prob )

	Ind = np.where( ice_labels == chr )[0]
	iced_chrm = nm[Ind,:]; iced_chrm = iced_chrm[:,Ind]

	Ind = ice_lenlabels[Ind]
	ice_em = em[Ind,:]; ice_em = ice_em[:,Ind]

	final_m = iced_chrm / ice_em
	for i in range(len(final_m)) :
        	final_m[i,i] = 10**-12

	cr_final_m =  np.corrcoef( final_m )
	(w,v) = np.linalg.eig(cr_final_m)
	

def  getGCcontent (GCcontentfile,resolution = 1000000, assembly = 'hg19' ) :
	GCarray = np.loadtxt(GCcontentfile,dtype={'names':('chr','start','end','AT','GC','x1','x2','x3','x4','x5'), 'formats': ('S5','i4','i4','f4','f4','f4','f4','f4','f4','f4') } )
	chrs = getChrs(assembly)

	GCcontentHash = {}
	for chr in chrs:
        	featureList = np.zeros( int( np.ceil( getChrLen(chr,assembly=assembly)/float(resolution) ) ) )
        	for i in range(len(GCarray)) :
                	thisChr = GCarray[i][0]; thisEnd = GCarray[i][2]; thisFeature = GCarray[i][4]
                	if thisChr == chr :
                        	featureList[ int( np.ceil( thisEnd/float(resolution) ) ) - 1 ] = thisFeature
        	GCcontentHash[chr] = featureList
	
	return GCcontentHash

def  getGeneCounts (Genecountfile,assembly ='hg19') :
	GCarray = np.loadtxt(Genecountfile,dtype={'names':('chr','start','end','gc'), 'formats': ('S5','i4','i4','i4') } )
	chrs = getChrs(assembly)

	GCHash = {}
	for chr in chrs:
        	Ind = np.where(GCarray['chr'] == chr )
        	GCHash[chr] = GCarray['gc'][Ind]
	
	return GCHash

def contactMatrixSpearman ( CM1, CM2, labels1 = None,labels2 = None ) :
	if labels1 is not None and labels2 is not None:
		np.labels1 = np.array(labels1); np.labels2 = np.array(labels2)
		intlabels = np.intersect1d(np.labels1,np.labels2)

		map1  = np.sort ( [ np.where( np.labels1 == query )[0][0] for query in intlabels  ] )
		map2  = np.sort ( [ np.where( np.labels2 == query )[0][0] for query in intlabels  ] )
		CM1 = CM1 [map1,:]; CM1 = CM1[:,map1 ]
		CM2 = CM2 [map2,:]; CM2 = CM2[:,map2 ]
	CMV1 = np.asarray(CM1).reshape(-1);
	CMV2 = np.asarray(CM2).reshape(-1);

	return spearmanr(CMV1,CMV2)[0]
	
def makeBedgraphFromVector( valueV, chr, resolution ) :
        chrLen = getChrLen(chr)
        bedgraph = []
        for i in range(chrLen/resolution+1) :
                s=i*resolution+1
                e=(i+1)*resolution
                if ( e > chrLen ) :
                        e = chrLen
                bedgraph.append( ( chr,s,e,valueV[i] ) )

        return(bedgraph)

def printBedgraph ( bedgraph, outfilename, trackname ='', printHeader = False ) :
        out = open(outfilename,'a')
        if trackname == '' :
                trackname = outfilename
        if printHeader :
                print >> out, "track type=bedGraph name=%s" % (trackname)
        for i in range(len(bedgraph)) :
                print >> out, "%s\t%d\t%d\t%f" % (bedgraph[i][0],bedgraph[i][1],bedgraph[i][2],bedgraph[i][3])
	
def RebinMatrix ( CM, orgRes, newRes ) :
	if newRes % orgRes != 0 :
		print 'New Res is not compatible with original Res'
		return
	newBinCount = int(np.ceil(float( np.shape(CM)[0] * orgRes ) / newRes))
	RCM = np.zeros((newBinCount,newBinCount))
	for i in range(np.shape(CM)[0]) :
		for k in range(i,np.shape(CM)[0]) :
			newI = int(np.ceil( float( (i+1)*orgRes ) / newRes ) - 1)
			newK = int(np.ceil( float( (k+1)*orgRes ) / newRes ) - 1)
			#print str(i) + ' ' + str(k) + ' ' + str(newI) + ' ' + str(newK)
			RCM[newI,newK] = RCM[newI,newK] + CM[i,k]
	RCM = RCM + np.transpose( np.triu(RCM,1) )
	return(RCM)


