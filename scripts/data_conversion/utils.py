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

