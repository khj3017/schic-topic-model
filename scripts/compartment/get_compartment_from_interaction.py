from utils import *
import matplotlib.pyplot as plt


def makeBedgraphFromVector( valueV, chr, resolution, assembly) :
        chrLen = getChrLen(chr, assembly = assembly)
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


interactionFile=sys.argv[1]
justFileName=interactionFile.split('/')[-1]
resolution=int(sys.argv[2])
outDir = sys.argv[3]

if not os.path.exists(outDir):
        os.makedirs(outDir)

print 'Parsing'

assembly = 'hg19'

( Mlist, labelList ) = ParseWholeInteractionBedByMatrix ( interactionFile, resolution = resolution, assembly = assembly )

if len(sys.argv) == 5 :
	contactProbFile=sys.argv[4]
	prob = np.loadtxt(contactProbFile)
else :
	print 'Calculating prob'
	prob = inferContactProbabilityFromInteraction ( Mlist, labelList, resolution = resolution, returnProb = True ) 

GCcontent = getGeneCounts('/net/noble/vol2/home/gurkan/proj/2015blood-differentiation/Visualizations/Data/GeneCountBeds/hg19.' + str(resolution/1000) + 'kb.Windows.PC.GeneCounts.bedgraph')
chrs = getChrs(assembly)
chrs = chrs[:22]

plotMatrices = []; plotEigens = []; bedGraphs = []

for i in range(len(chrs)) :
	chrCM = Mlist[i]
	chrLabels = labelList[i]
	chr = chrs[i]
	binStarts = np.array([ label.split(':')[1].split('-')[0] for label in chrLabels])
	( chrCM, chrLabels )= iceMatrix (chrCM,chrLabels,  percentile = 1 )
	( CorrMatrix, Eig1 ) = decomposeContactMatrix(chrCM,chrLabels,prob,resolution=resolution, chr = chr, validationV = GCcontent[chr] , EigI = 1, assembly = assembly)
	plotMatrices.append(CorrMatrix); plotEigens.append(Eig1)
	chrLen = getChrLen(chr, assembly = assembly)
	bedGraphs.append( makeBedgraphFromVector( Eig1, chr, resolution, assembly) )

plot=True
if plot :
	for  i in range(len(plotMatrices)) :
		plt.subplots(figsize=(6,9))
		ax3 = plt.subplot2grid((9,6),(0,0),rowspan=1,colspan=6)
		ax1 = plt.subplot2grid((9,6),(1,0),rowspan=6,colspan=6)
		ax2 = plt.subplot2grid((9,6),(7,0),rowspan=2,colspan=6)
		ref=plotHeatmap(plotMatrices[i],perc=100,axis=ax1 )
		coord = np.arange(0,len(plotEigens[i])*resolution,resolution)
		pos_idx = np.where(plotEigens[i] > 0)
		neg_idx = np.where(plotEigens[i] < 0)
		comp = coord[neg_idx]; comp2 = coord[pos_idx]
		colors1 = np.array([[1,0,0], [0,0,1]])
		ax2.bar( np.arange(0,len(plotEigens[i])*resolution,resolution), plotEigens[i], width=resolution  )
		ax2.set_xlim([0,len(plotEigens[i])*resolution]); 
		plt.colorbar(ref,cax=ax3,orientation='horizontal')
		plt.tight_layout()
		plt.savefig(os.path.join(outDir, chrs[i]+'.'+str(resolution)+'.compartments.png'))
		plt.close()

for i in range(len(bedGraphs)) :
	printBedgraph(bedGraphs[i], os.path.join(outDir, justFileName[:-7] + 'ICED.compartments.bedgraph'), trackname=justFileName[:-8],printHeader=(i==0) )

				
