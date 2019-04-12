#!/usr/bin/python
"""
run after assembleDistances.py
 
example:

 ~$ python runSpectralClustering.py -n drerio -a bp

"""

__author__ = "Adam Richards"

import sys,os,getopt,csv,cPickle,shutil
import matplotlib.pyplot as plt

import numpy as np
from htsint.stats import SpectralCluster
from htsint.database import GoTerm,db_connect

session,engine = db_connect()

## read in input file 
if len(sys.argv) < 3:
    print sys.argv[0] + "-n name -a aspect"
    sys.exit()

try:
    optlist, args = getopt.getopt(sys.argv[1:], 'a:n:')
except getopt.GetoptError:
    raise Exception(sys.argv[0] + "-n name -a aspect")
    sys.exit()

name,aspect = None,None
for o,a in optlist:
    if o == '-n':
        name = a 
    if o == '-a':
        aspect = a.lower()

if aspect not in ['bp','mf','cc']:
    raise Exception("Not a valid Gene Ontology aspect (-a [bp|mf|cc]")


def read_results(labelsPath):
    """
    read the results into an array
    """

    inFid = open(labelsPath,'r')
    reader = csv.reader(inFid)
    header = reader.next()
    genes,clusters = [],[]
    for linja in reader:
        genes.append(linja[0])
        clusters.append(int(linja[1]))

    inFid.close()

    return {'genes':np.array(genes),'clusters':np.array(clusters)}

def get_cluster_sizes(clusters):
    clusterSizes = []
    for k in np.arange(clusters.max()+1):
        clusterSizes.append(len(np.where(clusters==k)[0]))

    return np.array(clusterSizes)

## make sure we can find the labels
labelsPath = os.path.join("..","..","results","sc-labels-%s-%s.csv"%(name,aspect))
reclusteredPath = os.path.join("..","results","sc-labels-%s-%s-reclustered.csv"%(name,aspect))
if os.path.exists(labelsPath) == False:
    raise Exception("Cannot find labels")

clusterResults = read_results(labelsPath)
fig = plt.figure()
if os.path.exists(reclusteredPath):
    allResults = [read_results(rPath) for rPath in [labelsPath,reclusteredPath]]
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
else:
    allResults = [read_results(rPath) for rPath in [labelsPath]]
    ax1 = fig.add_subplot(111)
    ax2 = None

axisList = [ax1,ax2]
width = 0.6
BLUE = "#0066FF"
ORANGE = "#FF6600"
GREY = "#777777"
fontSize = 11
fontName = 'sans'

for rind,clusterResults in enumerate(allResults):
    genes = clusterResults['genes']
    clusters = clusterResults['clusters']
    clusterSizes = get_cluster_sizes(clusters)
    tooSmall = np.where(clusterSizes <= 2)[0]
    tooLarge = np.where(clusterSizes > 100)[0]
    excludedClusters = np.unique(np.hstack([tooSmall,tooLarge]))
    excludedGenes = clusterSizes[excludedClusters].sum()

    print("total gene: %s"%genes.size)
    print("total clusters: %s"%np.unique(clusters).size)
    print("total clusters <= 2: %s"%tooSmall.size) 
    print("total clusters > 100: %s"%tooLarge.size) 
    print("percent missing: %s/%s (%s)"%(excludedGenes,genes.size,round(float(excludedGenes)/genes.size*100,2)))

    ax = axisList[rind] 
    N = len(clusterSizes)
    indx = np.arange(N)
    rects = ax.bar(indx,clusterSizes,width,color=BLUE)

    ticks = np.array([int(round(i)) for i in np.linspace(0,len(indx),15)])
    ax.set_xticks(ticks+0.5*width)
    ax.set_xticklabels(ticks)

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(fontSize-2)
        tick.label.set_fontname(fontName)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(fontSize-2)
        tick.label.set_fontname(fontName)
    ax.set_ylabel("# Genes",fontsize=fontSize,fontname=fontName)
    ax.set_xlabel("Cluster ID",fontsize=fontSize,fontname=fontName)

    #leg1 = ax1.legend([rects1[0]], ['Original clustering'],loc='upper left')
    #ltext = leg1.get_texts()
    #for i in range(len(ltext)):
    #    plt.setp(ltext[i],fontsize=fontSize)

plt.savefig("sc-cluster-sizes-%s-%s.png"%(name,aspect),dpi=300)

## create a summary table and a csv file
columns = [75,16,16]
row = "+"
head = "+"
for col in columns:
    row += "-"*col+"-+"
    head += "="*col+"=+"

def show_contents(columns,items,withTrailing=True):
    if len(columns) != len(items):
        raise Exception("Dimension mismatch")
        
    toPrint = "| "

    for i,item in enumerate(items):
        item = str(item)

        if len(item) >= columns[i]:
            raise Exception("col %s not large enough min = %s"%(i,len(item)+2))

        toPrint += item+" "*(columns[i]-len(item))+"| "

    print(toPrint[:-1])
    
    if withTrailing:
        print(row)

print("\nGene sets\n_____________________")
print(row)
items = ['Cluster ID','Num. Genes','Unique Terms']
show_contents(columns,items,withTrailing=False)
print(head)

termsPath = os.path.join("..","..","results","go-terms-%s-%s.pickle"%(name,aspect))
summaryPath = os.path.join("summary-%s-%s.csv"%(name,aspect))
tmp = open(termsPath,'r')
gene2go,go2gene = cPickle.load(tmp)
tmp.close()

fidOut = open(summaryPath,'w')
writer = csv.writer(fidOut)
writer.writerow(["clusterId","clusterGenes","uniqueTerms","someTerms"])

## clean the results dir
resultsDir = os.path.realpath(os.path.join(".","%s-%s"%(name,aspect)))
if os.path.isdir(resultsDir):
    shutil.rmtree(resultsDir)
os.mkdir(resultsDir)

for clusterId in np.sort(np.unique(clusters)):
    clusterInds = np.where(clusters== clusterId)[0]
    clusterGenes = genes[clusterInds]

    if clusterGenes.size > 2:
        out = open(os.path.join(resultsDir,"%s-%s-%s.csv"%(name,aspect,clusterId)),'w')
        writeout = csv.writer(out)
        writeout.writerow(["geneId","terms"])

    ## get go terms for cluster 
    cterms = []
    for cgene in clusterGenes:
        cterms.extend(gene2go[cgene])
        if clusterGenes.size > 2:
            terms = [session.query(GoTerm).filter(GoTerm.go_id == key).first().name + " (%s)"%key for key in gene2go[cgene]]
            terms = ';'.join(map(str, terms))
            writeout.writerow([cgene,terms])

    ## rank terms in order of usage
    cterms = np.array(cterms)
    uniqueTerms = np.sort(np.unique(cterms))
    termCounts = [np.sum(cterms == ct) for ct in uniqueTerms]
    topFew = uniqueTerms[np.argsort(termCounts)[::-1]][:3]
    tfTerms = [session.query(GoTerm).filter(GoTerm.go_id == key).first().name for key in topFew]
    tfTerms = ';'.join(map(str, tfTerms))
    #topTerm = uniqueTerms[np.argsort(termCounts)[::-1]][:1]
    #tTerm = [session.query(GoTerm).filter(GoTerm.go_id == key).first().name for key in topTerm]
    #tTerm = ''.join(map(str, tTerm))

    ## truncate term if needed
    #if len(tTerm) > 70:
    #    tTerm = tTerm[:65] + "..."

    writer.writerow([clusterId,clusterGenes.size,uniqueTerms.size,tfTerms])
    if clusterGenes.size > 2:
        #print("%s,%s,%s,%s"%(clusterId,clusterGenes.size,uniqueTerms.size,tfTerms))
        linkedId = ":download:`%s <../../gene-sets/xeno-noiea-bp/xeno-noiea-bp-%s.csv>`"%(clusterId,clusterId)
        items = [linkedId,clusterGenes.size,uniqueTerms.size]
        show_contents(columns,items)
    

fidOut.close()
