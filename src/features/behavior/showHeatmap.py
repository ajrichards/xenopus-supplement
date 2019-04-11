#!/usr/bin/env python
"""
create the manuscript figure summarizing features

volcano plot - showing sig genes/isoforms
clustering heatmap - showing how genes/isoforms cluster


"""

import sys,os,getopt,re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from htsint.tools import read_matrix,read_de_results,Heatmap
from htsint.blast import get_blast_map
from htsint.database import db_connect, Taxon, Gene
from sqlalchemy.sql import select

## read input                                                                         
if len(sys.argv) < 2:
    raise Exception(sys.argv[0] + " -a assembly [-a] is 'dn' or 'gg' and transcript [-t] is 'genes' or 'isoforms'")
try:
    optlist, args = getopt.getopt(sys.argv[1:], 'a:t:')
except getopt.GetoptError:
    raise Exception(sys.argv[0] + " -a assembly [-a] is 'dn' or 'gg' and transcript [-t] is 'genes' or 'isoforms'")

assembly,transcript = None,None
for o,a in optlist:
    if o == '-a':
        assembly = a
    if o == '-t':
        transcript = a

#assembly = 'gg'
#transcript = 'genes'
homeDir = os.path.join(os.path.expanduser("~"),"sequencing","xenopus")

if transcript not in ['genes', 'isoforms']:
    raise Exception('Invalid transcript argument')
if assembly not in ['gg','dn']:
    raise Exception('Invalid assembly argument')

## variables
fontSize = 11
fontType = 'arial'
useColor = True
if useColor == True:
    myCmap = mpl.cm.gist_heat
else:
    myCmap = mpl.cm.gray

## get some dictionaries form htsint 
session,engine = db_connect()
conn = engine.connect()
s = select([Taxon.id,Taxon.ncbi_id,Taxon.name]).where(Taxon.ncbi_id.in_(['8364']))
_taxaQueries = conn.execute(s)
taxaQueries = _taxaQueries.fetchall()
gene2taxa,gene2desc,gene2sym = {},{},{}
for tquery in taxaQueries:
    s = select([Gene.taxa_id,Gene.ncbi_id,Gene.description,Gene.symbol],Gene.taxa_id==tquery['id'])
    _geneQueries = conn.execute(s)
    geneQueries = _geneQueries.fetchall()
    gene2taxa.update(dict([(str(r['ncbi_id']),str(r['taxa_id'])) for r in geneQueries]))
    gene2desc.update(dict([(str(r['ncbi_id']),str(r['description'])) for r in geneQueries]))
    gene2sym.update(dict([(str(r['ncbi_id']),str(r['symbol'])) for r in geneQueries]))

## load feature data
featuresDir = os.path.join(homeDir,"%s-trinity"%assembly,"features")
edgerResultsPath = os.path.join(featuresDir,"edger_%s_behavior_de.csv"%(transcript))
edgerIds, edgerColumns, edgerMat = read_de_results(edgerResultsPath,tool='edgeR')
deseqResultsPath = os.path.join(featuresDir,"deseq_%s_behavior_de.csv"%(transcript))
deseqIds, deseqColumns, deseqMat = read_de_results(deseqResultsPath,tool='DESeq')
dfeMatrixPath = os.path.join(featuresDir,"deseq_%s_behavior_de_samples.csv"%(transcript))
dfeIds,dfeColumns,dfeMat = read_matrix(dfeMatrixPath,mtype='float')

## load the blast map 
if transcript == 'genes':
    blastMap = get_blast_map(os.path.join("..","..","blast","summary_blast_%s.csv"%assembly),\
                             taxaList=["8364"],asGenes=True)
else:
    blastMap = get_blast_map(os.path.join("..","..","blast","summary_blast_%s.csv"%assembly),\
                             taxaList=["8364"],asGenes=False)

## setup filters for the transcripts
threshold = 0.1
padjInd = np.where(deseqColumns == 'padj')[0]
filter1Inds = np.where(deseqMat[:,padjInd] < threshold)[0]
fdrInd = np.where(edgerColumns=='FDR')[0]
filter2Inds = np.where(edgerMat[:,fdrInd] < threshold)[0]
print('filtering.......')
print('filter1 (deseq): %s transcripts'%filter1Inds.size)
print('filter2 (edger): %s transcripts'%filter2Inds.size)
includedIds = np.array(list(set(deseqIds[filter1Inds]).union(set(edgerIds[filter2Inds]))))
print('union on ids: %s'%len(includedIds))

## create dfe matrix based on filtered genes
mask = np.in1d(dfeIds,includedIds)
mat = dfeMat[mask,:]
matIds = dfeIds[mask]

## map any ids to gene names 
mappedIds = []
for transcriptId in matIds:
    if blastMap.has_key(transcriptId):
        ncbiId = blastMap[transcriptId][0]
        symbol = gene2sym[ncbiId]
        #print transcriptId,symbol
        mappedIds.append(symbol)
    else:
        mappedIds.append(transcriptId)
mappedIds = np.array(mappedIds)

## create the heatmap
hm = Heatmap(colLabels=np.array(["A","B","C","D","E","F","G","H"]),\
             rowLabels=mappedIds,hpad=0.11)
hm.cluster(mat,0)
hm.cluster(mat,1)
hm.draw_heatmap(mat,cmap='uy',clabels=True,rlabels=True,rowFont=5)
hm.save("heatmap-%s_%s.png"%(assembly,transcript),dpi=400)

sortedSamples = hm.colLabels[hm.indx['1']]
sortedTranscripts = matIds[hm.indx['0']]
sortedRemapped = hm.rowLabels[hm.indx['0']]

## print the significant transcritps
for t,transcriptId in enumerate(sortedTranscripts):
    ind1 = np.where(deseqIds==transcriptId)[0]
    sig1 = deseqMat[ind1,padjInd][0]
    ind2 = np.where(edgerIds==transcriptId)[0]
    sig2 = edgerMat[ind2,fdrInd][0]
    if sig1 <= 0.1 or sig2 <= 0.1:
        print t,transcriptId,sortedRemapped[t],sig1,sig2
    #ind = np.where(dfeIds==transcript)[0]
    #print "... samples", dfeMat[ind,:][0].tolist()
    #print "... samples", dfeMat[ind,:][0][hm.indx['1']].tolist()

hm.show()





sys.exit()
