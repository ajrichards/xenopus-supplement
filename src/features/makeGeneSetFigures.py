#!/usr/bin/env python

import os,sys,csv,re
import numpy as np
from htsint import GeneSet

## find a significant gene set that shares a bunch of genes
resultsFile = "../results/x-hs-noiea-bp-geneset-results.csv"
genesetFile = "../results/x-hs-noiea-bp-gg-genes.csv"
termsPath = "../results/go-terms-x-hs-noiea-bp.pickle"
termDistancesPath = os.path.join("..","results","distances-x-hs-noiea-bp.npy")

def read_results(resultsFile):
    """
    read a output from GSA file
    """
    
    fid = open(resultsFile,'r')
    reader = csv.reader(fid)
    header = reader.next()
    results = {}
    for item in header:
        results[item] = []

    for linja in reader:
        if linja[3] == '':
            continue

        for r,row in enumerate(linja):
            results[header[r]].append(row)

    for item in header:
        results[item] = np.array(results[item])

    fid.close()
    return results

results = read_results(resultsFile)
print results.keys()

## load distances
distMat = np.load(termDistancesPath)

## initialize the gene set
## layout: sprint, spectral
gsets = GeneSet(genesetFile,termsPath,distMat)
print gsets.geneset2gene.keys()

## select a gene set
#for g,gSet in enumerate(results['gene_set']):
#    print g,gSet
#    if len(gsets.geneset2gene[gSet]) > 2:
#        print gSet, len(gsets.geneset2gene[gSet]),len(gsets.geneset2transcript[gSet]),results['transcript'][g]

geneSet = 'gs-192'
#geneSet = 'gs-9-42'
gsets.draw_figure(geneSet,layout='spring',name='example.png',percentile=25)
#gsets.draw_figure(geneSet,layout='spring',name='Figure_Htsint_example.pdf',percentile=25)




#print geneInfo
#print gsets.gene2go['494851']

#print gsets.get_go_terms('494851')



print 'good'
