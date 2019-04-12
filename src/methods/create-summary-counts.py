#!/usr/bin/python
"""
read in the count file for each sample then create a summary csv file
"""

import re,os,sys,csv,cPickle
import numpy as np

## basic variables
seqDir = os.path.join(os.path.expanduser("~"),"sequencing","xenopus-data")
analysisDir = os.path.join(seqDir,"analyses")
sampleList = ["A","B","C","D","E","F","G","H"]

## read in jgi2gene
fidin = open(os.path.join("..","..","data","jgi2gen.pickle"),'r')
jgi2gene = cPickle.load(fidin)          
fidin.close()  

## read in the data and organize the counts
allCounts = {}
geneCounts = {}
uniqueGenes = set([])

for s,sample in enumerate(sampleList):
    countsFilePath = os.path.join(analysisDir,'counts_%s.txt'%sample.lower())
    fidin = open(countsFilePath,'rU')
    unmapped,total = 0,0

    for linja in fidin:

        linja = linja[:-1].split("\t")

        ## check for summary information
        if linja[0] in ["no_feature","ambiguous","too_low_aQual","not_aligned","alignment_not_unique"]:
            print(linja[0],linja[1])
            continue

        total += 1
        if not allCounts.has_key(linja[0]):
            allCounts[linja[0]] = [0]*8

        allCounts[linja[0]][s] = int(linja[1])

        if not jgi2gene.has_key(linja[0]):
            unmapped += 1
        else:
            mappedName = jgi2gene[linja[0]]['ncbi']
            uniqueGenes.update([mappedName])
            
            if not geneCounts.has_key(mappedName):
                geneCounts[mappedName] = [0]*8

            geneCounts[mappedName][s] = int(linja[1])

            ## debug
            if mappedName == 'apoe':
                print "...", "sample%s"%sample, linja[1]

    fidin.close()

    print("sample: %s"%sample)
    print("total transcripts: %s"%total)
    print("mapped transcripts: %s"%(total - unmapped))
    print("unique mapped transcripts %s"%len(uniqueGenes))

## write all transcripts to file
print("writing transcript-count-summary...")
fidout1 = open(os.path.join("..","..","data","transcript-count-summary.csv"),'w')
writer1 = csv.writer(fidout1)
writer1.writerow(["transcript","Sample-A","Sample-B","Sample-C","Sample-D","Sample-E","Sample-F","Sample-G","Sample-H"])

skipped = 0
for key in sorted(allCounts.keys()):
    ## debug
    if re.search("sample",key,flags=re.IGNORECASE):
        skipped += 1
        continue
    writer1.writerow([key]+allCounts[key])

fidout1.close()
print("skipped %s sample specific transcripts"%skipped)

## create a summary file of the counts with respect to the unique genes
print("writing gene-count-summary...")
fidout2 = open(os.path.join("..","..","data","gene-only-count-summary.csv"),'w')
writer2 = csv.writer(fidout2)
writer2.writerow(["geneid","Sample-A","Sample-B","Sample-C","Sample-D","Sample-E","Sample-F","Sample-G","Sample-H"])

for gene in sorted(list(uniqueGenes)):
    writer2.writerow([gene]+geneCounts[gene])

## wrapup
print("done.")
fidout2.close()
