#!/usr/bin/env python
"""
read the output of a DESeq2 analysis and summarize the results
create a new output file with associated gene names and go terms
"""

import os,sys,csv,cPickle,re
from htsint.database import db_connect,fetch_annotations,Gene,GoTerm,Taxon
session,engine = db_connect()

## functions
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

## read in jgi2gene 
fidin = open(os.path.join("..","..","data","jgi2gen.pickle"),'r')
jgi2gene = cPickle.load(fidin)
fidin.close()

## connect with results file
inPath = os.path.join(os.path.expanduser("~"),'documents','hts-integrate-docs','xenopus','data','deseq.csv')
if not os.path.exists(inPath):
    raise Exception("File path does not exist\n%s\nDid you run `Rscript standard-difexp.R`?"%inPath)

fidin = open(inPath,'rU')
reader = csv.reader(fidin)
header = reader.next()
print header

## create a summary table and a csv file 
outPath = os.path.join(os.path.expanduser("~"),'documents','hts-integrate-docs','xenopus','data','deseq-summary.csv')
fidout = open(outPath,'w')
writer = csv.writer(fidout)
writer.writerow(['transcript ID','mapped gene id','gene description','adjusted-pvalue'])

columns = [16,10,75,60]
row = "+"
head = "+"
for col in columns:
    row += "-"*col+"-+"
    head += "="*col+"=+"

print("\nGene sets\n_____________________")
print(row)
items = ['transcript ID','symbol','gene description','link']
show_contents(columns,items,withTrailing=False)
print(head)

## get xenopus trop. database id
query = session.query(Taxon).filter_by(ncbi_id='8364').first()
xenotropId = query.id

for linja in reader:
    if "NA" in linja:
        continue
    transcript = linja[0]
    baseMean = float(linja[1])
    log2FoldChange = float(linja[2])
    lfcSE = float(linja[3])
    stat = float(linja[4])
    pvalue = float(linja[5])
    padj = float(linja[6])

    ## get mapped name
    if jgi2gene.has_key(transcript):
        mappedGene = jgi2gene[transcript]
        mappedGeneId,mappedGeneDesc = mappedGene['ncbi'],mappedGene['desc']
    else:
        mappedGeneId,mappedGeneDesc,annotations = "NA","NA","NA"

    ## get terms
    #if mappedGeneId != "NA":
        #    geneQuery = session.query(Gene).filter_by(symbol=mappedGeneId).first()
        #    annotations = fetch_annotations([geneQuery.ncbi_id],session,idType='ncbi',useIea=True)#asTerms=False)
        #    associatedTerms = annotations[geneQuery.ncbi_id]
        
    ## write sig results output
    if padj <= 0.1:
        
        ## see if we can create a link
        geneQuery = session.query(Gene).filter_by(symbol=mappedGeneId).all()
        ncbiGeneId = None
        for gq in geneQuery:
            if gq.taxa_id == xenotropId:
                ncbiGeneId = gq.ncbi_id
                break

        if ncbiGeneId:
            ncbiUrl = "http://www.ncbi.nlm.nih.gov/gene/%s"%ncbiGeneId
        else:
            ncbiUrl = ""

        writer.writerow([transcript,mappedGeneId,mappedGeneDesc,padj])
        if ncbiGeneId:
            if len(mappedGeneDesc) > 60:
                mappedGeneDesc = mappedGeneDesc[:57] + "..."
            items = [transcript,mappedGeneId,mappedGeneDesc," `%s <%s>`_"%(ncbiGeneId,ncbiUrl)]
        else:
            items = [transcript,mappedGeneId,mappedGeneDesc]
        show_contents(columns,items)

fidout.close()
fidin.close()
print("done")
