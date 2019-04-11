#!/usr/bin/env python
"""
read the output of a DESeq2 analysis and summarize the results
create a new output file with associated gene names and go terms
"""

import os,sys,csv,getopt,re
import numpy as np
from htsint.database import db_connect,fetch_annotations,Gene,GoTerm,Taxon
from htsint.tools import read_matrix,read_de_results
from htsint.blast import get_blast_map
from sqlalchemy.sql import select
from htsint.database import db_connect,Taxon,Gene

session,engine = db_connect()

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

homeDir = os.path.join(os.path.expanduser("~"),"sequencing","xenopus")

if transcript not in ['genes', 'isoforms']:
    raise Exception('Invalid transcript argument')
if assembly not in ['gg','dn']:
    raise Exception('Invalid assembly argument')

homeDir = os.path.join(os.path.expanduser("~"),"sequencing","xenopus")
featuresDir = os.path.join(homeDir,"%s-trinity"%(assembly),"features")

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

evalue = 0.00001

edgerResultsPath = os.path.join(featuresDir,"edger_%s_behavior_de.csv"%(transcript))
deseqResultsPath = os.path.join(featuresDir,"deseq_%s_behavior_de.csv"%(transcript))
edgerMatIds, edgerMatColumns, edgerMat = read_de_results(edgerResultsPath,tool='edgeR')
deseqMatIds, deseqMatColumns, deseqMat = read_de_results(deseqResultsPath,tool='DESeq')

## create a summary table and a csv file 
outPath = os.path.join(featuresDir,'de-summary-%s-%s_behavior.csv'%(assembly,transcript))
fidout = open(outPath,'w')
writer = csv.writer(fidout)
writer.writerow(['transcript-id','ncbi-id','description','edgeR-FDR','edgeR-pvalue','DESeq-adjusted-pvalue','DESeq-pvalue'])

## prepare database connections
session,engine = db_connect()
conn = engine.connect()
s = select([Taxon.id,Taxon.ncbi_id,Taxon.name]).where(Taxon.ncbi_id.in_(['8364']))
_taxaQueries = conn.execute(s)
taxaQueries = _taxaQueries.fetchall()
totalQueries = set([])
filteredQueries = set([])
filteredHits = set([])
selectedTaxa = [str(tquery['id']) for tquery in taxaQueries]
taxa2name = dict([(str(tquery['id']),str(tquery['ncbi_id'])) for tquery in taxaQueries])

## create a gene2taxa dictionary
gene2taxa,gene2desc,gene2sym = {},{},{}
for tquery in taxaQueries:
    s = select([Gene.taxa_id,Gene.ncbi_id,Gene.description,Gene.symbol],Gene.taxa_id==tquery['id'])
    _geneQueries = conn.execute(s)
    geneQueries = _geneQueries.fetchall()
    gene2taxa.update(dict([(str(r['ncbi_id']),str(r['taxa_id'])) for r in geneQueries]))
    gene2desc.update(dict([(str(r['ncbi_id']),str(r['description'])) for r in geneQueries]))
    gene2sym.update(dict([(str(r['ncbi_id']),str(r['symbol'])) for r in geneQueries]))

if transcript == 'genes':
    blastMap = get_blast_map(os.path.join("..","..","blast","summary_blast_%s.csv"%assembly),\
                             evalue=evalue,taxaList=["8364"],asGenes=True)
else:
    blastMap = get_blast_map(os.path.join("..","..","blast","summary_blast_%s.csv"%assembly),\
                             evalue=evalue,taxaList=["8364"],asGenes=False)

## prepare supplment output
columns = [18,22,66,11,11]
row = "+"
head = "+"
for col in columns:
    row += "-"*col+"-+"
    head += "="*col+"=+"

print("\nGene sets\n_____________________")
print(row)
items = ['transcript ID','gene symbol','gene description','edgeR','DESeq']
show_contents(columns,items,withTrailing=False)
print(head)

## loop through transcripts based on FDR
rankedInds = np.argsort(edgerMat[:,3])
sphinxLinks = []
for rind in rankedInds:
    erAdjPvalue = edgerMat[rind,np.where(edgerMatColumns=='FDR')[0]][0]
    erPvalue = edgerMat[rind,np.where(edgerMatColumns=='PValue')[0]][0]
    transcriptId = edgerMatIds[rind]

    if blastMap.has_key(transcriptId):
        ncbiId = blastMap[transcriptId][0]
        symbol = gene2sym[ncbiId]
        ncbiUrl = "http://www.ncbi.nlm.nih.gov/gene/%s"%(ncbiId)
        ncbiLink =  "`%s`_"%symbol
        desc = gene2desc[ncbiId]
        desc = re.sub('"','',desc)        
    else:
        ncbiId,ncbiLink,desc,symbol = "-","-","-","-"

    deseqInd = np.where(deseqMatIds==transcriptId)[0]
    deseqFDR = deseqMat[deseqInd,np.where(deseqMatColumns=='padj')[0]][0]
    deseqPval = deseqMat[deseqInd,np.where(deseqMatColumns=='pvalue')[0]][0]
    writer.writerow([transcriptId,ncbiId,desc,erAdjPvalue,erPvalue,deseqFDR,deseqPval])

    ## print info for supplement
    if erAdjPvalue <= 0.1 or deseqFDR <= 0.1:
        if len(desc) >= 60:
            desc = desc[:60] + "..."

        if desc != '-':
            sphinxLinks.append(".. _%s: %s"%(symbol,ncbiUrl))
        items = [transcriptId,ncbiLink,desc,round(erAdjPvalue,7),round(deseqFDR,7)]
        show_contents(columns,items)
        
for sl in sphinxLinks:
    print sl
