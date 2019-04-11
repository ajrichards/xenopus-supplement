#!/usr/bin/python

import os,csv,cPickle
from htsint.database import fetch_annotations, db_connect


session,engine = db_connect()
homeDir = os.path.join(os.path.expanduser("~"),"sequencing","xenopus")
assembly = 'dn'

if assembly == 'ref':
    sourceDir = 'reference'
else:
    sourceDir = "%s-trinity"%(assembly)

## load feature data 
featuresDir = os.path.join(homeDir,"features")
inFilePath = os.path.join(featuresDir,"de-summary-%s.csv"%(tissue))
fidin = open(inFilePath,'r')
reader = csv.reader(fidin)
header = reader.next()

total = 0
hasBlast = 0 
hasGene = 0
geneList = set([])

for linja in reader:
    transcriptID = linja[0]
    hitId = linja[1]
    geneId = linja[2]
    hitSpecies = linja[3]
    eValue =  linja[4]
    pval = linja[5]
    adjPval = linja[6]

    total += 1 
    if hitId != '-':
        hasBlast += 1

    if geneId not in ['-','']:
        hasGene += 1
        geneList.update([geneId])
geneList = list(geneList)

## annotations
def save_annotations(annotations,fileName):
    go2gene = {}
    gene2go = {}
    for gene,terms in annotations.iteritems():
        gene2go[gene] = terms
        for term in terms:
            if go2gene.has_key(term) == False:
                go2gene[term] = set([])
            go2gene[term].update([gene])

    for term,genes in go2gene.iteritems():
        go2gene[term] = list(genes)
    
    tmp = open(fileName,'w')
    cPickle.dump([gene2go,go2gene],tmp)
    tmp.close()
    return gene2go, go2gene

annots1 = 0
print("fetching BP annotations for uniprot id...")
fileName1 = os.path.join(homeDir,'filtered','annotations-%s-%s.pickle'%(tissue,'noiea-bp'))
if os.path.exists(fileName1):
    tmp = open(fileName1,'r')
    gene2go1,go2gene1 = cPickle.load(tmp)
else:
    annotations1 = fetch_annotations(geneList,engine,idType='ncbi',useIea=False,verbose=True,aspect='biological_process')
    gene2go1, go2gene1 = save_annotations(annotations1,fileName1)

for key,item in gene2go1.iteritems():
    if len(item) > 0:
        annots1+=1

annots2 = 0
print("fetching MF annotations for uniprot id...")
fileName2 = os.path.join(homeDir,'filtered','annotations-%s-%s.pickle'%(tissue,'noiea-mf'))
if os.path.exists(fileName2):
    tmp = open(fileName2,'r')
    gene2go2,go2gene2 = cPickle.load(tmp)
else:
    annotations2 = fetch_annotations(geneList,engine,idType='ncbi',useIea=False,verbose=True,aspect='molecular_function')
    gene2go2,go2gene2 = save_annotations(annotations2,fileName2)

for key,item in gene2go2.iteritems():
    if len(item) > 0:
        annots2+=1

annots3 = 0
print("fetching MF annotations for uniprot id...")
fileName3 = os.path.join(homeDir,'filtered','annotations-%s-%s.pickle'%(tissue,'iea-bp'))
if os.path.exists(fileName3):
    tmp = open(fileName3,'r')
    gene2go3,go2gene3 = cPickle.load(tmp)
else:
    annotations3 = fetch_annotations(geneList,engine,idType='ncbi',useIea=True,verbose=True,aspect='biological_process')
    gene2go3,go2gene3 = save_annotations(annotations3,fileName3)

for key,item in gene2go3.iteritems():
    if len(item) > 0:
        annots3+=1

annots4 = 0
print("fetching MF annotations for uniprot id...")
fileName4 = os.path.join(homeDir,'filtered','annotations-%s-%s.pickle'%(tissue,'iea-mf'))
if os.path.exists(fileName4):
    tmp = open(fileName4,'r')
    gene2go4,go2gene4 = cPickle.load(tmp)
else:
    annotations4 = fetch_annotations(geneList,engine,idType='ncbi',useIea=True,verbose=True,aspect='molecular_function')
    gene2go4,go2gene4 = save_annotations(annotations4,fileName4)

for key,item in gene2go4.iteritems():
    if len(item) > 0:
        annots4+=1

print("-----------------------")
print(inFilePath)
print("has blast match: %s/%s"%(hasBlast,total))
print("has mapped gene: %s/%s"%(hasGene,total))
print("unique genes: %s"%(len(geneList)))
print("annotated with BP noiea: %s/%s"%(annots1,len(gene2go1.keys())))
print("annotated with MF noiea: %s/%s"%(annots2,len(gene2go2.keys())))
print("annotated with BP iea: %s/%s"%(annots3,len(gene2go3.keys())))
print("annotated with MF iea: %s/%s"%(annots4,len(gene2go4.keys())))
