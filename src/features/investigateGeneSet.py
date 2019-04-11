#!/usr/bin/env python
"""
create the manuscript figure summarizing features

volcano plot - showing sig genes/isoforms
clustering heatmap - showing how genes/isoforms cluster


"""

import sys,os,getopt,re,cPickle,csv
import numpy as np
from htsint import GeneSet
from htsint.blast import BlastMapper
from htsint.tools import read_de_results
from htsint.database import db_connect,GoTerm

session, engine = db_connect()
homeDir = os.path.join(os.path.expanduser("~"),"sequencing","xenopus")

def write_summary(name,aspect,transcript,assembly,geneset):
    ## load the go dictionaries 
    termsPath = os.path.join("..","results","go-terms-%s-%s.pickle"%(name,aspect))
    tmp = open(termsPath,'r')
    gene2go,go2gene = cPickle.load(tmp)
    tmp.close()

    ## load the blast map
    bm = BlastMapper()
    homeDir = os.path.join(os.path.expanduser("~"),"sequencing","xenopus")
    sizeMin, sizeMax = 5,100

    summaryFile = os.path.join(homeDir,"%s-trinity"%(assembly),'blast-%s-parsed_summary.csv'%assembly)
    if transcript == 'genes':
        bmap = bm.load_summary(summaryFile,trinityGene=True,best=False,taxaList=['8364','8355','9606'],evalue=0.0001)
    else:
        bmap = bm.load_summary(summaryFile,trinityGene=False,best=False,taxaList=['8364','8355','9606'],evalue=0.0001)

    ## get gene level differencial exp results
    featuresDir = os.path.join(homeDir,"%s-trinity"%assembly,"features")
    deseqResultsPath = os.path.join(featuresDir,"deseq_%s_de.csv"%(transcript))
    deseqIds, deseqColumns, deseqMat = read_de_results(deseqResultsPath,tool='DESeq')
    padjInd = np.where(deseqColumns == 'padj')[0]
    pvalInd = np.where(deseqColumns == 'pvalue')[0]

    ## input/output
    genesetSummaryFile = os.path.join("..","results","genesets",
                                      "%s-%s-%s-%s-%s.csv"%(name,aspect,transcript,assembly,re.sub("gs-","",geneset)))
    genesetFile = os.path.join("..","results","%s-%s-%s-%s.gmt"%(name,aspect,assembly,transcript))
    
    if not os.path.exists(genesetFile):
        raise Exception("cannot find gene set file")

    allGenesets = {}
    fid = open(genesetFile,'r')
    for linja in fid:
        linja = [re.sub("\s+","",l) for l in linja.split("\t")]
        allGenesets[linja[0]] = linja[2:]

    fid.close()

    gsTranscripts = allGenesets[geneset]


    ## map back to gene space and collect go terms
    transcript2genes = {}
    for t in gsTranscripts:
        transcript2genes[t] = {}
        species = list(set([hit[2] for hit in bmap[t]]))
 
        ## organize the hits by species
        for hit in bmap[t]:
            if not transcript2genes[t].has_key(hit[2]):
                transcript2genes[t][hit[2]] = []
            
            transcript2genes[t][hit[2]].append(hit[1])

    ## get inferred go terms for each transcript
    transcript2go = {}
    for t,hit in transcript2genes.iteritems():
        transcript2go[t] = []
        for genes in hit.itervalues():
            #gene = v[1]
            for gene in genes:
                if gene2go.has_key(gene):
                    transcript2go[t].extend(gene2go[gene])
        transcript2go[t] = list(set(transcript2go[t]))
        transcript2go[t].sort()
        
    ## write to file
    writer = csv.writer(open(genesetSummaryFile,'w'))
    writer.writerow(["transcript","p-value","genes","go-terms"])
    allTerms = []

    for ts in gsTranscripts:
        pvalue = deseqMat[np.where(deseqIds==ts)[0],pvalInd][0]
        reportedGenes = []
        for taxa, genes in transcript2genes[ts].iteritems():
            reportedGenes.extend(genes[:2])
        reportedGenes = list(set(reportedGenes))

        if len(reportedGenes) > 1:
            genes = ";".join(reportedGenes)
        else:
            genes = reportedGenes[0]
    
        terms = transcript2go[ts]

        if terms:
            allTerms.extend(terms)

        if not terms:
            terms = "None"
        elif len(terms) > 1: 
            terms = ";".join(terms)
        else:
            terms = terms[0]

        writer.writerow([ts,pvalue,genes,terms])

    writer.writerow(["--------"])
    ## write a summary of the go terms
    allTerms = np.array(list(set(allTerms)))
    allTermCounts = np.zeros(allTerms.size,)

    for t, term in enumerate(allTerms):
        for ts in gsTranscripts:
            allTermCounts[t] += np.where(np.array(transcript2go[ts]) == term)[0].size
    
    sortedTerms = allTerms[np.argsort(allTermCounts)[::-1]]    
    sortedCounts = allTermCounts[np.argsort(allTermCounts)[::-1]]
    writer.writerow(["ID","Counts","Description"])
    for t,term in enumerate(sortedTerms):
        desc = session.query(GoTerm).filter(GoTerm.go_id == term).first().name
        writer.writerow([term,sortedCounts[t],desc])

name = 'x-hs-noiea'
aspect = 'bp'
#transcript = 'genes'
#assembly = 'gg'
#geneset = 'gs-9'



if name == 'x-hs-noiea':

    ## gg, genes
    for gs in ['gs-33']:
        write_summary(name,aspect,"genes","gg",gs)

    ## gg, isoforms
    #for gs in ['gs-380','gs-31','gs-74','gs-395']:
    #    write_summary(name,aspect,"isoforms","gg",gs)

    #for gs in ['gs-301','gs-309','gs-367','gs-426','gs-183','gs-395']:
    #    write_summary(name,aspect,"genes","dn",gs)

#$negative
#     Gene_set Gene_set_name Score     p-value FDR
#[1,] "152"    "gs-483"      "-0.7464" "0"     "0"
#
#$positive
#     Gene_set Gene_set_name Score    p-value  FDR     
#[1,] "63"     "gs-205"      "1.5284" "0"      "0"     
#[2,] "86"     "gs-266"      "0.7075" "0"      "0"     
#[3,] "127"    "gs-395"      "1.205"  "0"      "0"     
#[4,] "31"     "gs-107"      "0.9474" "0.005"  "0.1925"
#[5,] "75"     "gs-234"      "0.1759" "0.0075" "0.231" 
#[6,] "35"     "gs-122"      "0.8115" "0.0125" "0.3208"
#[7,] "16"     "gs-59"       "0.723"  "0.02"   "0.44"  
#[8,] "13"     "gs-48"       "0.8178" "0.025"  "0.4812"



sys.exit()









for gene in gsGenes:
    geneTranscripts = []
    for tax, tscpts in gene2transcript[gene].iteritems():
        geneTranscripts.extend(tscpts)

    if len(geneTranscripts) > 1:
        pvalue = [deseqMat[np.where(deseqIds==t)[0],pvalInd][0] for t in geneTranscripts]
        pvalue = ";".join([str(round(pv,4)) for pv in pvalue])
    else:
        pvalue = deseqMat[np.where(deseqIds==geneTranscripts)[0],pvalInd][0]

    description = re.sub(",",";",gs.geneInfo[gene]['description'])
    taxa = gs.geneInfo[gene]['taxa']
    goTerms = gs.get_go_terms(gene)
    symbol = gs.geneInfo[gene]['symbol']
    
    ## get paralogs/orthologs
    similarGenes = set([])
    for gt in geneTranscripts:
        similarGenes.update([transcript2gene[gt]])

    similarGenes = list(similarGenes)

    ## make lists int string
    similarGenes.remove(gene)
    if len(similarGenes) > 1:
        similarGenes = ";".join(similarGenes[:3])
    elif len(similarGenes) == 1 :
        similarGenes = similarGenes[0]
    else:
        similarGenes = 'None'

    if len(geneTranscripts) > 1:
        geneTranscripts = ";".join(geneTranscripts)
    else:
        geneTranscripts = geneTranscripts[0]

    writer.writerow([gene,symbol,similarGenes,geneTranscripts,pvalue,description,taxa,goTerms])

sys.exit()












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
edgerResultsPath = os.path.join(featuresDir,"edger_%s_de.csv"%(transcript))
edgerIds, edgerColumns, edgerMat = read_de_results(edgerResultsPath,tool='edgeR')
deseqResultsPath = os.path.join(featuresDir,"deseq_%s_de.csv"%(transcript))
deseqIds, deseqColumns, deseqMat = read_de_results(deseqResultsPath,tool='DESeq')
#matrixFilePath = os.path.join(featuresDir,"Trinity_%s.counts.matrix"%source)
#countIds,countColumns,countMat = read_matrix(matrixFilePath,mtype='int')
dfeMatrixPath = os.path.join(featuresDir,"deseq_%s_de_samples.csv"%(transcript))
dfeIds,dfeColumns,dfeMat = read_matrix(dfeMatrixPath,mtype='float')

## load the blast map 
if transcript == 'genes':
    blastMap = get_blast_map(os.path.join("..","blast","summary_blast_%s.csv"%assembly),\
                             taxaList=["8364"],asGenes=True)
else:
    blastMap = get_blast_map(os.path.join("..","blast","summary_blast_%s.csv"%assembly),\
                             taxaList=["8364"],asGenes=False)

## setup filters for the transcripts
threshold = 0.5
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
hm.draw_heatmap(mat,cmap='uy',clabels=True,rlabels=True)
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
