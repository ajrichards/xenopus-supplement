#!/usr/bin/python
# tutorials and docs
# http://www-huber.embl.de/users/anders/HTSeq/doc/tour.html

import os,sys,csv
import HTSeq
from htsint import run_subprocess

## locations
seqDir = os.path.join(os.path.expanduser("~"),"sequencing","xenopus-data")
analysisDir = os.path.join(seqDir,"analyses")
if not os.path.isdir(analysisDir):
    os.mkdir(analysisDir)

def get_count_matrix(sample,clean=False):
    print("getting counts for sample %s"%sample)

    bamFilePath = os.path.join(seqDir,"Sample_%s"%sample,"tophat_remapping","accepted_hits.bam") 
    sbamFilePath = os.path.join(seqDir,"Sample_%s"%sample,"tophat_remapping","accepted_hits_sorted.bam") 
    samFilePath = os.path.join(seqDir,"Sample_%s"%sample,"tophat_remapping","accepted_hits.sam") 
    gtfFilePath = os.path.join(seqDir,"Sample_%s"%sample,"cufflinks_output","transcripts.gtf")
    countsFilePath = os.path.join(analysisDir,'counts_%s.txt'%sample.lower())
    geneCountsFilePath = os.path.join(analysisDir,'gene-counts_%s.txt'%sample.lower())

    ## clean 
    if clean == True:
        for filePath in [sbamFilePath,samFilePath,countsFilePath]:
            if os.path.exists(filePath):
                os.remove(filePath)

    ## error check
    if not os.path.exists(bamFilePath):
        raise Exception("Invalid file path\n...%s"%bamFilePath)
    if not os.path.exists(gtfFilePath):
        raise Exception("Invalid file path\n...%s"%gtfFilePath)

    ## sort bam file 
    if not os.path.exists(sbamFilePath):
        print('creating sorted bam file...')
        cmd = "/usr/bin/samtools sort -n %s %s"%(bamFilePath,sbamFilePath[:-4])
        run_subprocess(cmd)

    ## convert to sam file
    if not os.path.exists(samFilePath):
        print('converting bam to sam file...')
        cmd = "/usr/bin/samtools view -h -o %s %s"%(samFilePath,sbamFilePath)
        run_subprocess(cmd)

    ## save the count matrix
    if not os.path.exists(countsFilePath):
        print('running htseq-count')
        cmd = "/usr/bin/htseq-count %s %s > %s"%(samFilePath,gtfFilePath,countsFilePath)
        run_subprocess(cmd)
    
    ## create the exons
    gtf_file = HTSeq.GFF_Reader(gtfFilePath,end_included=True )
    sam_file = HTSeq.SAM_Reader(samFilePath)
    exons = HTSeq.GenomicArrayOfSets("auto",stranded=False)
    for feature in gtf_file:
        if feature.type == "exon":
            exons[feature.iv] += feature.name

    ## setup a count for each gene
    counts = {}
    for feature in gtf_file:
        if feature.type == "exon":
            counts[feature.name] = 0

    for alnmt in sam_file:
        if alnmt.aligned:
            iset = None
            for iv2, step_set in exons[ alnmt.iv ].steps():
                if iset is None:
                    iset = step_set.copy()
                else:
                    iset.intersection_update( step_set )
            if len( iset ) == 1:
                counts[ list(iset)[0] ] += 1
    
    ## save to file
    fidout = open(geneCountsFilePath,'w')
    writer = csv.writer(fidout)
    writer.writerow(["name","count"])
    for name in sorted(counts.keys()):
        writer.writerow([name,counts[name]])   

    fidout.close()
    print('done.')

if __name__ == "__main__":
    sampleList = ["A","B","C","D","E","F","G","H"]
    sampleList = ["A"]

    for sample in sampleList:
        get_count_matrix(sample,clean=False)
