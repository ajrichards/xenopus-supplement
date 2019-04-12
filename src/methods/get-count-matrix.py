#!/usr/bin/python
# tutorials and docs
# http://www-huber.embl.de/users/anders/HTSeq/doc/tour.html

import os,sys
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
    gtsFilePath = os.path.join(seqDir,"Sample_%s"%sample,"cufflinks_output","transcripts.gtf")
    countsFilePath = os.path.join(analysisDir,'counts_%s.txt'%sample.lower())
    warningsFilePath = os.path.join(analysisDir,'warnings_%s.txt'%sample.lower())

    ## clean 
    if clean == True:
        for filePath in [sbamFilePath,samFilePath,countsFilePath]:
            if os.path.exists(filePath):
                os.remove(filePath)

    ## error check
    if not os.path.exists(bamFilePath):
        raise Exception("Invalid file path\n...%s"%bamFilePath)
    if not os.path.exists(gtsFilePath):
        raise Exception("Invalid file path\n...%s"%gtsFilePath)

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
    print('running htseq-count')
    #cmd = "/usr/bin/htseq-count %s %s > %s"%(samFilePath,gtsFilePath,countsFilePath)
    cmd = "/usr/bin/htseq-count -m intersection-nonempty -s yes %s %s > %s 2> %s"%(samFilePath,gtsFilePath,countsFilePath,warningsFilePath)
    run_subprocess(cmd)


if __name__ == "__main__":
    sampleList = ["A","B","C","D","E","F","G","H"]

    for sample in sampleList:
        get_count_matrix(sample,clean=True)
