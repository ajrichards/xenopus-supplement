#!/usr/bin/env python 
"""From the docs... 

http://cufflinks.cbcb.umd.edu/tutorial.html

RNA-Seq is a powerful technology for gene and splice variant
discovery. You can use Cufflinks to help annotate a new genome or
find new genes and splice isoforms of known genes in even
well-annotated genomes. Annotating genomes is a complex and difficult
process, but we outline a basic workflow that should get you started
here. The workflow also excludes examples of the commands you'd run
to implement each step in the workflow. Suppose we have RNA-Seq reads
from human liver, brain, and heart.


Map the reads for each tissue to the reference genome

We recommend that you use TopHat to map your reads to the
reference genome. For this example, we'll assume you have
paired-end RNA-Seq data. You can map reads as follows:

tophat -r 50 -o tophat_brain /seqdata/indexes/hg19 brain_1.fq brain_2.fq tophat -r 50 -o 
tophat_liver /seqdata/indexes/hg19 liver_1.fq liver_2.fq tophat -r 50 -o tophat_heart
/seqdata/indexes/hg19 heart_1.fq heart_2.fq

The commands above are just examples of how to map reads with TopHat. Please see the TopHat manual for more details on RNA-Seq read mapping.

Run Cufflinks on each mapping file

The next step is to assemble each tissue sample independently using Cufflinks. Assemble each tissue like so:

cufflinks -o cufflinks_brain tophat_brain/accepted_hits.bam
cufflinks -o cufflinks_liver tophat_liver/accepted_hits.bam
cufflinks -o cufflinks_heart tophat_liver/accepted_hits.bam

"""

import os,sys,re,shutil
import HTSeq
from htsint import run_subprocess

## locations
homeDir = os.path.join(os.path.expanduser("~"),"sequencing","xenopus")
readsDir = os.path.join(homeDir,'reads')
genomePath = os.path.join("20100930","sequences","Xenopus_tropicalis.main_genome.scaffolds.fasta")
transcriptAnnot = os.path.join("Xentr7_2_Stable.gff3")
assemblyDir = os.path.join(homeDir,"reference")

if not os.path.isdir(assemblyDir):
    os.mkdir(assemblyDir)

def assemble_reads(sampleList):
    ## get the files
    def get_paired_files(sample,fullpath=True):
        leftFile,rightFile = None,None
        for fileName in os.listdir(readsDir):
            if re.search("unpaired",fileName):
                continue

            if re.search("^%s.*\.fq$"%sample,fileName) and re.search("left",fileName):
                if fullpath:
                    leftFile = os.path.realpath(os.path.join(readsDir,fileName))
                else:
                    leftFile = fileName
            elif re.search("^%s.*.fq$"%sample,fileName) and re.search("right",fileName):
                if fullpath:
                    rightFile = os.path.realpath(os.path.join(readsDir,fileName))
                else:
                    rightFile = fileName

        return leftFile,rightFile

    allLeft = []
    allRight = []
    for sample in sampleList:
        left,right = get_paired_files(sample)
        allLeft.append(left)
        allRight.append(right)
        print left,right

    ## check
    if len(allLeft) != len(sampleList):
        raise Exception("Invalid number of left sequences")
    if len(allRight) != len(sampleList):
        raise Exception("Invalid number of right sequences")

    print("\n\ncreate symbolic link\nln -s %s xtropicalis.fa"%genomePath)

    bowtieCmd = "bowtie2-build %s xtropicalis"%(genomePath)
    print("\n\nRun bowtie...")
    print(bowtieCmd)
    
    
    topHatCmd = "tophat -p 8 -G %s -o %s %s %s %s && "%(transcriptAnnot,os.path.join(assemblyDir,"A"),'xtropicalis',allLeft[0],allRight[0])+\
                "tophat -p 8 -G %s -o %s %s %s %s && "%(transcriptAnnot,os.path.join(assemblyDir,"B"),'xtropicalis',allLeft[1],allRight[1])+\
                "tophat -p 8 -G %s -o %s %s %s %s && "%(transcriptAnnot,os.path.join(assemblyDir,"C"),'xtropicalis',allLeft[2],allRight[2])+\
                "tophat -p 8 -G %s -o %s %s %s %s && "%(transcriptAnnot,os.path.join(assemblyDir,"D"),'xtropicalis',allLeft[3],allRight[3])+\
                "tophat -p 8 -G %s -o %s %s %s %s && "%(transcriptAnnot,os.path.join(assemblyDir,"E"),'xtropicalis',allLeft[4],allRight[4])+\
                "tophat -p 8 -G %s -o %s %s %s %s && "%(transcriptAnnot,os.path.join(assemblyDir,"F"),'xtropicalis',allLeft[5],allRight[5])+\
                "tophat -p 8 -G %s -o %s %s %s %s && "%(transcriptAnnot,os.path.join(assemblyDir,"G"),'xtropicalis',allLeft[6],allRight[6])+\
                "tophat -p 8 -G %s -o %s %s %s %s"%(transcriptAnnot,os.path.join(assemblyDir,"H"),'xtropicalis',allLeft[7],allRight[7])

    print("\nRun tophat...")
    print topHatCmd


              #"samtools index %s && "%(os.path.join(assemblyDir,"A","accepted_hits.sorted.bam"))+\
    sortCmd = "samtools sort -n %s %s && "%(os.path.join(assemblyDir,"A","accepted_hits.bam"),os.path.join(assemblyDir,"A","accepted_hits.sorted"))+\
              "samtools sort %s %s && "%(os.path.join(assemblyDir,"B","accepted_hits.bam"),os.path.join(assemblyDir,"B","accepted_hits_sorted"))+\
              "samtools sort %s %s && "%(os.path.join(assemblyDir,"C","accepted_hits.bam"),os.path.join(assemblyDir,"C","accepted_hits_sorted"))+\
              "samtools sort %s %s && "%(os.path.join(assemblyDir,"D","accepted_hits.bam"),os.path.join(assemblyDir,"D","accepted_hits_sorted"))+\
              "samtools sort %s %s && "%(os.path.join(assemblyDir,"E","accepted_hits.bam"),os.path.join(assemblyDir,"E","accepted_hits_sorted"))+\
              "samtools sort %s %s && "%(os.path.join(assemblyDir,"F","accepted_hits.bam"),os.path.join(assemblyDir,"F","accepted_hits_sorted"))+\
              "samtools sort %s %s && "%(os.path.join(assemblyDir,"G","accepted_hits.bam"),os.path.join(assemblyDir,"G","accepted_hits_sorted"))+\
              "samtools sort %s %s "%(os.path.join(assemblyDir,"H","accepted_hits.bam"),os.path.join(assemblyDir,"H","accepted_hits_sorted"))

    print("\nSort BAM files")
    print(sortCmd) 
    print("\n")
    
    cufflinksCmd = "cufflinks -p 8 -o %s %s && "%(os.path.join(assemblyDir,"A","cufflinks"),os.path.join(assemblyDir,"A","accepted_hits_sorted.bam"))+\
                   "cufflinks -p 8 -o %s %s && "%(os.path.join(assemblyDir,"B","cufflinks"),os.path.join(assemblyDir,"B","accepted_hits_sorted.bam"))+\
                   "cufflinks -p 8 -o %s %s && "%(os.path.join(assemblyDir,"C","cufflinks"),os.path.join(assemblyDir,"C","accepted_hits_sorted.bam"))+\
                   "cufflinks -p 8 -o %s %s && "%(os.path.join(assemblyDir,"D","cufflinks"),os.path.join(assemblyDir,"D","accepted_hits_sorted.bam"))+\
                   "cufflinks -p 8 -o %s %s && "%(os.path.join(assemblyDir,"E","cufflinks"),os.path.join(assemblyDir,"E","accepted_hits_sorted.bam"))+\
                   "cufflinks -p 8 -o %s %s && "%(os.path.join(assemblyDir,"F","cufflinks"),os.path.join(assemblyDir,"F","accepted_hits_sorted.bam"))+\
                   "cufflinks -p 8 -o %s %s && "%(os.path.join(assemblyDir,"G","cufflinks"),os.path.join(assemblyDir,"G","accepted_hits_sorted.bam"))+\
                   "cufflinks -p 8 -o %s %s "%(os.path.join(assemblyDir,"H","cufflinks"),os.path.join(assemblyDir,"H","accepted_hits_sorted.bam"))

    print("\nRun Cufflinks")
    print(cufflinksCmd)
    print("\n")

    ## create assembly_list.txt
    fid = open("assembly_list.txt","w")
    for sample in sampleList:
        fid.write("%s\n"%os.path.join(assemblyDir,sample,"cufflinks","transcripts.gtf"))
 
    fid.close()

    ## run cuff merge to put it all together


if __name__ == "__main__":

    sampleList = ["A", "B", "C", "D", "E", "F", "G", "H"]
    assemble_reads(sampleList)
