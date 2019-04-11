#!/usr/bin/python
"""
tutorials and docs


first run:
 ./align_and_estimate_abundance.pl --transcripts Trinity.fasta --est_method RSEM --aln_method bowtie --trinity_mode --prep_reference

then in parallel for each sample
./align_and_estimate_abundance.pl --transcripts Trinity.fasta --seqType fq --left reads_1.fq --right reads_2.fq --est_method RSEM --aln_method bowtie --trinity_mode

"""

import os,sys,re,shutil,getopt
from htsint import run_subprocess,RunSubprocess

## read input
if len(sys.argv) < 2:
    raise Exception(sys.argv[0] + " [-m] is 'dn' or 'gg' and [-c] specifies a cluster environment\n")

try:
    optlist, args = getopt.getopt(sys.argv[1:], 'm:c')
except getopt.GetoptError:
    raise Exception(sys.argv[0] + " [-m] is 'dn' or 'gg' and [-c] specifies a cluster environment\n")
    sys.exit()

cluster = False
method = None
for o,a in optlist:
    if o == '-c':
        cluster = True
    if o == '-m':
        method = a

## clean cluster dir
if cluster == True:
    clusterDir = os.path.join(os.path.realpath("."),'cluster')
    if os.path.isdir(clusterDir):
        shutil.rmtree(clusterDir)
    os.mkdir(clusterDir)

## specify the locations
homeDir = os.path.join(os.path.expanduser("~"),"sequencing","xenopus")
readsDir = os.path.join(homeDir,'reads')

if method == 'gg':
    assemblyDir = os.path.join(homeDir,"gg-trinity")
    transcriptsFilePath = os.path.join(assemblyDir,"Trinity-GG.fasta")
elif method == 'dn':
    assemblyDir = os.path.join(homeDir,"dn-trinity")
    transcriptsFilePath = os.path.join(assemblyDir,"Trinity.fasta")
    
featuresDir = os.path.join(assemblyDir,"features")
if not os.path.isdir(featuresDir):
    os.mkdir(featuresDir)

if not os.path.exists(transcriptsFilePath):
    raise Exception("Cannot find transcripts")

## more variables
email = "adam.richards@ecoex-moulis.cnrs.fr"
trinityDir = "/usr/src/trinityrnaseq-2.0.4"

def assemble_reads(sampleList):
    ## get the files
    def get_paired_files(sample):
        leftFile,rightFile = None,None
        for fileName in os.listdir(readsDir):
            if re.search("unpaired",fileName):
                continue

            if re.search("^%s.*\.fq$"%sample,fileName) and re.search("left",fileName):
                leftFile = os.path.realpath(os.path.join(readsDir,fileName))
            elif re.search("^%s.*.fq$"%sample,fileName) and re.search("right",fileName):
                rightFile = os.path.realpath(os.path.join(readsDir,fileName))
        return leftFile,rightFile

    allLeft = []
    allRight = []
    for sample in sampleList:
        left,right = get_paired_files(sample)
        allLeft.append(left)
        allRight.append(right)

    ## check
    if len(allLeft) != len(sampleList):
        raise Exception("Invalid number of left sequences")
    if len(allRight) != len(sampleList):
        raise Exception("Invalid number of right sequences")

    if None in allLeft or None in allRight:
        raise Exception("Were the sequences prepped?")


    ## first prepare the reference
    cmd = "export TRINITY_HOME=%s;\n"%(trinityDir)+\
          "$TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts %s "%(transcriptsFilePath)+\
          "--est_method RSEM --aln_method bowtie --trinity_mode --prep_reference --output_dir %s"%(featuresDir) 
    
    print('preping reference...')
    run_subprocess(cmd)
    print('reference prepreation complete.')
    
    for s,sample in enumerate(sampleList):
        cmd = "export TRINITY_HOME=%s;\n"%(trinityDir)+\
              "$TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts %s "%(transcriptsFilePath)+\
              "--seqType fq --left %s --right %s "%(allLeft[s],allRight[s])+\
              "--est_method RSEM --aln_method bowtie --trinity_mode "+\
              "--output_dir %s --output_prefix %s"%(featuresDir,sample)

        if cluster == True:
            submitFile = os.path.join(clusterDir,"%s-%s.sh"%(sample,s))
            submitLog =  os.path.join(clusterDir,"%s-%s.log"%(sample,s))

            f = open(submitFile, 'w')
            f.write("#!/bin/bash\n"+\
                    "#$ -S /bin/bash\n"+\
                    "#$ -j yes\n"+\
                    "#S -M %s\n"%email+\
                    "#$ -o %s\n"%submitLog+\
                    "export TRINITY_HOME=%s\n"%(trinityDir)+\
                    cmd)

            f.close()
            
            print("submitting %s"%submitFile)
            os.system("qsub " + submitFile)
            
        else:
            print("running...\n%s"%cmd)
            #run_subprocess(cmd)

if __name__ == "__main__":

    sampleList = ["A", "B", "C", "D", "E", "F", "G", "H"]
    assemble_reads(sampleList)
