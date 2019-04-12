#!/usr/bin/python
# tutorials and docs
# http://www-huber.embl.de/users/anders/HTSeq/doc/tour.html

# install notes
# sudo apt-get install cython samtools tabix libbam-dev python-pyrex
# sudo git clone https://github.com/pysam-developers/pysam.git
# cd pysam/
# sudo python setup.py build
# sudo python setup.py install
# sudo apt-get install python-htseq 

### 'coverage vector' - a 1D vector of the length of a chromosome, where 
### each element counts how many reads cover the corresponding base pair in their alignment. 

import os,sys
import numpy as np
import matplotlib.pyplot as plt
import HTSeq

seqDir = os.path.join(os.path.expanduser("~"),"sequencing","xenopus-data")

def make_coverage_plot(sample):
    """
    uses SAM Alignments 
    they have many uses: align.iv, align.read, align.read.qual etc

    """

    ## read in BAM file
    bamFilePath = os.path.join(seqDir,"Sample_%s"%sample,"tophat_remapping","accepted_hits.bam") 
    bam_reader = HTSeq.BAM_Reader(bamFilePath )

    ## get counts
    totalReads = 0
    alignedReads = 0
    alignedTo = set([])
    cvg = HTSeq.GenomicArray( "auto", stranded=True, typecode='i' )
    for alngt in bam_reader:
        totalReads += 1

        if alngt.aligned:
            alignedReads += 1
            alignedTo.update([alngt.iv])
            cvg[ alngt.iv ] += 1

    ## get average coverage for each region
    alignedTo = list(alignedTo)
    avgCoverage = []
    for at in alignedTo:
        avgCoverage.append(np.array(list(cvg[at])).mean())
    avgCoverage = np.array(avgCoverage)

    ## make figure
    fig = plt.figure()
    ax = fig.add_subplot(111)
    inds = np.argsort(avgCoverage)
    sortedAvgCoverage = avgCoverage[inds] 
    plt.plot(np.arange(avgCoverage.size), sortedAvgCoverage, 'y.-')
    ax.set_title("Sample %s coverage summary"%sample)
    ax.set_ylabel("Average coverage")
    ax.set_xlabel("Aligned reads")
    txtToAdd = ""
    txtToAdd += "total reads: %s"%totalReads
    txtToAdd += "\nregions: %s"%len(alignedTo) 
    txtToAdd += "\naligned reads: %s"%alignedReads

    ## calculate percents
    nSmall = np.where(avgCoverage < 30)[0].size
    percentSmall = round((float(nSmall) / float(avgCoverage.size)) * 100.0, 2)
    txtToAdd += "\n< 30: %s"%percentSmall + "%"
    txtToAdd += "\n median coverage: %s"%round(np.median(avgCoverage),3) 
    plt.figtext(0.2,0.8,txtToAdd,color='black',fontsize=12,
                ha="left", va="center",
                bbox = dict(boxstyle="round",facecolor='white',alpha=0.8)
            )

    plt.savefig('sample_%s_coverage.png'%sample.lower())

    ## write these files for viewing with a genome browser
    #cvg.write_bedgraph_file(os.path.join("..","results","%s_plus.wig"%(sample)), "+")
    #cvg.write_bedgraph_file(os.path.join("..","results","%s_minus.wig"%(sample)), "-")



## run script
for sample in ["A","B","C","D","E","F","G","H"]:
    print("creating coverage plot for Sample %s"%sample)
    make_coverage_plot(sample)


print('done')
