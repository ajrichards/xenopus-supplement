#!/usr/bin/python
"""

"""

import os,csv,sys,re,getopt
import numpy as np
import matplotlib.pyplot as plt
from htsint import run_subprocess
from htsint.tools import read_RSEM_counts_files
from htsint.blast import get_blast_map
from sqlalchemy.sql import select
from htsint.database import db_connect,Taxon,Gene


## specify the locations
homeDir = os.path.join(os.path.expanduser("~"),"sequencing","xenopus")
readsDir = os.path.join(homeDir,'reads')
    
## more variables
email = "adam.richards@ecoex-moulis.cnrs.fr"
trinityDir = "/usr/src/trinityrnaseq_r20140717/"


if __name__ == "__main__":

    for method in ['dn','gg']:
        featuresDir = os.path.join(homeDir,"%s-trinity"%method,"features")
        for source in ['genes','isoforms']:
            countMatrixPath = os.path.join(featuresDir,"Trinity_%s.counts.matrix"%source)

            ## run DESeq
            outputPath = os.path.join(featuresDir,"deseq_%s_%s_de.csv"%(source,'behavior'))
            cmd = "Rscript runDESeq.R %s %s"%(countMatrixPath,outputPath)
            print("running...\n%s"%cmd)
            run_subprocess(cmd)

            ## run edgeR
            outputPath = os.path.join(featuresDir,"edger_%s_behavior_de.csv"%source)
            cmd = "Rscript runEdgeR.R %s %s"%(countMatrixPath,outputPath)
            print("running...\n%s"%cmd)
            run_subprocess(cmd)

