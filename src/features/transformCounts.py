#!/usr/bin/python
"""
Use normalized count matrix provided by DESeq2 to get transformed values
"""

__author__ = "Adam Richards"

import sys,os,getopt,csv,time,re
import numpy as np

## declare variables
homeDir = os.path.join(os.path.expanduser("~"),"sequencing","xenopus")


def run_transformation(assembly,transcript):

    ## error check
    countMatrix = os.path.join(homeDir,"%s-trinity"%(assembly),'features','Trinity_%s.TMM.fpkm.matrix'%(transcript))
    transformedMatrix = os.path.join(homeDir,"%s-trinity"%(assembly),'features','Trinity_%s_counts_transformed_normalized.csv'%(transcript))

    if not os.path.exists(countMatrix):
        raise Exception('cannot find count matrix')

    ## infile
    fidin = open(countMatrix,'r')
    reader = csv.reader(fidin,delimiter="\t")
    header = reader.next()

    ## outfile
    fidout = open(transformedMatrix,'w')
    writer = csv.writer(fidout)
    writer.writerow(['transcript'] + header[1:])

    ## write file
    print("writing... %s"%transformedMatrix)
    lineCount = 0
    for linja in reader:
        lineCount += 1
        newLine = [linja[0]] + np.log2(1.0 + np.array([float(i) for i in linja[1:]])).tolist()
        writer.writerow(newLine)

run_transformation("gg","genes")
run_transformation("gg","isoforms")
run_transformation("dn","genes")
run_transformation("dn","isoforms")


print('complete.')
