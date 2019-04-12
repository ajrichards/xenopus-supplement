#!/usr/bin/env python
"""
probe the taxa in the list for annotation coverage
and summary information

"""

import sys,getopt,os

from htsint import TaxaSummary

#taxaFile = 'xenopus.tab'
#taxaFile = 'amphibia.tab'
#baseDir =  os.path.join("..",os.path.realpath(os.path.dirname(__file__)))
#taxaInput = os.path.join(baseDir,taxaFile)
#taxaInput = ["8355","8364"]
taxaInput = ["10090","9606","7227","7955"] 
tsum = TaxaSummary(taxaInput)

## print names summary
tsum.print_names_summary()

## print annotation summary
tsum.print_annotation_summary(useIea=False)
