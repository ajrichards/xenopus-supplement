#!/usr/bin/python
"""
Assembles the results from runBlast.py when done under a cluster envrionment
"""

import os,sys,getopt
from htsint.blast import ParseBlast

## read input
if len(sys.argv) < 2:
    raise Exception("\n\nInvalid syntax. use '%s'"%sys.argv[0] + " -s source\n")

try:
    optlist, args = getopt.getopt(sys.argv[1:], 's:')
except getopt.GetoptError:
        raise Exception("\n\nInvalid syntax. use '%s'"%sys.argv[0] + " -s source\n")

source = 'None'
for o,a in optlist:
    if o == '-s':
        source = a

if source not in ['dn','gg','reference']:
    raise Exception("bad source specified")

outfmt = '5'
blastFile = "blast-%s.%s"%(source, "outfmt"+outfmt)
outFile = "blast-%s-parsed.csv"%source

if not os.path.exists(blastFile):
    raise Exception("cannot find blast output file: %s"%blastFile)

## parse the resultsread the blast
parser = ParseBlast(blastFile,outFile=outFile)
parser.run()
