#!/usr/bin/python
"""
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Xenopus_Silurana_tropicalis/protein/protein.fa.gz
gunzip -c protein.fa.gz > Xenopus_Silurana_tropicalis_protein.fa
"""

import os,sys
from htsint.blast import ParallelBlast

filePath = os.path.realpath(os.path.join(".","Xenopus_Silurana_tropicalis_protein.fa"))

if os.path.exists(filePath) == False:
    raise Exception("Cannot find sequences")

BLASTDB = "/usr/local/share/htsint"
parallelBlast = ParallelBlast(filePath,'uniprot_sprot.fasta.db',cmd='blastp',BLASTDB=BLASTDB)
parallelBlast.evalue = 0.001

chunks = 29
parallelBlast.create_scripts(chunks,"adam.richards@ecoex-moulis.cnrs.fr")
parallelBlast.submit()
