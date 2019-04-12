#!/usr/bin/env python
"""
probe the taxa in the list for annotation coverage
and summary information

"""

import sys,getopt,os

from htsint import TaxaSummary
from htsint.database import db_connect

self.session,self.engine = db_connect()

## get how many genes code for proteins
codingQuery = self.session.query(Uniprot).filter_by(taxa_id=taxaQuery.id).all()
codingGenes = list(set([u.gene_id for u in uniprotQuery]))
remove_empty(codingGenes)

## get number of genes/proteins with at least one annotation                                  
annotatedGenes = list(set([a.gene_id for a in annotations]))
annotatedProts = list(set([a.uniprot_id for a in annotations]))
remove_empty(annotatedGenes)
remove_empty(annotatedProts)

apQuery = [self.session.query(Uniprot).filter_by(id=uid).first() for uid in annotatedProts]
#apQuery = self.session.query(Uniprot).filter(Uniprot.id.in_(annotatedProts)).all()
                                  
genesFromUniprot = list(set([uq.gene_id for uq in apQuery]))
remove_empty(genesFromUniprot)
genesCovered = list(set(annotatedGenes).union(set(genesFromUniprot)))










#taxaFile = 'foo.tab'
#baseDir =  os.path.join("..",os.path.realpath(os.path.dirname(__file__)))
#taxaFilePath = os.path.join(baseDir,taxaFile)

#tsum = TaxaSummary(taxaFilePath)

## print names summary
#tsum.print_names_summary()

## print annotation summary
tsum.print_annotation_summary()

#for taxon in tsum.taxaList:
#    print taxon
