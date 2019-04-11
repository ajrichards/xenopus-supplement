.. notes from the literature


Assembly
======================================

.. note:: 
  llumina HiSeq 2500, CASAVA 1.8.2, RTA 1.13, Run ID: 130530_SN982_0225_AC20RLACXX


Trinity assembly without a reference genome
-----------------------------------------------

We carried out the *de novo* transcriptome assembly with the software
pacakge `Trinity <http://trinityrnaseq.sourceforge.net>`_
[Grabherr11]_.  The Trinity software suite consists of three main
pieces:

   * **Inchworm** - assembles the RNA-seq data into the unique
     sequences of transcripts, often generating full-length
     transcripts for a dominant isoform, but then reports just the
     unique portions of alternative ly spliced transcripts.
   * **Chrysalis** - clusters the Inchworm contigs into clusters and
     constructs complete de Bruijn graphs for each cluster. Each
     cluster represents the full transcriptonal complexity for a given
     gene (or sets o f genes that share sequences in
     common). Chrysalis then partitions the full read set among these
     disjoint graphs.
   * **Butterfly** then processes the individual graphs in parallel,
     tracing the paths that reads and pairs of reads take within the
     graph, ultimately reporting full-length transcripts for
     alternatively spliced isoforms, and teasing apart transcripts
     that corresponds to paralogous genes.


Sequence preprocessing
-----------------------

A Python script was used to unzip, concatenate and trim the original reads via system calls (`subprocess <https://docs.python.org/2/library/subprocess.html>`_).
We used `Trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`_ [Bolger14]_ to trim reads based on quality with default settings (LEADING:5 TRAILING:5 MINLEN:36).

   * :download:`preprocessReads.py <../assembly/preprocessReads.py>`

There are four output files.  Two are for the *paired* output where both reads passed, and two are for corresponding *unpaired* output where only one read passed (see manual).


Running Trinity
-------------------

A Python script was used to generate the Trinity arguments.

   * :download:`runTrinity.py <../assembly/runTrinity.py>`


The generalized Trinity command is shown below.

   .. code-block:: bash

      ~$ export TRINITY_HOME="/usr/src/trinityrnaseq-2.0.4"
      ~$ $TRINITY_HOME/Trinity --seqType fq --output /path/to/out --trimmomatic --full_cleanup
         --SS_lib_type FR --max_memory 26G --CPU 29 --normalize_reads
         --left left1.fastq,left2.fastq,left3.fastq 
         --right right1.fastq,right2.fastq,right3.fastq 2>&1 | tee ./run-trinity-gg.log

This produces the output file `./<TRINITY_OUT>/Trinity.fasta` which we can run some basic statistics on.

Basic summary
--------------

   .. code-block:: bash

      ~$ export TRINITY_HOME="/usr/src/trinityrnaseq-2.0.4"  
      ~$ $TRINITY_HOME/util/TrinityStats.pl ~/sequencing/xenopus/dn-trinity/Trinity.fasta


      ################################
      ## Counts of transcripts, etc.
      ################################
      Total trinity 'genes':163981
      Total trinity transcripts:218541
      Percent GC: 42.25

      ########################################
      Stats based on ALL transcript contigs:
      ########################################

      Contig N10: 4381
      Contig N20: 2849
      Contig N30: 2102
      Contig N40: 1559
      Contig N50: 1106

      Median contig length: 378
      Average contig: 707.59
      Total assembled bases: 154637570


      #####################################################
      ## Stats based on ONLY LONGEST ISOFORM per 'GENE':
      #####################################################

      Contig N10: 3491
      Contig N20: 2113
      Contig N30: 1354
      Contig N40: 879
      Contig N50: 619

      Median contig length: 338
      Average contig: 547.26
      Total assembled bases: 89739624


Trinity groups transcripts into clusters that are loosely referred to
as a gene. The accession identifiers in the
./trinity_out_dir/Trinity.fasta encode gene and isoform
information. Per the documentation (see links below) if we have the
accession >c0_g1_i1 this refers to Trinity read cluster c0, gene g1
and isoform i1. The gene identifier in this case is c0_g1.

Map the reads using BLAST
-----------------------------

1. First, we BLAST the transcript against SwissProt (`-c` can be used to initiate cluster mode)

   * :download:`runBlast.py <../blast/runBlast.py>`

2. Then, these data were are parsed and summarized.

   * :download:`runBlastParse.py <../blast/runBlastParallelParse.py>`
   * :download:`runBlastSummarize.py <../blast/runBlastSummarize.py>`
   * :download:`showTaxaSummary.py <../blast/showTaxaSummary.py>`

All of SwissProt 
^^^^^^^^^^^^^^^^^^^^^^

The best hit by isoform.

   .. figure:: ../figures/dn-trinity-blast-pie-isoforms.png
      :scale: 70%
      :align: center
      :alt: de-novo trinity blast taxa by isoform
      :figclass: align-center

   .. code-block:: none

      SwissProt - isoforms
      transcripts: 81193
      genes: 22372

.. note:: Genes correspond to unique gene identifiers, but the number includes orthalogs 

Xenopus restricted SwissProt
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The best hits only considering *Xenopus sp.*

   .. figure:: ../figures/dn-trinity-blast-pie-frog.png 
      :scale: 70%
      :align: center
      :alt: de-novo trinity blast taxa by gene
      :figclass: align-center

   .. code-block:: none

      SwissProt [8355,8364] - isoforms
      transcripts: 36603
      genes: 4058

BLAST against Xenopus amino acid sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   .. code-block:: none

      ~$ wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Xenopus_Silurana_tropicalis/protein/protein.fa.gz
      ~$ gunzip -c protein.fa.gz > xtropicalis.fasta
      ~$ makeblastdb -in xtropicalis.fasta -dbtype 'prot' -out xtropicalis
      ~$ less xtropicalis.fasta | grep gi > foo.txt && wc -l foo.txt

There were 28495 amino acid sequences in the fasta file.

   .. code-block:: none

      X. tropicalis - isoforms
      transcripts: 93345
      genes: 16356

Links
-----------------


   * `Example pipeline <https://wiki.hpcc.msu.edu/display/Bioinfo/Pipeline+for+Illumina+Data>`_
   * `Xenbase files <ftp://ftp.xenbase.org/pub/Genomics/JGI/>`_
   * `Trinity Sourceforge page <http://trinityrnaseq.sourceforge.net>`_
   * `Trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`_
   * `Trinity output <http://trinityrnaseq.sourceforge.net/#trinity_output>`_
