.. notes from the literature


Preparation of the reference genome
======================================

In this section we download the reference genome for *X. tropicalis*, prepare the files, and map quality controlled reads to the genome.


Obtaining the necessary files
---------------------------------------------


1. obtain the genome and genome annotations

   .. code-block:: bash

      ~$ wget ftp://ftp.xenbase.org/pub/Genomics/JGI/Xentr7.1/xenopus_tropicalis_v7.1.tar.gz
      ~$ wget ftp://ftp.xenbase.org/pub/Genomics/JGI/Xentr7.1/Xentr7_2_Stable.gff3.gz
      ~$ tar xzf xenopus_tropicalis_v7.1.tar.gz 
      ~$ tar xzf Xentr7_2_Stable.gff3.gz  

2. convert gff3 to gtf file

   .. code-block:: bash

      ~$ gffread Xentr7_2_Stable.gff3 -T -o Xentr7_2_Stable.gtf
   
3. download xenbase gene mappings

   .. code-block:: bash

      ~$ wget ftp://ftp.xenbase.org/pub/GenePageReports/JgiToXenbaseGenePage_7.1.txt


Run STAR
------------------

STAR is a powerful and fast program to align RNA-Seq data [Dobin13]_.  A minimum of 30GB of RAM is recommended to run STAR.
Most of the following process is handled by a Python script and the generalized process is enumerated below.

   * :download:`runSTAR.py <../assembly/runSTAR.py>`

1. create a directory for the mapping and add necessary files

   .. code-block:: bash
 
      ~$ /usr/src/STAR/source/STAR --version
      STAR_2.4.0j
      ~$ mkdir ~/sequencing/xenopus/star
      ~$ cp ./20100930/sequences/Xenopus_tropicalis.main_genome.scaffolds.fasta ~/sequencing/xenopus/star/genome.fa
      ~$ cp Xentr7_2_Stable.gtf ~/sequencing/xenopus/star/Xentr7_2_Stable.gtf

2. index the genome

   .. code-block:: bash

      ~$ /usr/src/STAR/source/STAR --runMode genomeGenerate --runThreadN 24 --genomeDir ~/sequencing/xenopus/star 
         --genomeFastaFiles ~/sequencing/xenopus/star/genome.fa --sjdbGTFfile ~/sequencing/xenopus/star/Xentr7_2_Stable.gtf 
         --sjdbOverhang 100

3. align for each sample e.g.

   .. code-block:: bash

      ~$ /usr/src/STAR/source/STAR --runThreadN 24 --genomeDir --genomeDir ~/sequencing/xenopus/star
         --readFilesIn left_reads_sample_A.fq right_reads_sample_A.fq --outFileNamePrefix ~/sequencing/xenopus/star/A_

Convet the SAM files so that they are ready as input for genome assembly and feature analysis
------------------------------------------------------------------------------------------------

Each step need to be carried out on each sample. The following operations are carried out using SAMtools [Li09]_.

1. convert the SAM to BAM file

   .. code-block:: bash

      ~$ cd ~/sequencing/xenopus/star
      ~$ /usr/bin/samtools view -b -S A_Aligned.out.sam > A_aligned.bam

2. sort the BAM file

   .. code-block:: bash

      ~$  /usr/bin/samtools sort -n A_aligned.bam A_aligned_sorted

3. convert the sorted BAM into SAM

   .. code-block:: bash

      ~$ /usr/bin/samtools view -h A_aligned_sorted.bam > A_aligned_sorted.sam


From the sorted SAM files we can obtain a count matrix for feature analysis.  The sorted BAM files are concatenated and used as input into the genome guided assembly.

Concatenate the BAM files
----------------------------

1. merge the files

   .. code-block:: bash

      ~$ /usr/bin/samtools merge star_all_reads.bam A_aligned_sorted.bam B_aligned_sorted.bam ...

2. create a coordinate sorted bam file

   .. code-block:: bash

      ~$ /usr/bin/samtools sort star_all_reads.bam star_all_reads_sorted.bam

Links
----------

   * `STAR on GitHub <https://github.com/alexdobin/STAR>`_
   * `STAR Manual <https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf>`_
   * `SAMtools <http://samtools.sourceforge.net/>`_
