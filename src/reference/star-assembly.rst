.. notes from the literature


Assembly
======================================

.. note:: 
  llumina HiSeq 2500, CASAVA 1.8.2, RTA 1.13, Run ID: 130530_SN982_0225_AC20RLACXX


A Python script was used to unzip, concatenate and trim the original reads via system calls (`subprocess <https://docs.python.org/2/library/subprocess.html>`_).  We used `Trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`_ [Bolger14]_ to trim reads based on quality with default settings (LEADING:5 TRAILING:5 MINLEN:36).

   * :download:`preprocessReads.py <../assembly/preprocessReads.py>`

There are four output files.  Two are for the *paired* output where both reads passed, and two are for corresponding *unpaired* output where only one read passed (see manual).


Preparing the reference genome
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

See the section describing :doc:`the preparation of the reference genome </reference/reference-prep>`.


BLAST genome against SwissProt
---------------------------------

First, we obtain the known protein sequences.

.. code:: bash

   ~$ wget ftp://ftp.xenbase.org/pub/Genomics/JGI/Xentr7.1/Xentr7_2_Stable_Protein.fa.gz
   ~$ wget ftp://ftp.xenbase.org/pub/Genomics/JGI/Xentr7.1/Xentr7_2_Stable_Transcript.fa.gz
   ~$ gunzip -c Xentr7_2_Stable_Protein.fa.gz > Xentr7_2_Stable_Protein.fa
   ~$ gunzip -c Xentr7_2_Stable_Transcript.fa.gz > Xentr7_2_Stable_Transcript.fa

Next, we BLAST the transcript against SwissProt (`-c` can be used to initiate cluster mode)

   * :download:`runBlast.py <../blast/runBlast.py>`

Then, the data are parsed with one of the two following.

   * :download:`runBlastParallelParse.py <../blast/runBlastParallelParse.py>`
   * :download:`runBlastSummarize.py <../blast/runBlastSummarize.py>`


All of SwissProt 
^^^^^^^^^^^^^^^^^^^^^^

The best hit by transcript.

   .. figure:: ../figures/ref-trinity-blast-pie-isoforms.png
      :scale: 70%
      :align: center
      :alt: ref trinity blast taxa by isoform
      :figclass: align-center

   .. code-block:: none

      SwissProt
      transcripts: 38683
      genes: 17535
      
Xenopus restricted SwissProt
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The best hits only considering *Xenopus sp.*

   .. figure:: ../figures/ref-trinity-blast-pie-frog.png
      :scale: 70%
      :align: center
      :alt: ref trinity blast taxa by gene
      :figclass: align-center

   .. code-block:: none

      SwissProt [8355,8364]
      transcripts: 20638
      genes: 3750

.. note:: Genes correspond to unique gene identifiers, but the number includes orthalogs


BLAST against Xenopus amino acid sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   .. code-block:: none

      ~$ wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Xenopus_Silurana_tropicalis/protein/protein.fa.gz
      ~$ gunzip -c protein.fa.gz > xtropicalis.fasta
      ~$ makeblastdb -in xtropicalis.fasta -dbtype 'prot' -out xtropicalis
      ~$ less xtropicalis.fasta | grep gi > foo.txt && wc -l foo.txt

There were 28495 amino acid sequences in the fasta file.

   .. code-block:: none

      X. tropicalis
      transcripts: 41299
      genes: 20735
