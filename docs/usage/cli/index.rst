============
augur index
============

.. contents:: Table of Contents
   :local:

----

.. argparse::
    :module: augur
    :func: make_parser
    :prog: augur
    :path: index

Speed up filtering with a sequence index
========================================

As we describe in :doc:`the phylogenetic workflow tutorial <docs.nextstrain.org:tutorials/creating-a-phylogenetic-workflow>`, augur index precalculates the composition of the sequences (e.g., numbers of nucleotides, gaps, invalid characters, and total sequence length) prior to filtering.
The resulting sequence index speeds up subsequent filter steps especially in more complex workflows.

.. code-block:: bash

    mkdir -p results/
    augur index \
        --sequences data/sequences.fasta \
        --output results/sequence_index.tsv

The first lines in the sequence index look like this.

.. code-block:: bash

    strain	length	A	C	G	T	N	other_IUPAC	-	?	invalid_nucleotides
    PAN/CDC_259359_V1_V3/2015	10771	2952	2379	3142	2298	0	0	0	0	0
    COL/FLR_00024/2015	10659	2921	2344	3113	2281	0	0	0	0	0
    PRVABC59	10675	2923	2351	3115	2286	0	0	0	0	0
    COL/FLR_00008/2015	10659	2924	2344	3110	2281	0	0	0	0	0

We then provide the sequence index as an input to augur filter commands to speed up filtering on sequence-specific attributes.

.. code-block:: bash

    augur filter \
        --sequences data/sequences.fasta \
        --sequence-index results/sequence_index.tsv \
        --metadata data/metadata.tsv \
        --exclude config/dropped_strains.txt \
        --output results/filtered.fasta \
        --group-by country year month \
        --sequences-per-group 20 \
        --min-date 2012
