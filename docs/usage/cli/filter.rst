============
augur filter
============

.. contents:: Table of Contents
   :local:

----

.. argparse::
    :module: augur
    :func: make_parser
    :prog: augur
    :path: filter

Guides
======

Below are some examples of using ``augur filter`` to sample data.

Filtering
---------

The filter command allows you to select various subsets of your input data for different types of analysis.
A simple example use of this command would be

.. code-block:: bash

  augur filter \
    --sequences data/sequences.fasta \
    --metadata data/metadata.tsv \
    --min-date 2012 \
    --output-sequences filtered_sequences.fasta \
    --output-metadata filtered_metadata.tsv

This command will select all sequences with collection date in 2012 or later.
The filter command has a large number of options that allow flexible filtering for many common situations.
One such use-case is the exclusion of sequences that are known to be outliers (e.g. because of sequencing errors, cell-culture adaptation, ...).
These can be specified in a separate text file (e.g. ``exclude.txt``):

.. code-block::

  BRA/2016/FC_DQ75D1
  COL/FLR_00034/2015
  ...

To drop such strains, you can pass the filename to ``--exclude``:

.. code-block:: bash

  augur filter \
    --sequences data/sequences.fasta \
    --metadata data/metadata.tsv \
    --min-date 2012 \
    --exclude exclude.txt \
    --output-sequences filtered_sequences.fasta \
    --output-metadata filtered_metadata.tsv

Subsampling
-----------

Another common filtering operation is subsetting of data to a achieve a more even spatio-temporal distribution or to cut-down data set size to more manageable numbers.
The filter command allows you to select a specific number of sequences from specific groups, for example one sequence per month from each country:

.. code-block:: bash

  augur filter \
    --sequences data/sequences.fasta \
    --metadata data/metadata.tsv \
    --min-date 2012 \
    --exclude exclude.txt \
    --group-by country year month \
    --sequences-per-group 1 \
    --output-sequences subsampled_sequences.fasta \
    --output-metadata subsampled_metadata.tsv
