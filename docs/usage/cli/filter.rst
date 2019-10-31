============
augur filter
============

* `How we subsample sequences in the zika-tutoral <#how-we-subsample-sequences-in-the-zika-tutoral>`__

----

.. argparse::
    :module: augur
    :func: make_parser
    :prog: augur
    :path: filter
        
How we subsample sequences in the zika-tutoral
==============================================

As an example, we'll look that the ``filter`` command in greater detail using material form the `zika tutorial <../../tutorials/zika_tutorial.html>`__.
The filter command allows you to selected various subsets of your input data for different types of analysis.
A simple example use of this command would be

.. code-block:: bash

	augur filter --sequences data/sequences.fasta --metadata data/metadata.tsv --min-date 2012 --output filtered.fasta

This command will select all sequences with collection date in 2012 or later.
The filter command has a large number of options that allow flexible filtering for many common situations.
One such use-case is the exclusion of sequences that are known to be outliers (e.g. because of sequencing errors, cell-culture adaptation, ...).
These can be specified in a separate file:

.. code-block:: bash

    BRA/2016/FC_DQ75D1
    COL/FLR_00034/2015
    ...

To drop such strains, you can pass the name of this file to the augur filter command:

.. code-block:: bash

  augur filter --sequences data/sequences.fasta \
             --metadata data/metadata.tsv \
             --min-date 2012 \
             --exclude config/dropped_strains.txt \
             --output filtered.fasta

(To improve legibility, we have wrapped the command across multiple lines.)
If you run this command (you should be able to copy-paste this into your terminal) on the data provided in the `zika tutorial <zika_tutorial.html>`__, you should see that one of the sequences in the data set was dropped since its name was in the ``dropped_strains.txt`` file.

Another common filtering operation is subsetting of data to a achieve a more even spatio-temporal distribution or to cut-down data set size to more manageable numbers.
The filter command allows you to select a specific number of sequences from specific groups, for example one sequence per month from each country:

.. code-block:: bash

  augur filter \
    --sequences data/sequences.fasta \
    --metadata data/metadata.tsv \
    --min-date 2012 \
    --exclude config/dropped_strains.txt \
    --group-by country year month \
    --sequences-per-group 1 \
    --output filtered.fasta

This subsampling and filtering will reduce the number of sequences in the tutorial data set from 34 to 24.
