=======================
Augur command structure
=======================

Augur consists of a number of tools that allow the user to filter and align sequences, build trees, and integrate the phylogenetic analysis with meta data.
The different tools are meant to be composable and the output of one command will serve as the input of other commands.
All of Augur's commands are accessed through the ``augur`` program.
For example, to infer ancestral sequences from a tree, you'd run ``augur ancestral``.
Each command is documented below.
You can also run each command with the ``--help`` option, for example ``augur tree --help``, for more information at the command-line.


As an example, we'll look that the ``filter`` command in greater detail using the zika tutorial as template.
This command allows you to selected various subsets of your input data for different types of analysis.
A simple example use of this command would be

.. code-block:: bash

	augur filter --sequences data/sequences.fasta --metadata data/metadata.tsv --min-date 2012 --output filtered.fasta

This command will select all sequences with collection date in 2012 or later.
The filter command has a large number of options that allow flexible filtering for many common situations.
One such use-case is the exclusion of sequences that are known to be outliers (e.g.~because of sequencing errors, cell-culture adaptation, ...).
These can be specified in a separate file:
```
BRA/2016/FC_DQ75D1
COL/FLR_00034/2015
...
```
To drop such strains, you can pass the name of this file to the augur filter command:

.. code-block:: bash

  augur filter --sequences data/sequences.fasta \
             --metadata data/metadata.tsv \
             --min-date 2012 \
             --exclude config/dropped_strains.txt \
             --output filtered.fasta

(To improve legibility, we have wrapped the command across multiple lines.)
If you run this command (you should be able to copy-paste this into your terminal), you should see that one of the sequences in the data set was dropped since its name was in the ``dropped_strain.txt`` file.

.. code:: text

  BRA/2016/FC_DQ75D1
  COL/FLR_00034/2015
  ...

Another common filtering operation is subsetting of data to a achieve a more even spatio-temporal distribution or cut-down data set size to more manageable numbers.
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

This subsampling and filtering will reduce the number of sequences in this tutorial data set from 34 to 24.
