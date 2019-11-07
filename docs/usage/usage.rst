===========
Using Augur
===========

.. note:: We have just released version 6 of augur -- `check our upgrading guide <../releases/migrating-v5-v6.html>`__

Augur consists of a number of tools that allow the user to filter and align sequences, build trees, and integrate the phylogenetic analysis with meta data.
The different tools are meant to be composable and the output of one command will serve as the input of other commands.
All of Augur's commands are accessed through the ``augur`` program.
For example, to infer ancestral sequences from a tree, you'd run ``augur ancestral``.


Each command is documented below, and we are continually adding more examples to each command and providing real-life examples of their usage.
For instance, the documentation for `augur filter <./cli/filter.html>`__ shows how we use this to subsample the `nextstrain zika build <https://nextstrain.org/zika>`__.


`Note that you can also run each command with the` ``--help`` `option, for example` ``augur tree --help``, `for more information at the command-line.`


.. toctree::
   :maxdepth: 2

   cli/cli
   augur_snakemake
   json_format
   envvars
