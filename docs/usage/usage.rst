===========
Using Augur
===========

Augur consists of a number of tools that allow the user to filter and align sequences, infer trees, and integrate the phylogenetic analysis with meta data.
The different tools are meant to be composable and the output of one command will serve as the input of other commands.
All of Augur's commands are accessed through the ``augur`` program followed by the name of the command, e.g. ``augur ancestral`` to infer ancestral sequences.

Each command requires at least one input file and will produce one or more output files following the pattern

.. code-block:: bash

    augur ancestral --tree my_tree.nwk --alignment my_alignment.fasta --output-node-data ancestral_sequences.json

Most commands admit additional arguments to modify how the analysis is run.
Each command is documented below, and we are continually adding more examples to each command and providing real-life examples of their usage.
For instance, the documentation for :doc:`augur filter <./cli/filter>` shows how we use this to subsample the `nextstrain zika analysis <https://nextstrain.org/zika>`__.


`Note that you can also run each command with the` ``--help`` `option, for example` ``augur tree --help``, `for more information at the command-line.`


.. toctree::
   :maxdepth: 2

   cli/cli
   json_format
   envvars
