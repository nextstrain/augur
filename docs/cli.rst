==============
Augur commands
==============

All of Augur's commands are accessed through the ``augur`` program.
For example, to infer ancestral sequences from a tree, you'd run ``augur ancestral``.
Each command is documented below.
You can also run each command with the ``--help`` option, for example ``augur tree --help``, for more information at the command-line.

.. argparse::
   :ref: augur.make_parser
   :prog: augur
   :nodescription:
   :nosubcommands:

.. toctree::
   :maxdepth: 1
   :caption: Sub-commands
   :name: augur-commands
   :glob:

   cli/*
