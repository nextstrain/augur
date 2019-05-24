==================
Command-line usage
==================

All of Augur's commands are accessed through the ``augur`` program.
For example, to infer ancestral sequences from a tree, you'd run ``augur ancestral``.
Each command is documented below.
You can also run each command with the ``--help`` option, for example ``augur tree --help``, for more information at the command-line.

Augur command
=============

.. argparse::
    :module: augur
    :func: make_parser
    :prog: augur
    :nodescription:
