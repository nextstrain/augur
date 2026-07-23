===========================
augur refine
===========================


.. contents:: Table of Contents
   :local:

----

Command Line Options
====================

.. argparse::
    :module: augur
    :func: make_parser
    :prog: augur
    :path: refine

Configuration file options
==========================

``augur refine`` allows you to provide a configuration file (via ``--config``) containing values for any of its command-line options. 
Specifying an option on both the command line and in a config file will result in an error.

The following table shows the available configuration options and their descriptions. Note that boolean flags must be explicitly specified (e.g., ``true`` or ``false``).

.. argparse-config-table:: augur.refine register_parser

Guides
======

See :doc:`../../faq/refine`.
