============
augur export
============

.. note:: The behavior of ``augur export`` has changed in v6. Please see :doc:`here <../../releases/migrating-v5-v6>` for more details.

.. contents:: Table of Contents
   :local:
   :depth: 1

Command-line arguments
======================

.. argparse::
    :module: augur
    :func: make_parser
    :prog: augur
    :path: export

Configuration for v2
====================

Options can also be specified in a YAML configuration file supplied to
``--config`` with the following top-level keys:

.. schema-options-table:: augur/data/schema-export-v2-config.json
        
