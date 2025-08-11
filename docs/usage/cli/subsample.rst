===============
augur subsample
===============

.. contents:: Table of Contents
   :local:

Command line reference
======================

.. argparse::
   :module: augur
   :func: make_parser
   :prog: augur
   :path: subsample

Terminology
===========

.. glossary::

   sample

      This term can refer to either the process of creating a subset or the
      subset itself:

      1. **Process**: Selecting a subset of sequences from a dataset according
         to specific parameters for filtering and subsampling (e.g.
         minimum/maximum date, minimum/maximum sequence length, sample size).

         Example: *Run the focal sample …*

      2. **Resulting subset**: The set of sequences obtained from the process
         described in (1).

         Example: *The contextual sample consisted of …*

Configuration
=============

The ``--config`` option expects a YAML-formatted configuration file. This
section describes how the file should be structured.

.. tip::

    Use ``--config-section`` to read from a configuration file that puts these
    options under a specific section.

samples
-------

``samples`` must contain at least one sample.

.. code:: yaml

   samples:
     <sample 1>:
       # sample options (see below)
     <sample 2>:
       # sample options (see below)
     …

.. sample-options-schema-table:: augur/data/schema-subsample-config.json

Implementation details
======================

- Configurations containing a single :term:`sample` are run using a single call
  to :doc:`augur filter <./filter>`.

- Configurations containing multiple samples are run using multiple calls to
  :doc:`augur filter <./filter>`.

  Each sample has its own call to ``augur filter``, known as intermediate calls.
  These can run in parallel when ``--nthreads`` > 1.

  Each intermediate call uses ``--output-strains`` to write a text file
  containing the selected sequence ids for that sample.

  The output dataset is produced by a final ``augur filter`` call that uses the
  union of all sample id files to subset the input dataset.

- CLI and YAML config options map closely to augur filter options.

   The following table shows the mapping between ``augur subsample`` and ``augur
   filter`` CLI options.

   .. cli-option-table:: augur.subsample.GLOBAL_CLI_OPTIONS augur.subsample.FINAL_CLI_OPTIONS

   The following table shows the mapping between ``augur subsample`` sample
   configuration options and ``augur filter`` CLI options.

   .. yaml-option-table:: augur.subsample.SAMPLE_CONFIG

   Note that the following ``augur filter`` options are not supported:

   - ``--priority``
   - ``--output-log``
   - ``--output-group-by-sizes``
   - ``--output-strains``
   - ``--empty-output-reporting``
