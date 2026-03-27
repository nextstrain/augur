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

   filter sample
    
      The most common type of sample, one which is generated internally via ``augur filter``.
      Configured by specifying filtering parameters (date ranges, queries etc) and (optionally) a context sample to use as the input.

   proximal sample
      
      A specific type of sample where we compare a (small) focal sample against a (large) context sample and find the closest genetic matches.
      Uses ``augur proximity`` under the hood.
         
Configuration
=============

The ``--config`` option expects a YAML-formatted configuration file. This
section describes how the file should be structured.

.. code:: yaml

   defaults:
     # default sample options
   samples:
     <sample 1>:
       # sample options
     <sample 2>:
       # sample options
     …

.. tip::

    Use ``--config-section`` to read from a configuration file that puts these
    options under a specific section.


defaults
--------

The ``defaults`` section is optional and allows you to specify common options
that apply to all filter samples. This reduces repetition when multiple filter
samples share the same criteria.

Options specified in the ``defaults`` section can be overridden by individual
samples. If both defaults and a specific sample define the same option, the
sample-specific value takes precedence.

Note that some options are only available at the sample level and cannot be
specified in defaults.

.. schema-options-table:: augur/data/schema-subsample-config.json defaultProperties

samples
-------


``samples`` must contain at least one sample.

filter sample options
_____________________

These options override any values set in the ``defaults`` section.

.. schema-options-table:: augur/data/schema-subsample-config.json filterSampleProperties

proximal sample options
_______________________

.. schema-options-table:: augur/data/schema-subsample-config.json proximalSampleProperties


Implementation details
======================

- Configurations containing a single :term:`filter sample` are run using a single call
  to :doc:`augur filter <./filter>`.

- Configurations containing multiple samples are run using multiple "intermediate" calls to
  augur commands. Each :term:`filter sample` has its own call to ``augur filter`` and each
  :term:`proximal sample` has its own call to ``augur proximity``. Each intermediate call
  will write temporary files which may include a strains list, sequences FASTA and/or
  metadata TSV. The eventual output dataset is produced by a final ``augur filter`` call
  that uses the union of requested samples.
 
- Samples may be dropped from the final output (via the ``drop_sample`` config option.
  This is useful when we wish to use samples only for as an input for another sample.
  
- As samples may depend on other samples (via the ``focal_sample`` and ``context_sample``
  config values), internally we create a graph of samples which controls the order in
  which samples are evaluated.
  
- Multithreading (via ``--nthreads``) will allow samples to run as efficiently as possible,
  according to the sample dependency graph. Proximal samples always use all available threads
  which both allows them to run as efficiently as possible (they're often computationally expensive)
  as well as minimising the memory consumption of ``augur subsample``.

- CLI and YAML config options map closely to augur filter options.

   The following table shows the mapping between ``augur subsample`` and ``augur
   filter`` CLI options.

   .. cli-option-table:: augur.subsample.FILTER_GLOBAL_CLI_OPTIONS augur.subsample.FILTER_FINAL_CLI_OPTIONS

   The following table shows the mapping between ``augur subsample`` :term:`filter sample`
   configuration options and ``augur filter`` CLI options.

   .. yaml-option-table:: augur.subsample.FILTER_SAMPLE_CONFIG

   Note that the following ``augur filter`` options are not supported:

   - ``--priority``
   - ``--output-group-by-sizes``
   - ``--output-strains``
   - ``--empty-output-reporting``

   
   The following table shows the mapping between ``augur subsample`` :term:`proximal sample`
   configuration options and ``augur proximity`` CLI options.

   .. yaml-option-table:: augur.subsample.PROXIMAL_SAMPLE_CONFIG
