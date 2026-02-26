===============
augur proximity
===============

.. contents:: Table of Contents
   :local:

Command line reference
======================

.. argparse::
   :module: augur
   :func: make_parser
   :prog: augur
   :path: proximity

Overview
========

A common use case in outbreak investigation is to find the best set of related sequences
from All available sequences. ``augur proximity`` is designed to do just this, by finding
the nearest neighbor sequences for each sequence in a query set. 

Proximity calculations can be done as part of :doc:`augur subsample <./subsample>` which is
often more ergonomic for pipelines

.. note::
    Currently the only available method uses `Hamming distance <https://en.wikipedia.org/wiki/Hamming_distance>`__
    on nucleotide sequences (i.e. protein alignments are not currently supported).

.. note::
    Sequences must be aligned before using this tool

    
Example usage for outbreak tracking
-----------------------------------

If you have a list of outbreak strains, commonly injected into the analysis via our `multiple inputs <https://next.nextstrain.org/blog/2025-09-29-standardized-multiple-inputs>`__ support you can generate a set of samples via:

  1. Use ``augur filter`` to get the outbreak set
  2. Use ``augur proximity`` to generate a set of closely related strains, using the entire dataset as the context (background)
  3. Use ``augur filter`` to generate a (small) set of background sequences
  4. Use ``augur filter`` to merge the sets of strains produced from the above steps


  
Performance
-----------
  
  - Parallalisation is extremely efficient, use ``--nthreads auto`` (or a specific count) to process multiple query sequences in parallel.
  - Context sequences are loaded into a NumPy matrix for vectorised distance computation.
    Thus you will run out of memory is the number of context sequences is too large.
    Testing on ~500,000 influenza samples used only ~2GiB of memory.

SARS-CoV-2
----------
Due to the many millions of sequences available, this tool will not be able to search against all sequences (as they won't all fit into memory).
We suggest you downsample those using temporal filters, nextstrain clades or pango lineages to get a smaller contextual set before using this tool.



Missing data handling
---------------------

Non-ATCG characters in sequences are converted to N. The ``--missing-data``
option controls how positions with N are counted:

``none`` (default)
   N is treated as a regular base. Any position where two sequences
   differ (including N vs. a real base) counts toward the distance.

``all``
   Positions where *either* the query or context sequence has an N are ignored
   entirely (not counted as a mismatch).

``flanking``
   Runs of Ns at the start and end of each sequence are ignored. Interior Ns
   are still counted normally (i.e the same as ``none``).
