==================================
How do I specify ``refine`` rates?
==================================

How we use refine in the zika tutorial
======================================

In the Zika tutorial we used the following basic rule to run the :doc:`../usage/cli/refine` command:

.. code-block:: python

    rule refine:
        input:
            tree = rules.tree.output.tree,
            alignment = rules.align.output,
            metadata = "data/metadata.tsv"
        output:
            tree = "results/tree.nwk",
            node_data = "results/branch_lengths.json"
        shell:
            """
            augur refine \
                --tree {input.tree} \
                --alignment {input.alignment} \
                --metadata {input.metadata} \
                --timetree \
                --output-tree {output.tree} \
                --output-node-data {output.node_data}
            """


This rule will estimate the rate of the molecular clock, reroot the tree, and estimate a time tree.
The paragraphs below will detail how to exert more control on each of these steps through additional options the refine command.


Specify the evolutionary rate
=============================

By default ``augur`` (through ``treetime``) will estimate the rate of evolution from the data by regressing divergence vs sampling date.
In some scenarios, however, there is insufficient temporal signal to reliably estimate the rate and the analysis will be more robust and reproducible if one fixes this rate explicitly.
This can be done via the flag ``--clock-rate <value>`` where the implied units are substitutions per site and year.
In our zika example, this would look like this

.. code-block:: diff

    rule refine:
        input:
            tree = rules.tree.output.tree,
            alignment = rules.align.output,
            metadata = "data/metadata.tsv"
        output:
            tree = "results/tree.nwk",
            node_data = "results/branch_lengths.json"
    +    params:
    +    	clock_rate = 0.0008
        shell:
            """
            augur refine \
                --tree {input.tree} \
                --alignment {input.alignment} \
                --metadata {input.metadata} \
                --timetree \
    +           --clock-rate {params.clock_rate} \
                --output-tree {output.tree} \
                --output-node-data {output.node_data}
            """



Confidence intervals for divergence times
=========================================

Divergence time estimates are probabilistic and uncertain for multiple reasons, primarily because the accumulation of mutations is a probabilistic process and the rate estimate itself is not precise.
Augur/TreeTime will account for this uncertainty if the refine command is run with the flag ``--date-confidence`` and the standard deviation of the rate estimate is specified.

.. code-block:: diff

    rule refine:
        input:
            tree = rules.tree.output.tree,
            alignment = rules.align.output,
            metadata = "data/metadata.tsv"
        output:
            tree = "results/tree.nwk",
            node_data = "results/branch_lengths.json"
        params:
            clock_rate = 0.0008,
    +    	clock_std_dev = 0.0002
        shell:
            """
            augur refine \
                --tree {input.tree} \
                --alignment {input.alignment} \
                --metadata {input.metadata} \
                --timetree \
                --date-confidence \
    +            --clock-rate {params.clock_rate} \
    +            --clock-std-dev {params.clock_std_dev} \
                --output-tree {output.tree} \
                --output-node-data {output.node_data}
            """

If run with these parameters, augur will save an confidence interval (e.g. ``[2014.5,2014.7]``) for each node in the tree.

By default, augur runs TreeTime in a "covariance-aware" mode where the root-to-tip regression accounts for shared ancestry and covariance between terminal nodes.
This, however, is sometimes unstable when the temporal signal is low and can be switch off with the flag ``--no-covariance``.


Specifying the root of the tree
===============================

By default, augur/TreeTime reroots your input tree to optimize the temporal signal in the data. This is robust when there is robust temporal signal.
In other situations, you might want to specify the root explicitly, specify a rerooting mechanisms, or keep the root of the input tree.
The latter can be achieved by passing the argument ``--keep-root``.
To specify a particular strain (or the common ancestor of a group of strains), pass the name(s) of the(se) strain(s) like so:

.. code-block:: bash

    --root strain1 [strain2 strain3 ...]

Other available rooting mechanisms are

  * ``least-squares`` (default): minimize squared deviation of the root-to-tip regression
  * ``min-dev``: essentially midpoint rooting minimizing the variance in root-to-tip distance
  * ``oldest``: use the oldest strain as outgroup


Polytomy resolution
===================

if the data set contains many very similar sequences, their evolutionary relationship some times remains ambiguous resulting in zero-length branches or polytomies (that is internal nodes with more than 2 children).
Augur partially resolves those polytomies if such resolution helps the make the tree fit the temporal structure in the data.
If this is undesired, this can be switched-off using ``--keep-polytomies``.
