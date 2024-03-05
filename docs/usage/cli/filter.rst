============
augur filter
============

.. contents:: Table of Contents
   :local:

----

.. argparse::
    :module: augur
    :func: make_parser
    :prog: augur
    :path: filter

Guides
======

Below are some examples of using ``augur filter`` to sample data.

Filtering
---------

The filter command allows you to select various subsets of your input data for different types of analysis.
A simple example use of this command would be

.. code-block:: bash

  augur filter \
    --sequences data/sequences.fasta \
    --metadata data/metadata.tsv \
    --min-date 2012 \
    --output-sequences filtered_sequences.fasta \
    --output-metadata filtered_metadata.tsv

This command will select all sequences with collection date in 2012 or later.
The filter command has a large number of options that allow flexible filtering for many common situations.
One such use-case is the exclusion of sequences that are known to be outliers (e.g. because of sequencing errors, cell-culture adaptation, ...).
These can be specified in a separate text file (e.g. ``exclude.txt``):

.. code-block::

  BRA/2016/FC_DQ75D1
  COL/FLR_00034/2015
  ...

To drop such strains, you can pass the filename to ``--exclude``:

.. code-block:: bash

  augur filter \
    --sequences data/sequences.fasta \
    --metadata data/metadata.tsv \
    --min-date 2012 \
    --exclude exclude.txt \
    --output-sequences filtered_sequences.fasta \
    --output-metadata filtered_metadata.tsv

Subsampling within ``augur filter``
-----------------------------------

Another common filtering operation is subsetting of data to a achieve a more even spatio-temporal distribution or to cut-down data set size to more manageable numbers.
The filter command allows you to select a specific number of sequences from specific groups, for example one sequence per month from each country:

.. code-block:: bash

  augur filter \
    --sequences data/sequences.fasta \
    --metadata data/metadata.tsv \
    --min-date 2012 \
    --exclude exclude.txt \
    --group-by country year month \
    --sequences-per-group 1 \
    --output-sequences subsampled_sequences.fasta \
    --output-metadata subsampled_metadata.tsv

Subsampling using multiple ``augur filter`` commands
----------------------------------------------------

There are some subsampling strategies in which a single call to ``augur filter``
does not suffice. One such strategy is "tiered subsampling". In this strategy,
mutually exclusive sets of filters, each representing a "tier", are sampled with
different subsampling rules. This is commonly used to create geographic tiers.
Consider this subsampling scheme:

    Sample 100 sequences from Washington state and 50 sequences from the rest of the United States.

This cannot be done in a single call to ``augur filter``. Instead, it can be
decomposed into multiple schemes, each handled by a single call to ``augur
filter``. Additionally, there is an extra step to combine the intermediate
samples.

    1. Sample 100 sequences from Washington state.
    2. Sample 50 sequences from the rest of the United States.
    3. Combine the samples.

Calling ``augur filter`` multiple times
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A basic approach is to run the ``augur filter`` commands directly. This works
well for ad-hoc analyses.

.. code-block:: bash

  # 1. Sample 100 sequences from Washington state
  augur filter \
    --sequences sequences.fasta \
    --metadata metadata.tsv \
    --query "state == 'WA'" \
    --subsample-max-sequences 100 \
    --output-strains sample_strains_state.txt

  # 2. Sample 50 sequences from the rest of the United States
  augur filter \
    --sequences sequences.fasta \
    --metadata metadata.tsv \
    --query "state != 'WA' & country == 'USA'" \
    --subsample-max-sequences 50 \
    --output-strains sample_strains_country.txt

  # 3. Combine using augur filter
  augur filter \
    --sequences sequences.fasta \
    --metadata metadata.tsv \
    --exclude-all \
    --include sample_strains_state.txt \
              sample_strains_country.txt \
    --output-sequences subsampled_sequences.fasta \
    --output-metadata subsampled_metadata.tsv

Each intermediate sample is represented by a strain list file obtained from
``--output-strains``. The final step uses ``augur filter`` with ``--exclude-all``
and ``--include`` to sample the data based on the intermediate strain list
files. If the same strain appears in both files, ``augur filter`` will only
write it once in each of the final outputs.

Generalizing subsampling in a workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The approach above can be cumbersome with more intermediate samples. To
generalize this process and allow for more flexibility, a workflow management
system can be used. The following examples use `Snakemake`_.

1. Add a section in the `config file`_.

  .. code-block:: yaml

    subsampling:
      state: --query "state == 'WA'" --subsample-max-sequences 100
      country: --query "state != 'WA' & country == 'USA'" --subsample-max-sequences 50

2. Add two rules in a `Snakefile`_. If you are building a standard Nextstrain
   workflow, the output files should be used as input to sequence alignment. See
   :doc:`docs.nextstrain.org:learn/parts` to learn more about the placement of
   this step within a workflow.

  .. code-block:: python

    # 1. Sample 100 sequences from Washington state
    # 2. Sample 50 sequences from the rest of the United States
    rule intermediate_sample:
        input:
            metadata = "data/metadata.tsv",
        output:
            strains = "results/sample_strains_{sample_name}.txt",
        params:
            augur_filter_args = lambda wildcards: config.get("subsampling", {}).get(wildcards.sample_name, "")
        shell:
            """
            augur filter \
                --metadata {input.metadata} \
                {params.augur_filter_args} \
                --output-strains {output.strains}
            """

    # 3. Combine using augur filter
    rule combine_intermediate_samples:
        input:
            sequences = "data/sequences.fasta",
            metadata = "data/metadata.tsv",
            intermediate_sample_strains = expand("results/sample_strains_{sample_name}.txt", sample_name=list(config.get("subsampling", {}).keys()))
        output:
            sequences = "results/subsampled_sequences.fasta",
            metadata = "results/subsampled_metadata.tsv",
        shell:
            """
            augur filter \
                --sequences {input.sequences} \
                --metadata {input.metadata} \
                --exclude-all \
                --include {input.intermediate_sample_strains} \
                --output-sequences {output.sequences} \
                --output-metadata {output.metadata}
            """

3. Run Snakemake targeting the second rule.

  .. code-block:: bash

    snakemake combine_intermediate_samples

Explanation:

- The configuration section consists of one entry per intermediate sample in the
  format ``sample_name: <augur filter arguments>``.
- The first rule is run once per intermediate sample using `wildcards`_ and an
  `input function`_. The output of each run is the sampled strain list.
- The second rule uses `expand()`_ to define input as all the intermediate
  sampled strain lists, which are passed directly to ``--include`` as done in
  the previous example.

It is easy to add or remove intermediate samples. The configuration above can be
updated to add another tier in between state and country:

  .. code-block:: yaml

    subsampling:
      state: --query "state == 'WA'" --subsample-max-sequences 100
      neighboring_states: --query "state in {'CA', 'ID', 'OR', 'NV'}" --subsample-max-sequences 75
      country: --query "country == 'USA' & state not in {'WA', 'CA', 'ID', 'OR', 'NV'}" --subsample-max-sequences 50

.. _Snakemake: https://snakemake.readthedocs.io/en/stable/index.html
.. _config file: https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#snakefiles-standard-configuration
.. _Snakefile: https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html
.. _wildcards: https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#wildcards
.. _input function: https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#snakefiles-input-functions
.. _expand(): https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#the-expand-function
