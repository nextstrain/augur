==================
augur curate
==================

This suite of commands is intended to help users curate their metadata for Nextstrain analyses.
Each subcommand is designed to be tightly scoped to a single type of data transformation so curation pipelines can be easily customized.
All subcommands share the same input and output options so a curation pipeline can begin with any subcommand and the output can be directly piped to any other curate subcommand.

.. note::
    If you need to parse metadata from a FASTA headers, please continue to use :doc:`augur parse </usage/cli/parse>`.
    The output metadata TSV and FASTA files can then be used as inputs for any augur curate subcommand.

You'll find documentation for all augur curate subcommands below.
We will continue to add more subcommands as we identify other common data curation tasks.

.. toctree::
    :maxdepth: 1

    passthru
    normalize-strings
    format-dates
    titlecase
    apply-geolocation-rules
    apply-record-annotations
    abbreviate-authors
    parse-genbank-location
    transform-strain-name
    rename

