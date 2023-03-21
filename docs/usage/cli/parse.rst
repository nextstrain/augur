===========
augur parse
===========

.. contents:: Table of Contents
   :local:

----

.. argparse::
    :module: augur
    :func: make_parser
    :prog: augur
    :path: parse
        


Example: how to parse metadata from fasta-headers
=================================================

If you download sequence data from data bases like GISAID or fludb, there often is an option to include meta data such as dates into the header of fasta files.
This might for example look like this:

..

    >A/Canoas/LACENRS_1793/2015|A|H3N2|07/17/2015||Brazil|Human|KY925125
    ATG...
    >A/Canoas/LACENRS_773/2015|A|H3N2|05/06/2015||Brazil|Human|KY925599
    ATG...
    [...]

The fasta header contains information such as influenza lineage, dates (in an unpreferred format), country, etc...
To turn this metadata into a table, augur has a special command called `parse`.
A rule to parse the above file could look like this:

.. code-block:: python

    rule parse:
        input:
            sequences = "data/h3n2_ha.fasta"
        output:
            sequences = "results/sequences_h3n2_ha.fasta",
            metadata = "results/metadata_h3n2_ha.tsv"
        params:
            fields = "strain type subtype date season country host accession"
        shell:
            """
            augur parse \
                --sequences {input.sequences} \
                --fields {params.fields} \
                --output-sequences {output.sequences} \
                --output-metadata {output.metadata} \
                --fix-dates monthfirst
            """


Note the additional argument ``--fix-dates monthfirst``.
This triggers an attempt to parse these dates and turn them into ISO format assuming that the month preceeds the date in the input data.
Note that this is a brittle process that should be spot-checked.

