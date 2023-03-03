===============
augur translate
===============

.. contents::
    :local:

.. argparse::
    :module: augur
    :func: make_parser
    :prog: augur
    :path: translate

Example Node Data JSON
======================

Here's an example of the output node-data JSON where ``NODE_1`` has no
mutations compared to it's parent and ``NODE_2`` has multiple mutations in
multiple genes.

.. code-block:: json

    {
        "annotations": {
            "GENE_1": {
                "end": 1685,
                "seqid": "reference.gb",
                "start": 108,
                "strand": "+",
                "type": "CDS"
            },
            "GENE_2": {
                "end": 2705,
                "seqid": "reference.gb",
                "start": 1807,
                "strand": "+",
                "type": "CDS"
            },
        },
        "nodes": {
            "NODE_1": {
                "aa_muts": []
            },
            "NODE_2": {
                "aa_muts": [
                    "GENE_1": [
                        "S139N",
                        "R213K",
                        "R439G",
                        "V440A",
                        "D474N",
                        "S479W",
                        "S481T",
                        "P485L",
                        "R521K"
                    ],
                    "GENE_2": [
                        "P43S",
                        "D46N",
                        "C64R",
                        "R98K",
                        "D136G",
                        "M175V"
                    ]
                ]
            }
        },
        "reference": {
            "GENE_1": "MATLLRSLAL...",
            "GENE_2": "MAEEQARHVK..."
        }
    }
