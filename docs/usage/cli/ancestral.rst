===============
augur ancestral
===============

.. contents::
    :local:

.. argparse::
    :module: augur
    :func: make_parser
    :prog: augur
    :path: ancestral

Example Node Data JSON
======================

Here's an example of the output node-data JSON where ``NODE_1`` has no
mutations compared to it's parent and ``NODE_2`` has multiple mutations.

.. code-block:: json

    {
        "nodes": {
            "NODE_1": {
                "muts": [],
                "sequence": "TCCAAACAAAGT..."
            },
            "NODE_2": {
                "muts": [
                  "A4461G",
                  "A6591G",
                  "A9184C",
                  "A10385T",
                  "T15098C"
                ],
                "sequence": "TCCAAACAAAGT..."
            }
        }
    }
