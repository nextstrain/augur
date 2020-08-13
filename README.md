[![Build Status](https://travis-ci.com/nextstrain/augur.svg?branch=master)](https://travis-ci.com/nextstrain/augur)
[![PyPI version](https://badge.fury.io/py/nextstrain-augur.svg)](https://pypi.org/project/nextstrain-augur/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/augur/README.html)
[![Documentation Status](https://readthedocs.org/projects/nextstrain-augur/badge/?version=latest)](https://nextstrain-augur.readthedocs.io/en/stable/?badge=latest)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL%20v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)

## About Nextstrain

Nextstrain is an open-source project to harness the scientific and public health potential of pathogen genome data.
We provide a continually-updated view of publicly available data with powerful analytics and visualizations showing pathogen evolution and epidemic spread.
Our goal is to aid epidemiological understanding and improve outbreak response.

Resulting data and inferences are available live at the website [nextstrain.org](https://nextstrain.org).

## About Augur

*Definition: One held to foretell events by omens.*

Augur is the bioinformatics toolkit we use to track evolution from sequence and serological data.
It provides a collection of commands which are designed to be composable into larger processing pipelines.

The output of augur is a series of JSONs that can be used to visualize your results using [Auspice](https://github.com/nextstrain/auspice).

## Documentation

* [Overview of how Augur fits together with other Nextstrain tools](https://nextstrain.org/docs/getting-started/introduction#open-source-tools-for-the-community)
* [Overview of Augur usage](https://nextstrain.org/docs/bioinformatics/introduction-to-augur)
* [Technical documentation for Augur](https://nextstrain-augur.readthedocs.io/en/stable/installation/installation.html)
* [Contributor guide](https://github.com/nextstrain/.github/blob/master/CONTRIBUTING.md)
* [Project board with available issues](https://github.com/orgs/nextstrain/projects/6)
* [Developer docs for Augur](./DEV_DOCS.md)

## Quickstart

[Follow instructions to install augur](https://nextstrain-augur.readthedocs.io/en/stable/installation/installation.html).
Try out an analysis of real virus data by [completing the Zika tutorial](https://nextstrain.org/docs/tutorials/zika).

## Basic Usage

All of Augur's commands are accessed through the `augur` program.
For example, to infer ancestral sequences from a tree, you'd run `augur ancestral`.
If you've installed the `nextstrain-augur` package, you can just run `augur`.
Otherwise, you can run `./bin/augur` from a copy of the source code.

```
usage: augur [-h] {parse,filter,mask,align,tree,refine,ancestral,translate,clades,traits,sequence-traits,titers,export,validate,version} ...

Augur: A bioinformatics toolkit for phylogenetic analysis.

positional arguments:
  {parse,filter,mask,align,tree,refine,ancestral,translate,clades,traits,sequence-traits,titers,export,validate,version}
    parse               Parse delimited fields from FASTA sequence names into
                        a TSV and FASTA file.
    filter              Filter and subsample a sequence set.
    mask                Mask specified sites from a VCF file.
    align               Align multiple sequences from FASTA or VCF.
    tree                Build a tree using a variety of methods.
    refine              Refine an initial tree using sequence metadata.
    ancestral           Infer ancestral sequences based on a tree.
    translate           Translate gene regions from nucleotides to amino
                        acids.
    clades              Assign clades to nodes in a tree based on amino-acid
                        or nucleotide signatures.
    traits              Infer ancestral traits based on a tree.
    sequence-traits     Annotate sequences based on amino-acid or nucleotide
                        signatures.
    titers              Annotate a tree with actual and inferred titer
                        measurements.
    export              Export JSON files suitable for visualization with
                        auspice.
    validate            Validate a set of JSON files intended for
                        visualization in auspice.
    version             Print the version of augur.

optional arguments:
  -h, --help            show this help message and exit
```

For more information on a specific command, you can run it with the `--help` option, for example, `augur tree --help`.


## License and copyright

Copyright 2014-2019 Trevor Bedford and Richard Neher.

Source code to Nextstrain is made available under the terms of the [GNU Affero General Public License](LICENSE.txt) (AGPL). Nextstrain is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for more details.
