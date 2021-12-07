[![Build Status](https://travis-ci.com/nextstrain/augur.svg?branch=master)](https://travis-ci.com/nextstrain/augur)
[![PyPI version](https://badge.fury.io/py/nextstrain-augur.svg)](https://pypi.org/project/nextstrain-augur/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/augur/README.html)
[![Documentation Status](https://readthedocs.org/projects/nextstrain-augur/badge/?version=latest)](https://docs.nextstrain.org/projects/augur/en/stable/)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL%20v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02906/status.svg)](https://doi.org/10.21105/joss.02906)

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

## Quickstart

[Follow instructions to install Augur](https://docs.nextstrain.org/projects/augur/en/stable/installation/installation.html).
Try out an analysis of real virus data by [completing the Zika tutorial](https://nextstrain.org/docs/tutorials/zika).

## Documentation

* [Overview of how Augur fits together with other Nextstrain tools](https://nextstrain.org/docs/getting-started/introduction#open-source-tools-for-the-community)
* [Overview of Augur usage](https://nextstrain.org/docs/bioinformatics/introduction-to-augur)
* [Technical documentation for Augur](https://docs.nextstrain.org/projects/augur/en/stable/installation/installation.html)
* [Contributor guide](https://github.com/nextstrain/.github/blob/master/CONTRIBUTING.md)
* [Project board with available issues](https://github.com/orgs/nextstrain/projects/6)
* [Developer docs for Augur](./docs/contribute/DEV_DOCS.md)

## Citation

Huddleston J, Hadfield J, Sibley TR, Lee J, Fay K, Ilcisin M, Harkins E, Bedford T, Neher RA, Hodcroft EB, (2021). Augur: a bioinformatics toolkit for phylogenetic analyses of human pathogens. Journal of Open Source Software, 6(57), 2906, https://doi.org/10.21105/joss.02906

## License and copyright

Copyright 2014-2019 Trevor Bedford and Richard Neher.

Source code to Nextstrain is made available under the terms of the [GNU Affero General Public License](LICENSE.txt) (AGPL). Nextstrain is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for more details.
