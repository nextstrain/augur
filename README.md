## Introduction

Nextstrain is an open-source project to harness the scientific and public health potential of pathogen genome data. We provide a continually-updated view of publicly available data with powerful analytics and visualizations showing pathogen evolution and epidemic spread. Our goal is to aid epidemiological understanding and improve outbreak response.

Resulting data and inferences are available live at the website [nextstrain.org](https://nextstrain.org). Documentation is available at [nextstrain.org/docs](https://nextstrain.org/docs).

## Augur

*Definition: One held to foretell events by omens.*

Augur is the bioinformatic processing pipeline to track evolution from sequence and serological data. Documentation for augur is available at [nextstrain.org/docs/bioinformatics-pipeline](https://nextstrain.org/docs/bioinformatics-pipeline).

## Install

To install augur, clone the git repository and its submodules

```
git clone https://github.com/nextstrain/augur.git
cd augur
git submodule update --init
```

Augur has a number of python dependencies that are listed in `requirements.txt` and are best installed via a package manager like conda or pip.

```
pip install -r requirements.txt
```

Augur is written in Python 2.7 and requires Python 2.7 to run. Your version of Python can be confirmed by running `python --version`.

In addition, augur needs working installations of [mafft](https://mafft.cbrc.jp/alignment/software/) and one of the following tree builders

* DEFAULT: [RAxML](https://sco.h-its.org/exelixis/web/software/raxml/index.html). You'll probably need to create a symlink `raxml -> raxmlHPC` because `augur` expects an excutable named `raxml`.
* OPTIONAL: [FastTree](http://www.microbesonline.org/fasttree/)
* OPTIONAL: [IQ-TREE](http://www.iqtree.org/)


## Virus builds

Each virus build consists of a `prepare.py` and `process.py` file. Currently supported builds are listed in the [builds directory](builds/).

## License and copyright

Copyright 2014-2018 Trevor Bedford and Richard Neher.

Source code to Nextstrain is made available under the terms of the [GNU Affero General Public License](LICENSE.txt) (AGPL). Nextstrain is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for more details.
