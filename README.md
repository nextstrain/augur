[![Build Status](https://travis-ci.com/nextstrain/augur.svg?branch=master)](https://travis-ci.com/nextstrain/augur)
[![PyPI version](https://badge.fury.io/py/nextstrain-augur.svg)](https://pypi.org/project/nextstrain-augur/)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL%20v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)

# Introduction

Nextstrain is an open-source project to harness the scientific and public health potential of pathogen genome data.
We provide a continually-updated view of publicly available data with powerful analytics and visualizations showing pathogen evolution and epidemic spread.
Our goal is to aid epidemiological understanding and improve outbreak response.

Resulting data and inferences are available live at the website [nextstrain.org](https://nextstrain.org).
Documentation is available at [nextstrain.org/docs](https://nextstrain.org/docs).

# Augur

*Definition: One held to foretell events by omens.*

Augur is the bioinformatics toolkit we use to track evolution from sequence and serological data.
It provides a collection of commands which are designed to be composable into larger processing pipelines.
Documentation for augur is available at [nextstrain.org/docs/bioinformatics](https://nextstrain.org/docs/bioinformatics).

## Installation

Augur is written in Python 3 and requires at least Python 3.4.
It's published on [PyPi](https://pypi.org) as [nextstrain-augur](https://pypi.org/project/nextstrain-augur), so you can install it with `pip` (or `pip3`) like so:

    pip install nextstrain-augur

You can also install from a git clone or other copy of the source code by running:

    pip install .

If your system has both Python 2 and Python 3 installed side-by-side, you may need to use `pip3` or `python3 -m pip` instead of just `pip` (which often defaults to Python 2 when both Python versions are installed).

This install depends on a fairly minimal set of external Python libraries. There are some
functions in augur that require a larger set of dependencies. These can be installed via:

    pip install .[full]

Augur uses some common external bioinformatics programs which you'll need to install to have a fully functioning toolkit:

* `augur align` requires [mafft](https://mafft.cbrc.jp/alignment/software/)

* `augur tree` requires at least one of:
   - [IQ-TREE](http://www.iqtree.org/) (used by default)
   - [RAxML](https://sco.h-its.org/exelixis/web/software/raxml/) (optional alternative)
   - [FastTree](http://www.microbesonline.org/fasttree/) (optional alternative)

* Bacterial data (or any VCF usage) requires [vcftools](https://vcftools.github.io/)

Alternatively, all these dependencies (as well as augur itself) can be installed via Conda by running:

    conda env create -f environment.yml

Once installed, Conda the enviroment need to be activated whenever augur is to be used, by running:

    conda activate augur

## Usage

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


## Development

Development of `augur` happens at <https://github.com/nextstrain/augur>.

We currently target compatibility with Python 3.4 and higher.  This may be
increased to in the future.

Versions for this project from 3.0.0 onwards aim to follow the [Semantic
Versioning rules](https://semver.org).

### Running with local changes

From within a clone of the git repository you can run `./bin/augur` to test your local changes without installing them.
(Note that `./bin/augur` is not the script that gets installed by pip as `augur`; that script is generated by the `entry_points` configuration in `setup.py`.)

You can also install augur from source as an "editable" package so that your global `augur` command always uses your local source code copy:

    pip install -e .[dev]

This is not recommended if you want to be able to compare output from a stable version of augur to a development version (e.g. comparing output of `augur` installed with pip and `./bin/augur` from your local source code).

### Testing

Run doctests and unit tests for augur from Python 3 with pytest from the top-level of the augur repository.

    pytest -c pytest.python3.ini

Or, run tests for augur from Python 2.

    pytest -c pytest.python2.ini

As tests run on the development code in the augur repository, your environment should not have an existing augur installation that could cause a conflict in pytest.

### Releasing

New releases are tagged in git using a [_signed_ tag][].  The `release` branch
should always point to the latest release tag.  Source and wheel (binary)
distributions are uploaded to [the nextstrain-augur project on
PyPi](https://pypi.org/project/nextstrain-augur).

There is a `./devel/release` script which will prepare a new release from your
local repository.  It ends with instructions for you on how to push the release
commit/tag/branch and how to upload the built distributions to PyPi.  You'll
need [a PyPi account][] and [twine][] installed to do the latter.

[_signed_ tag]: https://git-scm.com/book/en/v2/Git-Tools-Signing-Your-Work
[a PyPi account]: https://pypi.org/account/register/
[twine]: https://pypi.org/project/twine

### Travis CI

Branches and PRs are tested by Travis CI jobs configured in `.travis.yml`.

New releases, via pushes to the `release` branch, trigger a new [docker-base][]
build to keep the Docker image up-to-date.

[docker-base]: https://github.com/nextstrain/docker-base

### Building documentation

[Documentation](https://nextstrain-augur.readthedocs.io) is built using [Sphinx](http://sphinx-doc.org/) and hosted on [Read The Docs](https://readthedocs.org/).
Versions of the documentation for each augur release and git branch are available and preserved.

Building the documentation locally is useful to test changes.
First, make sure you have the development dependencies of augur installed:

    pip install '.[dev]'

This installs packages listed in the `dev` section of `extras_require` in _setup.py_ (in addition to any normal augur dependencies if necessary).

Then build the HTML output format by running:

    make -C docs html

You can see other available formats by running:

    make -C docs help

To update the API documentation after adding or removing an augur submodule, autogenerate a new API file as follows.

    sphinx-apidoc -T -f -MeT -o docs augur

Sphinx caches built documentation by default, which is generally great, but can cause the sidebar of pages to be stale.  You can clean out the cache with:

    make -C docs clean

## License and copyright

Copyright 2014-2019 Trevor Bedford and Richard Neher.

Source code to Nextstrain is made available under the terms of the [GNU Affero General Public License](LICENSE.txt) (AGPL). Nextstrain is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for more details.
