## Introduction

Nextstrain is an open-source project to harness the scientific and public health potential of pathogen genome data. We provide a continually-updated view of publicly available data with powerful analytics and visualizations showing pathogen evolution and epidemic spread. Our goal is to aid epidemiological understanding and improve outbreak response.

Resulting data and inferences are available live at the website [nextstrain.org](https://nextstrain.org). Documentation is available at [nextstrain.org/docs](https://nextstrain.org/docs).

## Augur

*Definition: One held to foretell events by omens.*

Augur is the bioinformatic processing pipeline to track evolution from sequence and serological data. Documentation for augur is available at [nextstrain.org/docs/bioinformatics-pipeline](https://nextstrain.org/docs/bioinformatics-pipeline).

## Install

There are currently two parallel but separate versions of augur: __original
augur__ (or old augur) and __modular augur__ (or new augur).  They can be used
at the same time, but we're actively migrating pathogens to use the new,
modular augur and eventually old augur will be removed.  Currently, both
versions live inside this git repository, which can get pretty confusing.
These are the files associated with new, modular augur:

    augur/
      __init__.py
      align.py
      export.py
      ...
    bin/
      augur
    setup.py

and these are the files associated with the old augur:

    __init__.py
    base/
      auspice_export.py
      colorLogging.py
      config.py
      ...
    requirements.txt
    requirements-locked.txt
    scripts/
      beast-to-auspice-jsons-proof-of-principle.py
      json_tree_to_nexus.py
      plot_msa.py

New, modular augur is setup as a standard Python package with several modules
that are runnable from the command-line using the `augur` command.  Old augur
requires writing Python scripts which import functions and classes from the
`base/` directory.

Note that old augur requires Python 2.7 while new, modular augur requires
Python 3.

To install either version of augur, first clone the git repository:

```
git clone https://github.com/nextstrain/augur.git
cd augur
```

Both versions of augur have a number of Python dependencies.  The dependencies
of old augur are listed in `requirements.txt` and are best installed via a
package manager like conda or pip.  You may choose to use
`requirements-locked.txt` if you'd like a fixed set of known-good packages.

```
pip install -r requirements.txt
```

The dependencies of new, modular augur are specified in `setup.py`, which is
the standard location for Python packages.  You can install modular augur,
including the `augur` command, and all its dependencies using pip3:

```
pip3 install .
```

In addition, both versions of augur need working installations of
[mafft](https://mafft.cbrc.jp/alignment/software/) and one of the following
tree builders:

* DEFAULT: [RAxML](https://sco.h-its.org/exelixis/web/software/raxml/index.html)
* OPTIONAL: [FastTree](http://www.microbesonline.org/fasttree/)
* OPTIONAL: [IQ-TREE](http://www.iqtree.org/)

Additionally, bacterial builds require a working installation of [vcftools](https://vcftools.github.io/).

## Dev install

If you're going to be editing augur source code rather than just running augur from command line, it's recommended to install with

```
pip3 install -e .[dev]
```

to make module "editable" and install dev dependencies.


## Virus builds

Each virus build consists of a `prepare.py` and `process.py` file. Currently supported builds are listed in the [builds directory](builds/).


## Development

Development of `augur` happens at <https://github.com/nextstrain/augur>.

We currently target compatibility with Python 3.4 and higher.  This may be
increased to in the future.

Versions for this project from 3.0.0 onwards aim to follow the [Semantic
Versioning rules](https://semver.org).

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


## License and copyright

Copyright 2014-2018 Trevor Bedford and Richard Neher.

Source code to Nextstrain is made available under the terms of the [GNU Affero General Public License](LICENSE.txt) (AGPL). Nextstrain is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for more details.
