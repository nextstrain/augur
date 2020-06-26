# Installation

* [Using pip from PyPi](#using-pip-from-pypi)
* [Using Conda](#using-conda)
* [Install from source](#install-from-source)
* [Testing if it worked](#testing-if-it-worked)

---

## Using pip from PyPi

Augur is written in Python 3 and requires at least Python 3.6.
It's published on [PyPi](https://pypi.org) as [nextstrain-augur](https://pypi.org/project/nextstrain-augur), so you can install it with `pip` like so:

```bash
python3 -m pip install nextstrain-augur
```

Augur uses some common external bioinformatics programs which you'll need to install to have a fully functioning toolkit:

* `augur align` requires [mafft](https://mafft.cbrc.jp/alignment/software/)

* `augur tree` requires at least one of:
   - [IQ-TREE](http://www.iqtree.org/) (used by default)
   - [RAxML](https://sco.h-its.org/exelixis/web/software/raxml/) (optional alternative)
   - [FastTree](http://www.microbesonline.org/fasttree/) (optional alternative)

* Bacterial data (or any VCF usage) requires [vcftools](https://vcftools.github.io/)

On macOS, you can install these external programs using [Homebrew](https://brew.sh/) with:

    brew install mafft iqtree raxml fasttree vcftools

On Debian/Ubuntu, you can install them via:

    sudo apt install mafft iqtree raxml fasttree vcftools

Other Linux distributions will likely have the same packages available, although the names may differ slightly.

## Using Conda

Alternatively, augur itself and all of its dependencies can be installed into a [Conda](https://conda.io/miniconda.html) environment:

    conda env create -f environment.yml

> _By default this environment is named "augur" but you can change that by providing a name to the above command with `-n <your-env-name>`_

When that finishes, the enviroment needs to be activated whenever you want to use augur:

    conda activate augur

## Install from source

```bash
git clone https://github.com/nextstrain/augur.git
python3 -m pip install .
```

This install depends on a fairly minimal set of external Python libraries.
There are some functions in augur that require a larger set of dependencies.
These can be installed via:

```bash
python3 -m pip install '.[full]'
```

If you wish to also install the development dependencies, and install augur in an "editable" mode whereby changes to the source code are reflected in your version of `augur` then run:

```bash
python3 -m pip install -e '.[dev]'
```

[See above](#using-pip-from-pypi) for how to install the external bioinformatics programs which you'll need to have a fully functioning toolkit.


## Testing if it worked

If installation worked, you should be able to run `augur --help` and see
augur's primary help output.
