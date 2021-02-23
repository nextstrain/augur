# Installation

...there are quite a few ways to contract `augur`:

* [Using conda](#using-conda)
* [Using pip from PyPi](#using-pip-from-pypi)
* [Install from Debian repositories](#install-from-debian-repositories)
* [Install from source](#install-from-source)

Lastly... [testing if it worked](#testing-if-it-worked)

---

## Using conda

We recommend using `conda` because Python's `pip` does not protect dependencies between specific versions of Python packages.

[Install Miniconda with Python 3](https://docs.conda.io/en/latest/miniconda.html).
If you already have Miniconda installed with Python 2, download the latest Python 3 version and [follow conda's installation instructions](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).
If you already have an older Miniconda version installed with Python 3, you may need to update your installation prior to installing Nextstrain's tools with:

```bash
conda activate base
conda update conda
```

Create a conda environment to install augur into and activate that environment.

```bash
conda create -n nextstrain
conda activate nextstrain
```

Install augur and its dependencies into your environment. *(See variant below.)*

```bash
conda install -c conda-forge -c bioconda augur
```

For a much faster installation process, use [mamba](https://github.com/TheSnakePit/mamba) as a drop-in replacement for conda.

```bash
conda install -c conda-forge mamba
mamba install -c conda-forge -c bioconda augur
```

## Using pip from PyPi

Augur is written in Python 3 and requires at least Python 3.6.

It's published on [PyPi](https://pypi.org) as [nextstrain-augur](https://pypi.org/project/nextstrain-augur), so you can install it with `pip`:

```bash
python3 -m pip install nextstrain-augur
```

Augur uses some common external bioinformatics programs which you'll need to install to have a fully functioning toolkit:

* Nextstrain workflows and some tutorials require [Snakemake](https://snakemake.readthedocs.io)

* `augur align` requires [mafft](https://mafft.cbrc.jp/alignment/software/)

* `augur tree` requires at least one of:
   - [IQ-TREE](http://www.iqtree.org/) (used by default)
   - [RAxML](https://sco.h-its.org/exelixis/web/software/raxml/) (optional alternative)
   - [FastTree](http://www.microbesonline.org/fasttree/) (optional alternative)

* Bacterial data (or any VCF usage) requires [vcftools](https://vcftools.github.io/)

On __macOS__, you can install most of these external programs using [Homebrew](https://brew.sh/) with:

    brew tap brewsci/bio
    brew install mafft iqtree raxml fasttree vcftools

On __Debian/Ubuntu__, you can install them via:

    sudo apt install mafft iqtree raxml fasttree vcftools

Other Linux distributions will likely have the same packages available, although the names may differ slightly.
Follow [Snakemake's installation instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) for your operating system.


## Install from Debian repositories

If you're running a Debian-derived Linux, installing the `augur` package from the Ubuntu repositories would be the easiest way. It is provided in Ubuntu 20.10 and 21.04.

```bash
sudo apt install augur mafft
```

<!-- The missing mafft package was identified by running sudo apt install augur and checking whether mafft is included in the list of additional/new/suggested packages. --> 

(`mafft` does not get installed with `sudo apt install augur` alone, see note above in [Using pip from PyPi](#using-pip-from-pypi).)

If the step above does not work, you might have to add the appropriate Ubuntu repository to your `/etc/apt/sources.list` file:

```bash
sudo add-apt-repository http://ch.archive.ubuntu.com/ubuntu/ groovy main restricted
```
(For Ubuntu 21.04 - currently in development - replace `groovy` by `hirsute`.)

## Install from source

```bash
git clone https://github.com/nextstrain/augur.git
python3 -m pip install .
```

This install depends on a fairly minimal set of external Python libraries.
There are some functions in augur that require a larger set of dependencies.
These can be installed with:

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
