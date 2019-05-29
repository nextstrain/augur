# Installation

Augur is written in Python 3 and requires at least Python 3.4.
It's published on [PyPi](https://pypi.org) as [nextstrain-augur](https://pypi.org/project/nextstrain-augur), so you can install it with `pip` like so:

    python -m pip install nextstrain-augur

You can also install from a git clone or other copy of the source code by running:

    python -m pip install .

If your system has both Python 2 and Python 3 installed side-by-side, you may need to use `python3` instead of just `python` (which often defaults to Python 2 when both Python versions are installed).

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

## With Conda

Alternatively, augur itself and all of its dependencies can be installed into a [Conda](https://conda.io/miniconda.html) environment:

    conda env create -f environment.yml

When that finishes, the enviroment needs to be activated whenever you want to use augur:

    conda activate augur

## Testing if it worked

If installation worked, you should be able to run `augur --help` and see
augur's primary help output.
