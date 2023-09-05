============
Installation
============

.. note::
   This is an Augur-specific installation guide. If you wish to use Nextstrain as a whole, please refer to `the Nextstrain installation guide <https://docs.nextstrain.org/en/latest/install.html>`__.

.. contents::
   :local:

Installing dependencies
=======================

Augur uses some external bioinformatics programs:

- ``augur align`` requires `mafft <https://mafft.cbrc.jp/alignment/software/>`__

- ``augur tree`` requires at least one of:

   - `IQ-TREE <http://www.iqtree.org/>`__ (used by default)
   - `RAxML <https://sco.h-its.org/exelixis/web/software/raxml/>`__ (optional alternative)
   - `FastTree <http://www.microbesonline.org/fasttree/>`__ (optional alternative)

- Bacterial data (or any VCF usage) requires `vcftools <https://vcftools.github.io/>`__

If you use Conda or Mamba, you can install them in an active environment:

.. code:: bash

   conda install -c conda-forge -c bioconda mafft raxml fasttree iqtree vcftools --yes

On macOS using `Homebrew <https://brew.sh/>`__:

.. code:: bash

   brew tap brewsci/bio
   brew install mafft iqtree raxml fasttree vcftools

On Debian/Ubuntu:

.. code:: bash

   sudo apt install mafft iqtree raxml fasttree vcftools

Other Linux distributions will likely have the same packages available, although the names may differ slightly.

Install Augur as a user
=======================

Using Mamba
-----------

This assumes you have Conda installed and an environment active. If not, refer to instructions for ambient runtime setup on `the Nextstrain installation guide <https://docs.nextstrain.org/en/latest/install.html>`__.

.. code:: bash

   conda install -c conda-forge -c bioconda augur

If you encounter environment solving errors or want a faster installation process, use `mamba <https://github.com/TheSnakePit/mamba>`__ as a drop-in replacement for conda:

.. code:: bash

   mamba install -c conda-forge -c bioconda augur

Using pip from PyPi
-------------------

Augur is written in Python 3 and requires at least Python 3.8. It's published on `PyPi <https://pypi.org>`__ as `nextstrain-augur <https://pypi.org/project/nextstrain-augur>`__, so you can install it with ``pip`` like so:

.. code:: bash

   python3 -m pip install nextstrain-augur

From source
-----------

.. code:: bash

   git clone https://github.com/nextstrain/augur.git
   python3 -m pip install .

This installs Augur along with external Python dependencies.

Install Augur as a developer
============================

.. code:: bash

   python3 -m pip install -e '.[dev]'

This installs dependencies necessary for local development.

Testing if it worked
====================

If installation worked, you should be able to run ``augur --help`` and see augur's primary help output.
