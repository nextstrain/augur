Installation
============

-  `Using conda <#using-conda>`__
-  `Using pip from PyPi <#using-pip-from-pypi>`__
-  `Install from source <#install-from-source>`__
-  `Testing if it worked <#testing-if-it-worked>`__

--------------

Using conda
-----------

`Install Miniconda with Python 3 <https://docs.conda.io/en/latest/miniconda.html>`__. If you already have Miniconda installed with Python 2, download the latest Python 3 version and `follow conda’s installation instructions <https://conda.io/projects/conda/en/latest/user-guide/install/index.html>`__. If you already have an older Miniconda version installed with Python 3, you may need to update your installation prior to installing Nextstrain’s tools with:

.. code:: sh

   conda activate base
   conda update conda

Create a conda environment to install augur into and activate that environment.

.. code:: bash

   conda create -n nextstrain
   conda activate nextstrain

Install augur and its dependencies into your environment.

.. code:: bash

   conda install -c conda-forge -c bioconda augur

For a much faster installation process, use `mamba <https://github.com/TheSnakePit/mamba>`__ as a drop-in replacement for conda.

.. code:: bash

   conda install -c conda-forge mamba
   mamba install -c conda-forge -c bioconda augur

Using pip from PyPi
-------------------

Augur is written in Python 3 and requires at least Python 3.7. It’s published on `PyPi <https://pypi.org>`__ as `nextstrain-augur <https://pypi.org/project/nextstrain-augur>`__, so you can install it with ``pip`` like so:

.. code:: bash

   python3 -m pip install nextstrain-augur

Augur uses some common external bioinformatics programs which you’ll need to install to have a fully functioning toolkit:

-  Nextstrain workflows and some tutorials require `Snakemake <https://snakemake.readthedocs.io>`__

-  ``augur align`` requires `mafft <https://mafft.cbrc.jp/alignment/software/>`__

-  ``augur tree`` requires at least one of:

   -  `IQ-TREE <http://www.iqtree.org/>`__ (used by default)
   -  `RAxML <https://sco.h-its.org/exelixis/web/software/raxml/>`__ (optional alternative)
   -  `FastTree <http://www.microbesonline.org/fasttree/>`__ (optional alternative)

-  Bacterial data (or any VCF usage) requires `vcftools <https://vcftools.github.io/>`__

On macOS, you can install most of these external programs using `Homebrew <https://brew.sh/>`__ with:

::

   brew tap brewsci/bio
   brew install mafft iqtree raxml fasttree vcftools

On Debian/Ubuntu, you can install them via:

::

   sudo apt install mafft iqtree raxml fasttree vcftools

Other Linux distributions will likely have the same packages available, although the names may differ slightly. Follow `Snakemake’s installation instructions <https://snakemake.readthedocs.io/en/stable/getting_started/installation.html>`__ for your operating system.

Install from source
-------------------

.. code:: bash

   git clone https://github.com/nextstrain/augur.git
   python3 -m pip install .

This install depends on a fairly minimal set of external Python libraries.

If you wish to also install the development dependencies, and install augur in an “editable” mode whereby changes to the source code are reflected in your version of ``augur`` then run:

.. code:: bash

   python3 -m pip install -e '.[dev]'

`See above <#using-pip-from-pypi>`__ for how to install the external bioinformatics programs which you’ll need to have a fully functioning toolkit.

Testing if it worked
--------------------

If installation worked, you should be able to run ``augur --help`` and see augur’s primary help output.
