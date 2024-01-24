Augur uses some external bioinformatics programs that are not available on PyPI:

- ``augur align`` requires `mafft <https://mafft.cbrc.jp/alignment/software/>`__

- ``augur tree`` requires at least one of:

   - `IQ-TREE <http://www.iqtree.org/>`__ (used by default)
   - `RAxML <https://sco.h-its.org/exelixis/web/software/raxml/>`__ (optional alternative)
   - `FastTree <http://www.microbesonline.org/fasttree/>`__ (optional alternative)

- Bacterial data (or any VCF usage) requires `vcftools <https://vcftools.github.io/>`__

If you use Conda, you can install them in an active environment:

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
