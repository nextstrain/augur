Augur uses some external bioinformatics programs that are not available on PyPI:

- ``augur align`` requires `mafft <https://mafft.cbrc.jp/alignment/software/>`__

- ``augur tree`` requires at least one of:

   - `IQ-TREE <http://www.iqtree.org/>`__ (used by default)
   - `RAxML <https://cme.h-its.org/exelixis/web/software/raxml/>`__ (optional alternative)
   - `FastTree <http://www.microbesonline.org/fasttree/>`__ (optional alternative)

- ``augur merge`` requires:

    - ``sqlite3``, the `SQLite <https://sqlite.org>`__ CLI (version ≥3.39) for metadata
    - ``seqkit``, the `SeqKit program <https://bioinf.shenwei.me/seqkit/>`__, for sequences

- Bacterial data (or any VCF usage) requires `vcftools <https://vcftools.github.io/>`__

.. note::

   For local development, tests also require `tsv-utils <https://opensource.ebay.com/tsv-utils/>`__.

If you use Conda, you can install them in an active environment:

.. code:: bash

   conda install -c conda-forge -c bioconda mafft raxml fasttree iqtree vcftools sqlite seqkit --yes

On macOS using `Homebrew <https://brew.sh/>`__:

.. code:: bash

   brew tap brewsci/bio
   brew install mafft iqtree raxml fasttree vcftools sqlite seqkit

On Debian/Ubuntu:

.. code:: bash

   sudo apt install mafft iqtree raxml fasttree vcftools sqlite3

Other Linux distributions will likely have the same packages available, although the names may differ slightly.

The `SeqKit download page <https://bioinf.shenwei.me/seqkit/download/>`__ provides Linux binaries.
