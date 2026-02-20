============
Installation
============

.. note::
   This is an Augur-specific installation guide. If you wish to use Nextstrain as a whole, please refer to `the Nextstrain installation guide <https://docs.nextstrain.org/en/latest/install.html>`__.

.. contents::
   :local:

Install Augur
=============

.. warning::
   There are some package managers that have entries with a similar name.

   * `Augur on PyPI <https://pypi.org/project/Augur/>`__. This is a different project. Use `nextstrain-augur <https://pypi.org/project/nextstrain-augur/>`__ instead.
   * `augur on Debian <https://tracker.debian.org/pkg/augur>`__. This is an unofficial version not maintained by the Nextstrain core team. It may lag behind official Augur releases. We do not generally recommend its use.


There are several ways to install Augur, ordered from least to most complex.

.. tabs::

   .. group-tab:: Nextstrain

      Augur is part of the Nextstrain project and is available in all :term:`Nextstrain runtimes <docs.nextstrain.org:runtime>`.

      Continue by following the :doc:`Nextstrain installation guide <docs.nextstrain.org:install>`.

      Once installed, you can use :doc:`cli:commands/shell` to run ``augur`` directly.

   .. group-tab:: Conda

      Augur can be installed using Conda or another variant. This assumes you are familiar with how to `manage Conda environments <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`__.

      .. code:: bash

         conda install -c conda-forge -c bioconda augur

      This installs Augur along with all dependencies.

   .. group-tab:: PyPI

      .. warning::
         Installing other Python packages after Augur may cause dependency incompatibilities. If this happens, re-install Augur using the command in step 1.

      Augur is written in Python 3 and requires at least Python 3.10. It's published on `PyPI <https://pypi.org>`__ as `nextstrain-augur <https://pypi.org/project/nextstrain-augur>`__.

      1. Install Augur along with Python dependencies.

         .. code:: bash

            python3 -m pip install nextstrain-augur

      2. Install other dependencies.

         .. include:: non-python-dependencies.rst

   .. group-tab:: Source

      .. warning::
         Installing other Python packages after Augur may cause dependency incompatibilities. If this happens, re-install Augur using the command in step 1.

      Augur can be installed from source. This is useful if you want to use unreleased changes. To develop Augur locally, see instructions in the "Development" tab.

      1. Install Augur along with Python dependencies.

         .. code:: bash

            git clone https://github.com/nextstrain/augur.git
            cd augur
            python3 -m pip install .

      2. Install other dependencies.

         .. include:: non-python-dependencies.rst

   .. group-tab:: Development

      .. warning::
         Installing other Python packages after Augur may cause dependency incompatibilities. If this happens, re-create the environment from scratch.

      1. Install non-Python dependencies. The instructions below are for `conda <https://docs.conda.io/projects/conda>`__.

         .. code:: bash

            NAME=augur
            conda env create -n "$NAME" -f dev_env.yml

      2. Activate the new environment.

         .. code:: bash

            conda activate "$NAME"

      3. Install Augur along with Python dependencies.

         .. code:: bash

            git clone https://github.com/nextstrain/augur.git
            cd augur
            python3 -m pip install -e .'[dev]'

Testing if it worked
====================

If installation worked, you should be able to run ``augur --help`` and see augur's primary help output.
