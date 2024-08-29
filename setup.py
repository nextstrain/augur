from pathlib    import Path
import setuptools
import sys

py_min_version = (3, 8)  # minimal supported python version
since_augur_version = (23, 0)  # py_min_version is required since this augur version

if sys.version_info < py_min_version:
    error = """
Beginning with augur {0}, Python {1} or above is required.
You are using Python {2}.

This may be due to an out of date pip.

Make sure you have pip >= 9.0.1.
""".format('.'.join(str(n) for n in since_augur_version),
           '.'.join(str(n) for n in py_min_version),
           '.'.join(str(n) for n in sys.version_info[:3]))
    sys.exit(error)

base_dir = Path(__file__).parent.resolve()
version_file = base_dir / "augur/__version__.py"
readme_file = base_dir / "README.md"

# Eval the version file to get __version__; avoids importing our own package
with version_file.open() as f:
    exec(f.read())

# Get the long description from the README file
with readme_file.open(encoding = "utf-8") as f:
    long_description = f.read()



setuptools.setup(
    name = "nextstrain-augur",
    version = __version__,  # noqa: F821; This is imported from version_file.
    author = "Nextstrain developers",
    author_email = "trevor@bedford.io, richard.neher@unibas.ch",
    description = "A bioinformatics toolkit for phylogenetic analysis",
    long_description = long_description,
    long_description_content_type = "text/markdown",
    keywords = "nextstrain, molecular epidemiology",
    url = "https://github.com/nextstrain/augur",
    project_urls = {
        "Bug Reports": "https://github.com/nextstrain/augur/issues",
        "Change Log": "https://github.com/nextstrain/augur/blob/-/CHANGES.md#next",
        "Source": "https://github.com/nextstrain/augur",
    },
    packages = setuptools.find_packages(),
    package_data = {'augur': ['data/*']},
    python_requires = '>={}'.format('.'.join(str(n) for n in py_min_version)),
    install_requires = [
        "bcbio-gff >=0.7.1, ==0.7.*",
        # TODO: Remove biopython >= 1.80 pin if it is added to bcbio-gff: https://github.com/chapmanb/bcbb/issues/142
        "biopython >=1.80, ==1.*",
        "cvxopt >=1.1.9, ==1.*",
        "importlib_resources >=5.3.0; python_version < '3.11'",
        "isodate ==0.6.*",
        "jsonschema >=3.0.0, ==3.*",
        "networkx >= 2.5, <4",
        "numpy ==1.*",
        "packaging >=19.2",
        "pandas >=1.0.0, ==1.*",
        "phylo-treetime >=0.11.2, <0.12",
        "pyfastx >=1.0.0, <3.0",
        "python_calamine >=0.2.0",
        "scipy ==1.*",
        "xopen[zstd] >=1.7.0, <3" # TODO: Deprecated, remove v1 support around November 2024
    ],
    extras_require = {
        'dev': [
            "cram >=0.7",
            "deepdiff >=4.3.2, <8.0.0",
            "flake8",
            "freezegun >=0.3.15",
            "mypy",
            "nextstrain-sphinx-theme >=2022.5",
            "pandas-stubs >=1.0.0, ==1.*",
            "pylint >=1.7.6",
            "pytest >=5.4.1",
            "pytest-cov >=2.8.1",
            "pytest-mock >= 2.0.0",
            "recommonmark >=0.5.0",
            "Sphinx >=2.0.1",
            "sphinx-autobuild >=2021.3.14",
            "sphinx-argparse >=0.2.5, !=0.5.0",
            "sphinx-markdown-tables >= 0.0.9",
            "sphinx-rtd-theme >=0.4.3",
            "sphinx-autodoc-typehints >=1.21.4",
            "sphinx-tabs",
            "types-jsonschema >=3.0.0, ==3.*",
            "types-setuptools",
            "wheel >=0.32.3",
            "ipdb >=0.10.1"
        ]
    },
    classifiers = [
        "Development Status :: 5 - Production/Stable",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU Affero General Public License v3",

        # Python 3 only
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    # Install an "augur" program which calls augur.__main__.main()
    #   https://setuptools.readthedocs.io/en/latest/setuptools.html#automatic-script-creation
    entry_points = {
        "console_scripts": [
            "augur = augur.__main__:main",
        ]
    }
)
