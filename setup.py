from pathlib    import Path
import setuptools
import sys

min_version = (3, 6)

if sys.version_info < min_version:
    error = """
Beginning with augur 7.0.0, Python {0} or above is required.

This may be due to an out of date pip.

Make sure you have pip >= 9.0.1.
""".format('.'.join(str(n) for n in min_version)),
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
    version = __version__,
    author = "Nextstrain developers",
    author_email = "trevor@bedford.io, richard.neher@unibas.ch",
    description = "A bioinformatics toolkit for phylogenetic analysis",
    long_description = long_description,
    long_description_content_type = "text/markdown",
    keywords = "nextstrain, molecular epidemiology",
    url = "https://github.com/nextstrain/augur",
    project_urls = {
        "Bug Reports": "https://github.com/nextstrain/augur/issues",
        "Change Log": "https://github.com/nextstrain/augur/blob/master/CHANGES.md#next",
        "Source": "https://github.com/nextstrain/augur",
    },
    packages = setuptools.find_packages(),
    package_data = {'augur': ['data/*']},
    python_requires = '>={}'.format('.'.join(str(n) for n in min_version)),
    install_requires = [
        "bcbio-gff >=0.6.0, ==0.6.*",
        "biopython >=1.67, <=1.76",
        "jsonschema >=3.0.0, ==3.*",
        "packaging >=19.2",
        "pandas >=1.0.0, ==1.*",
        "phylo-treetime ==0.8.*",
        "xopen >=1.0.1, ==1.*"
    ],
    extras_require = {
        'full': [
            "cvxopt >=1.1.9, ==1.*",
            "matplotlib >=2.0, ==2.*",
            "seaborn >=0.9.0, ==0.9.*"
        ],
        'dev': [
            "cram >=0.7, ==0.*",
            "deepdiff >=4.3.2, ==4.3.*",
            "freezegun >=0.3.15, ==0.3.*",
            "nextstrain-sphinx-theme >=2020.3",
            "pylint >=1.7.6, ==1.7.*",
            "pytest >=5.4.1, ==5.4.*",
            "pytest-cov >=2.8.1, ==2.8.*",
            "pytest-mock >= 2.0.0, ==2.0.*",
            "recommonmark >=0.5.0, ==0.*",
            "snakemake >=5.4.0, <5.27",
            "Sphinx >=2.0.1, ==2.*",
            "sphinx-argparse >=0.2.5, ==0.*",
            "sphinx-markdown-tables >= 0.0.9",
            "sphinx-rtd-theme >=0.4.3, ==0.*",
            "wheel >=0.32.3, ==0.32.*",
            "ipdb >=0.10.1, ==0.*"
        ]
    },
    classifiers = [
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU Affero General Public License v3",

        # Python 3 only
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
    # Install an "augur" program which calls augur.__main__.main()
    #   https://setuptools.readthedocs.io/en/latest/setuptools.html#automatic-script-creation
    entry_points = {
        "console_scripts": [
            "augur = augur.__main__:main",
        ]
    }
)
