from setuptools import setup
from pathlib    import Path

base_dir     = Path(__file__).parent.resolve()
version_file = base_dir / "augur/__version__.py"
readme_file  = base_dir / "README.md"

# Eval the version file to get __version__; avoids importing our own package
with version_file.open() as f:
    exec(f.read())

# Get the long description from the README file
with readme_file.open(encoding = "utf-8") as f:
    long_description = f.read()

setup(
        name = "nextstrain-augur",
        version = __version__,
        author = "Nextstrain developers",
        author_email = "trevor@bedford.io, richard.neher@unibas.ch",
        description = "A bioinformatics toolkit for phylogenetic analysis.",
        long_description = long_description,
        long_description_content_type = "text/markdown",
        keywords = "nextstrain, molecular epidemiology",
        url = "https://github.com/nextstrain/augur",
        project_urls = {
            "Bug Reports": "https://github.com/nextstrain/augur/issues",
            "Change Log":  "https://github.com/nextstrain/augur/blob/master/CHANGES.md#next",
            "Source":      "https://github.com/nextstrain/augur",
        },
        packages=['augur'],
        package_data={'augur': ['data/*']},
        data_files = [("", ["LICENSE.txt"])],
        python_requires = '>=3.4',
        install_requires = [
            "bcbio-gff >=0.6.4, ==0.6.*",
            "biopython >=1.73, ==1.*",
            "boto >=2.38, ==2.*",
            "cvxopt >=1.1.9, ==1.1.*",
            "ipdb >=0.10.1",
            "jsonschema ==3.0.0a3",
            "matplotlib >=2.0, ==2.*",
            "pandas >=0.23.4, ==0.23.*",
            "phylo-treetime >=0.5.3, ==0.5.*",
            "seaborn >=0.9.0, ==0.9.*",
            "snakemake >=5.1.5, ==5.*"
        ],
        extras_require={
            'dev': [
                "pylint >=1.7.6, ==1.7.*",
                "pytest >=3.2.1, ==3.*",
                "wheel >=0.32.3, ==0.32.*"
            ]
        },
        classifiers=[
            "Development Status :: 3 - Alpha",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            "License :: OSI Approved :: GNU Affero General Public License v3",

            # Python 3 only; pathlib is >=3.4
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.4",
            "Programming Language :: Python :: 3.5",
            "Programming Language :: Python :: 3.6",
            ],

        # Install an "augur" program which calls augur.__main__.main()
        #   https://setuptools.readthedocs.io/en/latest/setuptools.html#automatic-script-creation
        entry_points = {
            "console_scripts": [
                "augur = augur.__main__:main",
            ],
        },
        )
