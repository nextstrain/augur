from setuptools import setup
from pathlib    import Path

base_dir     = Path(__file__).parent.resolve()
version_file = base_dir / "augur/__version__.py"

# Eval the version file to get __version__; avoids importing our own package
with version_file.open() as f:
    exec(f.read())

setup(
        name = "nextstrain-augur",
        version = __version__,
        author = "nextstrain developers",
        author_email = "trevor@bedford.io, richard.neher@unibas.ch",
        description = ("Pipelines for real-time phylogenetic analysis"),
        license = "MIT",
        keywords = "nextstrain, molecular epidemiology",
        url = "https://github.com/nextstrain/augur",
        packages=['augur'],
        package_data={'augur': ['data/*']},
        python_requires = '>=3.4',
        install_requires = [
            "bcbio-gff >=0.6.4, ==0.6.*",
            "biopython >=1.69, ==1.*",
            "boto >=2.38, ==2.*",
            "cvxopt >=1.1.8, ==1.1.*",
            "ipdb >=0.10.1",
            "jsonschema ==3.0.0a1",
            "matplotlib >=2.0, ==2.*",
            "pandas >=0.17.1",
            "phylo-treetime >=0.4.1, ==0.4.*",
            "seaborn >=0.6.0, ==0.6.*",
            "snakemake >=5.1.5, ==5.*"
        ],
        extras_require={
            'dev': [
                "pylint >=1.7.6, ==1.7.*",
                "pytest >=3.2.1, ==3.*",
            ]
        },
        classifiers=[
            "Development Status :: 3 - Alpha",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            "License :: OSI Approved :: MIT License",

            # Python 3 only; pathlib is >=3.4
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.4",
            "Programming Language :: Python :: 3.5",
            "Programming Language :: Python :: 3.6",
            ],
        scripts=['bin/augur']
        )
