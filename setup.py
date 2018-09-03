import os
from setuptools import setup

setup(
        name = "augur",
        version = "0.1.0",
        author = "nextstrain developers",
        author_email = "trevor@bedford.io, richard.neher@unibas.ch",
        description = ("Pipelines for real-time phylogenetic analysis"),
        license = "MIT",
        keywords = "nextstrain, molecular epidemiology",
        url = "https://github.com/nextstrain/augur",
        packages=['augur'],
        package_data={'augur': ['data/*']},
        install_requires = [
            "bcbio-gff >=0.6.4, ==0.6.*",
            "biopython >=1.69, ==1.*",
            "boto >=2.38, ==2.*",
            "cvxopt >=1.1.8, ==1.1.*",
            "ipdb >=0.10.1",
            "jsonschema ==3.0.0a1",
            "matplotlib >=2.0, ==2.*",
            "pandas >=0.17.1",
            "phylo-treetime >=0.5.1",
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
            "Topic :: Science",
            "License :: OSI Approved :: MIT License",
            ],
        scripts=['bin/augur']
        )
