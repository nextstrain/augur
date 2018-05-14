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
        classifiers=[
            "Development Status :: 3 - Alpha",
            "Topic :: Science",
            "License :: OSI Approved :: MIT License",
            ],
        scripts=['bin/augur']
        )
