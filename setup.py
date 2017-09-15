"""
Pipeline components for real-time virus analysis
Author: Richard Neher and Trevor Bedford
"""
import os
from setuptools import setup, find_packages

setup(
        name = "augur",
        version = "1.0.0",
        author = "Richard Neher and Trevor Bedford",
        author_email = "trevor@bedford.io",
        description = ("Pipeline components for real-time virus analysis"),
        license = "GNU Affero General Public License v3.0",
        keywords = "",
        url = "https://github.com/nextstrain/augur",
        packages=find_packages(exclude=['docs', 'tests']),
        classifiers=[
            "Development Status :: 5 - Production/Stable",
            "Topic :: Science",
            "License :: GNU Affero General Public License v3.0",
            ],
        )
