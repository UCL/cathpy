#!/usr/bin/env python

from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="cathpy",
    version="0.0.1",
    author="Ian Sillitoe",
    author_email="i.sillitoe@ucl.ac.uk",
<<<<<<< HEAD
    description="CathPy - Python Bioinformatics Toolkit (CATH Protein Classification)",
=======
    description="CathPy - Python Bioinformatics Toolkit for CATH resource (protein structure database).",
>>>>>>> 164c94286f2b3a3a59a5b686fa9b45e3cbcb355b
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/UCL/cathpy",
    packages=find_packages(include='cathpy.*'),
    test_suite="tests",
    scripts=[
        'scripts/cath-funfhmmer-api',
        'scripts/cath-align-scorecons',
        'scripts/cath-align-summary',
        'scripts/cath-sc-merge-alignment',
    ],
    install_requires=[
        'requests',
    ],
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
)
