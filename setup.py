#!/usr/bin/env python

from setuptools import setup

import cathpy

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="cathpy",
    version=cathpy.__version__,
    author="Ian Sillitoe",
    author_email="i.sillitoe@ucl.ac.uk",
    description="CathPy - Python Bioinformatics Toolkit for CATH (Protein Classification).",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/UCL/cathpy",
    packages=['cathpy'],
    test_suite="tests",
    package_dir={
        'cathpy': 'cathpy',
    },
    package_data={
        'cathpy': [
            'tools/GroupSim/*',
            'tools/*/scorecons',
            'tools/data/*',
        ],
    },
    scripts=[
        'scripts/cath-funfhmmer-api',
        'scripts/cath-align-scorecons',
        'scripts/cath-align-summary',
        'scripts/cath-sc-merge-alignment',
    ],
    install_requires=[
        'requests',
        'jsonpickle',
        'tqdm',
    ],
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
)
