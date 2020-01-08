#!/usr/bin/env python

from setuptools import setup, find_namespace_packages

import cathpy

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="cathpy",
    version='0.3.6',
    author="Ian Sillitoe",
    author_email="i.sillitoe@ucl.ac.uk",
    description="CathPy - Python Bioinformatics Toolkit for CATH (Protein Classification).",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/UCL/cathpy",
    packages=find_namespace_packages(include=['cathpy.*']),
    include_package_data=True,
    package_data={
        'cathpy.core': [
            'tools/GroupSim/*',
            'tools/*/cath-cluster',
            'tools/*/cath-resolve-hits',
            'tools/*/cath-ssap',
            'tools/*/cath-superpose',
            'tools/*/scorecons',
            'tools/data/*',
        ],
    },
    test_suite="tests",
    scripts=[
        'scripts/cath-funfhmmer-api',
        'scripts/cath-align-scorecons',
        'scripts/cath-align-summary',
        'scripts/cath-sc-merge-alignment',
    ],
    install_requires=[
        'Click',
        'requests',
        'jsonpickle',
        'tqdm',
        'dendropy',
        'redis',
        'celery',
        'pymongo',
    ],
    entry_points={
        'console_scripts': [
            'cath-cli = cathpy.core.scripts.cath_cli:cli',
        ],
    },
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
)
