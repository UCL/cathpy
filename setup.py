#!/usr/bin/env python

from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="cathpy",
    version="0.0.1",
    author="Ian Sillitoe",
    author_email="i.sillitoe@ucl.ac.uk",
    description="CathPy - Yet Another Python Bioinformatics Toolkit",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/UCL/cathpy",
    packages=find_packages(include='cathpy.*'),
    test_suite="tests.test_seqio",
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
)
