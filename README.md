# cathpy

[![Documentation Status](https://readthedocs.org/projects/cathpy/badge/?version=latest)](https://cathpy.readthedocs.io/en/latest/?badge=latest)
[![Build Status](https://travis-ci.com/UCL/cathpy.svg?branch=master)](https://travis-ci.com/UCL/cathpy)

`cathpy` is a Bioinformatics toolkit written in Python. It is developed and maintained by the
[Orengo Group](http://orengogroup.info/) at [UCL](https://www.ucl.ac.uk/orengo-group/) and is used for maintaining the [CATH protein structure database](http://www.cathdb.info) (and associated research).


## Getting Started

The easiest way to use this code is by installing the latest version into a virtual environment via `pip`:

```sh
$ python3 -m venv venv
$ source venv/bin/activate
$ pip install cathpy
```

If everything is installed and working properly then the following should work:

```sh
$ cath-align-summary -d tests/data/funfams/
file aln_len seq_count dops gap_per
tests/data/funfams/1.10.8.10-ff-14534.reduced.sto                          69     51  61.53  12.53
tests/data/funfams/1.10.8.10-ff-15516.reduced.sto                          66    429 100.00  13.04
tests/data/funfams/1.10.8.10-ff-5069.reduced.sto                           59     14   7.81   3.15
tests/data/funfams/1.10.8.10-ff-15593.reduced.sto                          63    203  95.88  17.70
```


## Development

If you are developing the code, then get access to the latest code and create a new branch:

```sh
$ git clone git@github.com:UCL/cathpy.git
$ cd cathpy
$ git checkout -b my-awesome-new-feature
```

Install the code (as editable package) into virtual environment

```sh
$ python3 -m venv venv
$ source venv/bin/activate
$ pip install -e .
```

Write your tests, make your changes then make sure all the tests still pass:

```sh
$ pytest
```

Then push your changes back to GitHub and raise a pull request.

```sh
$ git push
```

## FAQ

**What is cathpy?**

`cathpy` is a python package that contains bioinformatics tools and libraries
used in [CATH](http://cathdb.info) (protein structure classification resource at UCL).

**Hmmm.. that sounds like Yet Another Python Bioinformatics Toolkit?**

Well it is... sort of.

**Should I be using it?**

If you are looking for a general Bioinformatics toolkit, you should look at [BioPython](https://biopython.org/) first.

The `cathpy` project does contain some generic functionality that may overlap with BioPython,
however we are definitely not trying to rewrite that library. It has been published mainly for 
internal use (within CATH), however it has been released as open source in case others find the tools helpful.

## External software

This code base contains external tools that are not written and maintained by the authors
of this project. If you use the results of these tools, please reference the relevant papers.

*GroupSim*

Capra JA and Singh M. Characterization and Prediction of Residues Determining 
Protein Functional Specificity. Bioinformatics, 24(13): 1473-1480, 2008.

*Scorecons*

Valdar WSJ (2002). Scoring residue conservation. Proteins: Structure, Function, 
and Genetics. 43(2): 227-241.

## References

The most recent paper describing the CATH protein structure database:

Ian Sillitoe, et al (2018) CATH: expanding the horizons of structure-based functional annotations for genome sequences.
Nucleic Acids Research, Volume 47, Issue D1, 08 January 2019, Pages D280–D284, https://doi.org/10.1093/nar/gky1097
