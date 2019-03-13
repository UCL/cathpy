# cathpy

[![Documentation Status](https://readthedocs.org/projects/cathpy/badge/?version=latest)](https://cathpy.readthedocs.io/en/latest/?badge=latest)


## Getting Started

Get code

```sh
$ git clone git@github.com:UCL/cathpy.git
$ cd cathpy
```

Install code into virtual environment

```sh
$ python3 -m venv venv
$ source venv/bin/activate
$ pip install -e .
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

## FAQ

**What is cathpy?**

`cathpy` is a python package that contains bioinformatics tools and libraries
used in [CATH](http://cathdb.info) (protein structure classification resource at UCL).

**Hmmm.. that sounds like Yet Another Python Bioinformatics Toolkit?**

Well it is... sort of.

**Should I be using it?**

If you are looking for a general Bioinformatics toolkit, you probably want to look at [BioPython](https://biopython.org/) first.

This project does contain generic libraries and tools, however the code is new and the
API may move around. It has been published mainly for internal use (within CATH), however
it has been released as open source in case others find the tools helpful.

**So, why doesn't this use BioPython?**

We may well merge in some BioPython modules in the future. There are few features that
BioPython does not currently handle (eg regarding fairly low-level manipulation of
alignments). At some point, we may look into turning some of this code into suggestions
/ pull-requests.

## References

This code base contains external tools that are not written and maintained by the authors
of this project. If you use the results of these tools, please reference the relevant papers.

*GroupSim*

Capra JA and Singh M. Characterization and Prediction of Residues Determining 
Protein Functional Specificity. Bioinformatics, 24(13): 1473-1480, 2008.

*Scorecons*

Valdar WSJ (2002). Scoring residue conservation. Proteins: Structure, Function, 
and Genetics. 43(2): 227-241.
