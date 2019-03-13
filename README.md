# cathpy

## Getting Started

```sh
git clone git@github.com:UCL/cathpy.git
cd cathpy
```

## FAQ

**What is cathpy?**

`cathpy` is a python package that contains bioinformatics tools and libraries
used in [CATH](http://cathdb.info) (protein structure classification resource at UCL).

**Hmmm.. that sounds like Yet Another Python Bioinformatics Toolkit?**

Well it is... sort of.

**Should I be using it?**

If you are looking for a general Bioinformatics toolkit, you should look at BioPython first.

This project does contain generic libraries and tools, however the code is new and the
API may move around. It has been published mainly for internal use (within CATH) however
it has been release as open source in case others find the tools helpful.

**So, why doesn't this use BioPython?**

We may well merge in some BioPython modules in the future. There are few features that
BioPython does not currently handle (eg regarding fairly low-level manipulation of
alignments). At some point, we may look into turning some of this code into suggestions
/ pull-requests for BioPython.

## References

This code base contains external tools that are not written and maintained by the authors
of this project. If you use the results of these tools, please reference the relevant papers.

*GroupSim*

Capra JA and Singh M. Characterization and Prediction of Residues Determining 
Protein Functional Specificity. Bioinformatics, 24(13): 1473-1480, 2008.

*Scorecons*

Valdar WSJ (2002). Scoring residue conservation. Proteins: Structure, Function, 
and Genetics. 43(2): 227-241.
