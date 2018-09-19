# cathpy

## Getting Started

```sh
git clone git@github.com:UCL/cathpy.git
cd cathpy
```

## FAQ

**What is cathpy?**

`cathpy` is a python package that contains bioinformatics tools and libraries
used in CATH (protein structure classification at UCL).

**So this is Yet Another Python Bioinformatics Toolkit?**

Sort of. 

**Should I be using it?**

Probably not. You should probably be using BioPython.

This project does contain generic tools, however the code is new and the API 
will move around lots. It should be considered internal use (within CATH) until 
further notice.

**So, why doesn't this use BioPython?**

It may well merge in some BioPython modules in the future. There are few features that BioPython does not 
currently handle (eg regarding fairly low-level manipulation of alignments). At some point, I may look into 
turning some of this code into suggestions / pull-requests for BioPython.

