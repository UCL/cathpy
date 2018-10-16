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

If you are looking for a general Bioinformatics toolkit, you should look at BioPython.

This project does contain generic libraries and tools, however the code is new and the 
API will move around lots. It has been published mainly for internal use (within CATH) 
until further notice.

**So, why doesn't this use BioPython?**

We may well merge in some BioPython modules in the future. There are few features that 
BioPython does not currently handle (eg regarding fairly low-level manipulation of 
alignments). At some point, we may look into turning some of this code into suggestions 
/ pull-requests for BioPython.

