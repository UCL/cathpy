# Using GroupSim

##  What is GroupSim?

GroupSim can be used to predict specificity-determining residues between subgroups of sequences in a multiple-sequence alignment. The average similarity between each pair of amino acids within a subgroup is calculated using a similarity matrix to obtain a column score based on which SDS predictions are done.

The paper describing th method can be found [here](http://bioinformatics.oxfordjournals.org/cgi/content/abstract/24/13/1473) and the original code is available [here](http://compbio.cs.princeton.edu/specificity/).

### Script: `group_sim_sdp.py`

```
USAGE:
group_sim_sdp.py [options] alignfile group_id1 group_id2 ... group_idN

    -alignfile must be in fasta or clustal format.
    -group_id for each specificity group should be a unique string following
     the sequence name of member sequences in the alignment.  The start of the
     group_id should be indicated by a '|', e.g., seq1|group_id.

OPTIONS:
    -c [real in [0, 1)]
     column gap cutoff. Do not score columns that contain more than this
     fraction gaps. Default=.1

    -g [real in [0, 1)]
     group gap cutoff. Do not score columns that contain a group with more
     than this fraction gaps. Default=.3

    -h
     help. Print this message.

    -l [real in [0,1]]
     lambda for window heuristic linear combination. Default=.7

    -m [filename]
     similarity matrix file, e.g., blosum62.bla  Default=identity matrix

    -n
     map raw scores (pre-ConsWin) to [0,1]. This can be useful if the raw
     scores are on a different scale (e.g. due to a matrix).

    -o [filename]
     name of output file. Default=output to stdout

    -w [int]
     conservation window size. Number of residues on either side included in
     the window. Default=3

```
### Other variants of the script: group_sim_sdp_no_gap_cutoff.py

### Example alignment

Alignment: [example.aln](https://github.com/sayonidas03/cathpy/blob/master/tools/GroupSim/examples/example.aln)

The example alignment contains sequences from two subgroups - 1 and 2.

### Matrices

The identity matrix is used as default if no matrix is provided as an argument. In this case the groupsim scores range from 0-1.

The McLachlan (mc) matrix is useful for comparing chemical changes in AAs in MSA and is most suitable for identifying and ranking SDPs in a MSA. The GroupSim socres generated are within the range 0 - 5

### Example usage:

1. Most suitable when comparing 2 to 3 subgroups and inferring whether the groups are functionally similar or dissimilar (like in FunFHMMer). You can use either identity and McLachlan matrix.

```
python2 group_sim_sdp.py -c 0.3 -g 0.5 example.aln 1 2 > example_id_with_gap_cutoff.groupsim.txt
```
```
scored by GroupSim with identity matrix
# window params: len=3, lambda=0.70
# group gap cutoff: 0.50  column gap cutoff: 0.30
# col_num	score	column
0	None	RRRRRRRRRRRRRRRRRRRRR | --------
1	None	MMMMMMMMMMMMMMMMMMMMM | --------
2	None	FFFFFFFFFFFFFFFFFFFFF | --------
3	0.779126	VVVVVVVVVVVVVVVVVVVVV | GGGGGGGG
4	0.792105	NNNNNNNNNNNNNNNNNNNNN | LLLLLLLL
```

```
python2 group_sim_sdp.py -c 0.3 -g 0.5 -m mclachlan1972.aa.dat example.aln 1 2 > example_mc_with_gap_cutoff.groupsim.txt
```
```
scored by GroupSim with mclachlan1972.aa.dat matrix
# window params: len=3, lambda=0.70
# group gap cutoff: 0.30  column gap cutoff: 0.10
# col_num	score	column
0	None	RRRRRRRRRRRRRRRRRRRRR | --------
1	None	MMMMMMMMMMMMMMMMMMMMM | --------
2	None	FFFFFFFFFFFFFFFFFFFFF | --------
3	3.929126	VVVVVVVVVVVVVVVVVVVVV | GGGGGGGG
4	3.592105	NNNNNNNNNNNNNNNNNNNNN | LLLLLLLL
```

2. For identifying and ranking SDPs between more than 2 subgroups, it is best to use the McLachlan matrix so that the physicochemical characteristics of the AA changes are taken into consideration.

```
python2 group_sim_sdp.py -c 1 -g 1 -m mclachlan1972.aa.dat example.aln 1 2 > example_mc_without_gap_cutoff.groupsim.txt
```
```
scored by GroupSim with mclachlan1972.aa.dat matrix
# window params: len=3, lambda=0.70
# col_num	score	column
0	1.682266	RRRRRRRRRRRRRRRRRRRRR | --------
1	1.686611	MMMMMMMMMMMMMMMMMMMMM | --------
2	1.718003	FFFFFFFFFFFFFFFFFFFFF | --------
3	4.045361	VVVVVVVVVVVVVVVVVVVVV | GGGGGGGG
4	3.726776	NNNNNNNNNNNNNNNNNNNNN | LLLLLLLL
```
