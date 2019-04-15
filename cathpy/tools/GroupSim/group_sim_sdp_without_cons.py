#!/usr/bin/env python
"""
 group_sim_sdp.py - Copyright Tony Capra 2008 - Created: 05/13/08

 Change log: 
 09/11/09 - fixed alignment parser bug that required new line at end
            of alignment
 06/13/08 - changed format for group_id specificication to
            fix group_id prefix bug
 05/13/08 - created


 This program supports the paper: 
 Capra JA and Singh M. (2008) Characterization and Prediction of
 Residues Determining Protein Functional Specificity. Bioinformatics,
 24(13): 1473-1480, 2008.

 Please cite the paper if you use this code.

 See usage() and the README for more information and examples.

 Please contact the authors with any questions or bug reports.

 -----------------------------------------------------------------------------

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

"""

import sys
import math
import getopt


# Useful Definitions
amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M',
	       'F', 'P', 'S', 'T', 'W', 'Y', 'V', '-']
amino_acids_nogap = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K',
		     'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

iupac_alphabet = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M",
		  "N", "P", "Q", "R", "S", "T", "U", "V", "W", "Y", "Z", "X",
		  "*", "-"] 


ident = [] # identity matrix
for i in range(21):
    t = []
    for j in range(21):
	if i != j:
	    t.append(0.)
	else:
	    t.append(1.)
    ident.append(t)


def usage():
    usage_string = """
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
"""

    print >> sys.stderr, usage_string



################################################################################
# File IO functions
################################################################################

def read_scoring_matrix(sm_file):
    """ Read in a scoring matrix from a file, e.g., blosum62.bla, and return it
    as a list of lists. """
    aa_index = 0
    first_line = 1
    row = []
    list_sm = [] # hold the matrix in list form

    try:
	matrix_file = open(sm_file, 'r')

	for line in matrix_file:

	    if line[0] != '#' and first_line:
		first_line = 0

	    elif line[0] != '#' and first_line == 0:
		if len(line) > 1:
		    row = line.split()
		    list_sm.append(row)

    except IOError, e:
	#print >> sys.stderr, \
	    #"Could not load similarity matrix: %s. Using identity matrix...\n" \
	    #% sm_file
	return ident
	
    # if matrix is stored in lower tri form, copy to upper
    if len(list_sm[0]) < 20:
	for i in range(0,19):
	    for j in range(i+1, 20):
		list_sm[i].append(list_sm[j][i])

    for i in range(len(list_sm)):
	for j in range(len(list_sm[i])):
	    list_sm[i][j] = float(list_sm[i][j])

    return list_sm

def read_fasta_alignment(filename):
    """ Read in the alignment stored in the FASTA file, filename. Return two
    lists: the identifiers and sequences. """

    f = open(filename)

    names = []
    alignment = []
    cur_seq = ''

    for line in f:
        if line[-1].upper() not in iupac_alphabet: 
            line = line[:-1]
	if len(line) == 0: continue

	if line[0] == ';': continue
	if line[0] == '>':
	    names.append(line[1:])

	    if cur_seq != '':
		alignment.append(cur_seq.replace('B', 'D').replace('Z', 'Q')\
				 .replace('X', '-'))
		cur_seq = ''
	elif line[0] in iupac_alphabet:
	    cur_seq += line

    # add the last sequence
    alignment.append(cur_seq.replace('B', 'D').replace('Z', 'Q')\
		     .replace('X', '-'))

    return names, alignment
	
def read_clustal_alignment(filename):
    """ Read in the alignment stored in the CLUSTAL file, filename. Return
    two lists: the names and sequences. """

    names = []
    alignment = []

    f = open(filename)

    for line in f:
        if line[-1].upper() not in iupac_alphabet: 
            line = line[:-1]
	if len(line) == 0: continue
	if '*' in line: continue

	if 'CLUSTAL' in line: continue

	t = line.split()

	if len(t) == 2 and t[1][0] in iupac_alphabet:
	    if t[0] not in names:
		names.append(t[0])
		alignment.append(t[1].upper().replace('B', 'D')\
				 .replace('Z', 'Q').replace('X', '-'))
	    else:
		alignment[names.index(t[0])] += t[1].upper().replace('B', 'D')\
		.replace('Z', 'Q').replace('X','-')
		   
    return names, alignment



################################################################################
# Column Utilities
################################################################################

def get_column(col_num, alignment):
    """Return the col_num column of alignment as a list. alignment should
    be a list of lists/strings."""
    col = []
    for seq in alignment:
	if col_num < len(seq): col.append(seq[col_num])

    return col

def format_column(col_num, alignment, gid_to_seqs):
    """ Return a string with | dividing the residues of col_num column in
    different specificity groups as defined by the gid_to_seqs dictionary. """

    col = get_column(i, alignment)

    groups = []
    for gid in gids_to_seqs:
	groups.append('')
	for seq_num in gids_to_seqs[gid]:
	    groups[-1] += col[seq_num]

    return " | ".join(groups)

def gap_percentage(col):
    """Return the percentage of gaps in col."""
    num_gaps = 0.

    for aa in col:
	if aa == '-': num_gaps += 1

    return num_gaps / len(col)



################################################################################
# Functions for the ConsWin heuristic
################################################################################

def min_ignore(l, ignore_val):
    """ Return the min of a list ignoring a value."""
    m = max(l)
    for x in l:
        if x != ignore_val and x < m:
	    m = x

    return m

def norm_01(score_list):
    """ Map the values in score_list to the range [0,1]. Ignore None."""

    new_sl = []

    mini = min_ignore(score_list, None)
    maxi = max(score_list)

    range = float(maxi - mini)
    offset = mini

    for s in score_list:
	if s == None:
	    new_sl.append(None)
	elif s == None:
	    new_sl.append(None)
	else:
	    if range == 0:
		new_sl.append(1.)
	    else:
		new_sl.append((s - offset) / range)

    return new_sl

def freq_count(col, symbols, pc=.000001):
    """ Return the symbol frequency count for a column. """

    sym_num = 0
    if pc != 0:
	freq_counts = len(symbols) * [pc] # in aa order defined above
    else:
	freq_counts = len(symbols) * [0.]

    for sym in symbols:
	for j in range(len(col)):
	    if col[j] == sym:
		freq_counts[sym_num] += 1

	sym_num += 1

    num_counts = float(sum(freq_counts))
    if num_counts == 0: return len(symbols) * [0.]
    for j in range(len(freq_counts)):
	freq_counts[j] = freq_counts[j] / num_counts

    return freq_counts

def js_divergence(col1, col2=[]):
    """ Return the Jensen-Shannon Divergence for the column with the
    background distribution bd. If bd isn't specified, BLOSUM62 is
    used."""

    fc1 = freq_count(col1, amino_acids_nogap, .001)

    if col2 != []:
	fc2 = freq_count(col2, amino_acids_nogap, .001)
    else:
	# BLOSUM62 background distribution
	fc2 = [0.074, 0.052, 0.045, 0.054, 0.025, 0.034, 0.054, 0.074,
	       0.026, 0.068, 0.099, 0.058, 0.025, 0.047, 0.039, 0.057,
	       0.051, 0.013, 0.032, 0.073]

    if len(fc1) != len(fc2): return -1

    # make r distriubtion
    r = []
    for i in range(len(fc1)):
	r.append(.5 * fc1[i] + .5 * fc2[i])

    d = 0.
    for i in range(len(fc1)):
	if r[i] != 0.0:
	    if fc1[i] == 0.0:
		d += fc2[i] * math.log(fc2[i]/r[i], 2)
	    elif fc2[i] == 0.0:
		d += fc1[i] * math.log(fc1[i]/r[i], 2) 
	    else:
		d += (fc1[i] * math.log(fc1[i]/r[i], 2) + fc2[i] *
		      math.log(fc2[i]/r[i], 2))

    d /= 2 

    return (1 - gap_percentage(col1)) * d

def cons_window_score(sdp_scores, alignment, window_len, lam):
    """ Calculate the conservation (using JS-divergence) of each column in
    alignment and linearly combine the average conservation score over
    a window of size window_len on either side with the
    sdp_score. Weight the combination by lam. """
    
    win_scores = sdp_scores[:]
    cons_scores = []

    columns = []
    for i in range(len(alignment[0])):
	columns.append(get_column(i, alignment))

    for i, col in enumerate(columns):

	if sdp_scores[i] == None:
	    cons_scores.append(None)
	else:
	    col_string = "".join(col)
	    cons_scores.append(js_divergence(col_string))

    wincon_scores = []
    for i, sdp_score in enumerate(sdp_scores):
	if sdp_score == None: 
	    wincon_scores.append(None)
	    continue
	
	sum = 0.
	num_terms = 0.

	for j in range(max(i - window_len, 0),
		       min(i + window_len + 1, len(sdp_scores))):
	    if i != j and cons_scores[j] >= 0:
		num_terms += 1
		sum += cons_scores[j]
		
	if num_terms != 0:
	    wincon_scores.append(sum / num_terms)
	else:
	    wincon_scores.append(0.)

    wincon_scores = norm_01(wincon_scores)

    for i, sdp_score in enumerate(sdp_scores):
	if sdp_score != None and wincon_scores[i] != None:
	    win_scores[i] = (1 - lam) * wincon_scores[i] + lam * sdp_scores[i]
	    #print ("%d %.3f %.3f %.3f -- %.3f\t%s"
		#   % (i, sdp_score, cons_scores[i], wincon_scores[i],
		#      win_scores[i], "".join(columns[i])))

    return win_scores



################################################################################
# GroupSim function
################################################################################

def matrix_sum_pairs(groups, matrix=ident):
    """Give a score of + matrix(aa1, aa2) to every pair of aa in the same
    group and -matrix(aa1, aa2) for every pair of aa in different
    groups. Return the sum of the pair scores."""

    num_terms = 0.
    sum_within = 0.
    sum_between = 0.

    for group in groups:
	local_sum = 0.
	num_terms = 0.
	for i in range(len(group)):
	    for j in range(i):
		try:
		    i1 = amino_acids.index(group[i])
		    i2 = amino_acids.index(group[j])
		    local_sum += matrix[i1][i2]
		    num_terms += 1.
		except ValueError:
		    pass

        if num_terms > 0: 
            sum_within += (local_sum / num_terms)

    sum_within /= len(groups)

    num_pairs = 0.
    for g1_ind in range(len(groups)):
	for g2_ind in range(g1_ind):
	    num_terms = 0.
	    pair_sum = 0.
	    num_pairs += 1
	    for aa1 in groups[g1_ind]:
		for aa2 in groups[g2_ind]:
		    try:
			i1 = amino_acids.index(aa1)
			i2 = amino_acids.index(aa2)
			pair_sum += matrix[i1][i2]
			num_terms += 1.
		    except ValueError:
			# catch gaps
			pass

	    sum_between += (pair_sum / num_terms)

    sum_between /= num_pairs

    gr_str = ["".join(x) for x in groups]
    return (1. - gap_percentage("".join(gr_str))) * (sum_within - sum_between)
    # return math.sqrt( sum_within * (1.0 - sum_between))




################################################################################
# Begin Execution
################################################################################

# set parameter defaults - see usage()
group_gap_cutoff = .3
column_gap_cutoff = .1

cons_win_len = 3     # 0 = no ConsWin
lamb = .7     # for ConsWin linear combination

map_scores_to_01 = 0      # 1 = map raw scores to [0,1]

matrix_file = ''
outfile_name = ''

# parse options and args -- see usage()
try:
    opts, args = getopt.getopt(sys.argv[1:], "hno:m:w:l:c:g:")
except getopt.GetoptError:
    usage()
    sys.exit(1)

if len(args) < 3: 
    usage()
    sys.exit(2)


align_file = args[0]
group_ids = args[1:]

for opt, arg in opts:
    if opt == "-h":
	usage()
	sys.exit()

    elif opt == "-n":
	map_scores_to_01 = 1

    elif opt == "-o":
	outfile_name = arg

    elif opt == "-m":
	matrix_file = arg

    elif opt == "-w":
	try:
	    int(arg)
	    cons_win_len = int(arg)
	except ValueError:
	    print >> sys.stderr, \
		"ERROR: Window size must be an integer.",
	    print >> sys.stderr, "Using window_size %d...\n" % cons_win_len

    elif opt == "-l":
	try:
	    if not (0. <= float(arg) <= 1.): 
		raise ValueError
	    else:
		lamb = float(arg)
	except ValueError:
	    print >> sys.stderr, \
		"ERROR: ConsWin lambda must be a real in [0,1].",
	    print "Using lambda = %.1f...\n" % lamb

    elif opt == "-c":
	try:
	    if not (0. <= float(arg) < 1.): 
		raise ValueError
	    else:
		column_gap_cutoff = float(arg)
	except ValueError:
	    print >> sys.stderr, \
		"ERROR: Column gap cutoff must be a real in [0,1). ",
	    print >> sys.stderr, "Using a gap cutoff of %.1f...\n" \
		% column_gap_cutoff

    elif opt == "-g":
	try:
	    if not (0. <= float(arg) < 1.): 
		raise ValueError
	    else:
		group_gap_cutoff = float(arg)
	except ValueError:
	    print >> sys.stderr, \
		"ERROR: Group gap cutoff must be a real in [0,1).",
	    print "Using a gap cutoff of %.1f...\n" % group_gap_cutoff



# read matrix and alignment files
matrix = read_scoring_matrix(matrix_file)
if matrix == ident: matrix_file = 'identity'


try:
    names, alignment = read_clustal_alignment(align_file)
    if names == []:
	names, alignment = read_fasta_alignment(align_file)
except IOError, e:
    sys.exit("Problem reading %s. Exiting..." % align_file)

if len(alignment) != len(names) or alignment == []:
	sys.exit("Unable to parse alignment.\n")

#print >> sys.stderr, "Read %s.\n" % align_file



# make group to sequence number mapping
# group seqs may not be contiguous in alignment
gids_to_seqs = {}
for gid in group_ids:
    gids_to_seqs[gid] = []
    
for i, seq_name in enumerate(names):
    seq_matched = 0
    for gid in group_ids:
	if gid == seq_name.split('|')[-1]: 
	    gids_to_seqs[gid].append(i)
	    seq_matched += 1

    if seq_matched > 1:
	sys.exit("ERROR: %s matches multiple group names.\n" % seq_name)

for gid in gids_to_seqs:
    if len(gids_to_seqs[gid]) == 0:
        sys.exit("ERROR: No sequences found for specficity group \"%s\".\n" % gid)


# score columns with GroupSim
scores = []
for i in range(len(alignment[0])):
    col = get_column(i, alignment)

    groups = []
    for gid in gids_to_seqs:
	groups.append('')
	for seq_num in gids_to_seqs[gid]:
	    groups[-1] += col[seq_num]

    # filter out gappy columns on column and group gap criteria
    group_too_gappy = 0
    for group in groups:
	if gap_percentage(group) > group_gap_cutoff:
	    group_too_gappy = 1
	    break

    if group_too_gappy:
	col_score = None
    elif gap_percentage(col) > column_gap_cutoff:
	col_score = None
    else:
	col_score = matrix_sum_pairs(groups, matrix)


    scores.append(col_score)

# map to [0,1]
if map_scores_to_01:
    scores = norm_01(scores)

# apply ConsWin heuristic
#if cons_win_len > 0:
    #scores = cons_window_score(scores, alignment, cons_win_len, lamb)


#print >> sys.stderr, \
   # "Scored %s with %s matrix, window length=%d, and lambda=%.2f.\n" \
   # % (align_file, matrix_file, cons_win_len, lamb)


# print scores to stdout or file
if outfile_name == '':
    score_file = sys.stdout
else:
    score_file = open(outfile_name, 'w')


score_file.write('# %s scored by GroupSim with %s matrix\n'
		 % (align_file, matrix_file))
score_file.write('# window params: len=%d, lambda=%.2f\n'
		 % (cons_win_len, lamb))
score_file.write('# group gap cutoff: %.2f  column gap cutoff: %.2f\n'
		 % (group_gap_cutoff, column_gap_cutoff))
score_file.write('# col_num\tscore\tcolumn\n')

for i, col_score in enumerate(scores):
    if col_score != None:
	score_file.write('%d\t%f\t%s\n' 
			 % (i, col_score, 
			    format_column(i, alignment, gids_to_seqs)))
    else:
	score_file.write('%d\tNone\t%s\n' % 
			 (i, format_column(i, alignment, gids_to_seqs)))

score_file.close()

if outfile_name != '':
    print >> sys.stderr, "Wrote %s.\n" % outfile_name
