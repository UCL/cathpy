"""
General utility classes and functions
"""

# core
import logging
import io
import os
import platform
import re
import subprocess
from subprocess import Popen, PIPE, CalledProcessError
import tempfile

import pkg_resources

# local
from cathpy import error as err
from cathpy.tests import is_valid_domain_id
from cathpy.models import CathID, ClusterID, ClusterFile, FunfamID
from cathpy.align import Align, Correspondence
from cathpy.datafiles import ReleaseDir
from cathpy.version import CathVersion
from cathpy.error import ParseError, OutOfBoundsError

LOG = logging.getLogger(__name__)

TOOL_DIR = pkg_resources.resource_filename(__name__, "tools")
DATA_DIR = os.path.join(TOOL_DIR, "data")
PLATFORM_DIRNAME = "{}-{}".format(platform.system().lower(), platform.machine())
GROUPSIM_DIR = os.path.join(TOOL_DIR, 'GroupSim')
SCORECONS_EXE = os.path.join(TOOL_DIR, PLATFORM_DIRNAME, 'scorecons')
SCORECONS_MATRIX_FILE = os.path.join(DATA_DIR, 'PET91mod.mat2')

class GroupsimResult(object):
    """Represents the result from running the groupsim algorithm."""

    def __init__(self, *, scores = None):
        self.scores = scores

    @classmethod
    def new_from_file(cls, gs_file):
        """Create a new groupsim result from an output file."""
        with open(gs_file) as f:
            obj = cls.new_from_io(f)
            return obj

    @classmethod
    def new_from_io(cls, gs_io, *, maxscore=1):
        """Create a new groupsim result from an io source."""
        gs_scores = []
        for line in gs_io:
            line = line.strip()
            if line.startswith('#'):
                continue
            # aln_id, gs_score, col_data
            _, gs_score, _ = line.split(None, maxsplit=2)
            if gs_score == 'None':
                gs_score = None
            else:
                gs_score = float(gs_score) / float(maxscore)
            gs_scores.append( gs_score )
        return cls(scores=gs_scores)

    @property
    def count_positions(self):
        """Returns the number of positions in the groupsim result."""
        return len(self.scores)

    @property
    def to_string(self):
        # normalise 0-1 to 0-9
        return "".join([str(int(s*9)) if s else '-' for s in self.scores])


class GroupsimRunner(object):
    """Object that provides a wrapper around groupsim."""

    def __init__(self, *, groupsim_dir=GROUPSIM_DIR, python2path='python2', 
        column_gap=0.3, group_gap=0.5):

        self.groupsim_dir = groupsim_dir
        self.groupsim_path = os.path.join(self.groupsim_dir, 'group_sim_sdp.py')
        self.mclachlan_path = os.path.join(self.groupsim_dir, 'mclachlan1972.aa.dat')
        self.column_gap = column_gap
        self.group_gap = group_gap
        self.python2path = python2path

    def run_alignment(self, alignment, *, column_gap=None, group_gap=None, mclachlan=False):
        """Runs groupsim against a given alignment."""

        # mclachan max score is 5: normalise to 0-1 before storing  
        maxscore = 5 if mclachlan else 1

        fasta_tmp = tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".fa")
        fasta_tmp_filename = fasta_tmp.name

        if not column_gap:
            column_gap = self.column_gap 

        if not group_gap:
            group_gap = self.group_gap 

        column_gap = float(column_gap)
        group_gap = float(group_gap)

        assert(column_gap > 0 and column_gap < 1)
        assert(group_gap > 0 and group_gap < 1)

        # write out the alignment with funfam numbers appended to sequence ids
        # >1ebgB02/127-436|7431
        aln_copy = alignment.copy()
        for seq in aln_copy.seqs:
            if not seq.cluster_id:
                raise err.InvalidInputError("need to set_cluster_id() on alignment sequences before running groupsim: {}".format(
                    seq.__dict__
                ))
            seq._id = '|'.join([seq._id, str(seq.cluster_id)])

        source_ids = {s.cluster_id for s in aln_copy.seqs}

        # lower-case aa -> gaps
        # '.' -> '-'
        for s in aln_copy.seqs:
            s.set_all_gap_chars(gap_char='-')
            s.set_lower_case_to_gap(gap_char='-') 

        aln_copy.write_fasta(fasta_tmp_filename)

        groupsim_args = [self.python2path, self.groupsim_path, 
            '-c', str(self.column_gap), '-g', str(self.group_gap)]

        if mclachlan:
            groupsim_args.extend(['-m', self.mclachlan_path])

        groupsim_args.append(fasta_tmp_filename)
        groupsim_args.extend(source_ids)

        groupsim_args = [str(a) for a in groupsim_args]

        LOG.debug("running groupsim: sys: " + " ".join(groupsim_args))

        try:
            p = Popen(groupsim_args, stdout=PIPE, stderr=PIPE, universal_newlines=True)
            groupsim_out, groupsim_err = p.communicate()
        except CalledProcessError as e:
            LOG.error('CMD: {}\nCODE: {}\nOUTPUT: {}\nSTDERR: "{}"\nSTDOUT: "{}"\n'.format(
                e.cmd, e.returncode, e.output, e.stderr, e.stdout))
            raise e
        except:
            raise err.FileNotFoundError("Encountered unexpected error running GroupSim: `{}`".format(
                " ".join(groupsim_args)
            ))

        gs_io = io.StringIO(groupsim_out)

        res = GroupsimResult.new_from_io(gs_io, maxscore=maxscore)

        return res

class ScoreconsResult(object):
    """Represents the results from running scorecons."""

    def __init__(self, *, dops, scores):
        self.dops = dops
        self.scores = scores
    
    @property
    def to_string(self):
        """Returns the scorecons results as a string (one char per position)."""
        # normalise 0-1 to 0-9
        return "".join([str(int(s*9)) if s else '-' for s in self.scores])


class ScoreconsRunner(object):
    """Runs scorecons for a given alignment."""

    def __init__(self, *, scorecons_path=SCORECONS_EXE, matrix_path=SCORECONS_MATRIX_FILE):
        self.scorecons_path = scorecons_path
        self.matrix_path = matrix_path
        LOG.debug( "scorecons.init: bin:{} matrix:{}".format(
            self.scorecons_path, self.matrix_path))

    def run_alignment(self, alignment):
        """Runs scorecons on a given alignment."""

        fasta_tmp = tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".fa")
        fasta_tmp_filename = fasta_tmp.name

        aln_copy = alignment.copy()

        # lower-case aa -> gaps
        # '.' -> '-'
        for s in aln_copy.seqs:
            s.set_all_gap_chars(gap_char='-')
            s.set_lower_case_to_gap(gap_char='-') 

        LOG.info('Writing normalised alignment to %s', fasta_tmp_filename)
        aln_copy.write_fasta(fasta_tmp_filename)

        return self.run_fasta(fasta_tmp_filename)

    def run_stockholm(self, sto_file):
        """
        Returns scorecons data for the provided STOCKHOLM file.
        
        Returns:
            result (ScoreconsResult): scorecons result

        """

        fasta_tmp = tempfile.NamedTemporaryFile(mode='w+', delete=True, suffix=".fa")
        fasta_tmp_filename = fasta_tmp.name
        aln = Align.new_from_stockholm(sto_file)
        aln.write_fasta(fasta_tmp_filename)
        return self.run_fasta(fasta_tmp_filename)

    def run_fasta(self, fasta_file):
        """
        Returns scorecons data (ScoreconsResult) for the provided FASTA file.
        
        Returns:
            result (ScoreconsResult): scorecons result
        """

        tmp_scorecons_file = tempfile.NamedTemporaryFile(mode='w+', suffix='.scorecons', delete=True)

        scorecons_args = (self.scorecons_path, 
            '-a', fasta_file, 
            '-m', self.matrix_path, 
            '-o', tmp_scorecons_file.name)

        LOG.debug("running scorecons: sys: " + " ".join(scorecons_args))

        try:
            p = Popen(scorecons_args, stdout=PIPE, stderr=PIPE, universal_newlines=True)
            scorecons_out, scorecons_err = p.communicate()

        except subprocess.CalledProcessError as e:
            LOG.error('CMD: {}\nCODE: {}\nOUTPUT: {}\nSTDERR: "{}"\nSTDOUT: "{}"\n'.format(
                e.cmd, e.returncode, e.output, e.stderr, e.stdout))
            raise e
        except:
            raise err.FileNotFoundError("Encountered error running scorecons: `{}`".format(
                " ".join(scorecons_args)
            ))

        match = re.search(r'^DOPS score:\s+([0-9.]+)', scorecons_out, flags=re.MULTILINE)
        if not match:
            raise err.ParseError('Failed to find DOPS score in scorecons output: {} (STDERR: {})'.format(
                scorecons_out, scorecons_err))
        
        dops_score = float(match.group(1))

        sc_numbers = __class__.split_scorecons_file(tmp_scorecons_file.name)

        res = ScoreconsResult(dops = dops_score, scores = sc_numbers)

        return res

    @classmethod
    def split_scorecons_file(cls, scorecons_file):

        sc_numbers=[]
        with open(scorecons_file, 'r') as f:
            for line in f:
                sc_aln_num, _ = line.split(None, maxsplit=1)
                sc_numbers.append(float(sc_aln_num))

        return sc_numbers

class FunfamFileFinder(object):
    """Finds a Funfam alignment file within a directory."""

    grep_path = 'grep'

    def __init__(self, base_dir, *, ff_tmpl='__SFAM__-ff-__FF_NUM__.sto'):
        self.base_dir = base_dir
        if re.search(r'[^a-zA-Z_.\-*]', ff_tmpl):
            raise err.InvalidInputError("ff_tmpl contains unexpected characters")

        self.ff_tmpl = ff_tmpl

    def funfam_id_from_file(self, ff_file):
        """Extracts a FunfamID from the file name (based on the ff_tmpl)"""
        ff_re = self.ff_tmpl
        ff_re = ff_re.replace('__SFAM__', r'(?P<sfam_id>[0-9\.]+)')
        ff_re = ff_re.replace('__FF_NUM__', r'(?P<ff_num>[0-9]+)')
        filename = os.path.basename(ff_file)
        m = re.match(ff_re, filename)
        if not m:
            raise err.NoMatchesError("failed to match template '{}' against filename '{}'".format(
                ff_re, filename
            ))
        ff_id = FunfamID(m.group('sfam_id'), int(m.group('ff_num')))

        return ff_id

    def search_by_domain_id(self, domain_id):
        """Return the filename of the FunFam alignment containing the domain id."""
        if not is_valid_domain_id(domain_id):
            raise err.InvalidInputError('{} is not a valid domain id'.format(repr(domain_id)))

        # replace template placeholders with '*'
        glob_path = re.sub(r'__([A-Z_]+)__', '*', self.ff_tmpl)
        grep_args = (self.grep_path, '--include', glob_path, '-l', 
            '^' + domain_id, '-R', self.base_dir)
        LOG.debug("search_by_domain_id: sys: " + " ".join(grep_args))

        try:
            # note: this returns bytes (not strings)
            grep_out = subprocess.check_output(grep_args).decode('ascii')
        except subprocess.CalledProcessError as e:
            if e.returncode == 1:
                # grep telling us it didn't find any matches
                raise err.NoMatchesError('failed to find domain id {} with cmd {}'.format(domain_id, str(e.cmd)) )
            else:
                LOG.error('CMD: {}\nCODE: {}\nOUTPUT: {}\nSTDERR: "{}"\nSTDOUT: "{}"\n'.format(
                    e.cmd, e.returncode, e.output, e.stderr, e.stdout))
                raise
        except:
            raise err.FileNotFoundError("Encountered error trying to find domain_id '{}' (grep: `{}`)".format(
                domain_id, " ".join(grep_args)
            ))

        ff_files = grep_out.splitlines()

        if len(ff_files) == 0:
            raise err.FileNotFoundError("Failed to find FunFam alignment for domain_id '{}' (grep: `{}`)".format(
                domain_id, " ".join(grep_args)))
        elif len(ff_files) > 1:
            raise err.GeneralError("Found more than one FunFam file ({}) containing the domain id '{}' (grep: `{}`):\n{}\n".format(
                    len(ff_files), domain_id, " ".join(grep_args), "\n".join(ff_files),
                ))
        
        LOG.debug("search_by_domain_id: found funfam alignment {}".format(repr(ff_files[0])))

        return ff_files[0]

class StructuralClusterMerger(object):
    """
    Merges FunFams based on a structure-based alignment of representative sequences.
    
    Args:
        cath_version: version of CATH
        sc_file: structure-based alignment (`*.fa`) of funfam reps
        ff_dir: path of the funfam alignments (`*.sto`) to merge
        out_fasta: file to write merged alignment (FASTA)
        out_sto: file to write merged alignment (STOCKHOLM)
        ff_tmpl: template used to find the funfam alignment files
        add_groupsim: add groupsim data (default: True)
        add_scorecons: add scorecons data (default: True)
        cath_release: specify custom release data directory
    """

    def __init__(self, *, cath_version, sc_file, ff_dir, out_fasta=None, 
        out_sto=None, ff_tmpl="__SFAM__-ff-__FF_NUM__.sto", 
        add_groupsim=True, add_scorecons=True, cath_release=None):

        if type(cath_version) is str:
            cath_version = CathVersion.new_from_string(cath_version)

        self.cath_version = cath_version
        self.sc_file = sc_file
        self.out_fasta = out_fasta
        self.out_sto = out_sto
        self.ff_dir = ff_dir
        self.ff_tmpl = ff_tmpl
        self.wrap_width = None
        self.add_groupsim = add_groupsim
        self.add_scorecons = add_scorecons

        if not cath_release:
            cath_release = datafiles.ReleaseDir(self.cath_version)

        self.cath_release = cath_release

    def run(self):

        LOG.info("Running alignment merge")

        cath_release = self.cath_release

        # parse the structure-based alignment of representatives
        # eg /cath/data/v4_2_0/funfam/families/1.10.8.10/1.10.8.10__FF_SSG9__6.reps.fa
        sc_filename = os.path.basename(self.sc_file)
        sc_parts = re.match(r'(\d+\.\d+\.\d+\.\d+)__([A-Z0-9_]+)__(\d+)\b', sc_filename)
        
        if not sc_parts:
            raise Exception('failed to parse necessary meta info from sc_file name: ' + sc_filename)

        sfam_id, cluster_type, sc_num = sc_parts.group(1, 2, 3)

        LOG.info('Superfamily: ' + sfam_id)
        LOG.info('Cluster type: ' + cluster_type)
        LOG.info('Cluster number: ' + sc_num)

        LOG.info("Parsing structure-based alignment: ")
        sc_aln = Align.new_from_fasta(self.sc_file)
        LOG.info(" ... found {} representatives".format(sc_aln.count_sequences))

        cluster_id = '-'.join([sfam_id, cluster_type, sc_num])
        sc_aln.set_id(cluster_id)
        sc_aln.accession = cluster_id
        sc_aln.aln_type = cluster_type
        sc_aln.description = '{}, Structural Cluster ({}) {}'.format(sfam_id, cluster_type, sc_num)

        merge_count=1
        def next_merge_stage_file():
            nonlocal merge_count
            out_fasta = str(self.out_fasta)
            stage_file = re.sub(r'(\..*?)$', '.' + str(merge_count) + '\1', out_fasta)
            LOG.debug( "stage_file: merge_count={} out_fasta={} stage_file={}".format(
                merge_count, out_fasta, stage_file) )
            merge_count += 1
            return stage_file

        # create our funfam finder
        ff_finder = FunfamFileFinder(self.ff_dir, ff_tmpl=self.ff_tmpl)

        LOG.info("Searching for funfam files in dir: " + self.ff_dir)

        # for each representative in the structure-based alignment..
        sc_aln_orig = sc_aln.copy()
        for sc_rep_in_sc in sc_aln_orig.seqs:

            LOG.info('Working on SC rep: {}'.format(sc_rep_in_sc.accession) )

            sc_rep_acc = sc_rep_in_sc.accession
            
            # find the corresponding funfam alignment
            ff_aln_file = ff_finder.search_by_domain_id(sc_rep_acc)

            LOG.info('Reading FunFam alignment: {}'.format(ff_aln_file))

            # parse it into an alignment
            ff_aln = Align.new_from_stockholm(ff_aln_file)

            # we need the funfam_number for groupsim
            funfam_id = ff_finder.funfam_id_from_file(ff_aln_file)

            # find the sc_rep sequence within the funfam alignment
            sc_rep_in_ff = ff_aln.find_seq_by_accession(sc_rep_acc)
            if not sc_rep_in_ff:
                raise err.GeneralError('failed to find structural cluster representative {} in funfam {}'.format(
                    sc_rep_acc, ff_aln_file, ))

            LOG.debug('SC REP (SC): {}'.format(sc_rep_in_sc))
            LOG.debug('SC REP (FF): {}'.format(sc_rep_in_ff))

            # get the chain correspondence file 
            rep_chain_id = sc_rep_acc[:5]
            gcf_file = cath_release.get_file('chaingcf', rep_chain_id)

            chain_corr = Correspondence.new_from_gcf(gcf_file)
            
            # TODO: get a subset that only corresponds to the domain (not chain)
            seqres_segments = sc_rep_in_ff.segs
            LOG.warning("TODO: this code currently assumes that the start-stop information "
                "in the FunFam STOCKHOLM alignment matches the sequence and is based on SEQRES "
                "records (which needs to be double-checked)")

            if not seqres_segments:
                raise err.MissingSegmentsError(('need to have seqres segments defined in '
                    'structural cluster rep sequence (of funfam): {}').format(sc_rep_in_ff))

            LOG.info('applying segments to correspondence: {}'.format(repr(seqres_segments)))
            sc_rep_corr = chain_corr.apply_seqres_segments(seqres_segments)
            LOG.info('  ...correspondence changed from {} (first:{}, last:{}) to {} (first:{}, last:{})'.format(
                chain_corr.seqres_length, str(chain_corr.first_residue), str(chain_corr.last_residue),
                sc_rep_corr.seqres_length, str(sc_rep_corr.first_residue), str(sc_rep_corr.last_residue), ))

            # merge the funfam into the sc alignment
            sc_aln.merge_alignment(ff_aln, sc_rep_acc, sc_rep_corr, 
                cluster_label = funfam_id.cluster_num)

            merge_stage_file = next_merge_stage_file()
            #LOG.info("Writing tmp merge file to '{}'".format(merge_stage_file))
            #sc_aln.write_fasta(merge_stage_file, wrap_width=None)

        # add scorecons
        if self.add_scorecons:
            sc_aln.add_scorecons()

        # add groupsim
        if self.add_groupsim:
            sc_aln.add_groupsim()

        # write final merged alignment
        if self.out_fasta:
            LOG.info('Writing merged FASTA alignment: {}'.format(self.out_fasta))
            sc_aln.write_fasta(self.out_fasta, self.wrap_width)

        if self.out_sto:
            LOG.info('Writing merged STOCKHOLM alignment: {}'.format(self.out_sto))
            sc_aln.write_sto(self.out_sto)

        return sc_aln

class AlignmentSummary(object):
    """Stores summary information about an alignment."""

    def __init__(self, *, path, dops, aln_length, seq_count, gap_count, total_positions):
        self.path = path
        self.dops = float(dops) if dops is not None else None
        self.aln_length = int(aln_length)
        self.seq_count = int(seq_count)
        self.gap_count = int(gap_count)
        self.total_positions = int(total_positions)

    def __str__(self):
        return "{:<70s} {:>6d} {:>6d} {:>6s} {:>6.2f}".format(self.path,
            self.aln_length, self.seq_count, 
            "{:.2f}".format(self.dops) if self.dops is not None else 'NULL', 
            100 * self.gap_count / self.total_positions)

    @classmethod
    def headers(cls):
        return ['file', 'aln_len', 'seq_count', 'dops', 'gap_per']

    def to_string(self):
        return self.__str__()

class AlignmentSummaryRunner(object):
    """
    Provides a summary report for sequence alignment files.
    
    Args:
        aln_dir: input alignment directory
        aln_file: input alignment file 
        suffix: filter alignments by suffix
        skipempty: skip empty files    
    """

    def __init__(self, *, aln_dir=None, aln_file=None, suffix='.sto', skipempty=False):
        self.suffix = suffix
        self.skipempty = skipempty
        if aln_file:
            self._files = [aln_file]
        elif aln_dir:
            self._files = self._get_files(aln_dir, suffix=suffix)
        else:
            raise Exception("error: expected 'aln_file' or 'aln_dir'")
        
    def _get_files(self, input_dir, *, suffix):

        all_aln_files = []

        for dirpath, dirnames, filenames in os.walk(input_dir):
            aln_files = [os.path.join(dirpath, f) for f in filenames 
                if os.path.basename(f).endswith(self.suffix)]

            nonempty_files = []
            for f in aln_files:
                if os.stat(f).st_size == 0:
                    if self.skipempty:
                        LOG.warning('skipping empty file: {}'.format(f))
                    else:
                        raise err.FileEmptyError('file {} is empty (use skipempty to skip)'.format(f))
                else:
                    nonempty_files.append(f)
            aln_files = nonempty_files

            all_aln_files.extend(aln_files)
        
        return all_aln_files

    def run(self):

        print('{}'.format(" ".join(AlignmentSummary.headers())))

        for aln_file in self._files:
            try:
                aln = Align.new_from_stockholm(aln_file)
            except:
                raise Exception("Failed to parse STOCKHOLM alignment from file {}".format(aln_file))

            try:
                aln_sum = AlignmentSummary(path=aln_file, dops=aln.dops_score, aln_length=aln.aln_positions, 
                    seq_count=aln.count_sequences, gap_count=aln.total_gap_positions, 
                    total_positions=aln.total_positions)
            except:
                raise Exception("Failed to generate alignment summary from file: {}".format(aln_file))

            print(aln_sum.to_string())
