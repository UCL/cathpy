# core
import logging
import os
import re
import subprocess

# local
from cathpy import seqio, datafiles, error as err

logger = logging.getLogger(__name__)

def is_valid_domain_id(id_str: str) -> bool:
    """
    Returns whether the given input is a valid CATH domain identifier.    
    """
    return re.match('([0-9][a-zA-Z0-9]{3})([a-zA-Z])([0-9]{2})$', id_str)

class CathVersion(object):
    """Object that represents a CATH version."""

    def __init__(self, *args, **kwargs):
        """Creates a new instance of a CathVersion object.
        
        The following are all equivalent:

            cv = CathVersion('4.2')
            cv = CathVersion('v4_2_0')
            cv = CathVersion(4.2)
            cv = CathVersion(4, 2, 0)

        """
        if len(args) == 1:
            ver_parts = __class__._split_string(args[0])
        elif len(args) == 2:
            ver_parts = (args[0], args[1], 0)
        elif len(args) == 3:
            ver_parts = args
        else:
            raise Exception("expected either one")
        
        self.major = ver_parts[0]
        self.minor = ver_parts[1]
        self.trace = ver_parts[2]

    def join(self, join_char="."):
        return join_char.join([str(self.major), str(self.minor), str(self.trace)])

    @property
    def dirname(self):
        """Return the version represented as a directory name (eg 'v4_2_0')."""
        return "v" + self.join("_")

    @classmethod
    def new_from_string(cls, version_str):
        version_parts = cls._split_string(version_str)
        return cls(*version_parts)
    
    @classmethod
    def _split_string(cls, version_str):
        version_str = re.sub(r'^v', '', str(version_str))
        parts = re.split(r'[._]', version_str)
        if len(parts) == 2:
            return (parts[0], parts[1], 0)
        elif len(parts) == 3:
            return (parts[0], parts[1], parts[2])

    def __str__(self):
        return self.join(".")

class FunfamFileFinder(object):
    """Finds a Funfam alignment file within a directory."""

    grep_path = '/bin/grep'

    def __init__(self, base_dir, *, ff_tmpl='__SFAM__-ff-__FF_NUM__.sto'):
        self.base_dir = base_dir
        if re.search(r'[^a-zA-Z_.\-*]', ff_tmpl):
            raise err.InvalidInputError("ff_tmpl contains unexpected characters")

        self.ff_tmpl = ff_tmpl

    def search_by_domain_id(self, domain_id):
        """Return the filename of the FunFam alignment containing the domain id."""
        if not is_valid_domain_id(domain_id):
            raise err.InvalidInputError('{} is not a valid domain id'.format(repr(domain_id)))

        # replace template placeholders with '*'
        glob_path = re.sub(r'__([A-Z_]+)__', '*', self.ff_tmpl)
        grep_args = (self.grep_path, '--include', glob_path, '-l', '-P', '\\b' + domain_id + '\\b', '-R', self.base_dir)
        logger.debug("search_by_domain_id: sys: " + " ".join(grep_args))

        try:
            # note: this returns bytes (not strings)
            grep_out = subprocess.check_output(grep_args).decode('ascii')
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
        
        logger.debug("search_by_domain_id: found funfam alignment {}".format(repr(ff_files[0])))

        return ff_files[0]

class StructuralClusterMerger(object):
    """
    Merges FunFams based on a structure-based alignment of representative sequences.
    
    Args:
        cath_version: version of CATH
        sc_file: structure-based alignment (`*.fa`) of funfam reps
        ff_dir: path of the funfam alignments (`*.sto`) to merge
        out_file: file to write merged alignment
        ff_tmpl: template used to find the funfam alignment files

    """

    def __init__(self, *, cath_version, sc_file, ff_dir, out_file=None,
        ff_tmpl="__SFAM__-ff-__FF_NUM__.reduced.sto"):

        if type(cath_version) is str:
            cath_version = CathVersion.new_from_string(cath_version)

        if out_file is None:
            out_file = re.sub(r'\.sto$', '.merged.sto', sc_file)

        self.cath_version = cath_version
        self.sc_file = sc_file
        self.out_file = out_file
        self.ff_dir = ff_dir
        self.ff_tmpl = ff_tmpl

    def run(self):

        logger.info("Running alignment merge")

        cath_release = datafiles.ReleaseDir(self.cath_version)

        logger.info("Parsing structure-based alignment: ")

        # parse the structure-based alignment of representatives
        # eg /cath/data/v4_2_0/funfam/families/1.10.8.10/1.10.8.10__FF_SSG9__6.reps.fa
        sc_aln = seqio.Alignment.new_from_fasta(self.sc_file)

        logger.info(" ... found {} representatives".format(sc_aln.count_sequences))

        merge_count=1
        def next_merge_stage_file():
            nonlocal merge_count
            out_file = self.out_file
            logger.info( "stage_file: merge_count={} out_file={}".format(merge_count, str(self.out_file)) )
            stage_file = re.sub(r'(\..*?)$', '.' + str(merge_count) + '\1', self.out_file)
            merge_count += 1
            return stage_file

        # create our funfam finder
        ff_finder = FunfamFileFinder(self.ff_dir, ff_tmpl=self.ff_tmpl)

        logger.info("Searching for funfam files in dir: " + self.ff_dir)

        # for each representative in the structure-based alignment..
        sc_aln_orig = sc_aln.copy()
        for sc_rep_in_sc in sc_aln_orig.seqs:

            logger.info('Working on SC rep: {}'.format(sc_rep_in_sc.id) )

            sc_rep_id = sc_rep_in_sc.id
            
            # find the corresponding funfam alignment
            ff_aln_file = ff_finder.search_by_domain_id(sc_rep_id)

            # parse it into an alignment
            ff_aln = seqio.Alignment.new_from_stockholm(ff_aln_file)

            # find the sc_rep sequence within the funfam alignment
            sc_rep_in_ff = ff_aln.find_seq_by_id(sc_rep_id)

            logger.debug('SC REP (SC): {}'.format(sc_rep_in_sc))
            logger.debug('SC REP (FF): {}'.format(sc_rep_in_ff))

            # get the chain correspondence file 
            rep_chain_id = sc_rep_id[:5]
            gcf_file = cath_release.get_file('chaingcf', rep_chain_id)
            chain_corr = seqio.Correspondence.new_from_gcf(gcf_file)
            
            # TODO: get a subset that only corresponds to the domain (not chain)
            seqres_segments = sc_rep_in_ff.segs
            logger.warning("TODO: this code currently assumes that the start-stop information " 
                "in the FunFam STOCKHOLM alignment matches the sequence and is based on SEQRES "
                "records (which needs to be double-checked)")

            if not seqres_segments:
                raise err.NoSegmentsError(('need to have seqres segments defined in '
                    'structural cluster rep sequence (of funfam): {}').format(sc_rep_in_ff))

            logger.info('applying segments to correspondence: {}'.format(repr(seqres_segments)))
            sc_rep_corr = chain_corr.apply_seqres_segments(seqres_segments)
            logger.info('  ...correspondence changed from {} (first:{}, last:{}) to {} (first:{}, last:{})'.format(
                chain_corr.seqres_length, str(chain_corr.first_residue), str(chain_corr.last_residue),
                sc_rep_corr.seqres_length, str(sc_rep_corr.first_residue), str(sc_rep_corr.last_residue), ))

            # merge the funfam into the sc alignment
            sc_aln.merge_alignment(ff_aln, sc_rep_id, sc_rep_corr)

            merge_stage_file = next_merge_stage_file()
            logger.info("Writing tmp merge file to '{}'".format(merge_stage_file))
            sc_aln.write_fasta(merge_stage_file, wrap_width=None)

        return sc_aln