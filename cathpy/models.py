"""
Collection of classes used to model CATH data
"""

import logging
import os
import re
import json
import jsonpickle

import cathpy.error as err

LOG = logging.getLogger(__name__)


class AminoAcid(object):
    """Class representing an Amino Acid."""

    def __init__(self, one, three, word):
        self.one = one
        self.three = three
        self.word = word

    def __str__(self):
        """Return this AminoAcid as a single character string"""
        return self.one


class AminoAcids(object):
    """Provides access to recognised Amino Acids."""

    aa_by_id = {aa[2]: AminoAcid(aa[2], aa[1], aa[0]) for aa in (
        ('alanine', 'ala', 'A'),
        ('arginine', 'arg', 'R'),
        ('asparagine', 'asn', 'N'),
        ('aspartic acid', 'asp', 'D'),
        ('asparagine', 'asx', 'B'),
        ('cysteine', 'cys', 'C'),
        ('glutamic acid', 'glu', 'E'),
        ('glutamine', 'gln', 'Q'),
        ('glutamic acid', 'glx', 'Z'),
        ('glycine', 'gly', 'G'),
        ('histidine', 'his', 'H'),
        ('isoleucine', 'ile', 'I'),
        ('leucine', 'leu', 'L'),
        ('lysine', 'lys', 'K'),
        ('methionine', 'met', 'M'),
        ('phenylalanine', 'phe', 'F'),
        ('proline', 'pro', 'P'),
        ('serine', 'ser', 'S'),
        ('threonine', 'thr', 'T'),
        ('tryptophan', 'trp', 'W'),
        ('tyrosine', 'tyr', 'Y'),
        ('valine', 'val', 'V'),
    )}

    def __init__(self):
        pass

    @classmethod
    def is_valid_aa(cls, aa_letter):
        """Check if aa is a valid single character aa code."""
        return str(aa_letter).upper() in cls.aa_by_id

    @classmethod
    def get_by_id(cls, aa_letter):
        """Return the AminoAcid object by the given single character aa code."""
        return cls.aa_by_id[aa_letter.upper()]


class Residue(object):
    """Class to represent a protein residue."""

    def __init__(self, aa, seq_num=None, pdb_label=None, *, pdb_aa=None):
        assert isinstance(aa, str)
        if seq_num:
            seq_num = int(seq_num)

        self.aa = aa
        self.seq_num = seq_num
        self.pdb_label = pdb_label
        self.pdb_aa = pdb_aa if pdb_aa else aa

    def __str__(self):
        return self.aa

    def __repr__(self):
        return "res({},seq:{},pdb:{},pdb_aa:{})".format(
            self.aa, self.seq_num, self.pdb_label, self.pdb_aa)


class CathID(object):
    """Represents a CATH ID."""

    RE_CATH_ID = re.compile(r'^[1-9]+(\.[0-9]+){0,8}$')

    def __init__(self, cath_id):
        assert self.RE_CATH_ID.match(cath_id)
        self._cath_id_parts = cath_id.split('.')

    @property
    def depth(self):
        """Returns the depth of the CATH ID."""
        return len(self._cath_id_parts)

    @property
    def sfam_id(self):
        """Returns the superfamily id of the CATH ID."""

        if self.depth < 4:
            raise err.OutOfBoundsError("cannot get sfam_id for CATH ID '{}' (require depth >= 4, not {})".format(
                ".".join(self._cath_id_parts), len(self._cath_id_parts)
            ))
        else:
            return self.cath_id_to_depth(4)

    @property
    def cath_id(self):
        """Returns the CATH ID as a string."""
        return str(self.cath_id)

    def cath_id_to_depth(self, depth):
        """Returns the CATH ID as a string."""
        return ".".join(self._cath_id_parts[:depth])

    def __str__(self):
        return ".".join(self._cath_id_parts)


class ClusterID(object):
    """Represents a Cluster Identifier (FunFam, SC, etc)"""

    def __init__(self, sfam_id, cluster_type, cluster_num):
        self.sfam_id = sfam_id
        self.cluster_type = cluster_type
        self.cluster_num = cluster_num

    @classmethod
    def new_from_file(cls, file):
        """Parse a new :class:`ClusterID` from a filename."""
        cf = ClusterFile(file)
        cls(cf.sfam_id, cf.cluster_type, cf.cluster_num)

    def __str__(self):
        return "{}-{}-{}".format(self.sfam_id, self.cluster_type, self.cluster_num)


class FunfamID(ClusterID):
    """Object that represents a FunFam ID."""

    def __init__(self, sfam_id, cluster_num):
        super().__init__(sfam_id, 'FF', cluster_num)


class ClusterFile(object):
    """
    Object that represents a file relating to a CATH Cluster.

    Example:
        /path/to/1.10.8.10-ff-1234.sto

    """

    def __init__(self, path, *, dir=None, sfam_id=None, cluster_type=None, cluster_num=None,
                 join_char=None, desc=None, suffix=None):

        # explicitly declare attributes (to keep pylint happy at the very least)
        self.dir = None
        self.sfam_id = None
        self.cluster_type = None
        self.cluster_num = None
        self.join_char = None
        self.desc = None
        self.suffix = None

        # initialise from path (if given)
        attrs = ('dir', 'sfam_id', 'cluster_type',
                 'cluster_num', 'join_char', 'desc', 'suffix')
        if path:
            path_info = __class__.split_path(path)
            for attr in attrs:
                setattr(self, attr, path_info[attr])

        # allow other parts to be specified (or overridden) by extra args
        local_args = locals()
        for attr in attrs:
            if local_args[attr]:
                setattr(self, attr, local_args[attr])

    @property
    def cluster_id(self):
        """
        Returns the cluster id as a :class:`ClusterID` object
        """
        ClusterID(self.sfam_id, self.cluster_type, self.cluster_num)

    @classmethod
    def split_path(cls, path):
        """
        Returns information about a cluster based on the path (filename).
        """

        re_file = re.compile(r'(?P<sfam_id>[0-9.]+)'         # 1.10.8.10
                             # - or __
                             r'(?P<join_char>-|__)'
                             # ff
                             r'(?P<cluster_type>\w+)'
                             # - or __
                             r'(?:-|__)'
                             # 1234
                             r'(?P<cluster_num>[0-9]+)'
                             # .reduced
                             r'(?P<desc>.*?)'
                             r'(?P<suffix>\.\w+)$')                           # .sto

        m = re_file.match(os.path.basename(path))
        if m:
            info = m.groupdict()
            info['dir'] = os.path.dirname(path)
            return info

        raise err.NoMatchesError(
            'failed to parse cluster details from filename {}'.format(path))

    def to_string(self, join_char=None):
        """Represents the ClusterFile as a string (file path)."""
        if not join_char:
            join_char = self.join_char

        file_parts = [
            join_char.join(
                [self.sfam_id, self.cluster_type, self.cluster_num]),
            self.desc,
            self.suffix]

        # remove elements that are not defined (eg. desc, suffix)
        file_parts = [p for p in file_parts if p != None]

        filename = ''.join(file_parts)

        path = os.path.join(self.dir, filename) if self.dir else filename

        return path

    def __str__(self):
        return self.to_string()


class Segment(object):
    """Class to represent a protein segment."""

    def __init__(self, start: int, stop: int):
        self.start = int(start)
        self.stop = int(stop)

    def __str__(self):
        return "{}-{}".format(self.start, self.stop)

    def __repr__(self):
        return "Segment:{}-{}".format(self.start, self.stop)

    def __getitem__(self, idx):
        items = [self.start, self.stop]
        return items[idx]


class ScanHsp(object):
    """Object to store the High Scoring Pair (HSP) from a sequence scan."""

    def __init__(self, *, evalue, hit_start, hit_end, hit_string=None,
                 homology_string=None, length, query_start, query_end,
                 query_string=None, rank, score, **kwargs):
        """Create a new :class:`ScanHsp` object.

        Stores information about the High Scoring Pair
        within a given :class:`ScanHit`.

        Args:
            hit_start (int): position of the hsp start for the hit
            hit_end (int): position of the hsp end for the hit
            hit_string (str): hit residues in the matching hsp
            query_start (int): position of the hsp start for query
            query_end (int): position of the hsp end for query
            query_string
            evalue (float): evalue of the hsp
            homology_string (str): query  
            length (int): length of the hsp
            rank (int): the rank of this hsp
            score (float): the bit score 

        """

        self.evalue = evalue
        self.hit_start = hit_start
        self.hit_end = hit_end
        self.hit_string = hit_string
        self.homology_string = homology_string
        self.length = length
        self.query_start = query_start
        self.query_end = query_end
        self.query_string = query_string
        self.rank = rank
        self.score = score


class ScanHit(object):
    """Object to store a hit from a sequence scan.

    Each :class:`ScanResult` object (ie **query**),
    contains one or more :class:`ScanHit` objects (ie **match**)

    """

    def __init__(self, *, match_name, match_cath_id, match_description,
                 match_length, hsps, significance, data, **kwargs):
        """Create a new :class:`ScanHit` object.

        Args:
            match_name (str): name of the match
            match_cath_id (str): CATH superfamily id of the match
            match_description (str): description of the match
            match_length (str): number of residues in the match protein
            hsps (dict|ScanHsp): array of :class:`ScanHsp` objects 
            significance (float): significance score

        """
        self.match_name = match_name
        self.match_cath_id = match_cath_id
        self.match_description = match_description
        self.match_length = match_length
        self.hsps = [ScanHsp(**hsp) for hsp in hsps]
        self.data = data
        self.significance = significance


class ScanResult(object):
    """Object to store a result from a sequence :class:`Scan`."""

    def __init__(self, *, query_name, hits, **kwargs):
        self.query_name = query_name
        self.hits = [ScanHit(**hit) for hit in hits]


class Scan(object):
    """Object to store a sequence scan.

    A scan can contain one or more :class:`ScanResult` objects
    (eg **query**).

    """

    def __init__(self, *, results, **kwargs):
        self.results = [ScanResult(**res) for res in results]

    def as_json(self, *, pp=False):
        """Returns the Scan as JSON formatted string."""

        data = jsonpickle.encode(self)
        if pp:
            data = json.dumps(json.loads(data), indent=2, sort_keys=True)

        LOG.info("Serialized Scan as JSON string (length:%s)", len(data))
        return data

    def as_tsv(self, *, header=True):
        """Returns the Scan as CSV"""

        result = self.results[0]

        lines = []
        if header:
            headers = 'funfam members uniq_ec_terms query_region match_region evalue score description'.split()
            lines.append("\t".join(headers))

        for hit in result.hits:
            ec_term_count = hit.data['ec_term_count'] or 0
            for hsp in hit.hsps:
                line = '\t'.join([
                    '{}'.format(hit.match_name),
                    '{}'.format(hit.data['funfam_members']),
                    '{}'.format(ec_term_count),
                    '{}-{}'.format(hsp.query_start, hsp.query_end),
                    '{}-{}'.format(hsp.hit_start, hsp.hit_end),
                    '{:.1e}'.format(hsp.evalue),
                    '{:d}'.format(int(hsp.score)),
                    '"{}"'.format(hit.match_description),
                ])
                lines.append(line)
        return "".join([l + "\n" for l in lines])
