"""
Collection of classes used to model CATH data
"""

import logging
import os
import re
import json
import jsonpickle

import cathpy.core.error as err

LOG = logging.getLogger(__name__)


class AminoAcid(object):
    """Class representing an Amino Acid."""

    def __init__(self, one, three, word):
        self.one = str(one).upper()
        self.three = str(three).upper()
        self.word = word

    def __str__(self):
        """Return this AminoAcid as a single character string"""
        return self.one


class AminoAcids(object):
    """Provides access to recognised Amino Acids."""

    amino_acids = [AminoAcid(aa[2], aa[1], aa[0]) for aa in (
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
    )]

    _aa_by_one = {aa.one: aa for aa in amino_acids}
    _aa_by_three = {aa.three: aa for aa in amino_acids}

    def __init__(self):
        pass

    @classmethod
    def is_valid_aa(cls, aa_str):
        """Check if aa is a valid 1 or 3-letter AA code."""
        try:
            cls.get_by_id(aa_str)
            return True
        except Exception:
            return False

    @classmethod
    def get_by_id(cls, aa_str):
        """Return the AminoAcid object by the given single character aa code."""

        aa_str = str(aa_str)
        aa_obj = None
        if len(aa_str) == 1:
            aa_obj = cls._aa_by_one[aa_str.upper()]
        elif len(aa_str) == 3:
            aa_obj = cls._aa_by_three[aa_str.upper()]
        else:
            raise err.InvalidInputError(
                "expected either 1- or 3-character amino acid id (not: '{}')".format(aa_str))
        return aa_obj


class Residue(object):
    """Class to represent a protein residue."""

    def __init__(self, aa, seq_num=None, pdb_label=None, *, pdb_aa=None):
        assert isinstance(aa, str)
        if seq_num:
            seq_num = int(seq_num)

        self.pdb_residue_num = None
        self.pdb_insert_code = None
        if pdb_label:
            self.set_pdb_label(pdb_label)

        self.aa = aa
        self.seq_num = seq_num
        self.pdb_aa = pdb_aa if pdb_aa else aa

    def set_pdb_label(self, pdb_label):
        if pdb_label[-1].isalpha():
            self.pdb_residue_num = pdb_label[:-1]
            self.pdb_insert_code = pdb_label[-1]
        else:
            self.pdb_residue_num = pdb_label
            self.pdb_insert_code = None

    @property
    def pdb_label(self):
        if not self.pdb_residue_num:
            return None
        return '{}{}'.format(
            str(self.pdb_residue_num),
            self.pdb_insert_code if self.pdb_insert_code else '')

    @property
    def has_pdb_insert_code(self):
        return bool(self.pdb_insert_code)

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
        self._cath_id_parts = [int(n) for n in cath_id.split('.')]

    @property
    def depth(self):
        """Returns the depth of the CATH ID."""
        return len(self._cath_id_parts)

    @property
    def sfam_id(self):
        """Returns the superfamily id of the CATH ID."""

        if self.depth < 4:
            raise err.OutOfBoundsError("cannot get sfam_id for CATH ID '{}' (require depth >= 4, not {})".format(
                self.cath_id, len(self._cath_id_parts)
            ))
        else:
            return self.cath_id_to_depth(4)

    @property
    def cath_id(self):
        """Returns the CATH ID as a string."""
        return str(self)

    def cath_id_to_depth(self, depth):
        """Returns the CATH ID as a string."""
        return ".".join([str(c) for c in self._cath_id_parts[:depth]])

    def __lt__(self, other):
        """
        Compares this CATH ID against another one
        """
        if not isinstance(other, CathID):
            other = CathID(other)
        return self._cath_id_parts < other._cath_id_parts

    def __eq__(self, other):
        """
        Checks if this CATH ID is equal to another CATH ID (string or object)
        """
        if not isinstance(other, CathID):
            other = CathID(other)
        return self._cath_id_parts == other._cath_id_parts

    def __str__(self):
        return ".".join([str(c) for c in self._cath_id_parts])


class ClusterID(object):
    """
    Represents a Cluster Identifier (FunFam, SC, etc)

    Usage:

    ::

        # equivalent
        cluster_id = ClusterID('1.10.8.10-ff-1234')
        cluster_id = ClusterID.from_string('1.10.8.10-ff-1234')
        cluster_id = ClusterID(sfam_id='1.10.8.10', cluster_type='ff', cluster_num=1234)

        cluster_id.sfam_id        # '1.10.8.10'
        cluster_id.cath_id        # CathID('1.10.8.10')
        cluster_id.cluster_type   # 'ff'
        cluster_id.cluster_num    # 1234

        cluster_id.to_string()    # '1.10.8.10-ff-1234'

    """

    def __init__(self, *args, **kwargs):

        sfam_id = None
        cluster_type = None
        cluster_num = None

        if len(args) == 1:
            cluster_id = ClusterID.from_string(args[0])
            sfam_id = cluster_id.sfam_id
            cluster_type = cluster_id.cluster_type
            cluster_num = cluster_id.cluster_num

        if not sfam_id:
            sfam_id = kwargs['sfam_id']
        if not cluster_type:
            cluster_type = kwargs['cluster_type']
        if not cluster_num:
            cluster_num = kwargs['cluster_num']

        self.cath_id = CathID(sfam_id)
        self.cluster_type = cluster_type
        self.cluster_num = int(cluster_num)

    @property
    def sfam_id(self):
        return self.cath_id.sfam_id

    @classmethod
    def from_string(cls, file):
        """Parse a new :class:`ClusterID` from a filename."""
        cf = ClusterFile(file)
        return cls(sfam_id=cf.sfam_id, cluster_type=cf.cluster_type, cluster_num=cf.cluster_num)

    def to_string(self):
        """
        Returns this object as a string (ie filepath).
        """
        return "{}-{}-{}".format(self.sfam_id, self.cluster_type, self.cluster_num)

    def __str__(self):
        return self.to_string()


class FunfamID(ClusterID):
    """Object that represents a FunFam ID."""

    def __init__(self, *, sfam_id, cluster_num):
        super().__init__(sfam_id=sfam_id, cluster_type='FF', cluster_num=cluster_num)


class ClusterFile(object):
    """
    Object that represents a file relating to a CATH Cluster.

    ::

        cf = ClusterFile('/path/to/1.10.8.10-ff-1234.reduced.sto')
        cf.path          # '/path/to'
        cf.sfam_id       # '1.10.8.10'
        cf.cluster_type  # 'ff'
        cf.cluster_num   # '1234'
        cf.join_char     # '-'
        cf.desc          # '.reduced'
        cf.suffix        # '.sto'

    """

    def __init__(self, fullpath, *, path=None, sfam_id=None, cluster_type=None, cluster_num=None,
                 join_char=None, desc=None, suffix=None):

        # explicitly declare attributes (to keep pylint happy at the very least)
        self.path = None
        self.sfam_id = None
        self.cluster_type = None
        self.cluster_num = None
        self.join_char = None
        self.desc = None
        self.suffix = None

        # initialise from path (if given)
        attrs = ('path', 'sfam_id', 'cluster_type',
                 'cluster_num', 'join_char', 'desc', 'suffix')
        if fullpath:
            path_info = __class__.split_path(fullpath)
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
        ClusterID(sfam_id=self.sfam_id, cluster_type=self.cluster_type,
                  cluster_num=self.cluster_num)

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
                             # .reduced (optional)
                             r'(?P<desc>.*?)'
                             # .sto (optional)
                             r'(?P<suffix>\.\w+)?$')

        basename = os.path.basename(path)
        m = re_file.match(basename)
        if m:
            info = m.groupdict()
            info['path'] = os.path.dirname(path)
            return info

        raise err.NoMatchesError(
            'failed to parse cluster details from filename "{}"'.format(basename))

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
        file_parts = [p for p in file_parts if p is not None]

        filename = ''.join(file_parts)

        path = os.path.join(self.path, filename) if self.path else filename

        return path

    def __str__(self):
        return self.to_string()


class Segment(object):
    """
    Class to represent a protein segment.

    Args:
        start (int): numeric start position of the segment
        stop (int): numeric stop position of the segment

    """

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
    """
    Object to store the High Scoring Pair (HSP) from a sequence scan.

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

    def __init__(self, *, evalue, hit_start, hit_end, hit_string=None,
                 homology_string=None, length, query_start, query_end,
                 query_string=None, rank, score, **kwargs):
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

    Args:
        match_name (str): name of the match
        match_cath_id (str): CATH superfamily id of the match
        match_description (str): description of the match
        match_length (str): number of residues in the match protein
        hsps (dict|ScanHsp): array of :class:`ScanHsp` objects
        significance (float): significance score
        data (dict): additional meta data about the hit (optional)

    """

    def __init__(self, *, match_name, match_cath_id, match_description,
                 match_length, hsps, significance, data=None, **kwargs):
        self.match_name = match_name
        self.match_cath_id = match_cath_id
        self.match_description = match_description
        self.match_length = match_length
        self.hsps = [ScanHsp(**hsp) for hsp in hsps]
        if not data:
            data = {}
        self.data = data
        self.significance = significance


class ScanResult(object):
    """
    Object to store a result from a sequence :class:`Scan`.

    Args:
        query_name (str): name of the query sequence
        hits ([:class:`ScanHit`]): hits for the query sequence
    """

    def __init__(self, *, query_name, hits, **kwargs):
        self.query_name = query_name
        self.hits = [ScanHit(**hit) for hit in hits]


class Scan(object):
    """Object to store a sequence scan.

    A scan can contain one or more :class:`ScanResult` objects
    (eg **query**).

    Args:
        results ([:class:`ScanResult`]): results of this scan (one result per query sequence)

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

        lines = []
        if header:
            headers = 'funfam members uniq_ec_terms query_region match_region evalue score description'.split()
            lines.append("\t".join(headers))

        result = self.results[0]

        if result:
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


class SsapResult(object):
    """
    Represents a result from SSAP pairwise structure comparison

    Args:
        id1 (str): id of protein 1 
        id2 (str): id of protein 2
        len1 (int): length of protein 1
        len2 (int): length of protein 2
        seqid (int): number of aligned residues

    """

    def __init__(self, *,
                 id1=None, id2=None,
                 len1=None, len2=None,
                 seqid=None, alnov=None,
                 ssap_score=None, rmsd=None,
                 alignment=None):
        self.id1 = id1
        self.id2 = id2
        self.len1 = len1
        self.len2 = len2
        self.seqid = seqid
        self.alnov = alnov
        self.ssap_score = ssap_score
        self.rmsd = rmsd
        self.alignment = alignment

    @classmethod
    def from_string(cls, ssap_line, *, alignment=None):
        ssap_line = ssap_line.strip()
        id1, id2, len1, len2, ssap_score, seqid, alnov, rmsd = ssap_line.split()
        return cls(id1=id1, id2=id2,
                   len1=len1, len2=len2, seqid=seqid,
                   alnov=alnov, ssap_score=ssap_score,
                   rmsd=rmsd, alignment=alignment,)
