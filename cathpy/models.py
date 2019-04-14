"""
Collection of classes used to model CATH data
"""

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

class Segment(object):
    """Class to represent a protein segment."""

    def __init__(self, start: int, stop: int):
        self.start = int(start)
        self.stop = int(stop)

    def __str__(self):
        return "{}-{}".format(self.start, self.stop)

    def __repr__(self):
        return "Segment:{}-{}".format(self.start, self.stop)


class ScanHsp(object):
    """Object to store the High Scoring Pair (HSP) from a sequence scan."""
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
    """Object to store a hit from a sequence scan."""
    def __init__(self, *, match_name, match_cath_id, match_description,
                 match_length, hsps, significance, data, **kwargs):
        self.match_name = match_name
        self.match_cath_id = match_cath_id
        self.match_description = match_description
        self.match_length = match_length
        self.hsps = [ScanHsp(**hsp) for hsp in hsps]
        self.data = data
        self.significance = significance

class ScanResult(object):
    """Object to store a result from a sequence scan."""
    def __init__(self, *, query_name, hits, **kwargs):
        self.query_name = query_name
        self.hits = [ScanHit(**hit) for hit in hits]

class Scan(object):
    """Object to store a sequence scan."""
    def __init__(self, *, results, **kwargs):
        self.results = [ScanResult(**res) for res in results]
