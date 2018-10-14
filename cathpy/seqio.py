# core
import io
import logging
import re
import functools

# local
from cathpy import error as err, util

logger = logging.getLogger(__name__)

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

    aa_by_id = { aa[2]: AminoAcid(aa[2], aa[1], aa[0]) for aa in (
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
    ) }

    def __init__(self):
        pass

    @classmethod
    def is_valid_aa(cls, aa):
        """Check if aa is a valid single character aa code."""
        return str(aa).upper() in cls.aa_by_id

    @classmethod
    def get_by_id(cls, aa):
        """Return the AminoAcid object by the given single character aa code."""
        return cls.aa_by_id[aa.upper()]

class Residue(object):
    """Class to represent a protein residue."""

    def __init__(self, aa, seq_num=None, pdb_label=None):
        assert type(aa) is str
        if seq_num:
            seq_num = int(seq_num)
        
        self.aa = aa
        self.seq_num = seq_num
        self.pdb_label = pdb_label

    def __str__(self):
        return self.aa

    def __repr__(self):
        return "res({},seq:{},pdb:{})".format(
            self.aa, self.seq_num, self.pdb_label)

class Segment(object):
    """Class to represent a protein segment."""

    def __init__(self, start: int, stop: int):
        self.start = int(start)
        self.stop = int(stop)
    
    def __str__(self):
        return("{}-{}".format(self.start, self.stop))

    def __repr__(self):
        return("Segment:{}-{}".format(self.start, self.stop))

class Sequence(object):
    """Class to represent a protein sequence."""

    re_gap_chars = r'[.\-]'

    def __init__(self, hdr: str, seq: str, *, meta=None):
        self._hdr = hdr
        self._seq = seq
        try:
            hdr_info = Sequence.split_hdr(hdr)
        except:
            raise err.GeneralError('caught error while parsing sequence header: '+hdr)
        self.id = hdr_info['id']
        self.accession = hdr_info['accession']
        self.id_type = hdr_info['id_type']
        self.id_ver = hdr_info['id_ver']
        self.segs = hdr_info['segs']
        self.meta = hdr_info['meta']
        if meta:
            for key, val in meta.items():
                self.meta[key] = val
    
    def get_residues(self):
        """
        Returns an array of Residue objects based on this sequence.
        
        Note: if segment information has been specified then this
        will be used to calculate the `seq_num` attribute.

        Raises:
            OutOfBoundsError: problem mapping segment info to sequence
        """
        residues = []
        
        segs = self.segs
        if not segs:
            segs = [ Segment(1, len(self.seq_no_gaps)) ]
        
        current_seg_offset=0
        def next_seg():
            nonlocal current_seg_offset
            if current_seg_offset < len(segs):
                seg = segs[current_seg_offset]
                current_seg_offset += 1
                return seg
            else:
                return None
        
        # theoretical length according to segment info vs length according to sequence
        seg_length = 0
        for seg in segs:
            seg_length += seg.stop - seg.start + 1
        actual_length = len(self.seq_no_gaps)

        if seg_length != actual_length:
            # should this be a warning? (with 1-n numbering as fallback?)
            raise err.OutOfBoundsError(('segment information {} suggests that the sequence '
                'length should be {}, but the sequence has {} (non-gap) characters: {}').format(
                    repr(segs), seg_length, actual_length, self.seq) )

        current_seg = next_seg()

        seq_num = current_seg.start
        for offset, aa in enumerate(self.seq, 0):
            if current_seg and seq_num > current_seg.stop:
                current_seg = next_seg()
                if not current_seg:
                    if not Sequence.is_gap(aa):
                        raise err.OutOfBoundsError(('unable to map segment ({}) to sequence: '
                            'the final segment ends at {}, but the sequence has {} residues '
                            '(offset: {}, aa: {})').format(
                                repr(current_seg), seq_num-1, len(self.seq_no_gaps), offset, aa
                            ))
                    else:
                        seq_num = None
                else:
                    seq_num = current_seg.start
            
            if Sequence.is_gap(aa):
                res = Residue(aa)
            else:
                res = Residue(aa, seq_num)
                seq_num += 1
            
            residues.append(res)
        return residues

    def get_res_at_offset(self, offset):
        """Return the residue character at the given offset (includes gaps)."""
        try:
            res = self.seq[offset]
        except:
            raise err.SeqIOError("Error: failed to get residue at offset {} from sequence with length {}: '{}'".format(
                offset, self.length(), self.seq))

        return res

    def get_res_at_seq_position(self, seq_pos):
        """Return the residue character at the given sequence position (ignores gaps)."""

        seq_nogap = re.sub(Sequence.re_gap_chars, '', self.seq)
        try:
            res = seq_nogap[seq_pos-1]
        except:
            raise err.SeqIOError("Error: failed to get residue at position {} from sequence with {} non-gap sequence positions: '{}'".format(
                seq_pos, len(seq_nogap), self.seq))

        return res

    def get_seq_position_at_offset(self, offset):
        """Return the sequence position (ignores gaps) of the given residue offset (may include gaps)."""

        seq_to_offset = self.seq[:offset+1]
        if re.match(seq_to_offset[-1], Sequence.re_gap_chars):
            raise err.GapError("Cannot get sequence position at offset {} since this corresponds to a gap".format(offset))

        seq_nogap = re.sub(Sequence.re_gap_chars, '', seq_to_offset)
        return len(seq_nogap)

    def get_offset_at_seq_position(self, seq_pos):
        """Return the offset (with gaps) of the given sequence position (ignores gaps)."""

        current_seq_pos=0
        for offset in range(len(self.seq)):
            if not re.match(Sequence.re_gap_chars, self.seq[offset]):
                current_seq_pos += 1
            if current_seq_pos == seq_pos:
                return offset

        raise err.OutOfBoundsError("failed to find offset at sequence position {}".format(seq_pos))

    def length(self):
        """Return the length of the sequence."""
        return len(self.seq)
    
    @property
    def seq(self):
        """Return the amino acid sequence as a string."""
        return self._seq

    @property
    def seq_no_gaps(self):
        """Return the amino acid sequence as a string (after removing all gaps)."""
        seq = re.sub(self.re_gap_chars, '', self._seq)
        return seq

    def _set_sequence(self, seq):
        self._seq = seq

    def set_cluster_id(self, id):
        self.meta['CLUSTER_ID'] = id

    def get_cluster_id(self):
        return self.meta['CLUSTER_ID'] if 'CLUSTER_ID' in self.meta else None

    @classmethod
    def split_hdr(cls, hdr: str) -> dict:
        """
        Splits a sequence header into meta information.
        
        Args:
            hdr (str): header string (eg `'domain|1cukA01|4_2_0/3-23_56-123'`)

        Returns:
            info (dict): header info

            {
                'id': 'domain|1cukA01|4_2_0/3-23_56-123', 
                'accession': '1cukA01', 
                'id_type': 'domain', 
                'id_ver': '4_2_0', 
                'segs': [Segment(3, 23), Segment(56,123)], 
                'meta': {}
            }

        """

        id = None
        accession = None
        id_type = None
        id_ver = None
        segs = []
        meta = {}

        if not hdr:
            raise err.ParamError('hdr seems to be empty')

        # split meta features (after whitespace)
        hdr_parts = hdr.split(maxsplit=1)
        id_with_segs_str = hdr_parts[0]
        meta_str = hdr_parts[1] if len(hdr_parts) > 1 else None

        # split id / segments
        id_with_segs_parts = id_with_segs_str.split('/', maxsplit=1)
        id_str = id_with_segs_parts[0]
        segs_str = id_with_segs_parts[1] if len(id_with_segs_parts) > 1 else None

        # split id into type, id, version
        id_parts = id_str.split('|')

        # 1cukA01/23-123
        if (len(id_parts) == 1):
            accession = id_parts[0]
        # domain|1cukA01/23-123
        if (len(id_parts) == 2):
            (id_type, accession) = id_parts
        # domain|1cukA01|v4.2.0/23-123
        if (len(id_parts) == 3):
            (id_type, accession, id_ver) = id_parts

        # segments
        if segs_str:
            for seg_str in segs_str.split('_'):
                (start, stop) = seg_str.split('-')
                seg = Segment(int(start), int(stop))
                segs.append(seg)

        # features
        if meta_str:
            meta_parts = meta_str.split()
            for f in meta_parts.split('=', maxsplit=1):
                if len(f) == 2:
                    meta[f[0]] = f[1]
                else:
                    logger.warning("failed to parse meta feature from string {}".format(meta_str))

        return({'accession': accession, 'id': id_with_segs_str, 'id_type': id_type, 'id_ver': id_ver, 
            'segs': segs, 'meta': meta})

    def to_fasta(self, wrap_width=80):
        """Return a string for this Sequence in FASTA format."""

        str = ""
        str += '>' + self.id + '\n'
        if wrap_width:
            for line in Sequence._chunker(self.seq, wrap_width):
                str += line + '\n'
        else:
            str += self.seq + '\n'
        return(str)

    def copy(self):
        """Provide a deep copy of this sequence."""
        s = Sequence(self._hdr, self.seq, meta=self.meta)
        return s

    def insert_gap_at_offset(self, offset, gap_char="-"):
        """Insert a gap into the current sequence at a given offset."""
        new_seq = self.seq[:offset] + gap_char + self.seq[offset:]
        self._set_sequence( new_seq )

    def set_gap_char_at_offset(self, offset, gap_char):
        """If the residue at a given position is a gap, then override the gap char with the given character."""
        residues = list(self.seq)
        if Sequence.is_gap(residues[offset]) and residues[offset] != gap_char:
            residues[offset] = gap_char
            self._set_sequence("".join(residues))

    def lower_case_at_offset(self, start, end=None):
        """Lower case the residues in the given sequence window."""
        if end == None:
            end = start + 1
        
        old_seq = self.seq
        new_seq = old_seq[:start] + old_seq[start:end].lower() + old_seq[end:]
        self._set_sequence(new_seq)

    def set_all_gap_chars(self, gap_char='-'):
        """Sets all gap characters."""
        seqstr = re.sub(self.re_gap_chars, gap_char, self.seq)
        self._set_sequence(seqstr)

    def set_lower_case_to_gap(self, gap_char='-'):
        seqstr = re.sub(r'[a-z]', gap_char, self.seq)
        self._set_sequence(seqstr)

    def slice_seq(self, start, end=None):
        """Return a slice of this sequence."""
        return self.seq[start:end]

    @staticmethod
    def _chunker(str, width):
        return (str[pos:pos + width] for pos in range(0, len(str), width))

    @staticmethod
    def is_gap(res_char):
        """Test whether a character is considered a gap."""
        return res_char in ['-', '.']

    def __str__(self):
        """Represents this Sequence as a string."""
        return('{:<30} {}'.format(self.id, self.seq))


class Correspondence(object):
    """
    Provides a mapping between ATOM and SEQRES residues.

    A correspondence is a type of alignment that provides the equivalences
    between the residues in the protein sequence (eg ``SEQRES`` records) and 
    the residues actually observed in the structure (eg ``ATOM`` records).
    """

    GCF_GAP_CHAR = '*'
    FASTA_GAP_CHAR = '-'

    def __init__(self, id=None, residues=None, **kwargs):
        """Create a new Correspondence object."""

        self.id = id
        self.residues = residues if residues else []
        super().__init__(**kwargs)

    @classmethod
    def new_from_gcf(cls, gcf_io):
        """Create a new Correspondence object from a GCF io / filename / string."""

        """
        >gi|void|ref1
        A   1   5   A
        K   2   6   K
        G   3   7   G
        H   4   8   H
        P   5   9   P
        G   6  10   G
        P   7  10A  P
        K   8  10B  K
        A   9  11   A
        P  10   *   *
        G  11   *   * 
        """

        corr = Correspondence()

        if (type(gcf_io) == str):
            if (gcf_io[0] == '>'):
                gcf_io = io.StringIO(gcf_io)
            else:
                gcf_io = open(gcf_io)

        try:
            hdr = gcf_io.readline()
            id = hdr.strip().split('|')[-1]
            corr.id = id

        except AttributeError:
            # make a potentially confusing error slightly less so
            raise err.SeqIOError("encountered an error trying to readline() on GCF io ({})".format(gcf_io))

        line_no = 1
        for line in gcf_io:
            line_no += 1

            try:
                seqres_aa, seqres_num, pdb_label, pdb_aa = line.split()
                if pdb_aa is not seqres_aa and pdb_aa is not Correspondence.GCF_GAP_CHAR:
                    logger.warning("pdb_aa '{}' does not match seqres_aa '{}' (line: {})".format(pdb_aa, seqres_aa, line_no))
            except:
                raise err.SeqIOError("Error: failed to parse GCF '{}' ({}:{})".format(
                    line, str(gcf_io), line_no
                ))

            if pdb_label is Correspondence.GCF_GAP_CHAR:
                pdb_label = None
                pdb_aa = None

            res = Residue(seqres_aa, int(seqres_num), pdb_label)
            corr.residues.append(res)
        
        gcf_io.close()

        return corr

    @property
    def seqres_length(self) -> int:
        """Return the number of SEQRES residues"""
        return len(self.residues)

    @property
    def atom_length(self) -> int:
        """Return the number of ATOM residues"""
        atom_residues = [res for res in self.residues if res.pdb_label is not None]
        return len(atom_residues)

    def get_res_at_offset(self, offset: int):
        """Return the ``Residue`` at the given offset (zero-based)"""
        return self.residues[offset]

    def get_res_by_seq_num(self, seq_num: int):
        """Return the ``Residue`` with the given sequence number"""
        res = next((res for res in self.residues if res.seq_num == seq_num), None)
        return res

    def get_res_by_pdb_label(self, pdb_label):
        res = next((res for res in self.residues if res.pdb_label == pdb_label), None)
        return res

    def get_res_by_atom_pos(self, pos):
        """Return the Residue corresponding to the given position in the ATOM sequence (ignores gaps)."""
        assert pos >= 1
        atom_residues = [res for res in self.residues if res.pdb_label is not None]
        res = atom_residues[pos-1]
        return res

    def get_res_offset_by_atom_pos(self, pos):
        """Return the offset of the Residue at the given position in the ATOM sequence (ignores gaps)."""
        assert pos >= 1
        atom_pos=1
        for offset, res in enumerate(self.residues):
            # logger.debug("pos({}) -> res: offset: {}, res: {}, atom_pos: {}".format(
            #     pos, offset, repr(res), atom_pos))
            if atom_pos == pos:
                return offset
            if res.pdb_label is not None:
                atom_pos += 1

        atom_residues = [res for res in self.residues if res.pdb_label is not None]
        raise err.OutOfBoundsError('failed to find residue in atom pos {}, last atom residue is {} (position {})'.format(
            pos, repr(atom_residues[-1]), atom_pos,
        ))

    @property
    def first_residue(self):
        """Returns the first residue in the correspondence."""
        return self.get_res_at_offset(0)

    @property
    def last_residue(self):
        """Returns the last residue in the correspondence."""
        return self.get_res_at_offset(-1)

    @property
    def atom_sequence(self):
        """Returns a Sequence corresponding to the ATOM records."""

        id = "atom|{}".format(self.id)
        res = [res.aa if res.pdb_label else Correspondence.FASTA_GAP_CHAR for res in self.residues]
        return Sequence(id, "".join(res))

    @property
    def seqres_sequence(self):
        """Returns a Sequence corresponding to the SEQRES records."""

        id = "seqres|{}".format(self.id)
        res = [res.aa for res in self.residues]
        return Sequence(id, "".join(res))

    def apply_seqres_segments(self, segs):
        """Returns a new correspondence from just the residues within the segments."""

        current_seg_offset = 0
        def next_seg():
            nonlocal current_seg_offset
            # logger.debug("apply_seqres_segments.next_seg: current={} segs={}".format(
            #     current_seg_offset, repr(segs) ))
            if current_seg_offset < len(segs):
                seg = segs[current_seg_offset]
                current_seg_offset += 1
                return seg

        current_seg = next_seg()
        selected_residues = []

        for res in self.residues:
            # logger.debug('apply_seqres.res: [{}] {}-{} seq_num={}'.format(
            #     current_seg_offset, current_seg.start, current_seg.stop,
            #     res.seq_num))
            
            if res.seq_num >= current_seg.start and res.seq_num <= current_seg.stop:
                selected_residues.append(res)
            elif res.seq_num < current_seg.start:
                pass
            elif res.seq_num > current_seg.stop:
                current_seg = next_seg()
                if not current_seg:
                    break
            else:
                raise err.SeqIOError("unexpected error - shouldn't be able to reach this code")

        corr = __class__(id=self.id, residues=selected_residues)

        return corr

    def to_gcf(self):
        """Renders the current object as a GCF string."""

        """
        >gi|void|ref1
        A   1   5   A
        K   2   6   K
        G   3   7   G
        H   4   8   H
        P   5   9   P
        G   6  10   G
        P   7  10A  P
        K   8  10B  K
        A   9  11   A
        P  10   *   *
        G  11   *   * 
        """

        gcf_str = '>' + self.id + '\n'
        for res in self.residues:
            if res.pdb_label:
                vals = [res.aa, res.seq_num, res.pdb_label, res.aa]
            else:
                vals = [res.aa, res.seq_num, '*', '*']
            gcf_str += '{} {:>4} {:>4}  {}\n'.format(*vals)

        return gcf_str


    def to_sequences(self):
        seqs = (self.seqres_sequence, self.atom_sequence)
        return seqs

    def to_fasta(self, **kwargs):
        """Returns the Correspondence as a string (FASTA format)."""
        seqs = self.to_sequences()
        return seqs[0].to_fasta(**kwargs) + seqs[1].to_fasta(**kwargs)
        
    def to_aln(self):
        """Returns the Correspondence as an Align object."""
        seqs = self.to_sequences()
        return Align(seqs = seqs)
        
    def __str__(self):
        return self.to_fasta()

    def __repr__(self):
        return self.to_fasta()

class Align(object):
    """Object representing a protein sequence alignment."""

    REF_GAP_CHAR='-'
    MERGE_GAP_CHAR='.'

    STO_META_TO_ATTR=[
        ('ID', 'id'),
        ('DE', 'description'),
        ('AC', 'accession'),
        ('TP', 'aln_type'),
        ('DR', {
            'CATH': 'cath_version',
            'DOPS': 'dops_score', 
        }),
        ('TC', 'min_bitscore'),
        ('SQ', None),
    ]

    def __init__(self, seqs=None, *, id=None, accession=None, 
        description=None, aln_type=None, min_bitscore=None, 
        cath_version=None, dops_score=None):
        self.meta = {}
        self.id = id
        self.accession = accession
        self.description = description
        self.aln_type = aln_type
        self.min_bitscore = min_bitscore
        self.cath_version = cath_version
        self.dops_score = dops_score
        self.seqs = seqs if seqs else []
        self.__aln_positions = 0
        self._merge_counter = 0

    def _next_merge_id(self):
        self._merge_counter += 1
        return self._merge_counter

    @property
    def aln_positions(self):
        """Return the number of alignment positions."""
        return self.__aln_positions
    
    @aln_positions.setter
    def aln_positions(self, value):
        self.__aln_positions = value

    @property
    def count_sequences(self):
        """Return the number of sequences in the alignment."""
        return len(self.seqs)

    def find_first_seq_by_accession(self, acc):
        """Return the first Sequence with the given accession."""
        seqs_with_acc = [seq for seq in self.seqs if seq.accession == acc]
        return seqs_with_acc[0]

    def find_seq_by_id(self, id):
        """Return the Sequence corresponding to the provided id."""
        seqs_with_id = [seq for seq in self.seqs if seq.id == id]
        if len(seqs_with_id) > 1:
            raise err.SeqIOError("Found more than one ({}) sequence matching id '{}'".format(
                len(seqs_with_id), id))
        if len(seqs_with_id) == 0:
            raise err.NoMatchesError('failed to find sequence with id {} in alignment'.format(id))
        return seqs_with_id[0]

    def find_seq_by_accession(self, acc):
        """Return the Sequence corresponding to the provided id."""
        seqs_with_acc = [seq for seq in self.seqs if seq.accession == acc]
        if len(seqs_with_acc) > 1:
            raise err.TooManyMatchesError("Found more than one ({}) sequence matching accession '{}'".format(
                len(seqs_with_acc), acc), )
        if len(seqs_with_acc) == 0:
            raise err.NoMatchesError('failed to find sequence with accession {} in alignment'.format(acc))
        return seqs_with_acc[0]

    def get_seq_at_offset(self, offset):
        """Return the Sequence at the given offset (zero-based)."""
        return self.seqs[offset]

    @classmethod
    def new_from_fasta(cls, fasta_io):
        """Initialise an alignment object from a fasta file / string / io"""
        aln = Align()
        aln.read_sequences_from_fasta(fasta_io)
        return aln

    @staticmethod
    def _get_io_from_file_or_string(file_or_string):

        if isinstance(file_or_string, str):
            if file_or_string[0] == '>':
                _io = io.StringIO(file_or_string)
            else:
                _io = open(file_or_string)
        elif isinstance(file_or_string, io.IOBase):
            _io = file_or_string
        else:
            _io = file_or_string
            logger.warning("unexpected io type: " + repr(file_or_string))
        
        return _io

    @classmethod
    def new_from_stockholm(cls, sto_io):
        """Initialise an alignment object from a STOCKHOLM file / string / io""" 

        sto_io = __class__._get_io_from_file_or_string(sto_io)

        aln = cls()
        sto_header = sto_io.readline()
        
        assert sto_header.startswith('# STOCKHOLM 1.0')

        aln_meta = {}
        seq_meta_by_id = {}
        seq_aa_by_id = {}

        gc_meta_to_attr = { meta: attr for (meta, attr) in cls.STO_META_TO_ATTR }

        for line in sto_io:

            line = line.strip()

            if line.startswith('#=GF'):
                _, feature, per_file_ann = line.split(None, 2)

                if feature not in gc_meta_to_attr:
                    raise Exception('encountered unexpected GF tag {} in line "{}" (known tags: {})'.format(
                        feature, line, repr(gc_meta_to_attr)))

                attr = gc_meta_to_attr[feature]
                if type(attr) is dict:
                    key, value = per_file_ann.split(':')
                    key = key.strip()
                    value = value.strip()
                    if key not in attr:
                        raise Exception('encountered unexpected GF tag {}->{} in line "{}" (known tags: {})'.format(
                            feature, key, line, repr(attr)))
                    attr = attr[key]
                    per_file_ann = value
                
                if attr:
                    logger.debug('setting aln attr "{}" to "{}"'.format(attr, per_file_ann))
                    setattr(aln, attr, per_file_ann)

            elif line.startswith('#=GC'):
                _, feature, per_col_ann = line.split(None, 2)
                aln_meta[feature] = per_col_ann

            elif line.startswith('#=GS'):
                _, seq_id, feature, per_seq_ann = line.split(None, 3)
                if feature == 'DR':
                    dr_type, per_seq_ann = per_seq_ann.split(None, 1)
                    feature = feature + '_' + dr_type                
                if seq_id not in seq_meta_by_id:
                    seq_meta_by_id[seq_id] = {}
                seq_meta_by_id[seq_id][feature] = per_seq_ann

            elif line.startswith('#=GR'):
                _, seq_id, feature, per_res_ann = line.split(None, 3)
                seq_meta_by_id[seq_id][feature] = per_res_ann

            elif line.startswith('//'):
                pass
            
            else:
                seq_id, seq_aa = line.split()
                if seq_id not in seq_aa_by_id:
                    seq_aa_by_id[seq_id] = ''
                seq_aa_by_id[seq_id] += seq_aa

        for seq_id, seq_aa in seq_aa_by_id.items():
            seq_meta = seq_meta_by_id[seq_id] if seq_id in seq_meta_by_id else {}
            seq = Sequence(seq_id, seq_aa, meta=seq_meta)
            aln.add_sequence(seq)

        for key, val in aln_meta.items():
            aln.meta[key] = val

        sto_io.close()

        return aln

    def read_sequences_from_fasta(self, fasta_io):
        """Parse aligned sequences from FASTA (str, file, io) and adds them to the current
        Align object. Returns the number of sequences that are added."""

        fasta_io = __class__._get_io_from_file_or_string(fasta_io)

        re_seqstr = re.compile( r'^[a-zA-Z.\-]+$' )

        seq_added=0            
        current_hdr = None
        current_seq = ''
        line_count=0
        for line in fasta_io:
            line_count += 1
            if line == "":
                break
            
            line = line.rstrip()
            if line[0] == '>':
                if current_seq:
                    seq = Sequence(current_hdr, current_seq)
                    self.add_sequence(seq)
                    current_seq = ''
                    seq_added += 1
                current_hdr = line[1:]                
            else:
                if not re_seqstr.match(line):
                    raise err.SeqIOError(('encountered an error parsing FASTA: '
                        'string "{}" does not look like a sequence ({}, line {})').format(
                            line, str(fasta_io), line_count))

                if not current_hdr:
                    raise err.SeqIOError(('encountered an error parsing FASTA: '
                        'found sequence "{}" without a header ({}, line {})').format(
                            line, str(fasta_io), line_count))
                
                current_seq += str(line)

        fasta_io.close()

        if current_seq:
            seq = Sequence(current_hdr, current_seq)
            self.add_sequence(seq)
            seq_added += 1

        return seq_added

    def add_sequence(self, seq: Sequence):
        """Add a sequence to this alignment."""

        if self.aln_positions:
            if self.aln_positions != seq.length():
                raise err.SeqIOError( ("Error: cannot add a sequence (id:{}) "
                    "with {} positions to an alignment with {} positions.")
                    .format(seq.id, seq.length(), self.aln_positions) )
        else:
            self.__aln_positions = seq.length()

        self.seqs.append(seq)
        return seq        

    def remove_sequence_by_id(self, seq_id: str):
        """Removes a sequence from the alignment."""

        for idx, seq in enumerate(self.seqs):
            if seq.id == seq_id:
                logger.info("Removing sequence with '{}' from alignment".format(seq_id))
                del self.seqs[idx]
                return seq
        
        raise err.NoMatchesError('failed to find sequence with id {}'.format(seq_id))


    def remove_alignment_gaps(self):
        """Return a new alignment after removing alignment positions 
        that contain a gap for all sequences."""
        seqs = self.seqs
        seq_length = seqs[0].length()
        new_seq_strings = ["" for s in range(len(seqs))]
        for aln_offset in range(seq_length):
            total_gaps = 0
            for seq in seqs:
                if seq.seq[aln_offset] == '-' or seq.seq[aln_offset] == '.':
                    total_gaps += 1
            if total_gaps < len(seqs):
                for seq_pos in range(len(seqs)):
                    res = seqs[seq_pos].seq[aln_offset]
                    # print( "seq[{}:{}] pos:{} res:{}".format(aln_offset, seqs[seq_pos].id, seq_pos, res) )
                    new_seq_strings[seq_pos] += res
            else:
                logger.info("Removing complete gap from alignment offset: {}".format(aln_offset))

        new_aln = Align()
        for seq_pos in range(len(new_seq_strings)):
            hdr = seqs[seq_pos]._hdr
            seq_str = new_seq_strings[seq_pos]
            seq = Sequence(hdr, seq_str)
            new_aln.add_sequence(seq)

        return new_aln

    def insert_gap_at_offset(self, offset, gap_char='-'):
        """Insert a gap char at the given offset (zero-based)."""
        self.__aln_positions += 1
        for s in self.seqs:
            s.insert_gap_at_offset(offset, gap_char) 

    def set_gap_char_at_offset(self, offset, gap_char):
        """Override the gap char for all sequences at a given offset."""
        for s in self.seqs:
            s.set_gap_char_at_offset(offset, gap_char)

    def lower_case_at_offset(self, start, end=None):
        """Lower case all the residues in the given alignment window."""
        for s in self.seqs:
            s.lower_case_at_offset(start, end)

    def slice_seqs(self, start, end=None):
        """Return an array of Sequence objects from start to end."""
        return [Sequence(s._hdr, s.slice_seq(start, end)) for s in self.seqs]

    def merge_alignment(self, merge_aln, ref_seq_acc: str, ref_correspondence = None, *, 
            cluster_label=None):
        """
        Merges aligned sequences into the current object via a reference sequence.

        Sequences in ``merge_aln`` are brought into the current alignment using
        the equivalences identified in reference sequence `ref_seq_acc` (which
        must exist in both the ``self`` and ``merge_aln``).

        This function was originally written to merge FunFam alignments
        according to structural equivalences identified by CORA (a multiple
        structural alignment tool). Moving between structure and sequence 
        provides the added complication that
        sequences in the structural alignment (CORA) are based on ATOM records,
        whereas sequences in the merge alignment (FunFams) are based on SEQRES
        records. The ``ref_correspondence`` argument allows this mapping to be
        taken into account.

        Args: 
            merge_aln (Align): An Align containing the reference
                sequence and any additional sequences to merge. 
            ref_seq_acc (str): The accession that will be used to find the 
                reference sequence in the current alignment and merge_aln 
            ref_correspondence (Correspondence): An optional Correspondence 
                object that provides a mapping between the reference 
                sequence found in ``self`` (ATOM records) and reference 
                sequence as it appears in ``merge_aln`` (SEQRES records).

        Returns: 
            [Sequence]: Array of Sequences added to the current alignment.

        Raises: 
            MergeCorrespondenceError: problem mapping reference
                sequence between alignment and correspondence 

        """

        merge_aln = merge_aln.copy()

        if not cluster_label:
            cluster_label = self._next_merge_id()

        for seq in merge_aln.seqs:
            seq.set_cluster_id(cluster_label) 

        ref_seq_in_ref = self.find_seq_by_accession(ref_seq_acc)

        ref_seq_in_merge = merge_aln.find_seq_by_accession(ref_seq_acc)

        if ref_correspondence is None:
            # fake a 1:1 correspondence for internal use
            # ignore any residue that does not have a seq_num (ie gap)
            residues = [res for res in ref_seq_in_ref.get_residues() if res.seq_num]
            for r in residues:
                r.pdb_label = str(r.seq_num)
                # logger.debug("fake correspondence: residue={}".format(repr(r)))
            ref_correspondence = Correspondence(id=ref_seq_acc, residues=residues)

        # check: ref sequence (in self) must match the ATOM sequence in Correspondence
        ref_no_gaps = ref_seq_in_ref.seq_no_gaps
        corr_no_gaps = ref_correspondence.atom_sequence.seq_no_gaps
        if ref_no_gaps != corr_no_gaps:
            raise err.MergeCorrespondenceError(ref_seq_acc, 'current', 'ATOM', 
                ref_no_gaps, corr_no_gaps) 

        # check: ref sequence (in merge) must match the SEQRES sequence in Correspondence
        ref_no_gaps = ref_seq_in_merge.seq_no_gaps
        corr_no_gaps = ref_correspondence.seqres_sequence.seq_no_gaps
        if ref_no_gaps != corr_no_gaps:
            raise err.MergeCorrespondenceError(ref_seq_acc, 'merge', 'SEQRES', 
                ref_no_gaps, corr_no_gaps) 

        # clean up
        del ref_no_gaps
        del corr_no_gaps

        ref_aln_pos=0
        ref_corr_pos=0
        merge_aln_pos=0
        correspondence_length = ref_correspondence.seqres_length

        logger.debug("ref_alignment.positions: {}".format(self.aln_positions))
        logger.debug("merge_alignment.positions: {}".format(merge_aln.aln_positions))
        logger.debug("ref_seq_in_ref:   {}".format(str(ref_seq_in_ref)))
        logger.debug("ref_seq_in_merge: {}".format(str(ref_seq_in_merge)))

        while True:

            if merge_aln_pos >= merge_aln.aln_positions \
                and ref_aln_pos  >= self.aln_positions \
                and ref_corr_pos >= correspondence_length:
                break

            logger.debug("REF {}/{}; CORRESPONDENCE {}/{}; MERGE {}/{}".format(
                ref_aln_pos, self.aln_positions, ref_corr_pos, 
                correspondence_length, merge_aln_pos, merge_aln.aln_positions))

            # sort the gaps in the reference alignment
            if ref_aln_pos < self.aln_positions:

                for seq in self.slice_seqs(0, ref_aln_pos):
                    logger.debug( "{:<10} {}".format("REF", str(seq)) )

                ref_res_in_ref = ref_seq_in_ref.get_res_at_offset(ref_aln_pos)

                logger.debug("REF_POSITION   {:>3} of {:>3} => '{}'".format(
                    ref_aln_pos, self.aln_positions, ref_res_in_ref))

                # insert all the gaps in the reference alignment into the merge sequences
                # keep doing this until we don't have any more gaps
                if Sequence.is_gap(ref_res_in_ref):
                    logger.debug(("GAP '{}' in ref sequence in REF alignment [{}], "
                        "inserting gap '{}' at position [{}] in all merge sequences").format(
                            ref_res_in_ref, ref_aln_pos, ref_res_in_ref, merge_aln_pos))
                    merge_aln.insert_gap_at_offset(merge_aln_pos, gap_char=ref_res_in_ref)

                    # this is a gap: do NOT increment ref_corr_pos
                    ref_aln_pos   += 1
                    merge_aln_pos += 1
                    continue

            # sort the gaps in the merge alignment
            if merge_aln_pos < merge_aln.aln_positions:

                # for seq in merge_aln.slice_seqs(0, merge_aln_pos):
                #     logger.debug( "{:<10} {}".format("MERGE", str(seq)) )

                ref_res_in_merge = ref_seq_in_merge.get_res_at_offset(merge_aln_pos)

                logger.debug("MERGE_POSITION {:>3} of {:>3} => '{}'".format(
                    ref_aln_pos, self.aln_positions, ref_res_in_ref))

                # insert all the gaps in the merge alignment into the ref sequences
                # keep doing this until we don't have any more gaps
                if Sequence.is_gap(ref_res_in_merge):
                    logger.debug(("GAP '{}' in ref sequence in MERGE alignment [{}], "
                        "inserting gap '{}' at position [{}] in all ref sequences").format(
                            ref_res_in_merge, merge_aln_pos, Align.MERGE_GAP_CHAR, merge_aln_pos))
                    self.insert_gap_at_offset(ref_aln_pos, gap_char=Align.MERGE_GAP_CHAR)
                    merge_aln.lower_case_at_offset(merge_aln_pos)
                    merge_aln.set_gap_char_at_offset(merge_aln_pos, '.')

                    #ref_corr_pos  += 1
                    ref_aln_pos   += 1
                    merge_aln_pos += 1
                    continue

            # if there are gaps in the correspondence then we add gaps to the ref sequence here
            if ref_corr_pos < correspondence_length:

                for seq in ref_correspondence.to_sequences():
                    seq = seq.slice_seq(0, ref_corr_pos)
                    logger.debug( "{:<10} {}".format("CORR", str(seq)) )

                ref_res_in_corr = ref_correspondence.get_res_at_offset(ref_corr_pos)
                if ref_res_in_corr.pdb_label is None:

                    logger.debug(("GAP '{}' in ATOM records of correspondence [{}], "
                        "inserting gap '{}' at position [{}] in ref sequences").format(
                            '*', ref_corr_pos, Align.MERGE_GAP_CHAR, ref_aln_pos))

                    #merge_aln.insert_gap_at_offset(merge_aln_pos, gap_char=Align.MERGE_GAP_CHAR)
                    self.insert_gap_at_offset(ref_aln_pos, gap_char=Align.MERGE_GAP_CHAR)
                    merge_aln.lower_case_at_offset(merge_aln_pos)
                    merge_aln.set_gap_char_at_offset(merge_aln_pos, '.')

                    # IMPORTANT: do not increment merge_aln_pos
                    ref_corr_pos  += 1
                    ref_aln_pos   += 1
                    merge_aln_pos += 1
                    continue

            ref_corr_pos  += 1
            ref_aln_pos   += 1
            merge_aln_pos += 1

        logger.info("FINISHED MERGE")
        # for seq in ref_correspondence.to_sequences():
        #     seq = seq.slice_seq(0, ref_corr_pos)
        #     logger.debug( "{:<10} {}".format("CORR", str(seq)) )
        # for seq in self.seqs:
        #     logger.debug( "{:<10} {}".format("REF", str(seq)) )
        # for seq in merge_aln.seqs:
        #     logger.debug( "{:<10} {}".format("MERGE", str(seq)) )

        # add the merged sequences into this alignment
        for seq in merge_aln.seqs:
            self.add_sequence(seq)
        
        # remove the original reference sequence (from the structural sequence)
        self.remove_sequence_by_id(ref_seq_in_ref.id)

        # for seq in self.seqs:
        #     logger.debug( "{:<10} {}".format("MERGED", str(seq)) )

        # test the final, merged alignment
        # 1. get sequences that correspond to the input aln
        # 2. remove alignment positions where there's a gap in the reference sequence

        logger.debug("Checking merge results for {} ({}) ...".format(ref_seq_acc, repr(ref_seq_in_merge._hdr)))
        for original_seq in merge_aln.seqs:
            seq = self.find_seq_by_accession(original_seq.accession)

            # logger.debug('Working on sequence: {}'.format(str(original_seq)))

            # this provides the residues in the merge alignment with seqres numbering
            ref_merge_residues = ref_seq_in_merge.get_residues()
            # the lookup lets us go from the seq numbering to the sequence offset
            ref_merge_seqnum_to_seqpos = {}
            for seq_pos, res in enumerate([res for res in ref_merge_residues if res.seq_num], 1):
                ref_merge_seqnum_to_seqpos[res.seq_num] = seq_pos

            if not seq:
                raise err.SeqIOError("failed to find sequence with id '{}' in merge aln".format(seq.id))

            for aln_offset in range(self.aln_positions):

                ref_res = ref_seq_in_ref.get_res_at_offset(aln_offset)

                merged_res_at_aln_offset = seq.get_res_at_offset(aln_offset)

                if ref_res == self.MERGE_GAP_CHAR:
                    # everything else should be a '.' or a lowercase residue
                    assert merged_res_at_aln_offset == '.' or re.match(r'[a-z]', merged_res_at_aln_offset)
                elif ref_res == self.REF_GAP_CHAR:
                    # everything else should be a '-' or an uppercase residue
                    assert merged_res_at_aln_offset == '-' or re.match(r'[A-Z]', merged_res_at_aln_offset)
                else:
                    # find the sequence offset of this aln position in the ref sequence
                    ref_seq_pos_in_ref = ref_seq_in_ref.get_seq_position_at_offset(aln_offset)

                    # use the correspondence to find the equivalent reference residue in the merge alignment
                    ref_corr_res = ref_correspondence.get_res_by_atom_pos(ref_seq_pos_in_ref)
                    ref_seq_num_in_merge = ref_corr_res.seq_num

                    if ref_seq_num_in_merge is None:
                        raise err.GeneralError(('weird... found a residue without a seq_num in the correspondence record '
                            ' ref_seq_pos_in_ref: {}, res: {}, corr: {}').format(
                                ref_seq_pos_in_ref, repr(ref_corr_res), repr(ref_correspondence)))

                    if ref_seq_num_in_merge not in ref_merge_seqnum_to_seqpos:
                        raise err.OutOfBoundsError(('failed to find seq_num {} ({}) in seqnum/seqpos '
                            'lookup: {}\ncorrespondence (length: {})').format(
                                ref_seq_num_in_merge, repr(ref_corr_res), ref_merge_seqnum_to_seqpos, 
                                ref_correspondence.seqres_length, ))

                    # find out where this seq_num occurs in the merge sequence (account for segment numbering)
                    ref_seq_pos_in_merge = ref_merge_seqnum_to_seqpos[ref_seq_num_in_merge]

                    # find the aln offset for the equivalent position in the original merge alignment
                    ref_merge_offset = ref_seq_in_merge.get_offset_at_seq_position(ref_seq_pos_in_merge)

                    # logger.debug("ref_seq_pos (ref): {}, ref_seq_pos (merge): {}, correspondence_res: {}, ref_merge_offset: {}".format(
                    #     ref_seq_pos_in_ref, ref_seq_pos_in_merge, repr(ref_corr_res), ref_merge_offset
                    # ))

                    # find the residue at the equivalent position in the merge alignment
                    original_res = original_seq.get_res_at_offset(ref_merge_offset)

                    if merged_res_at_aln_offset != original_res:
                        raise err.MergeCheckError(("Expected the merged residue '{}' to "
                            "match the original residue '{}' at alignment "
                            "offset {} (sequence: '{}')\n\n"
                            "CORR_ATOM:        {}\n"
                            "CORR_SEQRES:      {}\n"
                            "\n\n"
                            "REF_SEQ_IN_REF:   {}\n"
                            "REF_SEQ_IN_MERGE: {}\n"
                            "ORIGINAL_SEQ:     {}\n"
                            "                  {aln_pointer:>{merge_pos}}\n"
                            "MERGED_SEQ:       {}\n"
                            "                  {aln_pointer:>{aln_pos}}\n"
                            "(aln_offset={}, seq_pos(ref)={}, seq_num(merge)={}, seq_pos(merge)={}, ref_merge_offset={})"
                            ).format(
                                merged_res_at_aln_offset, original_res, aln_offset, seq.id,
                                ref_correspondence.atom_sequence,
                                ref_correspondence.seqres_sequence,
                                ref_seq_in_ref.seq,
                                ref_seq_in_merge.seq,
                                original_seq.seq,
                                seq.seq,
                                aln_offset, ref_seq_pos_in_ref, ref_seq_num_in_merge, ref_seq_pos_in_merge, ref_merge_offset,
                                aln_pointer='^', aln_pos=(aln_offset+1), merge_pos=(ref_merge_offset+1)
                            ))

        logger.info("Finshed checking merge for {} ({})".format(ref_seq_acc, repr(ref_seq_in_merge._hdr)))

        return merge_aln.seqs

    def copy(self):
        """Return a deepcopy of this object."""
        new_aln = Align()
        new_seqs = [s.copy() for s in self.seqs]
        new_aln.seqs = new_seqs
        new_aln.aln_positions = new_aln.seqs[0].length()
        return new_aln

    def write_fasta(self, fasta_file, wrap_width=80):
        """Write the alignment to a file in FASTA format."""
        with open(fasta_file, 'w') as f:
            for seq in self.seqs:
                f.write( seq.to_fasta(wrap_width=wrap_width) )

    def add_scorecons(self):
        """Add scorecons annotation to this alignment."""
        scons = util.ScoreconsRunner()
        logger.info("Calculating scorecons / DOPS ...")
        # output alignment to tmp fasta file
        scons_result = scons.run_alignment(self)
        self.dops_score = scons_result.dops
        self.meta['scorecons'] = scons_result.to_string

    def add_groupsim(self):
        """Add groupsim annotation to this alignment."""
        gs = util.GroupsimRunner()
        logger.info("Calculating GroupSim ...")
        # output alignment to tmp fasta file
        gs_result = gs.run_alignment(self)
        self.meta['groupsim'] = gs_result.to_string

    def write_sto(self, sto_file, *, meta=None):
        """Write the alignment to a file in STOCKHOLM format."""

        # putting these here to separate the data from the formatting
        sto_format='1.0'

        # allow meta keys to be provided in args, otherwise fill with the 
        # appropriate alignment attributes
        aln_meta = {}
        if meta:
            for key, attr in self.STO_META_TO_ATTR:
                aln_meta[key] = meta.get(key, None)

        comment_pad=0
        for seq in self.seqs:
            comment_pad = max(comment_pad, len(seq.id) + 1)

        seq_pad = comment_pad + 8
        gc_pad  = seq_pad - 5

        # single data point about the file
        def _GF(f, key, val):
            f.write('#=GF {} {}\n'.format(key, val))

        # single data point about each sequence
        def _GS(f, seq_id, key, val):
            if key.startswith('DR_'):
                val = key[3:] + ' ' + val
                key = 'DR'
            f.write('#=GS {:<{comment_pad}} {} {}\n'.format(seq_id, key, val, comment_pad=comment_pad))

        # positional data about the file
        def _GC(f, key, per_pos_str):
            f.write('#=GC {:<{gc_pad}} {}\n'.format(key, per_pos_str, 
                gc_pad=gc_pad))

        # positional data about each sequence
        def _GR(f, seq_id, key, per_pos_str):
            f.write('#=GR {:<{comment_pad}} {} {}\n'.format(seq_id, key, per_pos_str, comment_pad=comment_pad))
        
        def _SEQ(f, seq):
            f.write('{:<{seq_pad}} {}\n'.format(seq.id, seq.seq, seq_pad=seq_pad))

        def _START(f):
            f.write('# STOCKHOLM {}\n'.format(sto_format))

        def _END(f):
            f.write('//\n')

        with open(sto_file, 'w') as f:
            _START(f)
            _GF(f, 'ID', aln_meta.get('ID', self.id) )
            _GF(f, 'DE', aln_meta.get('DE', self.description) )
            _GF(f, 'AC', aln_meta.get('AC', self.accession) )
            _GF(f, 'TP', aln_meta.get('TP', self.aln_type) )

            if self.cath_version:
                _GF(f, 'DR', 'CATH: ' + self.cath_version)

            if self.dops_score:
                _GF(f, 'DR', 'DOPS: {:.3f}'.format(float(self.dops_score)))

            for seq in self.seqs:
                for key, val in seq.meta.items():
                    _GS(f, seq.id, key, val)

            if self.min_bitscore:
                _GF(f, 'TC', self.min_bitscore)

            _GF(f, 'SQ', self.count_sequences)

            for seq in self.seqs:
                _SEQ(f, seq)

            for key, val in sorted(self.meta.items()):
                _GC(f, key, val)

            _END(f)

    def __str__(self):
        return "\n".join( [ str(seq) for seq in self.seqs ] )
