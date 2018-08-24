import re
import logging
import io

class SeqIOError(Exception):
    """General Exception class within the SeqIO module"""
    def __init__(self, str):
        super().__init__(str)

class SeqIOGapError(SeqIOError):
    """Exception raised when trying to find residue information about a gap position."""
    def __init__(self, str):
        super().__init__(str)

class SeqIOMergeCorrespondenceError(SeqIOError):
    """Exception raised when failing to match correspondence sequences during alignment merge."""

    def __init__(self, seq_id, aln_type, seq_type, ref_no_gaps, corr_no_gaps):
        message = ("Reference sequence '{}' in the {} alignment does not "
            "match the {} sequence provided in the Correspondence.\n"
            "{:>25}: {}\n"
            "{:>25}: {}\n").format(
                seq_id,
                aln_type,
                seq_type,
                'REF [{}]'.format(len(ref_no_gaps)), 
                ref_no_gaps,
                'CORRESPONDENCE [{}]'.format(len(corr_no_gaps)), 
                corr_no_gaps,
            )
        super().__init__(message)

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

    def __init__(self, aa, seq_num, pdb_label=None):
        assert type(aa) is str
        assert type(seq_num) is int
        self.aa = aa
        self.seq_num = seq_num
        self.pdb_label = pdb_label

    def __str__(self):
        return self.aa

    def __repr__(self):
        return "res: aa={}, seq_num={}, pdb_label={}".format(
            self.aa, self.seq_num, self.pdb_label)

class Segment(object):
    """Class to represent a protein segment."""

    def __init__(self, start: int, stop: int):
        self.start = start
        self.stop  = stop
    
    def __str__(self):
        return("seg: {}-{}".format(self.start, self.stop))

class Sequence(object):
    """Class to represent a protein sequence."""

    re_gap_chars = r'[.\-]'

    def __init__(self, hdr: str, seq: str):
        self.hdr = hdr
        self._seq = seq
        hdr_info = Sequence.split_hdr(hdr)
        self.id = hdr_info['id']
        self.id_type = hdr_info['id_type']
        self.id_ver = hdr_info['id_ver']
        self.segs = hdr_info['segs']
        self.features = hdr_info['features']
    
    def get_residues(self):
        """Return an array of Residue objects for this sequence"""
        residues = []
        for seq_num, aa in enumerate(self.seq, 1):
            res = Residue(aa, seq_num)
            residues.append(res)
        return residues

    def get_res_at_offset(self, offset):
        """Return the residue character at the given offset (includes gaps)."""
        try:
            res = self.seq[offset]
        except:
            raise SeqIOError("Error: failed to get residue at offset {} from sequence with length {}: '{}'".format(
                offset, self.length(), self.seq))

        return res

    def get_res_at_seq_position(self, seq_pos):
        """Return the residue character at the given sequence position (ignores gaps)."""
        seq_nogap = re.sub(Sequence.re_gap_chars, '', self.seq)
        try:
            res = seq_nogap[seq_pos-1]
        except:
            raise SeqIOError("Error: failed to get residue at position {} from sequence with {} non-gap sequence positions: '{}'".format(
                seq_pos, len(seq_nogap), self.seq))

        return res

    def get_seq_position_at_offset(self, offset):
        """Return the sequence position (ignores gaps) of the given residue offset (may include gaps)."""
        seq_to_offset = self.seq[:offset+1]
        if re.match(seq_to_offset[-1], Sequence.re_gap_chars):
            raise SeqIOGapError("Cannot get sequence position at offset {} since this corresponds to a gap".format(offset))

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

        raise SeqIOError("failed to find offset at sequence position {}".format(seq_pos))

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
        logging.debug("seq_no_gaps: BEFORE: " + self._seq)
        logging.debug("seq_no_gaps: AFTER:  " + seq)
        return seq

    def _set_sequence(self, seq):
        self._seq = seq

    @staticmethod
    def split_hdr(hdr: str) -> dict:
        id = None
        id_type = None
        id_ver = None
        segs = []
        features = {}

        # split meta features (after whitespace)
        hdr_parts = hdr.split(maxsplit=1)
        id_with_segs_str = hdr_parts[0]
        features_str = hdr_parts[1] if len(hdr_parts) > 1 else None

        # split id / segments
        id_with_segs_parts = id_with_segs_str.split('/', maxsplit=1)
        id_str = id_with_segs_parts[0]
        segs_str = id_with_segs_parts[1] if len(id_with_segs_parts) > 1 else None

        # split id into type, id, version
        id_parts = id_str.split('|')

        # 1cukA01/23-123
        if (len(id_parts) == 1):
            id = id_parts[0]
        # domain|1cukA01/23-123
        if (len(id_parts) == 2):
            (id_type, id) = id_parts
        # domain|1cukA01|v4.2.0/23-123
        if (len(id_parts) == 3):
            (id_type, id, id_ver) = id_parts

        # segments
        if segs_str:
            for seg_str in segs_str.split('_'):
                (start, stop) = seg_str.split('-')
                seg = Segment(int(start), int(stop))
                segs.append(seg)

        # features
        if features_str:
            features_parts = features_str.split()
            for f in features_parts.split('=', maxsplit=1):
                if len(f) == 2:
                    features[f[0]] = f[1]
                else:
                    logging.warning("failed to parse feature from string {}".format(features_str))

        return({'id': id, 'id_type': id_type, 'id_ver': id_ver, 
            'segs': segs, 'features': features})

    def to_fasta(self, wrap_width=80):
        """Return a string for this Sequence in FASTA format."""
        str = ""
        str += '>' + self.hdr + '\n'
        if wrap_width:
            for line in Sequence._chunker(self.seq, wrap_width):
                str += line + '\n'
        else:
            str += self.seq + '\n'
        return(str)

    def copy(self):
        """Provide a deep copy of this sequence."""

        # TODO: make this more efficient (possibly?)
        return Sequence(self.hdr, self.seq)

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

    def slice_seq(self, start, end=None):
        """Return a slice of this sequence."""
        return self.seq[start:end]

    @staticmethod
    def _chunker(str, width):
        return (str[pos:pos + width] for pos in range(0, len(str), width))

    @staticmethod
    def is_gap(res_char):
        return res_char in ['-', '.']

    def __str__(self):
        return('{:<30} {}'.format(self.hdr, self.seq))


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
        """Create a new Correspondence object from a GCF io / filename / string"""

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
            raise SeqIOError("encountered an error trying to readline() on GCF io ({})".format(gcf_io))

        line_no = 1
        for line in gcf_io:
            line_no += 1

            try:
                seqres_aa, seqres_num, pdb_label, pdb_aa = line.split()
                if pdb_aa is not seqres_aa and pdb_aa is not Correspondence.GCF_GAP_CHAR:
                    logging.warning("pdb_aa '{}' does not match seqres_aa '{}' (line: {})".format(pdb_aa, seqres_aa, line_no))
            except:
                raise SeqIOError("Error: failed to parse GCF '{}' ({}:{})".format(
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
        res = self.residues[seq_num - 1]
        assert res.seq_num is seq_num
        return res

    def get_res_by_pdb_label(self, pdb_label):
        res = next((res for res in self.residues if res.pdb_label == pdb_label), None)
        return res

    def get_res_by_atom_pos(self, pos):
        """Return the Residue corresponding to the given position in the ATOM sequence (ignores gaps)."""
        atom_residues = [res for res in self.residues if res.pdb_label is not None]
        res = atom_residues[pos-1]
        return res

    @property
    def atom_sequence(self):
        id = "atom|{}".format(self.id)
        res = [res.aa if res.pdb_label else Correspondence.FASTA_GAP_CHAR for res in self.residues]
        return Sequence(id, "".join(res))

    @property
    def seqres_sequence(self):
        id = "seqres|{}".format(self.id)
        res = [res.aa for res in self.residues]
        return Sequence(id, "".join(res))

    def to_sequences(self):
        seqs = (self.seqres_sequence, self.atom_sequence)
        return seqs

    def to_fasta(self, **kwargs):
        seqs = self.to_sequences()
        return seqs[0].to_fasta(**kwargs) + seqs[1].to_fasta(**kwargs)
        
    def to_aln(self):
        seqs = self.to_sequences()
        return Alignment(seqs = seqs)
        
    def __str__(self):
        return self.to_fasta()

    def __repr__(self):
        return self.to_fasta()

class Alignment(object):
    """Object storing a protein sequence alignment"""

    REF_GAP_CHAR='-'
    MERGE_GAP_CHAR='.'

    def __init__(self, seqs=None):
        self.seqs = seqs if seqs else []
        self.__aln_positions = 0

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

    def find_seq_by_id(self, id):
        """Return the Sequence corresponding to the provided id."""
        seqs_with_id = [seq for seq in self.seqs if seq.id == id]
        if len(seqs_with_id) > 1:
            raise SeqIOError("Found more than one sequence matching id '{}'".format(id))
        if len(seqs_with_id) == 0:
            return None
        return seqs_with_id[0]

    def get_seq_at_offset(self, offset):
        """Return the Sequence at the given offset (zero-based)."""
        return self.seqs[offset]

    @classmethod
    def new_from_fasta(cls, fasta_io):
        """Initialise an alignment object from a fasta file / string / io"""
        aln = Alignment()
        aln.read_sequences_from_fasta(fasta_io)
        return aln
    
    def read_sequences_from_fasta(self, fasta_io):
        """Parse aligned sequences from FASTA (str, file, io) and adds them to the current
        Alignment object. Returns the number of sequences that are added."""

        if (type(fasta_io) == str):
            if (fasta_io[0] == '>'):
                fasta_io = io.StringIO(fasta_io)
            else:
                fasta_io = open(fasta_io)
 
        seq_added=0
        current_hdr = None
        current_seq = ''
        while True:
            try:
                line = fasta_io.readline()
            except AttributeError:
                raise SeqIOError("encountered an error trying to readline() on fasta io ({})".format(fasta_io))

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
                raise SeqIOError( ("Error: cannot add a sequence (id:{}) "
                    "with {} positions to an alignment with {} positions.")
                    .format(seq.id, seq.length(), self.aln_positions) )
        else:
            self.__aln_positions = seq.length()

        self.seqs.append(seq)
        return seq        

    def remove_alignment_gaps(self):
        """Return a new alignment after removing alignment positions 
        that contain a gap for all sequences."""
        seqs = self.seqs
        seq_length = seqs[0].length()
        new_seq_strings = ["" for s in range(len(seqs))]
        for aln_pos in range(seq_length):
            total_gaps = 0
            for seq in seqs:
                if seq.seq[aln_pos] == '-' or seq.seq[aln_pos] == '.':
                    total_gaps += 1
            if total_gaps < len(seqs):
                for seq_pos in range(len(seqs)):
                    res = seqs[seq_pos].seq[aln_pos]
                    # print( "seq[{}:{}] pos:{} res:{}".format(aln_pos, seqs[seq_pos].id, seq_pos, res) )
                    new_seq_strings[ seq_pos ] += res
            else:
                logging.info("Removing complete gap from alignment position: {}".format(aln_pos))

        new_aln = Alignment()
        for seq_pos in range(len(new_seq_strings)):
            hdr = seqs[seq_pos].hdr
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
        return [Sequence(s.hdr, s.slice_seq(start, end)) for s in self.seqs]

    def merge_alignment(self, merge_aln, ref_seq_id: str, ref_correspondence = None):
        """
        Merges sequences from an alignment into the current object.

        Sequences in ``merge_aln`` are brought into the current alignment using
        the equivalences identified in reference sequence ``ref_seq_id`` (which
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
            merge_aln (Alignment): An Alignment containing the reference
                sequence and any additional sequences to merge. 
            ref_seq_id (str): A str that will be used to find the reference 
                sequence in the current alignment and merge_aln 
            ref_correspondence (Correspondence): An optional Correspondence 
                object that provides a mapping between the reference 
                sequence found in ``self`` (ATOM records) and reference 
                sequence as it appears in ``merge_aln`` (SEQRES records).

        Returns: 
            [Sequence]: Array of Sequences added to the current alignment.

        Raises: 
            SeqIOMergeCorrespondenceError: problem mapping reference
                sequence between alignment and correspondence 

        """

        merge_aln = merge_aln.copy()

        ref_seq_in_ref = self.find_seq_by_id(ref_seq_id)
        ref_seq_in_merge = merge_aln.find_seq_by_id(ref_seq_id)

        if ref_correspondence is None:
            # fake a 1:1 correspondence for internal use
            residues = ref_seq_in_ref.get_residues()
            for r in residues:
                r.pdb_label = str(r.seq_num)
                logging.info("fake correspondence: residue={}".format(repr(r)))
            ref_correspondence = Correspondence(id=ref_seq_id, residues=residues)

        logging.debug("merge_alignment: Correspondence: {}".format(repr(ref_correspondence)))

        # check: ref sequence (in self) must match the ATOM sequence in Correspondence
        ref_no_gaps = ref_seq_in_ref.seq_no_gaps
        corr_no_gaps = ref_correspondence.atom_sequence.seq_no_gaps
        if ref_no_gaps != corr_no_gaps:
            raise SeqIOMergeCorrespondenceError(ref_seq_id, 'current', 'ATOM', 
                ref_no_gaps, corr_no_gaps) 

        # check: ref sequence (in merge) must match the SEQRES sequence in Correspondence
        ref_no_gaps = ref_seq_in_merge.seq_no_gaps
        corr_no_gaps = ref_correspondence.seqres_sequence.seq_no_gaps
        if ref_no_gaps != corr_no_gaps:
            raise SeqIOMergeCorrespondenceError(ref_seq_id, 'merge', 'SEQRES', 
                ref_no_gaps, corr_no_gaps) 

        # clean up
        del ref_no_gaps
        del corr_no_gaps

        ref_aln_pos=0
        ref_corr_pos=0
        merge_aln_pos=0
        correspondence_length = ref_correspondence.seqres_length

        logging.debug("ref_alignment.positions: {}".format(self.aln_positions))
        logging.debug("merge_alignment.positions: {}".format(merge_aln.aln_positions))
        logging.debug("ref_seq_in_ref:   {}".format(str(ref_seq_in_ref)))
        logging.debug("ref_seq_in_merge: {}".format(str(ref_seq_in_merge)))

        while True:

            if merge_aln_pos >= merge_aln.aln_positions \
                and ref_aln_pos  >= self.aln_positions \
                and ref_corr_pos >= correspondence_length:
                break

            logging.debug("REF {}/{}; CORRESPONDENCE {}/{}; MERGE {}/{}".format(
                ref_aln_pos, self.aln_positions, ref_corr_pos, 
                correspondence_length, merge_aln_pos, merge_aln.aln_positions))

            # sort the gaps in the reference alignment
            if ref_aln_pos < self.aln_positions:

                for seq in self.slice_seqs(0, ref_aln_pos):
                    logging.debug( "{:<10} {}".format("REF", str(seq)) )

                ref_res_in_ref = ref_seq_in_ref.get_res_at_offset(ref_aln_pos)

                logging.debug("REF_POSITION   {:>3} of {:>3} => '{}'".format(
                    ref_aln_pos, self.aln_positions, ref_res_in_ref))

                # insert all the gaps in the reference alignment into the merge sequences
                # keep doing this until we don't have any more gaps
                if Sequence.is_gap(ref_res_in_ref):
                    logging.debug(("GAP '{}' in ref sequence in REF alignment [{}], "
                        "inserting gap '{}' at position [{}] in all merge sequences").format(
                            ref_res_in_ref, ref_aln_pos, ref_res_in_ref, merge_aln_pos))
                    merge_aln.insert_gap_at_offset(merge_aln_pos, gap_char=ref_res_in_ref)

                    # this is a gap: do NOT increment ref_corr_pos
                    ref_aln_pos   += 1
                    merge_aln_pos += 1
                    continue

            # sort the gaps in the merge alignment
            if merge_aln_pos < merge_aln.aln_positions:

                for seq in merge_aln.slice_seqs(0, merge_aln_pos):
                    logging.debug( "{:<10} {}".format("MERGE", str(seq)) )

                ref_res_in_merge = ref_seq_in_merge.get_res_at_offset(merge_aln_pos)

                logging.debug("MERGE_POSITION {:>3} of {:>3} => '{}'".format(
                    ref_aln_pos, self.aln_positions, ref_res_in_ref))

                # insert all the gaps in the merge alignment into the ref sequences
                # keep doing this until we don't have any more gaps
                if Sequence.is_gap(ref_res_in_merge):
                    logging.debug(("GAP '{}' in ref sequence in MERGE alignment [{}], "
                        "inserting gap '{}' at position [{}] in all ref sequences").format(
                            ref_res_in_merge, merge_aln_pos, Alignment.MERGE_GAP_CHAR, merge_aln_pos))
                    self.insert_gap_at_offset(ref_aln_pos, gap_char=Alignment.MERGE_GAP_CHAR)
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
                    logging.debug( "{:<10} {}".format("CORR", str(seq)) )

                ref_res_in_corr = ref_correspondence.get_res_at_offset(ref_corr_pos)
                if ref_res_in_corr.pdb_label is None:

                    logging.debug(("GAP '{}' in ATOM records of correspondence [{}], "
                        "inserting gap '{}' at position [{}] in ref sequences").format(
                            '*', ref_corr_pos, Alignment.MERGE_GAP_CHAR, ref_aln_pos))

                    #merge_aln.insert_gap_at_offset(merge_aln_pos, gap_char=Alignment.MERGE_GAP_CHAR)
                    self.insert_gap_at_offset(ref_aln_pos, gap_char=Alignment.MERGE_GAP_CHAR)
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

        logging.info("FINISHED MERGE")
        # for seq in ref_correspondence.to_sequences():
        #     seq = seq.slice_seq(0, ref_corr_pos)
        #     logging.debug( "{:<10} {}".format("CORR", str(seq)) )
        # for seq in self.seqs:
        #     logging.debug( "{:<10} {}".format("REF", str(seq)) )
        # for seq in merge_aln.seqs:
        #     logging.debug( "{:<10} {}".format("MERGE", str(seq)) )

        # add the merged sequences into this alignment
        for seq in merge_aln.seqs:
            if seq.id == ref_seq_id:
                continue
            self.add_sequence(seq)

        for seq in self.seqs:
            logging.debug( "{:<10} {}".format("MERGED", str(seq)) )

        # test the final, merged alignment
        # 1. get sequences that correspond to the input aln
        # 2. remove alignment positions where there's a gap in the reference sequence

        for original_seq in merge_aln.seqs:
            seq = self.find_seq_by_id(original_seq.id)

            if not seq:
                raise SeqIOError("failed to find sequence with id '{}' in merge aln".format(seq.id))

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
                    ref_seq_pos_in_merge = ref_corr_res.seq_num

                    # find the aln offset for the equivalent position in the original merge alignment
                    ref_merge_offset = ref_seq_in_merge.get_offset_at_seq_position(ref_seq_pos_in_merge)

                    # logging.debug("ref_seq_pos (ref): {}, ref_seq_pos (merge): {}, correspondence_res: {}, ref_merge_offset: {}".format(
                    #     ref_seq_pos_in_ref, ref_seq_pos_in_merge, repr(ref_corr_res), ref_merge_offset
                    # ))

                    # find the residue at the equivalent position in the merge alignment
                    original_res = original_seq.get_res_at_offset(ref_merge_offset)

                    if merged_res_at_aln_offset != original_res:
                        raise SeqIOError(("Expected the merged residue '{}' to "
                            "match the original residue '{}' at alignment "
                            "offset {} (sequence: '{}')\n\n"
                            "REF_SEQ_IN_REF:   {}\n"
                            "REF_SEQ_IN_MERGE: {}\n"
                            "ORIGINAL_SEQ:     {}\n"
                            "                  {aln_pointer:>{merge_pos}}\n"
                            "MERGED_SEQ:       {}\n"
                            "                  {aln_pointer:>{aln_pos}}\n"
                            "(aln_offset={}, ref_seq_pos(ref)={}, ref_seq_pos(merge)={}, ref_merge_pos={})"
                            ).format(
                                merged_res_at_aln_offset, original_res, aln_offset, seq.id,
                                ref_seq_in_ref.seq,
                                ref_seq_in_merge.seq,
                                original_seq.seq,
                                seq.seq,
                                aln_offset, ref_seq_pos_in_ref, ref_seq_pos_in_merge, ref_merge_offset,
                                aln_pointer='^', aln_pos=(aln_offset+1), merge_pos=(ref_merge_offset+1)
                            ))

        return merge_aln.seqs


    def copy(self):
        """Return a deepcopy of this object."""
        new_aln = Alignment()
        new_seqs = [s.copy() for s in self.seqs]
        new_aln.seqs = new_seqs
        new_aln.aln_positions = new_aln.seqs[0].length()
        return new_aln

    def write_fasta(self, fasta_file, wrap_width=80):
        """Write the sequences to a file in FASTA format."""
        with open(fasta_file, 'w') as f:
            for seq in self.seqs:
                f.write( seq.to_fasta(wrap_width=wrap_width) )

    def __str__(self):
        return "\n".join( [ str(seq) for seq in self.seqs ] )
