import re
import logging
import io

class SeqIOException(Exception):
    def __init__(self, str):
        super().__init__(str)

class SeqIOGapException(SeqIOException):
    def __init__(self, str):
        super().__init__(str)

class AminoAcid:
    """Class representing an Amino Acid."""

    def __init__(self, one, three, word):
        self.one = one
        self.three = three
        self.word = word

    def __str__(self):
        return self.one

class AminoAcids:

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
    def get_by_id(cls, aa):
        return cls.aa_by_id[aa.upper()]


class Residue:
    """Class to represent a protein residue."""

    def __init__(self, aa, seq_num, pdb_label):
        self.aa = aa
        self.seq_num = seq_num
        self.pdb_label = pdb_label
    
    def __str__(self):
        return "res: aa={}, seq_num={}, pdb_label={}".format(
            self.aa, self.seq_num, self.pdb_label)

class Segment:
    """Class to represent a protein segment."""

    def __init__(self, start: int, stop: int):
        self.start = start
        self.stop  = stop
    
    def __str__(self):
        return("seg: {}-{}".format(self.start, self.stop))

class Sequence:
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
    
    def get_res_at_offset(self, offset):
        """Return the residue character at the given sequence position (includes gaps)."""
        try:
            res = self.seq[offset]
        except:
            raise SeqIOException("Error: failed to get residue at offset {} from sequence with length {}: '{}'".format(
                offset, self.length(), self.seq))

        return res

    def get_res_at_seq_position(self, seq_pos):
        """Return the residue character at the given sequence position (ignores gaps)."""
        seq_nogap = re.sub(Sequence.re_gap_chars, '', self.seq)
        try:
            res = seq_nogap[seq_pos-1]
        except:
            raise SeqIOException("Error: failed to get residue at position {} from sequence with {} non-gap sequence positions: '{}'".format(
                offset, len(seq_nogap), self.seq))

        return res

    def get_seq_position_at_offset(self, offset):
        """Return the sequence position (ignores gaps) of the given residue offset (may include gaps)."""
        seq_to_offset = self.seq[:offset+1]
        if re.match(seq_to_offset[-1], Sequence.re_gap_chars):
            raise SeqIOGapException("Cannot get sequence position at offset {} since this corresponds to a gap".format(offset))

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

        raise SeqIOException("failed to find offset at sequence position {}".format(seq_pos))

    def length(self):
        """Return the length of the sequence."""
        return len(self.seq)
    
    @property
    def seq(self):
        """Return the amino acid sequence as a string."""
        return self._seq

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
                    logging.warning( "failed to parse feature from string {}".format(features_str) )

        return({'id': id, 'id_type': id_type, 'id_ver': id_ver, 
            'segs': segs, 'features': features})

    def to_fasta(self, wrap_width=80):
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
        self._set_sequence( new_seq )


    def slice_seq(self, start, end=None):
        """Return a slice of this sequence."""
        return self.seq[start:end]

    @staticmethod
    def _chunker(str, width):
        return (str[pos:pos + width] for pos in range(0, len(str), width))

    @staticmethod
    def is_gap(res_char):
        return res_char in [ '-', '.' ]

    def __str__(self):
        return('{:<30} {}'.format(self.hdr, self.seq))

class Alignment:
    """Object storing a protein sequence alignment"""

    REF_GAP_CHAR='-'
    MERGE_GAP_CHAR='.'

    def __init__(self):
        self.seqs = []
        self.__aln_positions = 0

    @property
    def aln_positions(self):
        return self.__aln_positions
    
    @aln_positions.setter
    def aln_positions(self, value):
        self.__aln_positions = value

    @property
    def count_sequences(self):
        return len(self.seqs)

    def find_seq_by_id(self, id):
        seqs_with_id = [seq for seq in self.seqs if seq.id == id]
        if len(seqs_with_id) > 1:
            raise SeqIOException("Found more than one sequence matching id '{}'".format(id))
        if len(seqs_with_id) == 0:
            return None
        return seqs_with_id[0]

    def get_seq_at_offset(self, offset):
        return self.seqs[offset]

    @classmethod
    def new_from_fasta(cls, fasta_io):
        """Initialise an alignment object from a fasta file / string / io"""
        aln = Alignment()
        aln.read_sequences_from_fasta(fasta_io)
        return aln
    
    def read_sequences_from_fasta(self, fasta_io):
        """Parse sequences from FASTA, return a list of Sequence objects."""

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
                raise SeqIOException("encountered an error trying to readline() on fasta io ({})".format(fasta_io))

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

    def add_sequence(self, seq):
        """Add a sequence to this alignment."""

        if self.aln_positions:
            if self.aln_positions != seq.length():
                raise SeqIOException( "Error: cannot add a sequence (id:{}) with {} positions to an alignment with {} positions.".format(
                    seq.id, seq.length(), self.aln_positions)
                )
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
                logging.info( "Removing complete gap from alignment position: {}".format(aln_pos) )

        new_aln = Alignment()
        for seq_pos in range(len(new_seq_strings)):
            hdr = seqs[seq_pos].hdr
            seq_str = new_seq_strings[seq_pos]
            seq = Sequence(hdr, seq_str)
            new_aln.add_sequence(seq)

        return new_aln

    def insert_gap_at_offset(self, offset, gap_char='-'):
        self.__aln_positions += 1
        for s in self.seqs:
            s.insert_gap_at_offset(offset, gap_char) 

    def set_gap_char_at_offset(self, offset, gap_char):
        for s in self.seqs:
            s.set_gap_char_at_offset(offset, gap_char)

    def lower_case_at_offset(self, start, end=None):
        """Make the residues in the given alignment window lower case."""
        for s in self.seqs:
            s.lower_case_at_offset(start, end)

    def slice_seqs(self, start, end=None):
        """Return an array of Sequence objects from start to end."""
        return [Sequence(s.hdr, s.slice_seq(start, end)) for s in self.seqs]

    def merge_alignment(self, merge_aln, ref_seq_id, ref_seq_alignment=None):
        """Merge the given alignment into the current alignment, using ref_seq 
        (reference sequence) to provide the equivalences used for the mapping.
        It is assumed that ref_seq can be found in both the current alignment and
        the reference alignment. If the reference sequence does not have a 1:1
        mapping between the residues in the reference alignment and the merge 
        alignment (eg ATOM vs SEQRES records) then this alignment can be provided
        in ref_seq_alignment."""
        
        # work through the ref sequence in the source alignment (self)
        #  - add gaps in the merge sequences whenever there's a gap in the reference sequence
        #  - 

        merge_aln = merge_aln.copy()

        ref_seq_in_ref = self.find_seq_by_id(ref_seq_id)
        ref_seq_in_merge = merge_aln.find_seq_by_id(ref_seq_id)
        ref_aln_pos=0
        merge_aln_pos=0
        logging.debug("ref_alignment.positions: {}".format(self.aln_positions))
        logging.debug("merge_alignment.positions: {}".format(merge_aln.aln_positions))
        logging.debug("ref_seq_in_ref:   {}".format(str(ref_seq_in_ref)))
        logging.debug("ref_seq_in_merge: {}".format(str(ref_seq_in_merge)))

        while True:

            if merge_aln_pos == merge_aln.aln_positions and ref_aln_pos == self.aln_positions:
                break

            logging.debug("REF {}/{}; MERGE {}/{}".format(
                ref_aln_pos, self.aln_positions, merge_aln_pos, merge_aln.aln_positions))

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

                    ref_aln_pos   += 1
                    merge_aln_pos += 1
                    continue

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

                    ref_aln_pos   += 1
                    merge_aln_pos += 1
                    continue

            ref_aln_pos   += 1
            merge_aln_pos += 1

        logging.info("FINISHED MERGE")
        for seq in self.seqs:
            logging.debug( "{:<10} {}".format("REF", str(seq)) )
        for seq in merge_aln.seqs:
            logging.debug( "{:<10} {}".format("MERGE", str(seq)) )

        # add the merged sequences into this alignment
        for seq in merge_aln.seqs:
            if seq.id == ref_seq_id:
                continue
            self.add_sequence(seq)

        # test the final, merged alignment
        # 1. get sequences that correspond to the input aln
        # 2. remove alignment positions where there's a gap in the reference sequence

        for original_seq in merge_aln.seqs:
            seq = self.find_seq_by_id(original_seq.id)

            if not seq:
                raise SeqIOException("failed to find sequence with id '{}' in merge aln".format(seq.id))

            for aln_offset in range(self.aln_positions):

                ref_res = ref_seq_in_ref.get_res_at_offset(aln_offset)

                merged_res_at_aln_offset = seq.get_res_at_offset(aln_offset)

                if ref_res == '.':
                    # everything else should be a '.' or a lowercase residue
                    assert( merged_res_at_aln_offset == '.' or re.match(r'[a-z]', merged_res_at_aln_offset) )
                elif ref_res == '-':
                    # everything else should be a '-' or an uppercase residue
                    assert( merged_res_at_aln_offset == '-' or re.match(r'[A-Z]', merged_res_at_aln_offset) )
                else:
                    # find the sequence offset of this aln position in the ref sequence
                    ref_seq_pos = ref_seq_in_ref.get_seq_position_at_offset(aln_offset)

                    # find the aln offset for the equivalent position in the original merge alignment
                    ref_merge_offset = ref_seq_in_merge.get_offset_at_seq_position(ref_seq_pos)

                    # find the residue at the equivalent position in the merge alignment
                    original_res = original_seq.get_res_at_offset(ref_merge_offset)

                    if merged_res_at_aln_offset != original_res:
                        raise SeqIOException(("Error: expected the merged residue '{}' to match the original residue '{}' at alignment offset {} (sequence: '{}')\n" +
                            "REF_SEQ_IN_REF:   {}\n" +
                            "REF_SEQ_IN_MERGE: {}\n" +
                            "ORIGINAL_SEQ:     {}\n" +
                            "MERGED_SEQ:       {}\n" +
                            "aln_offset: {}, ref_seq_pos: {}, ref_merge_offset: {}").format(
                            merged_res_at_aln_offset, original_res, aln_offset, seq.id,
                            str(ref_seq_in_ref),
                            str(ref_seq_in_merge),
                            str(original_seq),
                            str(seq),
                            aln_offset, ref_seq_pos, ref_merge_offset,
                        ))


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
