import logging
import os
import re
import tempfile
import unittest

from cathpy import seqio

from . import testutils

LOG = logging.getLogger(__name__)

class TestMergeAlignment(testutils.TestBase):

    def setUp(self):

        # merging funfam alignments into structural alignments (with atom/seqres mismatch)

        # note: ref1 is missing some residues that appear in the 'sequence' alignment
        self.aln_structure = '''
>ref1
---AKGHP--GPKA---AK--
>ref2
CGCAKGH-PKA--APGP--GT
'''[1:]

        # this is the alignment between ref1 as seen in the structural alignment (based on ATOM records)
        # and ref1 as seen in the sequence alignment being merged in (based on SEQRES records)
        self.gcf_ref1 = '''
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
P  12   *   *
A  13  12   A
K  14  13   K
'''[1:]

        # the equivalent alignment in fasta format
        self.gcf_ref1_as_fasta = '''
>seqres|ref1
AKGHPGPKAPGPAK
>atom|ref1
AKGHPGPKA---AK
'''[1:]

        # note: this does have the extra residues (from SEQRES records)
        self.aln_merge1 = '''
>ref1
-AK-GHPGP--KAPG--PAK
>src1.1
GAGGG-PGPGKKAPGG--AK
>src1.2
PAGGCCPGP--KAPGGSAA-
'''[1:]

        self.aln_merge2 = '''
>ref2
--CGC--AKGHPK-AAPGPGT-
>src2.1
--CGC-PAKHPPKGAA-GPGPA
>src2.2
GHCHCFSAKHPPK-AAHGPGPA
'''[1:]

        self.aln_after_merge1 = '''
>ref1
---.AK.GHP--GP..KA---.....AK--
>ref2
CGC.AK.GH-PKA-..-APGP.....--GT
>ref1_merge
---.AK.GHP--GP..KA---pg..pAK--
>src1.1
---gAGgG-P--GPgkKA---pgg..AK--
>src1.2
---pAGgCCP--GP..KA---pggsaA---
'''[1:]

        self.aln_after_merge2 = '''
>ref1
..---...AK.GHP--.GP..KA---.....AK--.
>ref2
..CGC...AK.GH-PK.A-..-APGP.....--GT.
>ref1_merge
..---...AK.GHP--.GP..KA---pg..pAK--.
>src1.1
..---g..AGgG-P--.GPgkKA---pgg..AK--.
>src1.2
..---p..AGgCCP--.GP..KA---pggsaA---.
>src2.1
..CGC..pAK.HP-PKgA-..-A-GP.....--GPa
>src2.2
ghCHC.fsAK.HP-PK.A-..-AHGP.....--GPa
'''[1:]
    
    def tearDown(self):
        pass

    @testutils.log_title
    def test_parse_gcf(self):
        gcf = seqio.Correspondence.new_from_gcf(self.gcf_ref1)
        self.assertEqual(gcf.seqres_length, 14)
        self.assertEqual(gcf.atom_length, 11)
        self.assertEqual(gcf.id, 'ref1')
        gcf_as_fasta = gcf.to_fasta()
        self.assertEqual(gcf_as_fasta, self.gcf_ref1_as_fasta)

    @testutils.log_title
#     @testutils.log_level('cathpy.seqio', 'DEBUG')
    def test_correspondence_apply_segments(self):

        corr_full = seqio.Correspondence.new_from_gcf(self.gcf_ref1)

        def _res(r):
            return (r.seq_num, r.pdb_label, r.aa)

        segs = [seqio.Segment(start=3, stop=11)]

        # [3-11]
        # 
        # G   3   7   G
        # H   4   8   H
        # P   5   9   P
        # G   6  10   G
        # P   7  10A  P
        # K   8  10B  K
        # A   9  11   A
        # P  10   *   *
        # G  11   *   * 

        corr = corr_full.apply_seqres_segments(segs)
        self.assertIsInstance(corr, seqio.Correspondence)
        self.assertEqual(corr.seqres_length, 9)
        self.assertEqual(corr.atom_length, 7)
        self.assertEqual(corr.seqres_sequence.seq, 'GHPGPKAPG')
        self.assertEqual(corr.atom_sequence.seq,   'GHPGPKA--')
        res = corr.get_res_at_offset(0)
        self.assertEqual(_res(res), (3, '7', 'G'))
        res = corr.get_res_at_offset(-1)
        self.assertEqual(_res(res), (11, None, 'G'))
        del segs, corr, res

        segs = [seqio.Segment(start=3, stop=7), seqio.Segment(start=9, stop=13)]

        # [3-7, 9-13]
        # 
        # G   3   7   G
        # H   4   8   H
        # P   5   9   P
        # G   6  10   G
        # P   7  10A  P
        # 
        # A   9  11   A
        # P  10   *   *
        # G  11   *   * 
        # P  12   *   *
        # A  13  12   A

        corr = corr_full.apply_seqres_segments(segs)
        self.assertIsInstance(corr, seqio.Correspondence)
        self.assertEqual(corr.seqres_length, 10)
        self.assertEqual(corr.atom_length, 7)
        self.assertEqual(corr.seqres_sequence.seq, 'GHPGPAPGPA')
        self.assertEqual(corr.atom_sequence.seq,    'GHPGPA---A')
        res = corr.get_res_at_offset(0)
        self.assertEqual(_res(res), (3, '7', 'G'))
        res = corr.get_res_at_offset(-1)
        self.assertEqual(_res(res), (13, '12', 'A'))
        del segs, corr, res

    @testutils.log_title
#     @testutils.log_level('cathpy.seqio', 'DEBUG')
    def test_merge_aln_with_correspondence(self):
        aln_ref = seqio.Align.new_from_fasta(self.aln_structure)
        self.assertEqual(aln_ref.count_sequences, 2)
        aln_merge1 = seqio.Align.new_from_fasta(self.aln_merge1)
        self.assertEqual(aln_merge1.count_sequences, 3)
        aln_merge2 = seqio.Align.new_from_fasta(self.aln_merge2)
        self.assertEqual(aln_merge2.count_sequences, 3)

        gcf = seqio.Correspondence.new_from_gcf(self.gcf_ref1)

        aln_ref.merge_alignment(aln_merge1, 'ref1', gcf)
        aln_after_merge1 = seqio.Align.new_from_fasta(self.aln_after_merge1)
        self.assertIn('ref1_merge', [s.id for s in aln_ref.seqs])
        #LOG.info("aln_after_merge1:\n%s", aln_ref.to_fasta())
        self.assertEqual(aln_ref.to_fasta(), aln_after_merge1.to_fasta())

        aln_ref.merge_alignment(aln_merge2, 'ref2')
        aln_after_merge2 = seqio.Align.new_from_fasta(self.aln_after_merge2)
        #LOG.info("aln_after_merge2:\n%s", aln_ref.to_fasta())
        self.assertEqual(aln_ref.to_fasta(), aln_after_merge2.to_fasta())

if __name__ == '__main__':
    unittest.main()
