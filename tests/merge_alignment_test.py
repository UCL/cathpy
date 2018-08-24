import logging
import os
import re
import tempfile
import unittest

from cathpy import seqio

logging.basicConfig(level=logging.DEBUG)

class TestMergeAlignment(unittest.TestCase):

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
---.AK.GHP--GP..KAPG..PAK--
>ref2
AKG.AK.GH-PKA-..-APG..P--GP
>src1.1
---gAGgG.P--GPgkKAPGg..AK--
>src1.2
---pAGgCcP--GP..KAPGGsAA.--
'''[1:]

        self.aln_after_merge2 = '''
>ref1  
..---...AK.GHP--.GP..KAPG..PAK--.
>ref2  
..CGC...AK.GH-PK.A-..-APG..P--GT.
>src1.1
..---g..AGgG-P--.GPgkKAPGg--AK--.
>src1.2
..---p..AGgCCP--.GP--KAPGgsAA---.
>src2.1
--CGC--pAK-HP-PKgA----A-G--P--GPa
>src2.2
ghCHC-fsAK-HP-PK-A----AHG--P--GPa
'''[1:]
    
    def tearDown(self):
        pass

    def test_parse_gcf(self):
        self.log_title("test_parse_gcf")
        gcf = seqio.Correspondence.new_from_gcf(self.gcf_ref1)
        self.assertEqual(gcf.seqres_length, 14)
        self.assertEqual(gcf.atom_length, 11)
        self.assertEqual(gcf.id, 'ref1')
        gcf_as_fasta = gcf.to_fasta()
        self.assertEqual(gcf_as_fasta, self.gcf_ref1_as_fasta)

    def test_merge_aln_with_correspondence(self):
        self.log_title("test_merge_aln")
        aln_ref = seqio.Alignment.new_from_fasta(self.aln_structure)
        self.assertEqual(aln_ref.count_sequences, 2)
        aln_merge1 = seqio.Alignment.new_from_fasta(self.aln_merge1)
        self.assertEqual(aln_merge1.count_sequences, 3)
        aln_merge2 = seqio.Alignment.new_from_fasta(self.aln_merge2)
        self.assertEqual(aln_merge2.count_sequences, 3)

        gcf = seqio.Correspondence.new_from_gcf(self.gcf_ref1)

        aln_ref.merge_alignment(aln_merge1, 'ref1', gcf)
        aln_after_merge1 = seqio.Alignment.new_from_fasta(self.aln_after_merge1)
        self.assertEqual(aln_after_merge1.count_sequences, 4)
        self.assertEqual(aln_ref.count_sequences, 4)

        aln_ref.merge_alignment(aln_merge2, 'ref2')
        aln_after_merge2 = seqio.Alignment.new_from_fasta(self.aln_after_merge2)
        self.assertEqual(aln_after_merge2.count_sequences, 6)
        self.assertEqual(aln_ref.count_sequences, 6)

    def log_title(self, title):
        hr = "=" * 80
        logging.info("")
        logging.info(hr)
        logging.info(" {} ".format(title))
        logging.info(hr)
        logging.info("")

if __name__ == '__main__':
    unittest.main()