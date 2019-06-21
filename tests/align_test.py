import logging
import os
import tempfile
import unittest

from cathpy.align import Segment, Sequence, Align

from . import testutils

LOG = logging.getLogger(__name__)


class TestSequence(testutils.TestBase):

    def test_create_sequence(self):
        seq = Sequence('id1/23-123', '---AKGHP--GPKAPGPAK--')
        self.assertEqual(seq.uid, 'id1/23-123')
        self.assertEqual(seq.accession, 'id1')
        self.assertEqual(len(seq.segs), 1)
        self.assertEqual(seq.segs[0].start, 23)
        self.assertEqual(seq.segs[0].stop, 123)
        self.assertEqual(seq.seq, '---AKGHP--GPKAPGPAK--')
        self.assertEqual(seq.get_res_at_offset(0), '-')
        self.assertEqual(seq.get_res_at_offset(3), 'A')
        seq.insert_gap_at_offset(5)
        self.assertEqual(seq.seq, '---AK-GHP--GPKAPGPAK--')
        seq.insert_gap_at_offset(-3, gap_char='.')
        self.assertEqual(seq.seq, '---AK-GHP--GPKAPGPA.K--')

    def test_sequence_methods(self):
        seq = Sequence('id1/23-123', '---AKGHP--GPKAPGPAK--')
        self.assertEqual(seq.get_offset_at_seq_position(1),
                         3)  # seq pos '1' 'A' -> offset '3'
        self.assertEqual(seq.get_res_at_seq_position(2), 'K')
        # offset 5 'G' -> seq pos 3 'AKG'
        self.assertEqual(seq.get_seq_position_at_offset(5), 3)
        self.assertEqual(seq.get_res_at_offset(5), 'G')

    def test_sequence_lower_case(self):
        seq = Sequence('id1/23-123', '---AKGHP--GPKAPGPAK--')
        seq.lower_case_at_offset(6)
        self.assertEqual(seq.seq, '---AKGhP--GPKAPGPAK--')

    def test_create_segment(self):
        seg = Segment(1, 10)
        self.assertEqual(seg.start, 1)
        self.assertEqual(seg.stop, 10)
        self.assertEqual(str(seg), '1-10')

    def test_split_hdr(self):
        hdr = Sequence.split_hdr('domain|1cukA01/12-134_178-234')
        self.assertEqual(hdr['id'], 'domain|1cukA01/12-134_178-234')
        self.assertEqual(hdr['accession'], '1cukA01')
        self.assertEqual(hdr['id_type'], 'domain')
        self.assertIsInstance(hdr['segs'][0], Segment)
        self.assertEqual(hdr['segs'][0].start, 12)
        self.assertEqual(hdr['segs'][0].stop, 134)
        self.assertIsInstance(hdr['segs'][1], Segment)
        self.assertEqual(str(hdr['segs'][1]), '178-234')
        self.assertEqual(hdr['id_ver'], None)

    def test_apply_segments(self):
        seq1str = 'AKGHP GPKAP GPAKK APHPP PAIIH PAPIL HADSA P'.replace(
            ' ', '')
        seq1 = Sequence('testid', seq1str)
        self.assertEqual(seq1.uid, 'testid')
        self.assertEqual(len(seq1.seq), len(seq1str))

        seq2 = seq1.apply_segments([Segment(3, 5)])
        self.assertEqual(seq1.uid, 'testid')
        self.assertEqual(len(seq1.seq), len(seq1str))
        self.assertEqual(seq2.uid, 'testid/3-5')
        self.assertEqual(seq2.seq, 'GHP')

        seq3 = seq1.apply_segments(
            [Segment(3, 10), Segment(15, 25)])

        seq3str = ''.join(['GHPGPKAP', 'KAPHPPPAIIH'])
        self.assertEqual(seq3.uid, 'testid/3-10_15-25')
        self.assertEqual(len(seq3.seq), 8+11)
        self.assertEqual(seq3.seq, seq3str)


class TestAlign(testutils.TestBase):

    def setUp(self):
        data_dir = os.path.join(os.path.dirname(__file__), 'data')
        self.stockholm_file = os.path.join(data_dir, 'test.sto')
        self.pfam_sto_file = os.path.join(data_dir, 'PF03770.sto')

        self.stockholm_gzip_file = self.stockholm_file + '.gz'

        self.fasta_file = tempfile.NamedTemporaryFile(mode='w+', delete=False)

        # removing gaps from alignments
        self.fasta_contents = '''
>id1
----AKGHP---GPKAPGPAK---
>id2
-AGPAK-HP--PGPKAPPPAK-G-
'''[1:]
        self.fasta_contents_without_gaps = '''
>id1
---AKGHP-GPKAPGPAK-
>id2
AGPAK-HPPGPKAPPPAKG
'''[1:]

        self.pir_aln_ref = '''
>P1;ref1
ref1 description
---AKGHP--GPKAPGPAK--*
>P1;ref2
ref2 description
CGCAKGH-PKA--APGP--GT*
'''[1:]

        # merging funfam alignments into structural alignments
        self.fasta_aln_ref = '''
>ref1
---AKGHP--GPKAPGPAK--
>ref2
CGCAKGH-PKA--APGP--GT
'''[1:]

        self.fasta_aln_merge1 = '''
>ref1
-AK-GHPGP--KAPG--PAK
>src1.1
GAGGG-PGPGKKAPGG--AK
>src1.2
PAGGCCPGP--KAPGGSAA-
'''[1:]

        self.fasta_aln_merge2 = '''
>ref2
--CGC--AKGHPK-AAPGPGT-
>src2.1
--CGC-PAKHPPKGAA-GPGPA
>src2.2
GHCHCFSAKHPPK-AAHGPGPA
'''[1:]

        self.fasta_aln_after_merge1 = '''
>ref1
---.AK.GHP--GP..KAPG..PAK--
>ref2
AKG.AK.GH-PKA-..-APG..P--GP
>src1.1
---gAGgG.P--GPgkKAPGg..AK--
>src1.2
---pAGgCcP--GP..KAPGGsAA.--
'''[1:]

        self.fasta_aln_after_merge2 = '''
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

        self.fasta_file.write(self.fasta_contents)
        self.fasta_file.seek(0)

    def tearDown(self):
        pass
        # os.remove(self.fasta_file.name)

    def test_read_fasta_filename(self):
        aln = Align.from_fasta(self.fasta_file.name)
        self.assertEqual(aln.count_sequences, 2)
        seqs = aln.seqs
        self.assertEqual(seqs[0].uid, 'id1')
        self.assertEqual(seqs[1].uid, 'id2')

    def test_read_fasta_fileio(self):
        self.fasta_file.seek(0)
        aln = Align.from_fasta(self.fasta_file)
        self.assertEqual(aln.count_sequences, 2)

    def test_read_stockholm_file(self):
        aln = Align.from_stockholm(self.stockholm_file)
        self.assertEqual(aln.count_sequences, 51)

    def test_read_stockholm_gzip_file(self):
        aln = Align.from_stockholm(self.stockholm_gzip_file)
        self.assertEqual(aln.count_sequences, 51)

    def test_read_fasta_str(self):
        aln = Align.from_fasta(self.fasta_contents)
        self.assertEqual(aln.count_sequences, 2)

    def test_remove_gaps(self):
        self.log_title('remove_gaps')
        self.fasta_file.seek(0)
        aln = Align.from_fasta(self.fasta_contents)
        self.assertEqual(aln.count_sequences, 2)
        new_aln = aln.remove_alignment_gaps()
        new_seqs = new_aln.seqs
        seqs_no_gap = "".join([s.to_fasta() for s in new_seqs])
        self.assertEqual(seqs_no_gap, self.fasta_contents_without_gaps)

    def test_copy_aln(self):
        self.log_title('copy_aln')
        aln_ref = Align.from_fasta(self.fasta_aln_ref)
        aln_copy = aln_ref.copy()
        self.assertNotEqual(aln_copy, aln_ref)
        self.assertEqual(str(aln_copy), str(aln_ref))

    def test_aln_add_gap(self):
        self.log_title('aln_add_gap')
        aln = Align.from_fasta(self.fasta_aln_ref)
        self.assertEqual(aln.seqs[0].seq, '---AKGHP--GPKAPGPAK--')
        self.assertEqual(aln.seqs[1].seq, 'CGCAKGH-PKA--APGP--GT')
        aln.insert_gap_at_offset(4)
        self.assertEqual(aln.seqs[0].seq, '---A-KGHP--GPKAPGPAK--')
        self.assertEqual(aln.seqs[1].seq, 'CGCA-KGH-PKA--APGP--GT')
        aln.insert_gap_at_offset(-3, gap_char='.')
        self.assertEqual(aln.seqs[0].seq, '---A-KGHP--GPKAPGPA.K--')
        self.assertEqual(aln.seqs[1].seq, 'CGCA-KGH-PKA--APGP-.-GT')

    def test_merge_aln(self):
        aln_ref = Align.from_fasta(self.fasta_aln_ref)
        self.assertEqual(aln_ref.count_sequences, 2)
        aln_merge1 = Align.from_fasta(self.fasta_aln_merge1)
        self.assertEqual(aln_merge1.count_sequences, 3)
        aln_merge2 = Align.from_fasta(self.fasta_aln_merge2)
        self.assertEqual(aln_merge2.count_sequences, 3)

        aln_ref.merge_alignment(aln_merge1, 'ref1')
        expected_aln_after_merge1 = Align.from_fasta(
            self.fasta_aln_after_merge1)
        self.assertEqual(expected_aln_after_merge1.count_sequences, 4)
        self.assertEqual([s.uid for s in aln_ref.seqs], [
            'ref1', 'ref2', 'src1.1', 'src1.2', ])

        aln_ref.merge_alignment(aln_merge2, 'ref2')
        expected_aln_after_merge2 = Align.from_fasta(
            self.fasta_aln_after_merge2)
        self.assertEqual(expected_aln_after_merge2.count_sequences, 6)
        self.assertEqual([s.uid for s in aln_ref.seqs], [
            'ref1', 'ref2', 'src1.1', 'src1.2', 'src2.1', 'src2.2', ])

        sto_tmp = tempfile.NamedTemporaryFile(
            mode='w+', delete=True, suffix='.sto')
        sto_out = sto_tmp.name

        aln_ref.add_groupsim()
        aln_ref.add_scorecons()
        aln_ref.write_sto(sto_out)

    def test_pir(self):
        aln_ref = Align.from_pir(self.pir_aln_ref)
        self.assertEqual(aln_ref.count_sequences, 2)
        pir_tmp = tempfile.NamedTemporaryFile(
            mode='w+', delete=True, suffix='.pir')
        aln_ref.write_pir(pir_tmp.name)

        with open(pir_tmp.name, 'r') as f:
            pir_expected = f.read()

        # self.assertMultiLineEqual(self.pir_aln_ref, pir_expected)
        self.assertEqual(self.pir_aln_ref, pir_expected)

    def test_write_sto(self):

        sto_tmp = tempfile.NamedTemporaryFile(mode='w+', delete=True)
        sto_out = sto_tmp.name

        aln = Align.from_stockholm(self.stockholm_file)

        self.assertEqual(aln.uid, '1.10.8.10/FF/14534')

        # make sure we have parsed the meta data okay
        first_seq = aln.get_seq_at_offset(0)
        seq_meta = first_seq.meta
        self.assertEqual(seq_meta['AC'], 'Q96CS3')
        self.assertEqual(seq_meta['OS'], 'Homo sapiens')
        self.assertEqual(seq_meta['DE'], 'FAS-associated factor 2')
        LOG.info('first seq: %s', repr(vars(first_seq)))

        LOG.info('Writing out tmp STOCKHOLM file to %s', sto_out)
        aln.write_sto(sto_out)

        sto_expected = ''
        with open(self.stockholm_file, 'r') as io:
            # some lines seem to have trailing spaces
            for line in io:
                sto_expected += line.strip() + '\n'

        with open(sto_out, 'r') as io:
            sto_got = io.read()

        self.maxDiff = None
        self.assertMultiLineEqual(sto_got, sto_expected)

    def test_parse_pfam_stockholm(self):
        aln = Align.from_stockholm(self.pfam_sto_file)
        self.assertEqual(aln.count_sequences, 107)

        # LOG.info("aln: %s", pprint.pformat(aln.__dict__))

        self.assertEqual(aln.cath_version, None)
        self.assertEqual(aln.dops_score, None)
        self.assertEqual(aln.accession, 'PF03770.16')
        self.assertRegex(aln.author, r'^Finn')

        with self.assertRaises(AttributeError):
            # make sure this is in meta, not object level attribute
            self.assertRegex(
                aln.source_seed, r'^Pfam-B')  # pylint: disable=no-member

        self.assertRegex(aln.meta['source_seed'], r'^Pfam-B')
        self.assertRegex(aln.meta['search_method'], r'^hmmsearch')

        self.assertEqual(list(aln.seq_meta.keys()), ['SS_cons', 'seq_cons'])

    def test_meta_summary(self):
        aln = Align.from_stockholm(self.stockholm_file)
        meta_info = aln.get_meta_summary()

        self.assertEqual(meta_info.ec_term_counts, None)
        self.assertEqual(meta_info.go_term_counts, {'GO:0005515': 2, 'GO:0005576': 2, 'GO:0005783': 3, 'GO:0005811': 3, 'GO:0030433': 4, 'GO:0030970': 3, 'GO:0031625': 3,
                                                    'GO:0034098': 3, 'GO:0034389': 3, 'GO:0035473': 3, 'GO:0035578': 2, 'GO:0043130': 3, 'GO:0043312': 2, 'GO:0055102': 3})
        self.assertEqual(meta_info.cath_domain_count, 1)
        self.assertEqual(meta_info.seq_count, 51)
        self.assertEqual(meta_info.dops_score, 61.529)
        self.assertEqual(meta_info.organism_newick, "((((((((((((((Homo,(Gorilla_gorilla)'Gorilla Gorilla gorilla (1)',Pan)'Homininae Homo (5)',(Pongo)'Ponginae Pongo (1)')'Hominidae Homininae (6)')'Hominoidea Hominidae (6)',(((Chlorocebus,Macaca,Papio)'Cercopithecinae Chlorocebus (3)')'Cercopithecidae Cercopithecinae (3)')'Cercopithecoidea Cercopithecidae (3)')'Catarrhini Hominoidea (9)')'Simiiformes Catarrhini (9)')'Haplorrhini Simiiformes (9)')'Primates Haplorrhini (9)',(((((Mus)'Mus Mus (1)',Rattus)'Murinae Mus (2)')'Muridae Murinae (2)')'Myomorpha Muridae (2)',((Heterocephalus,Fukomys)'Bathyergidae Heterocephalus (2)')'Hystricomorpha Bathyergidae (2)')'Rodentia Myomorpha (4)',((Oryctolagus)'Leporidae Oryctolagus (1)')'Lagomorpha Leporidae (1)')'Euarchontoglires Primates (14)',(((((Bos)'Bovinae Bos (1)')'Bovidae Bovinae (1)')'Pecora Bovidae (1)')'Ruminantia Pecora (1)',((Camelus)'Camelidae Camelus (1)')'Tylopoda Camelidae (1)',((((Pteropus)'Pteropodinae Pteropus (1)')'Pteropodidae Pteropodinae (1)')'Megachiroptera Pteropodidae (1)',((Myotis)'Vespertilionidae Myotis (1)')'Microchiroptera Vespertilionidae (1)')'Chiroptera Megachiroptera (2)',((((Felis)'Felinae Felis (1)')'Felidae Felinae (1)')'Feliformia Felidae (1)',(((Canis_lupus)'Canis Canis lupus (1)')'Canidae Canis (1)',(((Mustela_putorius)'Mustela Mustela putorius (1)')'Mustelinae Mustela (1)')'Mustelidae Mustelinae (1)')'Caniformia Canidae (2)')'Carnivora Feliformia (3)',(((Equus)'Equus Equus (2)')'Equidae Equus (2)')'Perissodactyla Equidae (2)')'Laurasiatheria Ruminantia (9)',(((Loxodonta)'Elephantidae Loxodonta (1)')'Proboscidea Elephantidae (1)')'Afrotheria Proboscidea (1)')'Mammalia Euarchontoglires (24)',(((((((Silurana)'Xenopus Silurana (4)')'Xenopodinae Xenopus (4)',(Hymenochirus)'Pipinae Hymenochirus (1)')'Pipidae Xenopodinae (5)')'Pipoidea Pipidae (5)')'Anura Pipoidea (5)')'Batrachia Anura (5)')'Amphibia Batrachia (5)',((((((Ophiophagus)'Elapinae Ophiophagus (1)')'Elapidae Elapinae (1)')'Colubroidea Elapidae (1)')'Serpentes Colubroidea (1)')'Squamata Serpentes (1)')'Lepidosauria Squamata (1)',(((((((((Poeciliopsis,Xiphophorus,Poecilia)'Poeciliinae Poeciliopsis (4)')'Poeciliidae Poeciliinae (4)',(Fundulus)'Fundulidae Fundulus (1)')'Cyprinodontoidei Poeciliidae (5)')'Cyprinodontiformes Cyprinodontoidei (5)',((((Oryzias)'Oryziinae Oryzias (1)')'Adrianichthyidae Oryziinae (1)')'Adrianichthyoidei Adrianichthyidae (1)')'Beloniformes Adrianichthyoidei (1)')'Atherinomorphae Cyprinodontiformes (6)',((((Astyanax)'Characidae Astyanax (1)')'Characoidei Characidae (1)')'Characiformes Characoidei (1)')'Characiphysae Characiformes (1)',((((Gasterosteus)'Gasterosteidae Gasterosteus (1)')'Gasterosteales Gasterosteidae (1)')'Cottioidei Gasterosteales (1)')'Perciformes Cottioidei (1)',((((Takifugu)'Tetraodontidae Takifugu (1)')'Tetradontoidea Tetraodontidae (1)')'Tetraodontoidei Tetradontoidea (1)')'Tetraodontiformes Tetraodontoidei (1)',((((Danio)'Cyprinidae Danio (2)')'Cyprinoidea Cyprinidae (2)')'Cypriniformes Cyprinoidea (2)')'Cypriniphysae Cypriniformes (2)',(((((Oreochromis)'Oreochromini Oreochromis (1)')'Pseudocrenilabrinae Oreochromini (1)')'Cichlidae Pseudocrenilabrinae (1)')'Cichliformes Cichlidae (1)')'Cichlomorphae Cichliformes (1)',(((Oncorhynchus)'Salmoninae Oncorhynchus (1)')'Salmonidae Salmoninae (1)')'Salmoniformes Salmonidae (1)')'Teleostei Atherinomorphae (13)',(((Lepisosteus)'Lepisosteidae Lepisosteus (1)')'Semionotiformes Lepisosteidae (1)')'Holostei Semionotiformes (1)')'Neopterygii Teleostei (14)')'Actinopteri Neopterygii (14)')'Actinopterygii Actinopteri (14)',(((((Gallus)'Phasianinae Gallus (1)')'Phasianidae Phasianinae (1)')'Galliformes Phasianidae (1)',((Picoides)'Picidae Picoides (1)')'Piciformes Picidae (1)',((((Taeniopygia)'Estrildinae Taeniopygia (1)')'Estrildidae Estrildinae (1)')'Passeroidea Estrildidae (1)',(Ficedula)'Muscicapidae Ficedula (1)')'Passeriformes Passeroidea (2)')'Neognathae Galliformes (4)')'Aves Neognathae (4)',((((Callorhinchus)'Callorhinchidae Callorhinchus (1)')'Chimaeriformes Callorhinchidae (1)')'Holocephali Chimaeriformes (1)')'Chondrichthyes Holocephali (1)',((Latimeria)'Coelacanthidae Latimeria (1)')'Coelacanthiformes Coelacanthidae (1)',(((Alligator)'Alligatorinae Alligator (1)')'Alligatoridae Alligatorinae (1)')'Crocodylia Alligatoridae (1)')'Craniata Mammalia (51)')'Chordata Craniata (51)')'Metazoa Chordata (51)')'Eukaryota Metazoa (51)')'ROOT (51)';")


if __name__ == '__main__':
    unittest.main()
