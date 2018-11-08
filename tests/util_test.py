import glob
import logging
import os
import re
import tempfile

from . import testutils

from cathpy import util, seqio, error as err

logger = logging.getLogger(__name__)

class TestUtil(testutils.TestBase):

    def setUp(self):
        self.cath_version = 'v4.2'
        self.data_dir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), 'data')
        self.sc_file = os.path.join(self.data_dir, 
            '1.10.8.10__FF_SSG9__6.aln_reps.cora.fa')
        self.ff_dir = os.path.join(self.data_dir, 'funfams')
        self.ff_file = os.path.join(self.ff_dir, '1.10.8.10-ff-14534.reduced.sto')
        self.ff_tmpl = '__SFAM__-ff-__FF_NUM__.reduced.sto'
        self.merge_sto_file = os.path.join(self.data_dir, 'merge.sto')
        self.example_fasta_file = self.sc_file

    @testutils.log_title
    def test_merge(self):
        tmp_fasta_file = tempfile.NamedTemporaryFile(mode='w+', suffix='.fa', delete=True)
        tmp_sto_file = tempfile.NamedTemporaryFile(mode='w+', suffix='.sto', delete=False)
        logger.info("Creating SC merger...")
        merger = util.StructuralClusterMerger(cath_version=self.cath_version, 
            sc_file=self.sc_file, 
            out_sto=tmp_sto_file.name, out_fasta=tmp_fasta_file.name,
            ff_dir=self.ff_dir, ff_tmpl=self.ff_tmpl)
        logger.info("Merging SC alignment {}".format(self.sc_file))
        merge_aln = merger.run()
        self.assertEqual(merge_aln.count_sequences, 701)

        with open(tmp_sto_file.name) as f:
            sto_got = f.read()
        with open(self.merge_sto_file) as f:
            sto_expected = f.read()

        logger.info("Checking {} versus {}".format(tmp_sto_file.name, self.merge_sto_file))
        self.assertMultiLineEqual(sto_got, sto_expected)

    @testutils.log_title
    def test_cath_version(self):
        cv = util.CathVersion('4.2.0')
        self.assertIsInstance(cv, util.CathVersion)
        self.assertEqual(cv.major, '4')
        self.assertEqual(cv.minor, '2')
        self.assertEqual(cv.trace, '0')
        self.assertEqual(cv.dirname, 'v4_2_0')
        self.assertEqual(cv.join('-'), '4-2-0')
        self.assertIsInstance(util.CathVersion('v4.2'), util.CathVersion)
        self.assertEqual(str(util.CathVersion(4.2)), '4.2.0')

    def test_funfam_file_finder(self):
        finder = util.FunfamFileFinder(base_dir=self.ff_dir, ff_tmpl='__SFAM__-ff-__FF_NUM__.reduced.sto')
        self.assertIsInstance(finder, util.FunfamFileFinder)
        ff_file = finder.search_by_domain_id('2damA00')
        self.assertEqual(os.path.basename(ff_file), '1.10.8.10-ff-14534.reduced.sto')

        with self.assertRaises( err.NoMatchesError ):
            finder.search_by_domain_id('1zzzA01')
        with self.assertRaises( err.InvalidInputError ):
            finder.search_by_domain_id('bingo')
        with self.assertRaises( err.InvalidInputError ):
            finder.search_by_domain_id(' file with &*! characters and spaces ')

    def test_ff_id_from_file(self):
        finder = util.FunfamFileFinder(base_dir=self.ff_dir, ff_tmpl='__SFAM__-ff-__FF_NUM__.reduced.sto')
        ff_file = finder.search_by_domain_id('2damA00')
        ff_id = finder.funfam_id_from_file(ff_file)
        self.assertEqual(ff_id.sfam_id, '1.10.8.10')
        self.assertEqual(ff_id.cluster_num, 14534)

    @testutils.log_level('cathpy.util', 'DEBUG')
    def test_scorecons(self):
        sc = util.ScoreconsRunner()
        aln = seqio.Align.new_from_fasta(self.example_fasta_file)
        
        sc_res = sc.run_fasta(self.example_fasta_file)
        self.assertEqual(sc_res.dops, 92.889)
        self.assertEqual(len(sc_res.scores), aln.aln_positions)

    def test_groupsim(self):
        gs = util.GroupsimRunner()
        aln = seqio.Align.new_from_fasta(self.example_fasta_file)

        seqs = aln.seqs

        for s in seqs[:2]:
            s.set_cluster_id('0001')
        for s in seqs[2:]:
            s.set_cluster_id('0002')

        gs_res = gs.run_alignment(aln)
        self.assertEqual(gs_res.count_positions, aln.aln_positions)
        print("GS: {}".format(repr(gs_res.__dict__)))

        # aln.write_sto('tmp.sto')
        # aln.add_scorecons()
        # aln.write_sto('tmp.plus_scorecons.sto')
        # aln.add_groupsim()
        # aln.write_sto('tmp.plus_groupsim.sto')

    def test_cluster_file(self):
        sc_path = os.path.abspath(self.sc_file)
        sc_dir = os.path.dirname(sc_path)
        sc_file = util.ClusterFile(sc_path)
        # 1.10.8.10__FF_SSG9__6.aln_reps.cora.fa
        self.assertDictEqual(sc_file.__dict__, { 
            'dir': sc_dir,
            'sfam_id': '1.10.8.10',
            'cluster_type': 'FF_SSG9',
            'cluster_num': '6',
            'desc': '.aln_reps.cora',
            'suffix': '.fa',
            'join_char': '__',
        })
        self.assertEqual(sc_file.to_string(), sc_path)

        ff_path = os.path.abspath(self.ff_file)
        ff_dir = os.path.dirname(ff_path)
        ff_file = util.ClusterFile(ff_path)
        # 1.10.8.10-ff-14534.reduced.sto
        self.assertDictEqual(ff_file.__dict__, { 
            'dir': ff_dir,
            'sfam_id': '1.10.8.10',
            'cluster_type': 'ff',
            'cluster_num': '14534',
            'desc': '.reduced',
            'suffix': '.sto',
            'join_char': '-',
        })
        self.assertEqual(ff_file.to_string(), ff_path)
