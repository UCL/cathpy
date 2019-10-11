import logging
import difflib
import os
import tempfile

from cathpy import error as err
from cathpy.util import (StructuralClusterMerger,
                         AlignmentSummaryRunner, GroupsimRunner, GroupsimResult, )
from cathpy.datafiles import ReleaseDir
from cathpy.align import Align
from cathpy import util

from .testutils import TestBase, log_title, log_level

LOG = logging.getLogger(__name__)


class TestUtil(TestBase):

    def setUp(self):
        self.cath_version = 'v4.2'
        self.data_dir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), 'data')
        self.sc_file = os.path.join(self.data_dir,
                                    '1.10.8.10__FF_SSG9__6.aln_reps.cora.fa')
        self.ff_dir = os.path.join(self.data_dir, 'funfams')
        self.ff_file = os.path.join(
            self.ff_dir, '1.10.8.10-ff-14534.reduced.sto')
        self.ff_tmpl = '__SFAM__-ff-__FF_NUM__.reduced.sto'
        self.merge_sto_file = os.path.join(self.data_dir, 'merge.sto')
        self.example_fasta_file = self.sc_file
        self.example_sto_file = self.ff_file
        self.cath_release = ReleaseDir(
            self.cath_version, base_dir=self.data_dir)

    def test_alignment_summary_file(self):

        runner = AlignmentSummaryRunner(
            aln_file=self.merge_sto_file)
        entries = runner.run()
        self.assertEqual(len(entries), 1)
        summary = entries[0]
        self.assertEqual(summary.aln_length, 92)
        self.assertEqual(summary.dops, 88.422)
        self.assertEqual(summary.gap_count, 25228)
        self.assertEqual(summary.total_positions, 64492)
        self.assertEqual(summary.seq_count, 701)
        self.assertEqual(round(summary.gap_per, 2), round(39.12, 2))

    def test_alignment_summary_dir(self):

        runner = AlignmentSummaryRunner(
            aln_dir=self.data_dir, suffix='.sto')
        entries = runner.run()
        self.assertEqual(len(entries), 3)

        runner = AlignmentSummaryRunner(
            aln_dir=self.data_dir, suffix='.sto', recursive=True)
        entries = runner.run()
        self.assertEqual(len(entries), 7)

    @log_title
    def test_merge(self):
        tmp_fasta_file = tempfile.NamedTemporaryFile(
            mode='w+', suffix='.fa', delete=True)
        tmp_sto_file = tempfile.NamedTemporaryFile(
            mode='w+', suffix='.sto', delete=False)
        LOG.info("Creating SC merger...")

        merger = StructuralClusterMerger(cath_version=self.cath_version,
                                         sc_file=self.sc_file,
                                         out_sto=tmp_sto_file.name,
                                         out_fasta=tmp_fasta_file.name,
                                         ff_dir=self.ff_dir,
                                         ff_tmpl=self.ff_tmpl,
                                         cath_release=self.cath_release)

        LOG.info("Merging SC alignment {}".format(self.sc_file))
        merge_aln = merger.run()
        self.assertEqual(merge_aln.count_sequences, 701)

        with open(tmp_sto_file.name) as f:
            sto_got = f.read()
        with open(self.merge_sto_file) as f:
            sto_expected = f.read()

        LOG.info("Checking {} versus {}".format(
            tmp_sto_file.name, self.merge_sto_file))
        self.assertMultiLineEqual(sto_got, sto_expected)

    def test_funfam_file_finder(self):
        finder = util.FunfamFileFinder(
            base_dir=self.ff_dir, ff_tmpl='__SFAM__-ff-__FF_NUM__.reduced.sto')
        self.assertIsInstance(finder, util.FunfamFileFinder)
        ff_file = finder.search_by_domain_id('2damA00')
        self.assertEqual(os.path.basename(ff_file),
                         '1.10.8.10-ff-14534.reduced.sto')

        with self.assertRaises(err.NoMatchesError):
            finder.search_by_domain_id('1zzzA01')
        with self.assertRaises(err.InvalidInputError):
            finder.search_by_domain_id('bingo')
        with self.assertRaises(err.InvalidInputError):
            finder.search_by_domain_id(' file with &*! characters and spaces ')

    def test_ff_id_from_file(self):
        finder = util.FunfamFileFinder(
            base_dir=self.ff_dir, ff_tmpl='__SFAM__-ff-__FF_NUM__.reduced.sto')
        ff_file = finder.search_by_domain_id('2damA00')
        ff_id = finder.funfam_id_from_file(ff_file)
        self.assertEqual(ff_id.sfam_id, '1.10.8.10')
        self.assertEqual(ff_id.cluster_num, 14534)

    @log_level('cathpy.util', 'DEBUG')
    def test_scorecons(self):
        sc = util.ScoreconsRunner()
        aln = Align.from_fasta(self.example_fasta_file)
        sc_res = sc.run_fasta(self.example_fasta_file)
        self.assertEqual(sc_res.dops, 92.889)
        self.assertEqual(len(sc_res.scores), aln.aln_positions)
        del aln

        aln = Align.from_stockholm(self.example_sto_file)
        sc_res = sc.run_stockholm(self.example_sto_file)
        self.assertEqual(sc_res.dops, 61.529)
        self.assertEqual(len(sc_res.scores), aln.aln_positions)

    def test_groupsim(self):
        gs = util.GroupsimRunner()
        aln = Align.from_fasta(self.example_fasta_file)

        seqs = aln.seqs

        for s in seqs[:2]:
            s.set_cluster_id('0001')
        for s in seqs[2:]:
            s.set_cluster_id('0002')

        gs_res = gs.run_alignment(aln)
        self.assertEqual(gs_res.count_positions, aln.aln_positions)
        LOG.info("GS: {}".format(repr(gs_res.__dict__)))

        sto_file = tempfile.NamedTemporaryFile(delete=False, suffix='.sto')
        sto_with_groupsim_file = tempfile.NamedTemporaryFile(delete=False,
                                                             suffix='.groupsim.sto')

        LOG.info("Writing STOCKHOLM file (without groupsim): %s", sto_file.name)
        aln.write_sto(sto_file.name)
        LOG.info("Adding groupsim data ... ")
        aln.add_groupsim()
        LOG.info("Writing STOCKHOLM file (with groupsim): %s",
                 sto_with_groupsim_file.name)
        aln.write_sto(sto_with_groupsim_file.name)

        with open(sto_file.name) as f1:
            with open(sto_with_groupsim_file.name) as f2:
                lines1 = f1.readlines()
                lines2 = f2.readlines()
                ndiff = difflib.ndiff(lines1, lines2)

        difflines = [l for l in ndiff if not l.startswith(' ')]
        LOG.info("DIFF: %s", ''.join(difflines))
        expected_groupsim = '#=GC groupsim                 --------------10014101040141141031--2151411010022021221001040000---0-1-10-----\n'
        self.assertEqual(''.join(difflines), '+ ' + expected_groupsim)

    def test_groupsim_runner(self):

        aln = Align.from_fasta(self.example_fasta_file)

        # need to set the cluster id on sequences
        runner = GroupsimRunner()
        with self.assertRaises(err.InvalidInputError):
            runner.run_alignment(aln)

        for seq_idx, seq in enumerate(aln.sequences):
            seq.set_cluster_id('cluster1' if seq_idx < 5 else 'cluster2')

        result = runner.run_alignment(aln)
        self.assertIsInstance(result, GroupsimResult)

    def test_cluster_file(self):
        sc_path = os.path.abspath(self.sc_file)
        sc_dir = os.path.dirname(sc_path)
        sc_file = util.ClusterFile(sc_path)
        # 1.10.8.10__FF_SSG9__6.aln_reps.cora.fa
        self.assertDictEqual(sc_file.__dict__, {
            'path': sc_dir,
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
            'path': ff_dir,
            'sfam_id': '1.10.8.10',
            'cluster_type': 'ff',
            'cluster_num': '14534',
            'desc': '.reduced',
            'suffix': '.sto',
            'join_char': '-',
        })
        self.assertEqual(ff_file.to_string(), ff_path)
