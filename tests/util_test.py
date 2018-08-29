
import glob
import logging
import os
import re
import tempfile

from . import testutils

from cathpy import util, error as err

logger = logging.getLogger(__name__)

class TestUtil(testutils.TestBase):

    def setUp(self):
        self.cath_version = 'v4.2'
        self.data_dir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), 'data')
        self.sc_file = os.path.join(self.data_dir, 
            '1.10.8.10__FF_SSG9__6.aln_reps.cora.fa')
        self.ff_dir = os.path.join(self.data_dir, 'funfams')
        self.ff_tmpl = '__SFAM__-ff-__FF_NUM__.reduced.sto'

    @testutils.log_title
#    @testutils.log_level('cathpy.seqio', 'DEBUG')
    def test_merge(self):
        tmp_file = tempfile.NamedTemporaryFile(mode='w+', delete=False)
        logger.info("Creating SC merger...")
        merger = util.StructuralClusterMerger(cath_version=self.cath_version, 
            sc_file=self.sc_file, out_file=tmp_file.name,
            ff_dir=self.ff_dir, ff_tmpl=self.ff_tmpl)
        logger.info("Merging SC alignment {} to {}".format(self.sc_file, tmp_file.name))
        merge_aln = merger.run()
        self.assertEqual(merge_aln.count_sequences, 697)

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
        finder = util.FunfamFileFinder(base_dir=self.ff_dir, ff_tmpl='*.sto')
        self.assertIsInstance(finder, util.FunfamFileFinder)
        ff_file = finder.search_by_domain_id('2damA00')
        self.assertEqual(os.path.basename(ff_file), '1.10.8.10-ff-14534.reduced.sto')

        with self.assertRaises( err.FileNotFoundError ):
            finder.search_by_domain_id('1zzzA01')
        with self.assertRaises( err.InvalidInputError ):
            finder.search_by_domain_id('bingo')
        with self.assertRaises( err.InvalidInputError ):
            finder.search_by_domain_id(' file with &*! characters and spaces ')
