# core
from .testutils import TestBase, log_title, log_level
import logging
from os import path
import tempfile
from unittest.mock import patch

# local
from cathpy.core.ssap import (SsapRunner, SsapResult, SsapStore,
                              SsapStorageFactory, MongoSsapStore, RedisSsapStore, FileSsapStore)


LOG = logging.getLogger(__name__)


class TestSsap(TestBase):

    def setUp(self):
        self.datadir = path.join(path.dirname(__file__), 'data')
        self.pdbdir = path.join(self.datadir, 'pdb')
        self.cath_version = 'current'
        self.example_ssap_result = '4b0hB00  3hzaA01  143  127  86.65  114   79   26   1.60'

    def test_file_ds(self):
        tmpfile = tempfile.NamedTemporaryFile(suffix='.ssaplist')
        dsn = f'file:{tmpfile.name}'

        runner = SsapRunner(cath_version=self.cath_version,
                            pdb_path=self.pdbdir)

        id1, id2, *_ = self.example_ssap_result.split()

        # run the ssap
        ssap_result = runner.run(id1=id1, id2=id2, datastore_dsn=dsn)
        self.assertIsInstance(ssap_result, SsapResult)

        # check the result looks okay
        self.assertEqual(ssap_result.to_string(),
                         self.example_ssap_result, 'expected ssap result')

        self.assertTrue(ssap_result.simax_score >
                        2 and ssap_result.simax_score < 2.5)

        # get the results from the datastore (file)
        ds = SsapStorageFactory.get(dsn)
        self.assertIsInstance(ds, SsapStore)
        ds_ssap_result = ds.get_ssap_result(
            id1=id1, id2=id2, cath_version=self.cath_version)

        # check the result looks okay
        self.assertIsInstance(ds_ssap_result, SsapResult)
        self.assertEqual(ds_ssap_result.to_string(),
                         self.example_ssap_result)

    # @patch('pymongo')
    # def test_mongo_ds(self, pymongo):
    #     ds = SsapStorageFactory.get('mongo://localhost')
    #     self.assertIsInstance(ds, SsapStore)

    # @patch('redis')
    # def test_redis_ds(self, redis):
    #     ds = SsapStorageFactory.get('redis://localhost')
