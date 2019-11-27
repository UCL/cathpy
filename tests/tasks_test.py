import logging
from os import path
import tempfile

import redis
from cathpy.core import ssap, tasks
from celery import Celery
from . import testutils

LOG = logging.getLogger(__name__)


REDIS_TEST_DATABASE = 3


class TestSsap(testutils.TestBase):

    def setUp(self):
        self.id1 = '4b0hB00'
        self.id2 = '3hzaA01'
        self.expected_ssap_line = '4b0hB00  3hzaA01  143  127  86.65  114   79   26   1.60'
        self.expected_ssap_key = ssap.SsapResult.mk_key(
            cath_version='current', id1=self.id1, id2=self.id2)
        self.test_data_dir = path.join(path.dirname(__file__), 'data')
        self.ssap_pairs_file = path.join(self.test_data_dir, 'ssap_pairs.txt')
        self.ssap_results_file = path.join(
            self.test_data_dir, 'ssap_results.txt')
        self.cath_version = 'current'
        self.pdb_path = path.join(self.test_data_dir, 'pdb')
        self.redis_ds = ssap.SsapStorageFactory.get(
            'redis://localhost', db=REDIS_TEST_DATABASE)

        pairs = []
        with open(self.ssap_pairs_file) as fh:
            for line in fh:
                id1, id2 = line.strip().split()[:2]
                pairs.extend([[id1, id2]])
        self.ssap_pairs = pairs

    def test_redis_datastore(self):
        ds = self.redis_ds
        key = self.expected_ssap_key
        ssap_res = ssap.SsapResult.from_string(
            self.expected_ssap_line, cath_version=self.cath_version)
        ds.redis.delete(key)
        self.assertFalse(ds.redis.get(key),
                         f'redis key {key} is empty')
        ssap_res.to_redis(ds.redis)
        self.assertTrue(ds.redis.get(key), f'redis key {key} has been set')
        got_ssap_res = ssap.SsapResult.from_redis(
            ds.redis, id1=ssap_res.prot1, id2=ssap_res.prot2, cath_version=ssap_res.cath_version)
        self.assertIsInstance(got_ssap_res, ssap.SsapResult,
                              'got SsapResult instance from Redis')

    def test_ssap_sync(self):
        ds = self.redis_ds
        res = tasks.run_ssap(
            self.id1, self.id2, cath_version=self.cath_version, datastores=['redis://localhost'])

    def test_ssap_pairs_sync(self):
        ds = self.redis_ds
        res = tasks.run_ssap_pairs(
            pairs=self.ssap_pairs, cath_version=self.cath_version, datastores=['redis://localhost'])

    # def test_ssap_async(self):
    #     ds = ssap.SsapStorageFactory.get('redis://localhost')
    #     res = tasks.run_ssap.delay(
    #         self.id1, self.id2, cath_version=self.cath_version, datastores=[ds])
    #     res.get(timeout=5)
