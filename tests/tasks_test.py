import logging
from os import path
import tempfile

import pytest
from cathpy.core import ssap, tasks
from celery import Celery
from . import testutils

import redis
import pymongo

LOG = logging.getLogger(__name__)


REDIS_TEST_DATABASE = 3
MONGO_TEST_DATABASE = 'test'


def redis_is_available():
    LOG.info('Testing redis connection...')
    try:
        redis.Redis(db=REDIS_TEST_DATABASE).ping()
        return True
    except redis.RedisError:
        LOG.info('Failed to connect to redis instance (skipping tests ...)')
        return False


def mongo_is_available():
    LOG.info('Testing mongodb connection...')
    try:
        pymongo.MongoClient()
        return True
    except Exception as err:
        LOG.info(
            f'Failed to connect to mongo instance (ERR:{err}) (skipping tests ...)')
        return False


@pytest.mark.skipif(not redis_is_available(), reason="Failed to connect to local redis instance (skipping tests)")
@pytest.mark.skipif(not mongo_is_available(), reason="Failed to connect to local mongo instance (skipping tests)")
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
        self.mongo_ds = ssap.SsapStorageFactory.get(
            'mongodb://localhost', db=MONGO_TEST_DATABASE)

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
        ds.set_ssap_result(ssap_res)
        self.assertTrue(ds.redis.get(key), f'redis key {key} has been set')
        got_ssap_res = ds.get_ssap_result(
            id1=ssap_res.prot1, id2=ssap_res.prot2, cath_version=ssap_res.cath_version)
        self.assertIsInstance(got_ssap_res, ssap.SsapResult,
                              'got SsapResult instance from Redis')

    def test_mongo_datastore(self):
        ds = self.mongo_ds
        key = self.expected_ssap_key
        ssap_res = ssap.SsapResult.from_string(
            self.expected_ssap_line, cath_version=self.cath_version)
        ds._collection.delete_one({'_id': key})
        self.assertFalse(ds._collection.find_one({'_id': key}),
                         f'redis key {key} is empty')
        ds.set_ssap_result(ssap_res)
        self.assertTrue(ds._collection.find_one({'_id': key}),
                        f'redis key {key} has been set')
        got_ssap_res = ds.get_ssap_result(
            id1=ssap_res.prot1, id2=ssap_res.prot2, cath_version=ssap_res.cath_version)
        self.assertIsInstance(got_ssap_res, ssap.SsapResult,
                              'got SsapResult instance from Redis')

    def test_ssap_redis_sync(self):
        ds = self.redis_ds
        res = tasks.run_ssap(
            self.id1, self.id2, cath_version=self.cath_version, datastore_dsn='redis://localhost')

    def test_ssap_redis_pairs_sync(self):
        ds = self.redis_ds
        res = tasks.run_ssap_pairs(
            pairs=self.ssap_pairs, cath_version=self.cath_version, datastore_dsn='redis://localhost')

    def test_ssap_mongo_sync(self):
        dsn = f"{self.mongo_ds}"
        res = tasks.run_ssap(
            self.id1, self.id2, cath_version=self.cath_version, datastore_dsn=dsn)

    def test_ssap_mongo_pairs_sync(self):
        dsn = f"{self.mongo_ds}"
        res = tasks.run_ssap_pairs(
            pairs=self.ssap_pairs, cath_version=self.cath_version, datastore_dsn=dsn)

    # def test_ssap_async(self):
    #     ds = ssap.SsapStorageFactory.get('redis://localhost')
    #     res = tasks.run_ssap.delay(
    #         self.id1, self.id2, cath_version=self.cath_version, datastores=[ds])
    #     res.get(timeout=5)
