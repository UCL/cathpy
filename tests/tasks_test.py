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
        LOG.info('   ... redis available')
        return True
    except redis.RedisError:
        LOG.info('Failed to connect to redis instance (skipping tests ...)')
        return False


def mongo_is_available():
    LOG.info('Testing mongodb connection...')
    try:
        client = pymongo.MongoClient(serverSelectionTimeoutMS=1000)
        client.server_info()
        LOG.info('   ... mongo available')
        return True
    except Exception as err:
        LOG.info(
            f'Failed to connect to mongo instance (ERR:{err}) (skipping tests ...)')
        return False


class TestObject:
    def __init__(self):
        self.id1 = '4b0hB00'
        self.id2 = '3hzaA01'
        self.expected_ssap_line = '4b0hB00  3hzaA01  143  127  86.65  114   79   26   1.60'
        self.expected_ssap_key = ssap.SsapResult.mk_key(
            cath_version='current', id1=self.id1, id2=self.id2)
        self.test_data_dir = path.join(path.dirname(__file__), 'data')
        self.ssap_pairs_file = path.join(
            self.test_data_dir, 'ssap_pairs.txt')
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


requires_redis = pytest.mark.skipif(not redis_is_available(),
                                    reason="Failed to connect to local redis instance (skipping tests)")
requires_mongo = pytest.mark.skipif(not mongo_is_available(),
                                    reason="Failed to connect to local mongo instance (skipping tests)")


@pytest.fixture
def test_object():
    return TestObject()


@requires_redis
@requires_mongo
def test_redis_datastore(test_object):
    ds = test_object.redis_ds
    key = test_object.expected_ssap_key
    ssap_res = ssap.SsapResult.from_string(
        test_object.expected_ssap_line,
        cath_version=test_object.cath_version)

    ds.redis.delete(key)
    assert ds.redis.get(key) is False

    ds.set_ssap_result(ssap_res)
    assert ds.redis.get(key) is True

    got_ssap_res = ds.get_ssap_result(
        id1=ssap_res.prot1,
        id2=ssap_res.prot2,
        pdb_path=test_object.pdb_path,
        cath_version=ssap_res.cath_version)
    assert isinstance(got_ssap_res, ssap.SsapResult)


@requires_mongo
def test_mongo_datastore(test_object):
    ds = test_object.mongo_ds
    key = test_object.expected_ssap_key
    ssap_res = ssap.SsapResult.from_string(
        test_object.expected_ssap_line,
        cath_version=test_object.cath_version)
    ds._collection.delete_one({'_id': key})
    assert ds._collection.find_one({'_id': key}) is False

    ds.set_ssap_result(ssap_res)
    assert ds._collection.find_one({'_id': key}) is True

    got_ssap_res = ds.get_ssap_result(
        id1=ssap_res.prot1,
        id2=ssap_res.prot2,
        cath_version=ssap_res.cath_version,
        pdb_path=test_object.pdb_path,)
    assert isinstance(got_ssap_res, ssap.SsapResult)


@requires_redis
def test_ssap_redis_sync(test_object):
    ds = test_object.redis_ds
    res = tasks.run_ssap(
        test_object.id1,
        test_object.id2,
        cath_version=test_object.cath_version,
        pdb_path=test_object.pdb_path,
        datastore_dsn='redis://localhost')


@requires_redis
def test_ssap_redis_pairs_sync(test_object):
    ds = test_object.redis_ds
    res = tasks.run_ssap_pairs(
        pairs=test_object.ssap_pairs,
        cath_version=test_object.cath_version,
        pdb_path=test_object.pdb_path,
        datastore_dsn='redis://localhost')


@requires_mongo
def test_ssap_mongo_sync(test_object):
    dsn = f"{test_object.mongo_ds}"
    res = tasks.run_ssap(
        test_object.id1,
        test_object.id2,
        cath_version=test_object.cath_version,
        pdb_path=test_object.pdb_path,
        datastore_dsn=dsn)


@requires_mongo
def test_ssap_mongo_pairs_sync(test_object):
    dsn = f"{test_object.mongo_ds}"
    res = tasks.run_ssap_pairs(
        pairs=test_object.ssap_pairs,
        cath_version=test_object.cath_version,
        pdb_path=test_object.pdb_path,
        datastore_dsn=dsn)

# def test_ssap_async(self):
#     ds = ssap.SsapStorageFactory.get('redis://localhost')
#     res = tasks.run_ssap.delay(
#         self.id1, self.id2, cath_version=self.cath_version, datastores=[ds])
#     res.get(timeout=5)
