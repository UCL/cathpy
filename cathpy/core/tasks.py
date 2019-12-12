import logging

from cathpy.core import models, ssap
from celery import Celery
import redis

DEFAULT_REDIS_HOST = 'localhost'
DEFAULT_REDIS_DB = '1'
DEFAULT_REDIS_PORT = 6379
DEFAULT_PDB_PATH = '/cath/data/current/pdb'

DEFAULT_CELERY_BROKER = 'redis://localhost'
DEFAULT_CELERY_BACKEND = 'redis://localhost'

app = Celery('cathpy.core.tasks',
             backend=DEFAULT_CELERY_BACKEND,
             broker=DEFAULT_CELERY_BROKER)


LOG = logging.getLogger(__name__)


@app.task
def run_ssap(id1, id2, *, cath_version, pdb_path=DEFAULT_PDB_PATH, datastore_dsn=None):

    ssap_runner = ssap.SsapRunner(
        pdb_path=pdb_path, cath_version=cath_version)

    ssap_runner.run(id1, id2, datastore_dsn=datastore_dsn)


@app.task
def run_ssap_pairs(pairs, *, cath_version, pdb_path=DEFAULT_PDB_PATH, datastore_dsn=None, max_consecutive_errors=5):

    ssap_runner = ssap.SsapRunner(
        pdb_path=pdb_path, cath_version=cath_version)

    consecutive_error_count = 0

    for idx, pair in enumerate(pairs):
        id1, id2 = pair
        try:
            ssap_runner.run(id1, id2, datastore_dsn=datastore_dsn)
            consecutive_error_count = 0
        except Exception as e:
            consecutive_error_count += 1
            LOG.error(
                f"Caught exception when running ssap {id1} {id2} (skipping) (err:{e})")

        if consecutive_error_count > max_consecutive_errors:
            msg = f"Encounted max consecutive errors ({consecutive_error_count}) - bailing out"
            LOG.error(msg)
            raise Exception(msg)
