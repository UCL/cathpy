#!/usr/bin/env python3

import logging
import pprint as pp
import time

import click
import click_log
from cathpy.core import tasks, ssap, version

logging.basicConfig(
    level='DEBUG', format='%(asctime)s %(levelname)6s | %(message)s', datefmt='%Y/%m/%d %H:%M:%S')
LOG = logging.getLogger(__name__)
# click_log.basic_config(LOG)

DEFAULT_PDB_PATH = '/cath/data/{cath_version}/pdb'
DEFAULT_REDIS_DSN = 'redis://localhost'
DEFAULT_CHUNK_SIZE = 1000


class BatchContext(object):
    def __init__(self, *, pdb_path, redis, cath_version, chunk_size, ssap_batch, force=False, debug=False):
        self.pdb_path = pdb_path
        self.redis = redis
        self.ssap_batch = ssap_batch
        self.chunk_size = chunk_size
        self.cath_version = str(cath_version)
        self.force = force
        self.debug = debug


@click.group()
@click.help_option('-h', '--help')
@click.version_option('0.1', '-v', '--version', message='%(prog)s %(version)s')
@click_log.simple_verbosity_option(LOG, '--verbosity')
def cli():
    pass


def validate_cath_version(ctx, param, value):
    try:
        return version.CathVersion(value)
    except ValueError as e:
        raise click.BadParameter(
            f'failed to parse cath version "{value}" (eg "current") err:{e}')


def validate_redis_ds(ctx, param, value):
    try:
        ds = ssap.SsapStorageFactory.get(value)
        return value
    except ValueError as e:
        raise click.BadParameter(
            f'failed to parse dsn "{value}" (eg "redis://localhost") err:{e}')


@click.group()
@click.option('--cath_version', type=str, callback=validate_cath_version, required=True)
@click.option('--pairs', type=click.Path(exists=True, file_okay=True), required=True)
@click.option('--redis', type=str, callback=validate_redis_ds, default=DEFAULT_REDIS_DSN)
@click.option('--pdb_path', type=str, default=DEFAULT_PDB_PATH)
@click.option('--chunk_size', type=int, default=DEFAULT_CHUNK_SIZE)
@click.option('--force/--no-force', default=False)
@click.option('--debug/--no-debug', default=False, envvar='CATHPY_DEBUG')
@click.pass_context
def batch(ctx, redis, chunk_size, cath_version, pairs, pdb_path, force, debug):
    '''
    Groups all the options used to manage batches of SSAP jobs
    '''

    ssap_batch = ssap.SsapBatch(cath_version=cath_version,
                                datastore=redis,
                                pairs_file=pairs,)
    opts = {
        'redis': redis,
        'chunk_size': chunk_size,
        'cath_version': cath_version,
        'debug': debug,
        'pdb_path': pdb_path,
        'ssap_batch': ssap_batch,
        'force': force,
    }

    ctx.obj = BatchContext(**opts)


@click.command()
@click.pass_obj
def batch_load(batch_ctx):
    """
    Load batches of SSAP pairs onto a queue
    """
    LOG.info('Loading missing SSAP pairs onto the queue')

    batch = batch_ctx.ssap_batch
    chunk_size = batch_ctx.chunk_size
    redis = batch_ctx.redis
    pdb_path = batch_ctx.pdb_path
    cath_version = batch_ctx.cath_version
    force = batch_ctx.force

    batch_count = 1
    total_processed = 0
    celery_tasks = []
    for pairs_missing, processed in batch.read_pairs(chunk_size=chunk_size, force=force):
        if not pairs_missing:
            LOG.info("All pairs present, not submitting any ssap tasks")
            break

        LOG.info(
            f'Submitting {len(pairs_missing)} missing records from datastore ({batch.datastore})')
        LOG.debug(
            f'cath_version:{cath_version} pdb_path:{pdb_path} datastores={[redis]}')
        result = tasks.run_ssap_pairs.delay(
            pairs_missing, cath_version=str(cath_version), pdb_path=pdb_path, datastores=[redis])
        celery_tasks.extend([result])

    LOG.info(
        f"Created {len(celery_tasks)} batch tasks")

    while True:
        status_types = ('PENDING', 'STARTED', 'RETRY', 'FAILURE', 'SUCCESS')
        status_count = {}
        for status_type in status_types:
            status_count[status_type] = len(
                [t for t in celery_tasks if t.status == status_type])

        LOG.info("Batch status at {} ({}):".format(
            time.strftime('%Y-%m-%d %H:%M:%S'), redis))
        LOG.info(
            ' '.join(['{:15s}'.format(f'{t}:{status_count[t]:<3}') for t in status_types]))
        LOG.info("")

        if len(celery_tasks) == len([t for t in celery_tasks if t.status not in ('PENDING', 'STARTED', 'RETRY')]):
            LOG.info("All tasks finished")
            break

        time.sleep(5)


@click.command()
@click.pass_obj
def batch_info(batch_ctx):
    """Report information on this batch"""

    LOG.info('SSAP pairs in ')
    batch = batch_ctx.ssap_batch

    total_processed = 0
    total_missing = 0
    total_found = 0
    for pairs_missing, processed in batch.read_pairs():
        LOG.info(
            f'Processed {processed} records, found {len(pairs_missing)} missing from datastore ({batch.datastore})')

        total_missing += len(pairs_missing)
        total_found += processed - len(pairs_missing)
        total_processed += processed

    LOG.info(f'TOTAL_RECORDS: {total_processed}')
    LOG.info(f'TOTAL_MISSING: {total_missing}')
    LOG.info(f'TOTAL_FOUND:   {total_found}')


@click.command()
def batch_worker(redis, cath_version, pairs):
    """Start a worker to process SSAP pairs"""
    LOG.info('Starts worker(s) to process SSAP pairs from queue')


batch.add_command(batch_load, name='load')
batch.add_command(batch_info, name='info')
batch.add_command(batch_worker, name='worker')


cli.add_command(batch)


if __name__ == '__main__':
    cli()
