#!/usr/bin/env python3

import logging
import json
import pprint as pp
import time
from urllib.parse import urlparse

import click
import redis

from cathpy.core import tasks, ssap, version

#import click_log
# click_log.basic_config(LOG)

logging.basicConfig(
    level='INFO', format='%(asctime)s %(levelname)7s | %(message)s', datefmt='%Y/%m/%d %H:%M:%S')
LOG = logging.getLogger(__name__)

DEFAULT_PDB_PATH = '/cath/data/{cath_version}/pdb'
DEFAULT_DATASTORE_DSN = 'mongodb://localhost/cath'
DEFAULT_CHUNK_SIZE = 1000
DEFAULT_DUMP_MAX_ENTRIES = 1000


class SsapBatchContext(object):
    def __init__(self, *, pdb_path, datastore_dsn, cath_version, chunk_size, ssap_batch, force=False, inline=False, debug=False):
        self.pdb_path = pdb_path
        self.datastore_dsn = datastore_dsn
        self.ssap_batch = ssap_batch
        self.chunk_size = chunk_size
        self.cath_version = str(cath_version)
        self.force = force
        self.inline = inline
        self.debug = debug


@click.group()
@click.help_option('-h', '--help')
@click.version_option('0.1', '-v', '--version', message='%(prog)s %(version)s')
# @click_log.simple_verbosity_option(LOG, '--verbosity')
def cli():
    pass


def validate_cath_version(ctx, param, value):
    if value is None:
        return None
    try:
        return version.CathVersion(value)
    except ValueError as e:
        raise click.BadParameter(
            f'failed to parse cath version "{value}" (eg "current") err:{e}')


def validate_datastore_dsn(ctx, param, value):
    if value is None:
        return None
    try:
        ds = ssap.SsapStorageFactory.get(value)
        return value
    except ValueError as e:
        raise click.BadParameter(
            f'failed to parse dsn "{value}" (eg "mongodb://localhost/cath", "redis://localhost/0") err:{e}')


@click.command()
@click.option('--cath_version', type=str, callback=validate_cath_version, required=True)
@click.option('--datastore_dsn', type=str, callback=validate_datastore_dsn, default=DEFAULT_DATASTORE_DSN)
@click.option('--max_entries', type=int, default=DEFAULT_DUMP_MAX_ENTRIES)
def ssap_batch_dump(datastore_dsn, cath_version, max_entries):
    '''
    dump data from SSAP datastore
    '''
    ds = ssap.SsapStorageFactory.get(datastore_dsn)

    if not isinstance(ds, ssap.MongoSsapStore):
        raise Exception(
            'Sorry, currently only able to dump ssaps from mongodb datastore')

    # this was here for the redis datastore
    match_key = ssap.SsapResult.mk_key(
        cath_version=cath_version, id1='ID1', id2='ID2').replace('ID1-ID2', '*')

    LOG.info(f"Dumping all datastore records from CATH v'{cath_version}'")

    report_chunk_size = int(max_entries / 10) if max_entries else 10000

    total_entries_dumped = 0
    for idx, ssap_data in enumerate(ds.find({'cath_version': str(cath_version)}), 1):
        # key = key.decode('utf-8')
        # ssap_data = json.loads(ds.get(key))
        key = ssap_data.pop('_id', None)
        ssap_str = json.dumps({'id': key, **ssap_data})
        print(ssap_str)
        if report_chunk_size and idx % report_chunk_size == 0:
            LOG.info(f"Dumped {idx} records")
        total_entries_dumped += 1
        if max_entries and idx >= max_entries:
            break

    LOG.info(f"Dumped {total_entries_dumped} entries")

    LOG.info("DONE")


@click.group()
@click.option('--cath_version', type=str, callback=validate_cath_version, required=True)
@click.option('--pairs', type=click.Path(exists=True, file_okay=True), required=True)
@click.option('--datastore_dsn', type=str, callback=validate_datastore_dsn, default=DEFAULT_DATASTORE_DSN)
@click.option('--pdb_path', type=str, default=DEFAULT_PDB_PATH)
@click.option('--chunk_size', type=int, default=DEFAULT_CHUNK_SIZE)
@click.option('--force/--no-force', default=False)
@click.option('--inline', is_flag=True, default=False)
@click.option('--debug/--no-debug', default=False, envvar='CATHPY_DEBUG')
@click.pass_context
def ssap_batch_group(ctx, datastore_dsn, chunk_size, cath_version, pairs, pdb_path, force, inline, debug):
    '''
    manage SSAP jobs
    '''

    pdb_path = pdb_path.format(cath_version=cath_version.dirname)

    ssap_batch = ssap.SsapBatch(cath_version=cath_version,
                                datastore_dsn=datastore_dsn,
                                pairs_file=pairs,)
    opts = {
        'datastore_dsn': datastore_dsn,
        'chunk_size': chunk_size,
        'cath_version': cath_version,
        'debug': debug,
        'pdb_path': pdb_path,
        'ssap_batch': ssap_batch,
        'force': force,
        'inline': inline,
    }

    ctx.obj = SsapBatchContext(**opts)


@click.command()
@click.pass_obj
def ssap_batch_submit(batch_ctx):
    """
    Submits batches of SSAP pairs onto a queue
    """
    LOG.info('Loading missing SSAP pairs onto the queue')

    batch = batch_ctx.ssap_batch
    chunk_size = batch_ctx.chunk_size
    datastore_dsn = batch_ctx.datastore_dsn
    pdb_path = batch_ctx.pdb_path
    cath_version = batch_ctx.cath_version
    inline = batch_ctx.inline
    force = batch_ctx.force

    batch_count = 1
    task_results = []
    pairs_submitted = 0

    ssap_func = tasks.run_ssap_pairs if inline else tasks.run_ssap_pairs.delay

    for pairs_missing, processed in batch.read_pairs(chunk_size=chunk_size, force=force):
        if not pairs_missing:
            LOG.info("All pairs present, not submitting any ssap tasks")
            break

        LOG.info(
            f'Submitting {len(pairs_missing)} missing records from datastore ({datastore_dsn})')
        LOG.debug(
            f'cath_version:{cath_version} pdb_path:{pdb_path} datastore_dsn={datastore_dsn}')

        pairs_submitted += len(pairs_missing)

        result = ssap_func(pairs_missing, cath_version=str(
            cath_version), pdb_path=pdb_path, datastore_dsn=datastore_dsn)
        task_results.extend([result])

    if inline:
        LOG.info(
            f"Ran {len(task_results)} batch tasks ({pairs_submitted} SSAP pairs)")

    else:
        LOG.info(
            f"Created {len(task_results)} batch tasks")

        while True:
            status_types = ('PENDING', 'STARTED',
                            'RETRY', 'FAILURE', 'SUCCESS')
            status_count = {}
            for status_type in status_types:
                status_count[status_type] = len(
                    [t for t in task_results if t.status == status_type])

            LOG.info("Batch status at {} ({}):".format(
                time.strftime('%Y-%m-%d %H:%M:%S'), datastore_dsn))
            LOG.info(
                ' '.join(['{:15s}'.format(f'{t}:{status_count[t]:<3}') for t in status_types]))
            LOG.info("")

            if len(task_results) == len([t for t in task_results if t.status not in ('PENDING', 'STARTED', 'RETRY')]):
                LOG.info("All tasks finished")
                break

            time.sleep(5)


@click.command()
@click.pass_obj
def ssap_batch_info(batch_ctx):
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
def ssap_batch_worker(datastore_dsn, cath_version, pairs):
    """Start a worker to process SSAP pairs"""
    LOG.info('Starts worker(s) to process SSAP pairs from queue')


ssap_batch_group.add_command(ssap_batch_submit, name='submit')
ssap_batch_group.add_command(ssap_batch_info, name='info')
ssap_batch_group.add_command(ssap_batch_worker, name='worker')
#ssap_batch_group.add_command(ssap_batch_dump, name='dump')

cli.add_command(ssap_batch_dump, name='ssap-batch-dump')
cli.add_command(ssap_batch_group, name='ssap-batch')


if __name__ == '__main__':
    cli()
