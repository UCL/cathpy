# core
import logging
from os import path
import tempfile
from unittest.mock import patch
import pytest

# deps
import celery
from click.testing import CliRunner

# local
from cathpy.core.scripts.cath_cli import cli

from .testutils import TestBase, log_title, log_level

LOG = logging.getLogger(__name__)


class TestCli(TestBase):

    def test_version(self):
        run = CliRunner()
        res = run.invoke(cli, ['--version'])
        self.assertEqual(res.exit_code, 0, '--version exits okay')

    def test_help(self):
        run = CliRunner()
        res = run.invoke(cli, ['--help'])
        self.assertEqual(res.exit_code, 0, '--help exits okay')
        self.assertRegex(
            res.output, r'--help\s+Show this message and exit', '--help output okay')

    def test_list(self):
        run = CliRunner()
        res = run.invoke(cli, [])
        self.assertEqual(res.exit_code, 0, 'list options okay')
        self.assertRegex(
            res.output, r'Commands:', 'list options okay')


class TestSsapBatch(TestBase):

    def setUp(self):
        self.cmd = 'ssap-batch'
        self.datadir = path.join(path.dirname(__file__), 'data')

    def cli(self, extra_args=None):
        if not extra_args:
            extra_args = []
        run = CliRunner()
        args = [self.cmd, *extra_args]
        LOG.info(f"CMD: {' '.join(args)}")
        res = run.invoke(cli, args)
        self.assertEqual(res.exit_code, 0)
        return res

    def test_help(self):
        res = self.cli()
        self.assertRegex(
            res.output, r'Commands:', 'list options okay')

    def test_no_alignment(self):
        tmpfile = tempfile.NamedTemporaryFile(suffix='.ssaplist')
        res = self.cli([
            '--cath_version', '4.3',
            '--pairs', path.join(self.datadir, 'ssap_pairs.no_align.txt'),
            '--pdb_path', path.join(self.datadir, 'pdb'),
            '--datastore_dsn', tmpfile.name,
            '--inline',
            'submit'
        ])
        LOG.info(f"output: {res.output}, tmpfile: {tmpfile.name}")
