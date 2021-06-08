import glob
import logging
import os
import re
import tempfile

from .testutils import TestBase, log_title, log_level

from cathpy.core.funfhmmer import Client, ResultResponse
import cathpy.core.error as err

LOG = logging.getLogger(__name__)


class TestFunfhmmer(TestBase):

    def setUp(self):
        self.test_data_dir = os.path.join(os.path.dirname(__file__), 'data')
        self.test_fasta_file1 = os.path.join(self.test_data_dir, '2damA00.fa')
        self.test_fasta_file2 = os.path.join(
            self.test_data_dir, 'A0A0Q0Y989_9BACI.fa')

    def test_basic(self):
        client = Client()
        response = client.search_fasta(fasta_file=self.test_fasta_file1)
        self.assertIsInstance(response, ResultResponse)
        re_expected_id = r'1.10.8.10[\-/]FF[\-/][0-9]+'
        self.assertRegex(response.as_json(pp=True), re_expected_id)

        self.assertRegex(response.funfam_scan.as_json(), re_expected_id)
        self.assertRegex(response.funfam_scan.as_tsv(), re_expected_id)
        self.assertRegex(
            response.funfam_resolved_scan.as_tsv(), re_expected_id)

    def test_bad_url(self):
        client = Client(base_url='http://invalid.base_url.that.does.not.exist')
        with self.assertRaises(err.HttpError):
            client.search_fasta(fasta_file=self.test_fasta_file1)

    def test_no_results(self):
        client = Client()
        with self.assertRaises(err.NoMatchesError):
            client.search_fasta(
                fasta='>test\nAAAAAAGGGGGAAAAAGGGGAAAAAAGGGGGAAAAAGGGGGAAAAGGGGGAAAAA')

    def test_resolved_results(self):
        client = Client()
        response = client.search_fasta(fasta_file=self.test_fasta_file2)

        LOG.info("funfam_resolved_scan.as_tsv: %s",
                 response.funfam_resolved_scan.as_tsv())

        self.assertEqual(len(response.funfam_resolved_scan.results),
                         1, 'resolved_scan has one result')
        first_result = response.funfam_resolved_scan.results[0]
        self.assertEqual(len(first_result.hits), 2,
                         'first resolved_scan result has correct number of hits')
