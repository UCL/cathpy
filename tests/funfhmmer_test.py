import glob
import logging
import os
import re
import tempfile

from .testutils import TestBase, log_title, log_level

from cathpy.funfhmmer import Client, ResultResponse

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
        re_expected_id = r'1.10.8.10/FF/14534'
        self.assertRegex(response.as_json(), re_expected_id)
        self.assertRegex(response.funfam_scan.as_json(), re_expected_id)
        self.assertRegex(response.funfam_scan.as_tsv(), re_expected_id)
        self.assertRegex(
            response.funfam_resolved_scan.as_tsv(), re_expected_id)

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
