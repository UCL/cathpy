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
        self.test_fasta_file = os.path.join(self.test_data_dir, '2damA00.fa')

    def test_basic(self):
        client = Client()
        response = client.search_fasta(fasta_file=self.test_fasta_file)
        self.assertIsInstance(response, ResultResponse)
        self.assertRegex(response.as_json(), r'1.10.8.10/FF/14534')
        self.assertRegex(response.as_csv(), r'1.10.8.10/FF/14534')
