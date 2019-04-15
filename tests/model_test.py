import logging

from cathpy.models import AminoAcid, AminoAcids

from . import testutils

logger = logging.getLogger(__name__)


class TestAminoAcids(testutils.TestBase):

    def test_aa(self):
        ala = AminoAcids.get_by_id('A')
        self.assertEqual(ala.one, 'A')
        self.assertEqual(ala.three, 'ala')
        self.assertEqual(ala.word, 'alanine')
