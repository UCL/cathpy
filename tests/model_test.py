import logging

from cathpy.models import AminoAcid, AminoAcids, CathID
from cathpy.error import OutOfBoundsError
from . import testutils

LOG = logging.getLogger(__name__)


class TestAminoAcids(testutils.TestBase):

    def test_aa(self):

        ala = AminoAcids.get_by_id('A')
        self.assertEqual(ala.one, 'A')
        self.assertEqual(ala.three, 'ALA')
        self.assertEqual(ala.word, 'alanine')

        ala = AminoAcids.get_by_id('ALA')
        self.assertEqual(ala.one, 'A')
        self.assertEqual(ala.three, 'ALA')
        self.assertEqual(ala.word, 'alanine')

        val = AminoAcids.get_by_id('val')
        self.assertEqual(val.one, 'V')


class TestCathID(testutils.TestBase):

    def test_cath_id(self):
        self.assertEqual(str(CathID("1")), "1")
        self.assertEqual(str(CathID("1.10.8")), "1.10.8")
        self.assertEqual(str(CathID("1.10.8.10.1.1.1.2.3")),
                         "1.10.8.10.1.1.1.2.3")

        self.assertEqual(CathID("1.10.8.10").sfam_id, "1.10.8.10")
        self.assertEqual(CathID("1.10.8.10.1").sfam_id, "1.10.8.10")
        with self.assertRaises(OutOfBoundsError) as err:
            cath_id = CathID("1.10.8").sfam_id
            self.assertRegex(err.exception, r'require depth',
                             'sfam_id fails when depth < 4')
