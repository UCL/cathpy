import glob
import logging
import os
import re
import tempfile

from .testutils import TestBase, log_title, log_level

from cathpy.version import CathVersion
import cathpy.error as err

LOG = logging.getLogger(__name__)


class TestVersion(TestBase):

    def test_init(self):
        self.assertEqual(str(CathVersion('4.2.0')), '4.2.0')
        self.assertEqual(str(CathVersion('4.2')), '4.2.0')
        with self.assertRaises(err.ParseError):
            CathVersion('4')
        self.assertEqual(str(CathVersion(4, 1, 1)), '4.1.1')
        self.assertEqual(str(CathVersion(4, 1)), '4.1.0')

    def test_v4_2_0(self):
        cv = CathVersion('4.2.0')
        self.assertIsInstance(cv, CathVersion)
        self.assertEqual(cv.major, '4')
        self.assertEqual(cv.minor, '2')
        self.assertEqual(cv.trace, '0')
        self.assertEqual(cv.dirname, 'v4_2_0')
        self.assertEqual(cv.join('-'), '4-2-0')
        self.assertFalse(cv.is_current)

    def test_current(self):
        cv_current = CathVersion.from_string('current')
        self.assertEqual(cv_current.pg_dbname, 'cathdb_current')
        self.assertEqual(cv_current.dirname, 'current')
        self.assertEqual(cv_current.join('-'), 'current')
        self.assertEqual(cv_current.major, 'current')
        self.assertEqual(cv_current.minor, None)
        self.assertEqual(cv_current.trace, None)
        self.assertTrue(cv_current.is_current)
