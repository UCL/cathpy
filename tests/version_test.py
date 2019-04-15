import glob
import logging
import os
import re
import tempfile

from .testutils import TestBase, log_title, log_level

from cathpy.version import CathVersion

logger = logging.getLogger(__name__)


class TestVersion(TestBase):

    def test_cath_version(self):
        cv = CathVersion('4.2.0')
        self.assertIsInstance(cv, CathVersion)
        self.assertEqual(cv.major, '4')
        self.assertEqual(cv.minor, '2')
        self.assertEqual(cv.trace, '0')
        self.assertEqual(cv.dirname, 'v4_2_0')
        self.assertEqual(cv.join('-'), '4-2-0')
        self.assertIsInstance(CathVersion('v4.2'), CathVersion)
        self.assertEqual(str(CathVersion(4.2)), '4.2.0')

        cv_current = CathVersion.new_from_string('current')
        self.assertEqual(cv_current.pg_dbname, 'cathdb_current')
