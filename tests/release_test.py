import difflib
import filecmp
import logging
from random import shuffle
import os
import tempfile

from cathpy.release import (
    CathDomainList, CathNamesList, CathDomall,
    CathDomainListEntry, CathDomallEntry, )

from . import testutils

LOG = logging.getLogger(__name__)


def cmp_file_contents(f1, f2, rstrip=False, max_diff=50):

    with open(f1) as fh1:
        with open(f2) as fh2:
            lines1 = [l for l in fh1 if not l.startswith('#')]
            lines2 = [l for l in fh2 if not l.startswith('#')]

    if rstrip:
        lines1 = [l.rstrip() for l in lines1]
        lines2 = [l.rstrip() for l in lines2]

    diff = difflib.unified_diff(
        lines1, lines2, fromfile=f1, tofile=f2)

    diff_lines = list(diff)

    if diff_lines:
        LOG.info("DIFF: %s %s", f1, f2)
        for idx, d in enumerate(diff_lines):
            LOG.info("%s", d.strip())
            if idx > max_diff:
                break

    return len(diff_lines)


class TestDomainList(testutils.TestBase):

    def setUp(self):
        self.domainlist_file = os.path.join(os.path.dirname(
            __file__), 'data', 'release', 'CathDomainList')
        self.domainlist = CathDomainList.from_file(self.domainlist_file)

    def test_domainlist(self):
        tmplistfile = tempfile.NamedTemporaryFile(mode='wt')

        domainlist = self.domainlist
        self.assertEqual(len(domainlist), 984)

        domainlist.to_file(tmplistfile.name)
        cmp_file_contents(self.domainlist_file, tmplistfile.name)
        self.assertEqual(cmp_file_contents(
            self.domainlist_file, tmplistfile.name), 0)

        domentry = domainlist[3]
        self.assertEqual(domentry.domain_id, '3friA01')
        self.assertEqual(domentry.cath_id, '1.10.8.10.2.1.1.1.2')
        self.assertEqual(domentry.sfam_id, '1.10.8.10')
        self.assertEqual([d.domain_id for d in domainlist[2:4]], [
                         '3frhA01', '3friA01'])

        self.assertEqual(domainlist['3friA01'].domain_id, '3friA01')
        self.assertEqual(domainlist['3friA01'].cath_id, '1.10.8.10.2.1.1.1.2')

    def test_sort(self):
        entries = [e for e in self.domainlist.entries]
        listcopy = CathDomainList(entries=entries)
        self.assertEqual(self.domainlist.entries, listcopy.entries)
        shuffle(listcopy.entries)
        self.assertNotEqual(self.domainlist.entries, listcopy.entries)
        sorted_entries = sorted(listcopy.entries)
        listsorted = CathDomainList(entries=sorted_entries)
        self.assertEqual(self.domainlist.entries, listsorted.entries)
        # sort in place doesn't work yet
        listcopy.sort()
        self.assertEqual(self.domainlist.entries, listcopy.entries)

    def test_filter_cath_id(self):
        filteredlist = self.domainlist.filter_cath_id('1.10.8.20')
        expected_domains = {
            '3hx3A01': '1.10.8.20.1.1.1.1.1',
            '4cj6A01': '1.10.8.20.1.1.1.2.1',
            '3hy5A01': '1.10.8.20.5.1.1.1.1',
            '4cizA01': '1.10.8.20.5.1.1.1.2',
        }
        entries = filteredlist.entries
        self.assertEqual(
            {c.domain_id: c.cath_id for c in entries[:2] + entries[-2:]}, expected_domains)

    def test_filter_reps(self):
        domainlist = self.domainlist
        hreps = domainlist.filter_reps(4)
        sreps = domainlist.filter_reps(5)
        self.assertIsInstance(hreps, CathDomainList)
        expected_hreps = ['1.10.8.10', '1.10.8.20',
                          '1.10.8.40', '1.10.8.50', '1.10.8.60']
        self.assertEqual([c.cath_id_to_depth(4)
                          for c in hreps.entries], expected_hreps)

        self.assertIsInstance(sreps[0], CathDomainListEntry)
        self.assertEqual(sreps[0].cath_id, '1.10.8.10.1.1.1.1.1')
        self.assertEqual(len(sreps[:3]), 3)
        self.assertEqual(len(sreps[-3:]), 3)

        expected_sreps = {
            '1oaiA00': '1.10.8.10.1',
            '3frhA01': '1.10.8.10.2',
            '4g3oA00': '1.10.8.10.3',
            '1c1kA01': '1.10.8.60.10',
            '4lrtC02': '1.10.8.60.11',
            '2fnaA02': '1.10.8.60.12', }

        self.assertEqual({c.domain_id: c.cath_id_to_depth(5)
                          for c in sreps[:3] + sreps[-3:]}, expected_sreps)


class TestNamesList(testutils.TestBase):

    def setUp(self):
        self.namelist_file = os.path.join(os.path.dirname(
            __file__), 'data', 'release', 'CathNames')

    def test_nameslist(self):
        tmplistfile = tempfile.NamedTemporaryFile(mode='wt')

        namelist = CathNamesList.from_file(self.namelist_file)
        self.assertEqual(len(namelist), 984)

        namelist.to_file(tmplistfile.name)
        self.assertEqual(cmp_file_contents(
            self.namelist_file, tmplistfile.name, rstrip=True), 0)

        self.assertEqual(namelist[0].cath_id, '1')
        self.assertEqual([n.cath_id for n in namelist[5:7]], ['1.20', '1.25'])
        self.assertEqual([n.name for n in namelist[5:7]], [
                         'Up-down Bundle', 'Alpha Horseshoe'])


class TestDomallList(testutils.TestBase):

    def setUp(self):
        self.domall_file = os.path.join(os.path.dirname(
            __file__), 'data', 'release', 'CathDomall')

    def test_line(self):
        entry_strings = (
            '10gsA D02 F01  2  A    2 - A   78 -  A  187 - A  208 -  1  A   79 - A  186 -  A  209 - A  209 - (1)',
            '1adiB D03 F00  2  B    1 - B  100 -  B  201 - B  265 -  1  B  101 - B  200 -  1  B  266 - B  431 -',
        )

        entry = CathDomallEntry.from_string(entry_strings[0])
        self.assertEqual(entry.chain_id, '10gsA')
        self.assertEqual(len(entry.domains), 2)
        self.assertEqual(len(entry.fragments), 1)

        for entry_string in entry_strings:
            entry = CathDomallEntry.from_string(entry_string)
            self.assertEqual(entry.to_string(), entry_string)

    def test_domall(self):
        tmplistfile = tempfile.NamedTemporaryFile(mode='wt')

        domall = CathDomall.from_file(self.domall_file)
        self.assertEqual(len(domall), 982)
        self.assertEqual([d.chain_id for d in domall[3:5]], ['103lA', '103mA'])

        domall.to_file(tmplistfile.name)
        self.assertEqual(cmp_file_contents(
            self.domall_file, tmplistfile.name), 0)

        domall.to_file(tmplistfile.name)
        newlist = CathDomall.from_file(tmplistfile.name)
        self.assertEqual(len(newlist), 982)
