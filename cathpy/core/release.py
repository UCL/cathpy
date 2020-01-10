"""
Provides access to read/write/manipulate CATH Release files
"""

# core
import logging
import re

import cathpy.core.error as err
from cathpy.core.models import CathID

LOG = logging.getLogger(__name__)


class HasCathIDMixin(object):
    """
    Mixin for classes that contain a :class:`cathpy.core.models.CathID`

    Usage:

    .. code:: python

        class MyClass(HasCathIDMixin, object):
            def __init__(self, *, cath_id, **kwargs):
                super(MyClass, self).super(cath_id=cath_id)

    Provides:

    .. code:: python

        self.cath_id
        self.cath_id_depth
        self.sfam_id
        self.cath_id_to_depth(n)

    """

    def __init__(self, *, cath_id, **kwargs):
        super().__init__(**kwargs)
        self._cath_id = CathID(cath_id)

    @property
    def cath_id(self):
        """
        Returns the full CATH ID
        """
        return self._cath_id.cath_id

    @property
    def cath_id_depth(self):
        """
        Returns the depth of the CATH ID (eg '1.10.8' = 3)
        """
        return self._cath_id.depth

    @property
    def sfam_id(self):
        """
        Returns the Superfamily ID of the CATH ID 
        """
        return self._cath_id.cath_id_to_depth(4)

    def cath_id_to_depth(self, depth):
        """
        Returns the CATH ID to a given depth
        """
        return self._cath_id.cath_id_to_depth(depth)

    def __lt__(self, other):
        if isinstance(other, str):
            other = CathID(other)
        return self._cath_id < other._cath_id


class HasEntriesWithCathIDMixin(object):
    """
    Mixin for container classes that have entries with :class:`HasCathIDMixin`
    """

    def filter_cath_id(self, cath_id):
        """
        Returns a new container after filtering to only include entries within a given CATH ID
        """
        if not isinstance(cath_id, CathID):
            cath_id = CathID(cath_id)
        depth = cath_id.depth
        filtered_entries = [
            c for c in self.entries if c.cath_id_to_depth(depth) == cath_id]
        return self.__class__(entries=filtered_entries)

    def filter_reps(self, depth):
        """
        Returns a new container after filtering to only include one rep at a given depth 
        """

        # sort by cath_id
        sorted_entries = sorted(self.entries)

        # take first entry at given depth
        reps = {}
        for entry in sorted_entries:
            rep_id = entry.cath_id_to_depth(depth)
            if rep_id not in reps:
                reps[rep_id] = entry

        # return the entries
        return self.__class__(entries=list(reps.values()))


class BaseReleaseFileList(list):
    """
    Base class for CATH release lists
    """

    entry_cls = None
    pk_key = None

    def __init__(self, *, entries=None, **kwargs):
        super().__init__(**kwargs)
        if not isinstance(entries, list):
            entries = list(entries)
        self._entries = entries

    def __getitem__(self, key):
        if isinstance(key, str) and self.pk_key:
            entries = [e for e in self._entries if getattr(
                e, self.pk_key) == key]
            if len(entries) == 1:
                return entries[0]
            elif entries == 0:
                return None
            else:
                raise err.TooManyMatchesError(
                    'found more than one match where {}={} (matches={})'.format(
                        self.pk_key, key, len(entries)))
        else:
            return self._entries[key]

    def __setitem__(self, key, value):
        self._entries[key] = value

    def __iter__(self):
        return iter(self._entries)

    @property
    def entries(self):
        """
        Returns the entries
        """
        return self._entries

    @classmethod
    def from_file(cls, fileio):
        """
        Creates a new instance of this object from a file
        """

        if isinstance(fileio, str):
            fileio = open(fileio, 'rt')

        entries = []
        try:
            for line in fileio:
                if line.startswith('#'):
                    continue
                entry = cls.entry_cls.from_string(line)
                entries.extend([entry])
        finally:
            fileio.close()

        return cls(entries=entries)

    def to_file(self, fileio):
        """
        Writes the entries out to file
        """

        if isinstance(fileio, str):
            fileio = open(fileio, 'wt')

        try:
            for entry in self.entries:
                fileio.write(entry.to_string() + "\n")
        finally:
            fileio.close()

    def sort(self, key=None):
        """
        Sort the entries inplace 
        """
        self._entries = sorted(self.entries, key=key)

    def __len__(self):
        """
        Returns the number of entries in the list
        """
        return len(self._entries)


class SegFrag(object):

    def __init__(self, *, chain_code, start_pdb, start_insert, end_pdb, end_insert, atom_length=None):
        self.chain_code = chain_code
        self.start_pdb = int(start_pdb)
        self.start_insert = start_insert
        self.end_pdb = int(end_pdb)
        self.end_insert = end_insert
        self.atom_length = atom_length


class Segment(SegFrag):
    pass


class Fragment(SegFrag):
    pass


class Domain(object):
    def __init__(self, *, segments):
        self.segments = segments


class CathDomallEntry(object):
    """
    Class representing a CATH Domall entry (chopping)
    """

    def __init__(self, *, chain_id, domains, fragments=None):
        self.chain_id = chain_id
        self.domains = domains
        if not fragments:
            fragments = []
        self.fragments = fragments

    def to_string(self):
        """
        Returns the entry as a Domall line
        """

        domall_line = "{} D{:02d} F{:02d} ".format(
            self.chain_id, len(self.domains), len(self.fragments))

        def insert_code(ins_code):
            return ins_code if ins_code else '-'

        for d in self.domains:
            domall_line += " {} ".format(len(d.segments))
            for s in d.segments:
                domall_line += " {} {:>4d} {} {} {:>4d} {} ".format(
                    s.chain_code, s.start_pdb, insert_code(s.start_insert),
                    s.chain_code, s.end_pdb, insert_code(s.end_insert),
                )

        for f in self.fragments:
            domall_line += " {} {:>4d} {} {} {:>4d} {} ({}) ".format(
                f.chain_code, f.start_pdb, insert_code(f.start_insert),
                f.chain_code, f.end_pdb, insert_code(f.end_insert),
                f.atom_length,
            )

        domall_line = domall_line.strip()

        return domall_line

    @classmethod
    def from_string(cls, domall_line):
        """
        Create a new instance from a Domall string

        Usage:

        ::

            domall_str = '10gsA D02 F01  2  A    2 - A   78 -  A  187 - A  208 -  1  A   79 - A  186 -  A  209 - A  209 - (1)'
            domall = Domall.from_string(domall_str)

            domall.domains[0].segments[0].chain_code    # 'A'
            domall.domains[0].segments[0].start_pdb     # 2
            domall.domains[0].segments[0].start_insert  # None
            domall.domains[0].segments[0].end_pdb       # 78
            domall.domains[0].segments[0].end_insert    # None

        """
        domall_line = domall_line.strip()
        cols = domall_line.split()
        chain_id, dom_count, frag_count = cols[0:3]

        dom_count = int(dom_count[1:])
        frag_count = int(frag_count[1:])

        domains = []
        fragments = []

        idx = 3
        dom_idx = 0
        while dom_idx < dom_count:
            seg_count = int(cols[idx])
            idx += 1
            segments = []
            # LOG.info("dom[%s] seg_count=%s", dom_idx, seg_count)
            for seg_idx in range(seg_count):
                start_chain, start_pdb, start_ins, end_chain, end_pdb, end_ins = cols[idx:idx+6]
                if start_ins == '-':
                    start_ins = None
                if end_ins == '-':
                    end_ins = None
                idx += 6
                seg = Segment(chain_code=start_chain, start_pdb=start_pdb, start_insert=start_ins,
                              end_pdb=end_pdb, end_insert=end_ins)
                # LOG.info("seg[%s]: %s", seg_idx, seg.__dict__)
                segments.extend([seg])
            dom = Domain(segments=segments)
            domains.extend([dom])
            dom_idx += 1

        frag_idx = 0
        while frag_idx < frag_count:
            start_chain, start_pdb, start_ins, end_chain, end_pdb, end_ins, frag_len = cols[
                idx:idx+7]
            idx += 7
            frag_idx += 1
            frag_len_match = re.match(r'^\((\d+)\)', frag_len)
            if not frag_len_match:
                raise err.ParseError(
                    'failed to parse frag len from "{}": {} (idx={})'.format(frag_len, domall_line, idx))
            atom_length = frag_len_match.group(1)
            frag = Fragment(chain_code=start_chain, start_pdb=start_pdb, start_insert=start_ins,
                            end_pdb=end_pdb, end_insert=end_ins, atom_length=atom_length)
            fragments.extend([frag])

        if idx != len(cols):
            raise err.ParseError('col index is {}, but there are {} columns: {}'.format(
                idx, len(cols), domall_line,
            ))

        return cls(
            chain_id=str(chain_id),
            domains=domains,
            fragments=fragments,
        )


class CathDomall(HasEntriesWithCathIDMixin, BaseReleaseFileList):
    """
    Class representing a CathDomall release file - domain boundaries

    Inherits:
     - :class:`BaseReleaseFileList`
     - :class:`HasEntriesWithCathIDMixin`

    """

    entry_cls = CathDomallEntry
    pk_key = 'chain_id'


class CathNamesEntry(HasCathIDMixin, object):
    """
    Class representing a CATH Names entry - name of node in CATH hierarchy
    """

    def __init__(self, *, cath_id, example_domain_id, name,):
        super(CathNamesEntry, self).__init__(cath_id=cath_id)
        self.example_domain_id = example_domain_id
        self.name = name

    def to_string(self):
        return "{}    {}    :{}".format(
            self.cath_id,
            self.example_domain_id,
            self.name,
        )

    @classmethod
    def from_string(cls, nameline):
        nameline = nameline.strip()
        cols = nameline.split(sep=None, maxsplit=2)
        if len(cols) != 3:
            raise err.ParseError("expected 3 cols in CathNames line '{}' (found {})".format(
                nameline, len(cols)))

        name = cols[2]
        if not name.startswith(':'):
            raise err.ParseError(
                "expected name '{}' to start with ':'".format(nameline))
        name = name[1:]

        return cls(
            cath_id=str(cols[0]),
            example_domain_id=str(cols[1]),
            name=str(name)
        )


class CathNamesList(HasEntriesWithCathIDMixin, BaseReleaseFileList):
    """
    Represents a CathNamesList - file containing names for nodes in the CATH hiearchy

    Inherits:
     - :class:`BaseReleaseFileList`
     - :class:`HasEntriesWithCathIDMixin`
    """

    entry_cls = CathNamesEntry
    pk_key = 'cath_id'


class CathDomainListEntry(HasCathIDMixin, object):
    def __init__(self, *, domain_id,
                 class_code, arch_code, top_code, homol_code,
                 s35_code, s60_code, s95_code, s100_code, domain_code,
                 atom_length, resolution,):

        codes = [str(code) for code in [class_code, arch_code, top_code,
                                        homol_code, s35_code, s60_code,
                                        s95_code, s100_code, domain_code]]
        cath_id = '.'.join(codes)
        super(CathDomainListEntry, self).__init__(cath_id=cath_id)

        self.domain_id = domain_id
        self.class_code = class_code
        self.arch_code = arch_code
        self.top_code = top_code
        self.homol_code = homol_code
        self.s35_code = s35_code
        self.s60_code = s60_code
        self.s95_code = s95_code
        self.s100_code = s100_code
        self.domain_code = domain_code
        self.atom_length = atom_length
        self.resolution = resolution

    def to_string(self):
        """
        Returns the entry as a CathDomainList string
        """
        return "{} {:5d} {:5d} {:5d} {:5d} {:5d} {:5d} {:5d} {:5d} {:5d} {:5d} {:.3f}".format(
            self.domain_id,
            self.class_code, self.arch_code, self.top_code,
            self.homol_code, self.s35_code, self.s60_code,
            self.s95_code, self.s100_code, self.domain_code,
            self.atom_length,
            self.resolution,
        )

    @classmethod
    def from_string(cls, domainlist):
        """
        Creates a new entry from a CathDomainList string
        """
        domainlist = domainlist.strip()
        cols = domainlist.split()
        if len(cols) != 12:
            raise err.ParseError("expected 12 cols in CathDomainList line '{}' (found {})".format(
                domainlist, len(cols)))

        return cls(
            domain_id=str(cols[0]),
            class_code=int(cols[1]),
            arch_code=int(cols[2]),
            top_code=int(cols[3]),
            homol_code=int(cols[4]),
            s35_code=int(cols[5]),
            s60_code=int(cols[6]),
            s95_code=int(cols[7]),
            s100_code=int(cols[8]),
            domain_code=int(cols[9]),
            atom_length=int(cols[10]),
            resolution=float(cols[11]),
        )

    def __repr__(self):
        return self.to_string()


class CathDomainList(HasEntriesWithCathIDMixin, BaseReleaseFileList):
    """
    Represents a CathDomainList - file containing classification of CATH domains

    Inherits:
     - :class:`BaseReleaseFileList`
     - :class:`HasEntriesWithCathIDMixin`
    """

    entry_cls = CathDomainListEntry
    pk_key = 'domain_id'
