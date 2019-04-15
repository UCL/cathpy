"""
cathpy.version - manipulating database / release versions
"""

# core
import logging
import re

from cathpy.error import ParseError

LOG = logging.getLogger(__name__)


class CathVersion(object):
    """Object that represents a CATH version."""

    CATH_VERSION_CURRENT = 'current'

    def __init__(self, *args, **kwargs):
        """Creates a new instance of a CathVersion object.

        The following are all equivalent:

            cv = CathVersion('4.2')
            cv = CathVersion('v4_2_0')
            cv = CathVersion(4.2)
            cv = CathVersion(4, 2, 0)

        """
        if len(args) == 1:
            ver_parts = self._split_string(args[0])
        elif len(args) == 2:
            ver_parts = (args[0], args[1], 0)
        elif len(args) == 3:
            ver_parts = args
        else:
            raise Exception("expected 1, 2, or 3 args (not {}): {}".format(
                len(args), " ".join([*args])))

        self.major = ver_parts[0]
        self.minor = ver_parts[1]
        self.trace = ver_parts[2]

    def join(self, join_char="."):
        """Returns the version string (with an optional join_char)."""
        if not self.is_current:
            return join_char.join([str(self.major), str(self.minor), str(self.trace)])
        else:
            return str(self.major)

    @property
    def is_current(self):
        """Returns whether the version corresponds to 'current' (eg HEAD)"""
        return self.major == self.CATH_VERSION_CURRENT

    @property
    def pg_dbname(self):
        """Return the version represented as a postgresql database (eg 'cathdb_v4_2_0')."""
        if self.is_current:
            return "cathdb_current"
        else:
            return "cathdb_v" + self.join("_")

    @property
    def dirname(self):
        """Return the version represented as a directory name (eg 'v4_2_0')."""
        if self.is_current:
            return "current"
        else:
            return "v" + self.join("_")

    @classmethod
    def new_from_string(cls, version_str):
        """Create a new CathVersion object from a string."""
        version_parts = cls._split_string(version_str)
        return cls(*version_parts)

    @classmethod
    def _split_string(cls, version_str):
        version_str = re.sub(r'^v', '', str(version_str))
        parts = re.split(r'[._]', version_str)
        if len(parts) == 1:
            if parts[0].lower() == cls.CATH_VERSION_CURRENT.lower():
                parts = (cls.CATH_VERSION_CURRENT, None, None)
            else:
                raise ParseError(
                    "failed to parse version '{}': expected '<major>.<minor>' or 'current'".format(version_str))
        if len(parts) == 2:
            parts = (parts[0], parts[1], 0)
        elif len(parts) == 3:
            parts = (parts[0], parts[1], parts[2])

        return parts

    def __str__(self):
        return self.join(".")

    def to_string(self):
        """Returns the CATH version in string form."""
        return str(self)
