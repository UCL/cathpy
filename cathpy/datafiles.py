"""
Access data files
"""

# core
import logging
import os
import re

from cathpy.version import CathVersion

LOG = logging.getLogger(__name__)


class GenericFileType(object):
    """Represents a type of CATH Data file."""

    type_id = None
    path_parts = []
    suffix = ''


class GcfFileType(GenericFileType):
    """Represents a GCF file type (registered as 'chaingcf')."""

    type_id = 'chaingcf'
    path_parts = ['chaingcf']
    suffix = '.gcf'


class AtomFastaFileType(GenericFileType):
    """Represents a FASTA file type (registered as 'atomfasta')."""

    type_id = 'atomfasta'
    path_parts = ['fasta']
    suffix = '.ATOM.fa'


class CombsFastaFileType(GenericFileType):
    """Represents a FASTA file type (registered as 'combsfasta')."""

    type_id = 'combsfasta'
    path_parts = ['fasta']
    suffix = '.COMBS.fa'


class FileTypes(object):

    _types = {}
    for cls in GenericFileType.__subclasses__():
        _types[cls.type_id] = cls

    @classmethod
    def get(cls, file_type):
        return cls._types.get(file_type)


class ReleaseDir(object):
    """
    Provides access to files relating to an official release of CATH.

    Args:
        cath_version: version of CATH (eg 'v4_2_0')
        base_dir: root directory for all data files (default: '/cath/data')
    """

    def __init__(self, cath_version, *, base_dir="/cath/data"):
        self.cath_version = CathVersion(cath_version)
        self.base_dir = base_dir

    def get_file(self, file_type, entity_id):
        """
        Returns the path for the given file type and identifier.

        Args:
            file_type: type of file (eg 'chaingcf')
            entity_id: identifier of the CATH entity (eg '1cukA')

        """

        if not isinstance(file_type, GenericFileType):
            file_type = FileTypes.get(file_type)

        path = os.path.join(self.base_dir, self.cath_version.dirname,
                            *file_type.path_parts, entity_id + file_type.suffix)

        LOG.debug("get_file(%s, %s): %s", file_type, entity_id, path)

        return path
