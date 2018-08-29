# core
import logging
import os
import re

logger = logging.getLogger(__name__)

class GenericFileType(object):
    """Represents a type of CATH Data file."""

    type_id = None
    path_parts = []
    suffix = ''

class GcfFileType(GenericFileType):

    type_id = 'chaingcf'
    path_parts = ['chaingcf']
    suffix = '.gcf'

class FileTypes(object):

    _types = {}
    for cls in GenericFileType.__subclasses__():
        _types[cls.type_id] = cls

    @classmethod
    def get(cls, file_type):
        return cls._types.get(file_type)

class ReleaseDir(object):
    """Provides access to files relating to an official release of CATH."""

    def __init__(self, cath_version, *, base_dir="/cath/data"):
        """TODO"""
        self.cath_version = cath_version
        self.base_dir = base_dir

    def get_file(self, file_type, id):
        
        if not isinstance(file_type, GenericFileType):
            file_type = FileTypes.get(file_type)
        
        return os.path.join(self.base_dir, self.cath_version.dirname, 
            *file_type.path_parts, id + file_type.suffix)
