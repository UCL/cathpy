import io
import logging
import re

logger = logging.getLogger(__name__)

class GeneralError(Exception):
    """General Exception class within the cathpy package."""
    pass

class ParamError(GeneralError):
    """Incorrect parameters."""
    pass

class InvalidInputError(GeneralError):
    """Exception raised when an error is encountered due to incorrect input."""
    pass

class OutOfBoundsError(GeneralError):
    """Exception raised when code has moved outside expected boundaries."""
    pass

class MergeCheckError(GeneralError):
    """Exception raised when an error is encountered when checking the merge."""
    pass

class MissingSegmentsError(GeneralError):
    """Exception raised when segment information is missing."""
    pass

class FileNotFoundError(GeneralError):
    """File not found."""
    pass

class SeqIOError(GeneralError):
    """General Exception class within the SeqIO module"""
    pass

class GapError(SeqIOError):
    """Exception raised when trying to find residue information about a gap position."""
    pass


class MergeCorrespondenceError(SeqIOError):
    """Exception raised when failing to match correspondence sequences during alignment merge."""

    def __init__(self, seq_id, aln_type, seq_type, ref_no_gaps, corr_no_gaps):
        message = ("Reference sequence '{}' in the {} alignment does not "
            "match the {} sequence provided in the Correspondence.\n"
            "{:>25}: {}\n"
            "{:>25}: {}\n").format(
                seq_id,
                aln_type,
                seq_type,
                'REF [{}]'.format(len(ref_no_gaps)), 
                ref_no_gaps,
                'CORRESPONDENCE [{}]'.format(len(corr_no_gaps)), 
                corr_no_gaps,
            )
        super().__init__(message)
