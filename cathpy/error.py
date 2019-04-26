"""
CATH exception classes
"""

import logging

LOG = logging.getLogger(__name__)


class GeneralError(Exception):
    """General Exception class within the cathpy package."""


class DuplicateSequenceError(GeneralError):
    """More than one sequence in an alignment has the same id"""


class MissingExecutableError(GeneralError):
    """Missing an external executable."""


class MissingScoreconsError(MissingExecutableError):
    """Failed to find scorecons executable"""


class MissingGroupsimError(MissingExecutableError):
    """Failed to find groupsim executable"""


class JsonError(GeneralError):
    """Problem parsing JSON"""


class HttpError(GeneralError):
    """Problem getting/sending data over HTTP"""


class ParamError(GeneralError):
    """Incorrect parameters."""


class NoMatchesError(GeneralError):
    """No matches."""


class ParseError(GeneralError):
    """Failed to parse information."""


class TooManyMatchesError(GeneralError):
    """Found more matches than expected."""


class InvalidInputError(GeneralError):
    """Exception raised when an error is encountered due to incorrect input."""


class OutOfBoundsError(GeneralError):
    """Exception raised when code has moved outside expected boundaries."""


class MergeCheckError(GeneralError):
    """Exception raised when an error is encountered when checking the merge."""


class MissingSegmentsError(GeneralError):
    """Exception raised when segment information is missing."""


class FileEmptyError(GeneralError):
    """File is empty."""


class SeqIOError(GeneralError):
    """General Exception class within the SeqIO module"""


class GapError(SeqIOError):
    """Exception raised when trying to find residue information about a gap position."""


class MergeCorrespondenceError(SeqIOError):
    """Exception raised when failing to match correspondence sequences during alignment merge."""

    def __init__(self, *, seq_id, aln_type, seq_type, ref_no_gaps, corr_no_gaps):
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
