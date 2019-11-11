"""
Simple validation checks
"""

import logging
import re

LOG = logging.getLogger(__name__)


def is_valid_domain_id(id_str: str) -> bool:
    """
    Returns whether the given input is a valid CATH domain identifier.
    """
    return re.match('([0-9][a-zA-Z0-9]{3})([a-zA-Z0-9])([0-9]{2})$', id_str)
