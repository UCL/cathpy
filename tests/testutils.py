
import logging
import unittest
from functools import wraps


def log_title(f):
    @wraps(f)
    def wrapper(*args, **kwargs):
        hr = "=" * 80
        log = logging.getLogger(f.__qualname__)
        log.info("")
        log.info(hr)
        log.info(" {} ".format(f.__name__))
        log.info(hr)
        log.info("")
        return f(*args, **kwargs)
    return wrapper

def log_level(name, level):
    """Decorator to temporarily change the log level of 'name' within this function."""
    def wrap(f):
        @wraps(f)
        def decorator(*args, **kwargs):
            log = logging.getLogger(name)
            # logging.info("LOG_LEVEL: getting logger for name '{}'".format(name))
            original_level = log.getEffectiveLevel()
            # logging.info("LOG_LEVEL: original level '{}'".format(original_level))
            log.setLevel(level)
            # logging.info("LOG_LEVEL: setting level to '{}' and running function ...".format(level))
            try:
                return_args = f(*args, **kwargs)
            except:
                raise
            finally:
                # logging.info("LOG_LEVEL: setting level back to '{}'".format(original_level))
                log.setLevel(original_level)
                # logging.info("LOG_LEVEL: log: ".format(repr(log)))

            return return_args
        return decorator
    return wrap

class TestBase(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.log = logging.getLogger(self.__class__.__name__)

    def log_title(self, title):
        hr = "=" * 80
        self.log.info("")
        self.log.info(hr)
        self.log.info(" {} ".format(title))
        self.log.info(hr)
        self.log.info("")
