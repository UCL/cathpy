"""
Search protein sequence remotely against CATH FunFams
"""

import functools
import json
import logging
import time

import jsonpickle
import requests

from tqdm import trange

import cathpy.error as err
from cathpy.models import Scan

LOG = logging.getLogger(__name__)

DEFAULT_HMMSCAN_QUEUE = 'hmmscan_api'


def log_progress(_func=None, *, msg=None):
    """Provides a generic method to log progress"""

    def decorator_log_success(func):
        @functools.wraps(func)
        def wrapper_decorator(*args, **kwargs):
            title = "{:<70s} ... ".format(msg if msg else func.__name__)
            LOG.info(title)
            try:
                value = func(*args, **kwargs)
                LOG.info("%s %s", title, '[OK]')
            except:
                LOG.info("%s %s", title, '[FAIL]')
                raise
            return value
        return wrapper_decorator
    if _func is None:
        return decorator_log_success
    else:
        return decorator_log_success(_func)


class ApiClientBase(object):
    """Base class implementing default local behaviour of an API client."""

    def __init__(self, base_url, *, default_accept='application/json'):
        self.base_url = base_url
        self.default_accept = default_accept

    def get(self, url, *, accept=None):
        """Performs a GET request"""

        if not accept:
            accept = self.default_accept
        headers = {'accept': accept}
        req = requests.get(url, headers=headers)
        return req

    def post(self, url, *, accept=None):
        """Performs a POST request"""

        if not accept:
            accept = self.default_accept
        headers = {'accept': accept}
        req = requests.post(url, headers=headers)
        return req


class ResponseBase(object):
    """Base class that represents the HTTP response."""

    def __init__(self, **kwargs):
        pass


class SubmitResponse(ResponseBase):
    """Class that represents the response from FunFHMMER SUBMIT request."""

    def __init__(self, **kwargs):
        try:
            self.task_id = kwargs.pop('task_id')
        except KeyError:
            raise TypeError("Missing 'task_id' parameter")

        super().__init__(**kwargs)


class CheckResponse(ResponseBase):
    """Class that represents the response from FunFHMMER STATUS request."""

    def __init__(self, *, data, message, success, **kwargs):
        self.data = data
        self.message = message
        self.success = success

        super().__init__(**kwargs)


class ResultResponse(ResponseBase):
    """Class that represents the response from FunFHMMER RESULTS request."""

    def __init__(self, *, query_fasta, funfam_scan, funfam_resolved_scan, cath_version, **kwargs):
        self.query_fasta = query_fasta
        self.funfam_scan = Scan(**funfam_scan)
        self.funfam_resolved_scan = Scan(**funfam_resolved_scan)
        self.cath_version = cath_version

        super().__init__(**kwargs)

    def as_json(self, *, pp=False):
        """Returns the Scan as JSON formatted string."""

        data = jsonpickle.encode(self)
        if pp:
            data = json.dumps(json.loads(data), indent=2, sort_keys=True)

        LOG.info("Serialized ResultResponse as JSON string (length:%s)", len(data))
        return data


class Client(ApiClientBase):
    """
    Client for the CATH FunFhmmer API (protein sequence search server).

    The CATH FunFhmmer server allows users to locate matching CATH Functional
    Families (FunFams) in their protein sequence.
    """

    def __init__(self, *, base_url='http://www.cathdb.info', queue=DEFAULT_HMMSCAN_QUEUE,
                 sleep=2, retries=50, log=None):
        """Create a new API client.

        Args:
            base_url (str): override the default base URL
            queue (str): override the default queue
            sleep (int): number of seconds to wait between checks
            retries (int): number of times to check before quitting

        """

        super().__init__(base_url)
        self.sleep = sleep
        self.submit_url = base_url + '/search/by_funfhmmer'
        self.check_url = base_url + '/search/by_funfhmmer/check/:task_id'
        self.results_url = base_url + '/search/by_funfhmmer/results/:task_id'
        self.headers = {'accept': 'application/json'}
        self.retries = retries
        self.queue = queue
        if not log:
            log = LOG
        self.log = log

    def search_fasta(self, fasta=None, fasta_file=None):
        """Submits a sequence search and retrieves results."""

        log = self.log

        if not fasta and not fasta_file:
            raise Exception('expected either "fasta" or "fasta_file"')

        if not fasta:
            log.info('Reading sequence from file: %s ...', fasta_file)
            with open(fasta_file, 'r') as f:
                fasta = f.read()

        submit_response = self.submit(fasta, queue=self.queue)

        task_id = submit_response.task_id

        @log_progress(msg="Waiting for job to complete ...")
        def wait_for_task(task_id, retries, sleep):
            retry_count = 0
            with trange(retries) as progress:
                progress.set_description("...")
                while True:
                    if retry_count > retries:
                        raise Exception(
                            'failed to get task within {} retries'.format(retries))
                    r = self.check(task_id)
                    if r.message == 'done':
                        progress.update(retries)
                        break
                    retry_count += 1
                    progress.set_description(r.message)
                    progress.update(retry_count)
                    time.sleep(sleep)

        wait_for_task(task_id, self.retries, self.sleep)

        response = self.results(task_id)

        return response

    @log_progress(msg="Submitting sequence search")
    def submit(self, fasta, *, queue=None):
        """Submits a protein sequence to be searched; returns :class:`SubmitResponse`."""

        if not queue:
            queue = self.queue

        try:
            self.log.info("submit.POST: %s", self.submit_url)
            data = {'fasta': fasta, 'queue': queue}
            r = requests.post(self.submit_url, data=data, headers=self.headers)
        except:
            raise err.HttpError(
                'http request failed: POST {}'.format(self.submit_url))

        try:
            data = r.json()
        except:
            raise err.JsonError(
                'failed to parse json from http response: {}'.format(r.text))

        self.log.info("submit.POST.results: %s", str(data))

        return SubmitResponse(**data)

    def check(self, task_id):
        """Checks the status of an existing search, returns :class:`CheckResponse`"""

        url = self.check_url.replace(':task_id', task_id)

        self.log.info("check.GET: %s", url)
        req = requests.get(url, headers=self.headers)
        try:
            data = req.json()
        except Exception as e:
            LOG.error(
                "encountered error when trying to convert check request to JSON: (%s) %s", type(e), e)
            raise
        return CheckResponse(**data)

    @log_progress(msg="Retrieving results")
    def results(self, task_id):
        """Retrieves the results of a search, returns :class:`ResultResponse`"""

        url = self.results_url.replace(':task_id', task_id)

        self.log.info("results.GET: %s", url)
        r = requests.get(url, headers=self.headers)
        if not r.content:
            LOG.warning(
                "funfhmmer results returned empty content, assuming this means no results found")
            raise err.NoMatchesError
        else:
            data = r.json()
            return ResultResponse(**data)
