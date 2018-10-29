"""
CATH FunFHMMER - sequence search tool
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

def log_progress(_func=None, *, msg=None):
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
        if not accept:
            accept = self.default_accept
        headers = {'accept': accept}
        req = requests.get(url, headers=headers)
        return req

    def post(self, url, *, accept=None):
        if not accept:
            accept = self.default_accept
        headers = {'accept': accept}
        req = requests.post(url, headers=headers)
        return req

class SubmitResponse(object):
    def __init__(self, *, task_id):
        """Class storing the response from FunFHMMER submit."""
        self.task_id = task_id

class CheckResponse(object):
    def __init__(self, *, data, message, success):
        """Class storing the response from FunFHMMER status."""
        self.data = data
        self.message = message
        self.success = success

class ResultResponse(object):
    def __init__(self, *, query_fasta, funfam_scan, cath_version):
        """Class storing the response from FunFHMMER results."""
        self.query_fasta = query_fasta
        self.funfam_scan = Scan(**funfam_scan)
        self.cath_version = cath_version

    def as_json(self):
        data = jsonpickle.encode(self)
        return json.dumps(json.load(data), indent=4, sort_keys=True)

    def as_csv(self):
        """Returns the result as CSV"""
        result = self.funfam_scan.results[0]
        out = ( '# cath_version: {}\n'
                '# funfam members uniq_ec_terms query_region match_region evalue score description\n'
                ).format(self.cath_version)

        for hit in result.hits:
            ec_term_count = hit.data['ec_term_count'] or 0

            for hsp in hit.hsps:
                out += ' '.join([
                    '{:<20}'.format(hit.match_name),
                    '{:<5}'.format(hit.data['funfam_members']),
                    '{:<5}'.format(ec_term_count),
                    '{:>4}-{:<4}'.format(hsp.query_start, hsp.query_end),
                    '{:>4}-{:<4}'.format(hsp.hit_start, hsp.hit_end),
                    '{:<7.1e}'.format(hsp.evalue),
                    '{:<5d}'.format(int(hsp.score)),
                    '"{}"'.format(hit.match_description),
                ]) + '\n'
        return out


class Client(ApiClientBase):
    """
    Client for the CATH FunFhmmer API (protein sequence search server).

    The CATH FunFhmmer server allows users to locate matching CATH Functional
    Families (FunFams) in their protein sequence.
    """

    def __init__(self, *, base_url='http://www.cathdb.info', sleep=2, retries=50, log=None):
        super().__init__(base_url)
        self.sleep = sleep
        self.submit_url  = base_url + '/search/by_funfhmmer'
        self.check_url   = base_url + '/search/by_funfhmmer/check/:task_id'
        self.results_url = base_url + '/search/by_funfhmmer/results/:task_id'
        self.headers     = { 'accept': 'application/json' }
        self.retries     = retries
        if not log:
            log = logging.getLogger(__name__)
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
        
        task_id = self.submit(fasta).task_id

        @log_progress(msg="Waiting for job to complete ...")
        def wait_for_task(task_id, retries, sleep):
            retry_count = 0
            with trange(retries) as progress:
                progress.set_description("...")
                while True:
                    if retry_count > retries:
                        raise Exception('failed to get task within {} retries'.format(retries))
                    r = self.check(task_id)
                    if r.message == 'done':
                        progress.update(retries)
                        break
                    retry_count += 1
                    progress.set_description(r.message)
                    progress.update(retry_count)
                    time.sleep(sleep)

        wait_for_task(task_id, self.retries, self.sleep)

        r = self.results(task_id)
        
        return r
    
    @log_progress(msg="Submitting sequence search")
    def submit(self, fasta):
        """Submits a protein sequence to be searched and returns a task_id."""

        try:
            r = requests.post(self.submit_url, data={'fasta': fasta}, headers=self.headers)
        except:
            raise err.HttpError('http request failed: POST {}'.format(self.submit_url))

        try:
            data = r.json()
        except:
            raise err.JsonError('failed to parse json from http response: {}'.format(r.text))

        return SubmitResponse(**data)

    def check(self, task_id):
        """Checks the status of an existing search."""

        url = self.check_url.replace(':task_id', task_id)
        req = requests.get(url, headers=self.headers)
        data = req.json()
        return CheckResponse(**data)

    @log_progress(msg="Retrieving results")
    def results(self, task_id):
        """Retrieves the results of a search."""

        url = self.results_url.replace(':task_id', task_id)
        r = requests.get(url, headers=self.headers)
        data = r.json()
        return ResultResponse(**data)
