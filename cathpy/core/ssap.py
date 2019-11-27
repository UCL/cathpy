import json
import logging
import os
import re
import subprocess
import tempfile
from urllib.parse import urlparse

from cathpy.core import version, tasks

LOG = logging.getLogger()

DEFAULT_CHUNK_SIZE = 1000


class SsapBatch:
    """
    Provides management functions for performing a batch of SSAPs
    """

    def __init__(self, *, pairs_file, datastore, cath_version):
        self.pairs_file = pairs_file
        self.datastore = SsapStorageFactory.get(datastore)
        self.cath_version = version.CathVersion(cath_version)
        self._pairs_fh = None
        self._pairs_line_count = None

    def read_pairs(self, *, chunk_size=DEFAULT_CHUNK_SIZE, require_alignment=False, force=False):
        """
        Read through SSAP pairs, return results missing from the datastore
        """

        processed = 0
        pairs_missing = []

        with open(self.pairs_file) as fh:
            for line_number, line in enumerate(fh, 1):

                id1, id2 = line.strip().split()[0:2]

                key = SsapResult.mk_key(
                    id1=id1, id2=id2, cath_version=self.cath_version)

                found_ssap = True if key in self.datastore else False

                processed += 1

                # LOG.debug("[{:3d}] id1:{} id2:{}  [{}]".format(
                #     line_number, id1, id2, found_ssap))

                if not found_ssap or force:
                    pairs_missing.extend([[id1, id2]])

                if require_alignment:
                    raise Exception(
                        'need to add extra code to check for alignment in datastore')

                if len(pairs_missing) >= chunk_size:
                    yield (pairs_missing, processed)
                    pairs_missing = []
                    processed = 0

        yield (pairs_missing, processed)


class SsapResult:
    """
    Represents a SSAP Result

    ::

        1gaxA02  1qu3A04  227  150  79.71  137   60   20   3.06

    Args:
        prot1 (str): protein 1 identifier
        prot2 (str): protein 2 identifier
        length1 (int): number of residues in protein 1
        length2 (int): number of residues in protein 2
        ssap_score (float): SSAP score (out of 100)
        num_equivs (int): number of equivalent residues
        overlap_pc (int): percentage of overlapping residues
        seq_id_pc (int): percentage sequence identity
        rmsd (float): RMSD
        cath_version (str): CATH version
        alignment_text (str): SSAP alignment

    """

    def __init__(self, *, prot1, prot2, length1, length2, ssap_score, num_equivs, overlap_pc, seq_id_pc,
                 rmsd, cath_version, alignment_text=None):
        self.prot1 = prot1
        self.prot2 = prot2
        self.length1 = int(length1)
        self.length2 = int(length2)
        self.ssap_score = float(ssap_score)
        self.num_equivs = int(num_equivs)
        self.overlap_pc = int(overlap_pc)
        self.seq_id_pc = int(seq_id_pc)
        self.rmsd = float(rmsd)
        self.alignment_text = alignment_text
        self.cath_version = str(cath_version)

    @classmethod
    def from_string(cls, ssap_line, *, cath_version):
        """
        Create a new SsapResult object from text output of `cath-ssap`
        """
        prot1, prot2, length1, length2, ssap_score, num_equivs, overlap_pc, seq_id_pc, rmsd = ssap_line.split()
        return cls(prot1=prot1, prot2=prot2, length1=length1, length2=length2, ssap_score=ssap_score,
                   num_equivs=num_equivs, overlap_pc=overlap_pc, seq_id_pc=seq_id_pc, rmsd=rmsd,
                   cath_version=str(cath_version))

    def to_string(self):
        """
        Render this result as `cath-ssap` text output  
        """
        return '{:6s}  {:6s} {:4d} {:4d} {:6.2f} {:4d} {:4d} {:4d} {:6.2f}'.format(
            self.prot1, self.prot2, self.length1, self.length2, self.ssap_score,
            self.num_equivs, self.overlap_pc, self.seq_id_pc, self.rmsd)

    @classmethod
    def mk_key(cls, *, cath_version, id1, id2):
        return f'ssap-{cath_version}-{id1}-{id2}'

    def unique_key(self):
        return SsapResult.mk_key(cath_version=self.cath_version, id1=self.prot1, id2=self.prot2)

    @property
    def simax_score(self):
        simax = None
        if self.num_equivs > 0:
            simax = max(self.length1, self.length2) * \
                (self.rmsd / self.num_equivs)
        return simax

    def to_dict(self, *, ignore_empty=False):
        d = self.__dict__
        d['simax_score'] = self.simax_score
        if ignore_empty:
            d = {k: v for k, v in d.items() if v}
        return d

    def to_redis(self, redis, key=None):
        """
        Sets the SSAP result in the provided Redis datastore

        Args:
            redis (Redis): redis connection
        """
        if not key:
            key = self.unique_key()
        v = self.to_dict(ignore_empty=True)
        j = json.dumps(v)
        # LOG.info(f'Setting redis key {key}={v[:10]}...')
        try:
            redis.set(key, j)
        except:
            LOG.error(f"Encountered error when setting redis key '{key}'")
            raise

    @classmethod
    def from_redis(cls, redis, id1, id2, cath_version):
        """
        Attempts to create a SSAP result from the provided Redis datastore

        Args:
            redis (Redis): redis connection
        """
        key = cls.mk_key(cath_version=cath_version, id1=id1, id2=id2)
        j = redis.get(key)
        ssap_args = json.loads(j)
        LOG.info(f"Creating SsapResult from redis data={ssap_args}")
        del ssap_args['simax_score']
        ssap_result = cls(**ssap_args)
        return ssap_result


class SsapRunner:
    """
    Provides a wrapper around `cath-ssap` jobs

    Args:
    * results_file (str):
    * ssap_exe (str): path to `cath-ssap` executable
    * pdb_path (str): path to find PDB files
    """

    def __init__(self, *, ssap_exe='cath-ssap', pdb_path=None, cath_version):
        self.ssap_exe = ssap_exe
        self.pdb_path = pdb_path
        self.cath_version = cath_version
        self.tmp_aln_dir = tempfile.TemporaryDirectory()

    def run(self, id1, id2, *, lock=None, datastores=None):
        cath_version = self.cath_version

        ssap_key = SsapResult.mk_key(
            cath_version=cath_version, id1=id1, id2=id2)

        redis_stores = [
            ds for ds in datastores if isinstance(ds, RedisSsapStore)]
        if redis_stores:
            redis_ds = SsapStorageFactory.get(redis_stores[0])
            ssap_result = None
            if ssap_key in redis_ds:
                LOG.info(f"Found SSAP: [{ssap_key}]")
                try:
                    ssap_result = SsapResult.from_redis(
                        redis=redis_ds.redis, id1=id1, id2=id2, cath_version=cath_version)
                except Exception as e:
                    LOG.warning(
                        f"Failed to create SsapResult from existing data in redis ({e}), running again...")

            if ssap_result:
                return ssap_result

        ssap_args = [self.ssap_exe, '--aligndir', self.tmp_aln_dir.name]
        if self.pdb_path:
            ssap_args.extend(['--pdb-path', self.pdb_path])
        ssap_args.extend([id1, id2])
        ssap_result = None
        aln_file = os.path.join(self.tmp_aln_dir.name,
                                "{}{}.list".format(id1, id2))

        LOG.info("Running SSAP: %s %s", id1, id2)

        try:
            proc = subprocess.run(ssap_args,
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                  encoding='utf-8', check=True)
            ssap_result_text = proc.stdout.strip()

            aln_text = None
            ssap_result = SsapResult.from_string(
                ssap_result_text, cath_version=cath_version)
            with open(aln_file) as aln_fh:
                aln_text = aln_fh.read()
            ssap_result.alignment_text = aln_text

            if datastores:
                for datastore in datastores:
                    try:
                        ds = SsapStorageFactory.get(datastore)
                    except TypeError as err:
                        LOG.error(
                            "Failed to get ssap datastore with definition '{datastore}'")
                        raise

                    LOG.info(
                        f'Saving SSAP results [{ds.__class__.__name__}:{ssap_result.unique_key()}]')
                    ds.set_ssap_result(
                        ssap_result, lock=lock, alnfile=aln_file)

        except subprocess.CalledProcessError as error:
            LOG.error(("caught exception when running 'cath-ssap':\n",
                       "ARGS: %s\n",
                       "CODE: %s\n",
                       "STDOUT: %s\n",
                       "STDERR: %s\n"),
                      " ".join(ssap_args), error.returncode, error.stdout, error.stderr)
            # raise

        try:
            os.unlink(aln_file)
        except:
            pass

        return ssap_result


class SsapStore:
    pass


class FileSsapStore(SsapStore):
    def __init__(self, filepath):
        self.filepath = filepath
        self._already_processed = set()

    def set_ssap_result(self, ssap_result, *, lock, **kwargs):
        filepath = self.filepath
        lock.acquire(timeout=DEFAULT_SSAP_OUTPUT_TIMEOUT)
        try:
            with open(filepath, 'a') as ssap_logfh:
                ssap_logfh.write(f'{ssap_result.to_string()}\n')
        finally:
            lock.release()

        self.filepath

    def _process_existing_results(self):
        results_file = self.filepath
        try:
            LOG.info(
                "Building list of SSAP results we have already processed  ... ")
            with open(results_file, 'rt') as ssap_io:
                for line_count, ssapline in enumerate(ssap_io, 1):
                    if ssapline.startswith('#'):
                        continue
                    try:
                        ssap_result = SsapResult.from_string(ssapline)
                    except Exception as err:
                        raise Exception(
                            'failed to parse SSAP file {}, line {}: "{}": {}'.format(
                                results_file, line_count, ssapline, err))

                    ssap_key = ssap_result.unique_key()
                    self._already_processed.add(ssap_key)
            LOG.info("  ... found %d existing records",
                     len(self.already_processed))
        except FileNotFoundError:
            pass

    def __contains__(self, key):
        return key in self._already_processed


class TarSsapStore(SsapStore):
    def __init__(self, tarfile):
        self.tarfile = tarfile

    def set_ssap_result(self, ssap_result, *, alnfile, lock, **kwargs):
        tarfile = self.tarfile
        if not lock:
            raise MissingLockError(
                f'require lock when writing SSAP alignment to tarfile "{tarfile}"')

        LOG.debug("Locking thread to write to SSAP alignment tarfile ... ")
        with lock:
            try:
                LOG.info("Adding SSAP alignment '%s' to '%s' ... ",
                         alnfile, tarfile)
                with tarfile.open(tarfile, mode='a') as tar:
                    tar.add(alnfile)
                LOG.info("   ... added SSAP alignment: %s", alnfile)
            except Exception as error:
                LOG.error("caught error %s: %s",
                          error.__class__, str(error))
                raise

    def __contains__(self, key):
        raise Exception(
            'not going to try and search archive for key (choose another datastore)')
        return None


class RedisSsapStore(SsapStore):
    def __init__(self, url=None, *, db=0, host='localhost', port=6379):
        redis_args = {'host': host, 'port': port, 'db': db}
        if url:
            u = urlparse(url)
            if u.hostname:
                redis_args['host'] = u.hostname
            if u.port:
                redis_args['port'] = u.port

        self.host = redis_args['host']
        self.port = redis_args['port']
        self.db = redis_args['db']

        LOG.debug(f"Creating Redis connection: {redis_args}")
        redis = __import__('redis')
        self._redis = redis.Redis(**redis_args)

    @property
    def redis(self):
        return self._redis

    def set_ssap_result(self, ssap_result, *args, **kwargs):
        ssap_result.to_redis(self.redis)

    def __contains__(self, key):
        key_exists = self._redis.exists(key)
#        LOG.debug(f'redis.exists[{key}] {key_exists}')
        return key_exists

    def __str__(self):
        return f'redis://{self.host}:{self.port}/{self.db}'


class SsapStorageFactory:

    initialisers = {
        'file': lambda arg_str: FileSsapStore(arg_str),
        'tar': lambda arg_str: TarSsapStore(arg_str),
        'redis': lambda arg_str: RedisSsapStore(arg_str),
        #            **{k: v for kv in arg_str.split(';') for k, v in kv.split('=')}),
    }

    @classmethod
    def get(cls, dsn, *args, **kwargs):

        if isinstance(dsn, SsapStore):
            return dsn

        storage = None
        if not re.match(r'(\w+):', dsn):
            dsn = f'file:{dsn}'

        m = re.match(r'(\w+):(.*)', dsn)

        if not m:
            raise Exception(f'badly formatted data storage dsn "{dsn}"')

        storage_type, storage_args_str = m.groups()
        storage_type = storage_type.lower()
        if storage_type not in cls.initialisers:
            raise Exception(f'unknown storage type "{storage_type}"')

        storage = cls.initialisers[storage_type](storage_args_str)
        return storage
