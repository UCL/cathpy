"""CATH Models used"""

class ScanHsp(object):
    """Object to store the High Scoring Pair (HSP) from a sequence scan."""
    def __init__(self, *, evalue, hit_start, hit_end, hit_string=None,
                 homology_string=None, length, query_start, query_end,
                 query_string=None, rank, score, **kwargs):
        self.evalue = evalue
        self.hit_start = hit_start
        self.hit_end = hit_end
        self.hit_string = hit_string
        self.homology_string = homology_string
        self.length = length
        self.query_start = query_start
        self.query_end = query_end
        self.query_string = query_string
        self.rank = rank
        self.score = score

class ScanHit(object):
    """Object to store a hit from a sequence scan."""
    def __init__(self, *, match_name, match_cath_id, match_description,
                 match_length, hsps, significance, data, **kwargs):
        self.match_name = match_name
        self.match_cath_id = match_cath_id
        self.match_description = match_description
        self.match_length = match_length
        self.hsps = [ScanHsp(**hsp) for hsp in hsps]
        self.data = data
        self.significance = significance

class ScanResult(object):
    """Object to store a result from a sequence scan."""
    def __init__(self, *, query_name, hits, **kwargs):
        self.query_name = query_name
        self.hits = [ScanHit(**hit) for hit in hits]

class Scan(object):
    """Object to store a sequence scan."""
    def __init__(self, *, results, **kwargs):
        self.results = [ScanResult(**res) for res in results]
