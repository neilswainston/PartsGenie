"""
PartsGenie (c) GeneGenie Bioinformatics 2025

@author: neilswainston
"""


def pairwise(iterable):
    """s -> (s0,s1), (s1,s2), (s2, s3), ..."""
    return [(iterable[i], iterable[i + 1]) for i in range(len(iterable) - 1)]
