"""
PartsGenie (c) GeneGenie Bioinformatics 2025

@author: neilswainston
"""
# pylint: disable=too-few-public-methods
import os

from partsgenie import client


def test_submit_simple():
    """Tests submit method with simple query."""
    _test_submit('simple_query.json')


def test_submit_complex():
    """Tests submit method with complex query."""
    _test_submit('complex_query.json')


def test_submit_promoter():
    """Tests submit method with promoter query."""
    _test_submit('promoter_query.json')


def test_submit_multiple():
    """Tests submit method with simple query."""
    _test_submit('multiple_query.json')


def _test_submit(filename):
    """Tests submit method."""
    directory = os.path.dirname(os.path.realpath(__file__))
    filepath = os.path.join(directory, filename)

    try:
        assert client.submit(filepath)
    except ValueError as err:
        if 'unfinished' in str(err):
            return

        raise err
