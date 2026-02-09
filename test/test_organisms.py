"""
PartsGenie (c) GeneGenie Bioinformatics 2025

@author: neilswainston
"""
import organisms


def test_get_organisms():
    """Tests get_organisms first time, generating cached file."""
    _test_get_organisms(2)


def test_get_organisms_cached():
    """Tests get_organisms second time, using cached file."""
    _test_get_organisms(2)


def _test_get_organisms(parent_id):
    """Tests submit method."""
    orgs = organisms.get_organisms(parent_id)

    assert len(orgs) == 32863
    assert 'Escherichia coli' in orgs
    assert isinstance(orgs['Escherichia coli'], str)
