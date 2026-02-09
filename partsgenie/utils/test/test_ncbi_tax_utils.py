"""
PartsGenie (c) GeneGenie Bioinformatics 2025

@author: neilswainston
"""
from partsgenie.utils.ncbi_tax_utils import TaxonomyFactory


def test_get_child_ids():
    """Tests get_child_ids method."""
    factory = TaxonomyFactory()

    assert factory.get_child_ids('208962') == ['502347', '550692', '550693',
                                               '909209', '910238', '1115511',
                                               '1440052']


def test_get_names():
    """Tests get_names method."""
    factory = TaxonomyFactory()

    assert 'Escherichia coli' in factory.get_names('562')
