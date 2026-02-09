"""
PartsGenie (c) GeneGenie Bioinformatics 2025

@author: neilswainston
"""
from partsgenie.utils import uniprot_utils


def test_get_uniprot_values():
    """Tests get_uniprot_values method."""
    result = uniprot_utils.get_uniprot_values(['P19367', 'P42212'],
                                              ['organism_id',
                                               'protein_name'], 1)

    expected = {'P19367': {'Entry': 'P19367',
                           'Organism (ID)': '9606',
                           'Protein names': ['Hexokinase-1',
                                             'EC 2.7.1.1',
                                             'Brain form hexokinase',
                                             'Hexokinase type I',
                                             'HK I',
                                             'Hexokinase-A']},
                'P42212': {'Entry': 'P42212',
                           'Organism (ID)': '6100',
                           'Protein names': ['Green fluorescent protein']}}

    assert result == expected
