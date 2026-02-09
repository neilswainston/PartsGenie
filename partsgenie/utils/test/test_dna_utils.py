"""
PartsGenie (c) GeneGenie Bioinformatics 2025

@author: neilswainston
"""
# pylint: disable=invalid-name
# pylint: disable=no-name-in-module
import json
import locale
import os

from Bio.Restriction import HgaI, MlyI

from partsgenie.utils import dna_utils


def test_copy():
    """Tests copy method."""
    directory = os.path.dirname(os.path.realpath(__file__))
    dna1 = _read_json(os.path.join(directory, 'sbol.json'))
    dna2 = dna1.copy()
    assert dna1 == dna2


def test_json():
    """Tests json roundtrip."""
    directory = os.path.dirname(os.path.realpath(__file__))
    dna1 = _read_json(os.path.join(directory, 'sbol.json'))
    dna2 = json.loads(json.dumps(dna1))
    assert dna1 == dna2


def test_app_restrict_site_linear():
    """Tests apply_restriction_site method."""
    _, dnas = _get_apply_restrict_site_dnas(MlyI, False)
    assert [len(dna['seq']) for dna in dnas] == [25, 831, 25]


def test_app_restrict_site_circular():
    """Tests apply_restriction_site method."""
    _, dnas = _get_apply_restrict_site_dnas('MlyI', True)
    assert [len(dna['seq']) for dna in dnas] == [50, 831]


def test_app_restrict_site_nomatch():
    """Tests aplly_restriction_site method."""
    parent, dnas = _get_apply_restrict_site_dnas(HgaI, False)
    assert len(dnas) == 1
    assert parent == dnas[0]


def _get_apply_restrict_site_dnas(restr, circ):
    """Tests apply_restriction_site method."""
    directory = os.path.dirname(os.path.realpath(__file__))
    par = _read_json(os.path.join(directory, 'restrict.json'))
    return par, dna_utils.apply_restricts(par, [restr], circ)


def _read_json(filename):
    """Parses json."""
    with open(filename, encoding=locale.getpreferredencoding()) as fle:
        return json.load(fle)
