"""
PartsGenie (c) GeneGenie Bioinformatics 2025

@author: neilswainston
"""
import pytest

from partsgenie import vienna_utils
from partsgenie.rbs_calculator import RbsCalculator


def test_calc_kinetic_score():
    """Tests calc_kinetic_score method."""
    r_rna = 'acctcctta'
    calc = RbsCalculator(r_rna, vienna_utils)

    m_rna = 'TTCTAGAGGGGGGATCTCCCCCCAAAAAATAAGAGGTACACATGACTAAAACTTTCA' + \
        'AAGGCTCAGTATTCCCACTGAG'

    start_pos = 41

    assert calc.calc_kinetic_score(m_rna, start_pos) == \
        pytest.approx(0.528571428571428)


def test_get_calc_dgs():
    """Tests calc_dgs method."""
    r_rna = 'acctcctta'
    calc = RbsCalculator(r_rna, vienna_utils)

    m_rna = 'TTCTAGAGGGGGGATCTCCCCCCAAAAAATAAGAGGTACACATGACTAAAACTTTCA' + \
        'AAGGCTCAGTATTCCCACTGAG'

    dgs = calc.calc_dgs(m_rna)
    assert list(dgs.keys()) == [41, 74]
    assert dgs[41][0] == pytest.approx(-6.088674036389431)
    assert dgs[74][0] == pytest.approx(5.793940143051147)

    dgs = calc.calc_dgs(m_rna)
    assert list(dgs.keys()) == [41, 74]
    assert dgs[41][0] == pytest.approx(-6.088674036389431)
    assert dgs[74][0] == pytest.approx(5.793940143051147)


def test_mfe_fail():
    """Tests mfe method."""
    m_rna = 'GCGGGAATTACACATGGCATGGACGAACTTTATAAATGA'

    energies, bp_xs, bp_ys = vienna_utils.run(
        'mfe', [m_rna], temp=37.0, dangles='none')

    assert energies == [0.0]
    assert not bp_xs[0]
    assert not bp_ys[0]


def test_subopt():
    """Tests subopt method."""
    r_rna = 'ACCTCCTTA'
    m_rna = 'AACCTAATTGATAGCGGCCTAGGACCCCCATCAAC'

    _, _, bp_ys = vienna_utils.run(
        'subopt', [m_rna, r_rna], temp=37.0,  dangles='all', energy_gap=3.0)

    nt_in_r_rna = False

    for bp_y in bp_ys:
        for nt_y in bp_y:
            if nt_y > len(m_rna):
                nt_in_r_rna = True

    assert nt_in_r_rna


def test_subopt_fail():
    """Tests subopt method."""
    r_rna = 'CCC'
    m_rna = 'CCC'

    energies, bp_xs, bp_ys = vienna_utils.run(
        'subopt', [m_rna, r_rna], temp=37.0, dangles='all', energy_gap=3.0)

    assert not energies
    assert not bp_xs
    assert not bp_ys
