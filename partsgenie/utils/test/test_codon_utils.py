"""
PartsGenie (c) GeneGenie Bioinformatics 2025

@author: neilswainston
"""
import random

import pytest

from partsgenie.utils import codon_utils, seq_utils


def test_get_codon_optim_seq():
    """Tests get_codon_optim_seq method."""
    cod_opt = codon_utils.CodonOptimiser('9606')
    aa_codes = codon_utils.AA_CODES
    aa_codes.pop('Stop')

    aa_seq = ''.join([random.choice(list(aa_codes.values()))
                      for _ in range(random.randint(100, 2500))])

    max_repeat_nuc = 5
    restr_enzyms = ['BsaI']
    dna_seq = cod_opt.get_codon_optim_seq(aa_seq,
                                          max_repeat_nuc=max_repeat_nuc,
                                          restr_enzyms=restr_enzyms)

    assert not seq_utils.is_invalid(dna_seq, max_repeat_nuc, restr_enzyms)


def test_get_cai():
    """Tests get_cai method."""
    cod_opt = codon_utils.CodonOptimiser('83333')
    assert cod_opt.get_cai('AAACCC') == pytest.approx(0.6203738317757009)


def test_get_cai_human():
    """Tests get_cai method."""
    cod_opt = codon_utils.CodonOptimiser('9606')
    assert cod_opt.get_cai('AAACCC') == pytest.approx(0.8620689655172413)


def test_get_codon_usage_organisms_normal():
    """Tests get_random_codon method."""
    organisms = codon_utils.get_codon_usage_organisms()
    assert 'Escherichia coli' in organisms


def test_get_codon_usage_organisms_expand():
    """Tests get_random_codon method."""
    organisms = codon_utils.get_codon_usage_organisms(expand=True)
    assert 'Escherichia coli' in organisms
