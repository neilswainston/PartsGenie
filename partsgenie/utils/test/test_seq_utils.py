"""
PartsGenie (c) GeneGenie Bioinformatics 2025

@author: neilswainston
"""
import random

from Bio.Seq import Seq

from partsgenie.utils import seq_utils


def _get_seqs():
    """Get sequences."""
    nucl_seq = \
        'TACACACGAATAAAAGATAACAAAGATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGT' \
        'CCCAATTCTTGTTGAATTAGATGGCGATGTTAATGGGCAAAAATTCTCTGTCAGTGGAGA' \
        'GGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGGAA' \
        'GCTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTC' \
        'AAGATACCCAGATCATATGAAACAGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTA' \
        'TGTACAGGAAAGAACTATATTTTACAAAGATGACGGGAACTACAAGACACGTGCTGAAGT' \
        'CAAGTTTGAAGGTGATACCCTTGTTAATAGAATCGAGTTAAAAGGTATTGATTTTAAAGA' \
        'AGATGGAAACATTCTTGGACACAAAATGGAATACAACTATAACTCACATAATGTATACAT' \
        'CATGGCAGACAAACCAAAGAATGGAATCAAAGTTAACTTCAAAATTAGACACAACATTAA' \
        'AGATGGAAGCGTTCAATTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCC' \
        'TGTCCTTTTACCAGACAACCATTACCTGTCCACACAATCTGCCCTTTCCAAAGATCCCAA' \
        'CGAAAAGAGAGATCACATGATCCTTCTTGAGTTTGTAACAGCTGCTGGGATTACACATGG' \
        'CATGGATGAACTATACAAATAAATGTCCAGACTTCCAATTGACACTAAAGTGTCCGAACA' \
        'ATTACTAAATTCTCAGGGTTCCTGGTTAAATTCAGGCTGAGACTTTATTTATATATTTAT' \
        'AGATTCATTAAAATTTTATGAATAATTTATTGATGTTATTAATAGGGGCTATTTTCTTAT' \
        'TAAATAGGCTACTGGAGTGTAT'

    cds_seq = \
        'ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTT' \
        'GAATTAGATGGCGATGTTAATGGGCAAAAATTCTCTGTCAGTGGAGAGGGTGAAGGTGAT' \
        'GCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGGAAGCTACCTGTTCCA' \
        'TGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCAAGATACCCAGAT' \
        'CATATGAAACAGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAAAGA' \
        'ACTATATTTTACAAAGATGACGGGAACTACAAGACACGTGCTGAAGTCAAGTTTGAAGGT' \
        'GATACCCTTGTTAATAGAATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATT' \
        'CTTGGACACAAAATGGAATACAACTATAACTCACATAATGTATACATCATGGCAGACAAA' \
        'CCAAAGAATGGAATCAAAGTTAACTTCAAAATTAGACACAACATTAAAGATGGAAGCGTT' \
        'CAATTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCA' \
        'GACAACCATTACCTGTCCACACAATCTGCCCTTTCCAAAGATCCCAACGAAAAGAGAGAT' \
        'CACATGATCCTTCTTGAGTTTGTAACAGCTGCTGGGATTACACATGGCATGGATGAACTA' \
        'TACAAATAA'

    aa_seq = \
        'MSKGEELFTGVVPILVELDGDVNGQKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTL' \
        'VTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFYKDDGNYKTRAEVKFEGDTLV' \
        'NRIELKGIDFKEDGNILGHKMEYNYNSHNVYIMADKPKNGIKVNFKIRHNIKDGSVQLAD' \
        'HYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMILLEFVTAAGITHGMDELYK*'

    return nucl_seq, cds_seq, aa_seq


def test_get_frame_forward():
    """Tests get_frame method."""
    nucl_seq, cds_seq, aa_seq = _get_seqs()
    translations = seq_utils.translate(nucl_seq)

    expected = {'frame': 2,
                'start': 25,
                'end': 742,
                'nucl_seq': cds_seq,
                'aa_seq': aa_seq}

    assert translations[0] == expected


def test_get_frame_reverse():
    """Tests get_frame method (reverse)."""
    nucl_seq, cds_seq, aa_seq = _get_seqs()
    nucl_seq = Seq(nucl_seq).reverse_complement()
    translations = seq_utils.translate(nucl_seq)

    expected = {'frame': -2,
                'start': 180,
                'end': 897,
                'nucl_seq': str(Seq(cds_seq).reverse_complement()),
                'aa_seq': aa_seq}

    assert translations[0] == expected


def test_get_all_rev_trans():
    """Tests get_all_rev_trans method."""
    aa_seq = 'LS'
    rev_trans = seq_utils.get_all_rev_trans(aa_seq)
    assert len(rev_trans) == 36


def test_find_invalid():
    """Tests find_invalid method."""
    seq = 'ggtctaaaaatttttttaaaaaccagagtttttt'
    assert seq_utils.find_invalid(seq, 5, ['BsaI']) == [10, 11, 28]


def test_is_invalid():
    """Tests is_invalid method."""
    seq_inv = 'ggtctaaaaatttttttaaaaaccagagtttttt'
    assert seq_utils.is_invalid(seq_inv, 5, ['BsaI'])

    seq_val = 'ggtctaaaa'
    assert not seq_utils.is_invalid(seq_val, 5, ['BsaI'])


def test_get_random_dna():
    """Tests get_random_dna method."""
    lngth = random.randint(10, 100)

    assert lngth == len(seq_utils.get_random_dna(lngth, 4, ['BsaI']))


def test_get_seq_by_melt_temp():
    """Tests get_seq_by_melt_temp method."""
    seq, _ = seq_utils.get_seq_by_melt_temp('agcgtgcgaagcgtgcgatcctcc', 70)
    assert seq == 'agcgtgcgaagcgtgcgatc'


def test_get_rand_seq_by_melt_temp():
    """Tests get_rand_seq_by_melt_temp method."""
    target_temp = random.randint(50, 100)
    _, temp = seq_utils.get_rand_seq_by_melt_temp(target_temp, 4, ['BsaI'])
    assert abs(target_temp - temp) / target_temp < 0.025
