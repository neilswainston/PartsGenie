"""
PartsGenie (c) GeneGenie Bioinformatics 2025

@author: neilswainston
"""
# pylint: disable=no-member
# pylint: disable=no-name-in-module
# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals
import re
import uuid

from Bio.Restriction import RestrictionBatch
from Bio.Seq import Seq


SO_GENE = 'http://identifiers.org/so/SO:0000704'
SO_PROM = 'http://identifiers.org/so/SO:0000167'
SO_RBS = 'http://identifiers.org/so/SO:0000139'
SO_CDS = 'http://identifiers.org/so/SO:0000316'
SO_PART = 'http://identifiers.org/so/SO:0000804'
SO_ASS_COMP = 'http://identifiers.org/so/SO:0000143'
SO_PLASMID = 'http://identifiers.org/so/SO:0000637'
SO_DESIGNED = 'http://identifiers.org/so/SO:0000546'
SO_RANDOM = 'http://identifiers.org/so/SO:0000449'


def get_disp_id():
    """Get disp_id."""
    return '_' + str(uuid.uuid4()).replace('-', '_')


def apply_restricts(dna, restricts, circular=False):
    """Applies restriction site cleavage to forward and reverse strands."""
    out_dnas = [dna]

    for restrict in restricts:
        batch = RestrictionBatch()
        batch.add(str(restrict))
        restrict = batch.get(str(restrict))
        out_dnas = _apply_restrict_to_dnas(out_dnas, restrict, circular)

    return out_dnas


def add(dna1, dna2):
    """Adds two DNA objects together."""
    # Add names, etc.
    dna1.disp_id = get_disp_id()
    dna1['name'] = _concat([dna1['name'], dna2['name']])
    dna1['desc'] = _concat([dna1['desc'], dna2['desc']])

    # Add sequences:
    orig_seq_len = len(dna1['seq'])
    dna1['seq'] += dna2['seq']
    dna1['end'] = len(dna1['seq'])

    # Update parameters:
    for key, value in dna2['parameters'].items():
        param = dna1['parameters'].get(key, None)

        if param is None:
            dna1['parameters'][key] = value
        elif isinstance(param, list):
            param.append(value)
            dna1['parameters'][key] = param
        else:
            dna1['parameters'][key] = [param, value]

    # Update features:
    for feature in dna2.get('features', []):
        feature = feature.copy()
        feature['start'] += orig_seq_len
        feature['end'] += orig_seq_len
        dna1['features'].append(feature)

    return dna1


def expand(dna):
    """Add missing fields to dna dictionary."""
    dna['seq'] = ''.join([feat['seq'] for feat in dna['features']])
    dna['start'] = 1
    dna['end'] = dna['start'] + len(dna['seq'])

    start = 1

    for feat in dna['features']:
        feat['start'] = start
        feat['end'] = feat['start'] + len(feat['seq'])
        start = feat['end']

    return dna


def _concat(strs):
    """Concatenates non-None strings."""
    return ' - '.join([string for string in strs
                       if string is not None and len(string)])


def _apply_restrict_to_dnas(dnas, restrict, circular):
    """Applies restriction site cleavage to forward and reverse strands."""
    out_dnas = []

    for dna in dnas:
        restrict.catalyse(Seq(dna['seq']), linear=not circular)
        start = 0
        seqs = []

        for result in restrict.results:
            end = result - 1
            seqs.append((dna['seq'][start:end], start, end))
            start = end

        seqs.append((dna['seq'][start:], start, len(dna['seq'])))

        if circular:
            seqs = [(seqs[-1][0] + seqs[0][0], seqs[-1][1], seqs[0][2])] + \
                seqs[1:-1]

        for seq in seqs:
            out_dnas.append(_get_concat_dna(dna, seq[0], seq[1], seq[2]))

    return out_dnas


def _apply_restrict_to_seq(seq, restrict):
    """Applies restriction site cleavage to a sequence."""
    sub_seqs = [(match.group(0), match.start())
                for match in re.finditer(restrict, seq)]
    end = sub_seqs[0][1] if sub_seqs else len(seq)
    return [(seq[:end], 0)] + sub_seqs


def _get_concat_dna(parent_dna, seq, start, end):
    """Returns a DNA object from the supplied subsequence from a parent DNA
    object."""
    if parent_dna['seq'] == seq:
        return parent_dna

    # else:
    disp_id = get_disp_id()
    frag_str = ' [' + str(start) + ':' + str(end) + ']'

    dna = {
        'disp_id': disp_id,
        'name': parent_dna['name'] + frag_str,
        'desc': parent_dna['desc'] + frag_str,
        'seq': seq
    }

    for feature in parent_dna['features']:
        if feature['start'] >= start and feature['end'] <= end:
            copy_feature = feature.copy()
            copy_feature['start'] -= start
            copy_feature['end'] -= start
            dna['features'].append(copy_feature)

    return dna
