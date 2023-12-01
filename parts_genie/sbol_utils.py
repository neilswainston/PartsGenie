'''
PartsGenie (c) University of Liverpool 2020

All rights reserved.

@author:  neilswainston
'''
from sbml2sbol.sbol import Document, SO_CDS, SO_RBS

from liv_utils import dna_utils, uniprot_utils
import pandas as pd
import requests
from io import StringIO

def to_query(filename, taxonomy_id):
    '''Convert SBOL documents to PartsGenie query.'''
    doc = Document()
    doc.read(filename)
    return _to_query(doc, taxonomy_id)


def _to_query(doc, taxonomy_id):
    '''Get query.'''
    query = {}
    query['app'] = 'PartsGenie'

    query['designs'] = []

    for comp_def in [c for c in doc.componentDefinitions
                     if dna_utils.SO_GENE in c.roles]:
        query['designs'].append(_get_design(doc, comp_def))

    query['filters'] = {
        'max_repeats': 5,
        'gc_min': 0.25,
        'gc_max': 0.65,
        'local_gc_window': 50,
        'local_gc_min': 0.15,
        'local_gc_max': 0.8,
        'restr_enzs': [],
        'excl_codons': []
    }

    query['organism'] = {
        'taxonomy_id': taxonomy_id,
        'r_rna': 'acctccttt'
    }

    return query


def _get_design(doc, comp_def):
    '''Get design.'''
    design = {}
    design['name'] = comp_def.displayId
    design['desc'] = comp_def.identity
    design['features'] = []

    for sub_comp_def in [doc.getComponentDefinition(c.definition)
                         for c in comp_def.components]:
        design['features'].append(_get_feature(sub_comp_def))

    return design


def _get_feature(comp_def):
    '''Get feature.'''
    if dna_utils.SO_ASS_COMP in comp_def.roles:
        # Assembly component:
        return {
            'typ': dna_utils.SO_ASS_COMP,
            'name': comp_def.identity,
            'seq': '',
            'parameters': {
                'Tm target': 70
            },
            'temp_params': {
                'fixed': True,
                'required': ['name', 'tm'],
                'valid': True,
                'id': comp_def.displayId
            }
        }

    if SO_RBS in comp_def.roles:
        # RBS:
        return {
            'typ': dna_utils.SO_RBS,
            'name': comp_def.identity,
            'end': 60,
            'parameters': {
                'TIR target': float(comp_def.displayId.split('_')[1])
            },
            'temp_params': {
                'fixed': False,
                'required': ['name', 'tir'],
                'min_end': 35,
                'max_end': 10000,
                'valid': True,
                'id': comp_def.displayId
            }
        }

    if SO_CDS in comp_def.roles:
        # CDS:
        uniprot_id = comp_def.displayId.split('_')[0]
#        uniprot_vals = uniprot_utils.get_uniprot_values(
#            [uniprot_id], ['sequence'])
        query = ("https://rest.uniprot.org/uniprotkb/search?query="+
                 uniprot_id+
                 "&format=tsv&fields=id,sequence")
        tab = pd.read_table(StringIO(requests.get(query).text))
        
        uniprot_vals = { uniprot_id: {'Sequence': tab.loc[0,'Sequence']} }

        if uniprot_id not in uniprot_vals:
            raise ValueError('Uniprot id not found: %s' % uniprot_id)

        return {
            'typ': dna_utils.SO_CDS,
            'name': comp_def.identity,
            'temp_params': {
                'fixed': False,
                'required': ['name', 'prot'],
                'valid': True,
                'id': comp_def.displayId,
                'aa_seq': uniprot_vals[uniprot_id]['Sequence'],
                'orig_seq': uniprot_vals[uniprot_id]['Sequence']
            },
            'desc': '',
            'links': ['http://identifiers.org/uniprot/%s' % uniprot_id]
        }

    raise ValueError('Invalid roles in component definition: %s' % comp_def)
