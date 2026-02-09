"""
PartsGenie (c) GeneGenie Bioinformatics 2025

@author: neilswainston
"""
# pylint: disable=broad-except
# pylint: disable=too-many-arguments
import re
import urllib

import requests


def get_uniprot_values(uniprot_ids, fields, batch_size=128, verbose=False):
    """Gets dictionary of ids to values from Uniprot."""
    values = []

    for i in range(0, len(uniprot_ids), batch_size):
        _get_uniprot_batch(uniprot_ids, i, batch_size, fields, values,
                           verbose)

    return {value['Entry']: value for value in values}


def search_uniprot(query, fields, limit=128):
    """Gets dictionary of ids to values from Uniprot."""
    values = []

    url = 'https://rest.uniprot.org/uniprotkb/search?query=' + \
        urllib.parse.quote(query) + \
        '&size=' + str(limit) + \
        '&format=tsv' + \
        '&fields=id,' + ','.join(
            [urllib.parse.quote(field) for field in fields])

    _parse_uniprot_data(url, values)

    return values


def _get_uniprot_batch(uniprot_ids, i, batch_size, fields, values, verbose):
    """Get batch of Uniprot data."""
    if verbose:
        print('seq_utils: getting Uniprot values ' + str(i) + ' - ' +
              str(min(i + batch_size, len(uniprot_ids))) + ' / ' +
              str(len(uniprot_ids)))

    batch = uniprot_ids[i:min(i + batch_size, len(uniprot_ids))]
    query = '+or+'.join(['accession:' + uniprot_id for uniprot_id in batch])

    url = 'https://rest.uniprot.org/uniprotkb/search?query=' + \
        query + \
        '&format=tsv' + \
        '&fields=accession,' + ','.join(
            [urllib.parse.quote(field) for field in fields])

    _parse_uniprot_data(url, values)


def _parse_uniprot_data(url, values):
    """Parses Uniprot data."""
    headers = None

    try:
        resp = requests.get(url, allow_redirects=True, timeout=3600)

        for line in resp.iter_lines():
            line = line.decode('utf-8')
            tokens = line.strip().split('\t')

            if headers is None:
                headers = tokens
            else:
                resp = dict(zip(headers, tokens))

                if 'Protein names' in resp:
                    regexp = re.compile(r'(?<=\()[^)]*(?=\))|^[^(][^()]*')
                    names = regexp.findall(resp.pop('Protein names'))
                    resp['Protein names'] = [nme.strip() for nme in names]

                for key in resp:
                    if key.startswith('Cross-reference'):
                        resp[key] = resp[key].split(';')

                values.append(resp)
    except Exception as err:
        print(err)
