"""
PartsGenie (c) GeneGenie Bioinformatics 2025

@author: neilswainston
"""
# pylint: disable=invalid-name
# pylint: disable=no-member
# pylint: disable=wrong-import-order
import json
import os
import sys
import tempfile
import traceback
import uuid
import zipfile

import pandas as pd
from Bio import Restriction
from flask import Flask, jsonify, request, Response

import manager
from partsgenie.utils import uniprot_utils
import organisms


# Configuration:
SECRET_KEY = str(uuid.uuid4())

# Create application:
_STATIC_FOLDER = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                              'static')

app = Flask(__name__, static_folder=_STATIC_FOLDER)
app.config.from_object(__name__)
app.config['MAX_CONTENT_LENGTH'] = 1024 * 1024 * 1024
app.config['UPLOAD_FOLDER'] = tempfile.gettempdir()


_MANAGER = manager.Manager()

_ORG_PARENT_IDS = {'2157': 'Archaea',
                   '2': 'Bacteria',
                   '2759': 'Eukaryota'}

_ORGANISMS = {parent_id: organisms.get_organisms(parent_id)
              for parent_id in _ORG_PARENT_IDS}

DEBUG = False
TESTING = False


@app.route('/')
def home():
    """Renders homepage."""
    return app.send_static_file('index.html')


@app.route('/<path:path>')
def get_path(path):
    """Renders homepage."""
    return_path = path if path.startswith('export') else 'index.html'
    return app.send_static_file(return_path)


@app.route('/submit', methods=['POST'])
def submit():
    """Responds to submission."""
    return json.dumps({'job_ids': _MANAGER.submit(request.data)})


@app.route('/progress/<job_id>')
def progress(job_id):
    """Returns progress of job."""
    return Response(_MANAGER.get_progress(job_id),
                    mimetype='text/event-stream')


@app.route('/cancel/<job_id>')
def cancel(job_id):
    """Cancels job."""
    return _MANAGER.cancel(job_id)


@app.route('/organism_parents/')
def get_organism_parents():
    """Get organism parents."""
    return json.dumps(_ORG_PARENT_IDS)


@app.route('/organisms/', methods=['POST'])
def get_organisms():
    """Gets organisms from search term.
    Updated to assume r_rna corresponds to most prevalent
    See https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6107228/."""
    query = json.loads(request.data)

    data = [{'taxonomy_id': taxonomy_id,
             'name': name,
             'r_rna': 'acctccttt'}
            for name, taxonomy_id in _ORGANISMS[query['parent_id']].items()
            if query['term'].lower() in name.lower()]

    return json.dumps(data)


@app.route('/restr_enzymes')
def get_restr_enzymes():
    """Gets supported restriction enzymes."""
    return json.dumps([str(enz) for enz in Restriction.AllEnzymes])


@app.route('/uniprot/<query>')
def search_uniprot(query):
    """Search Uniprot."""
    fields = ['accession', 'protein_name', 'sequence', 'ec', 'organism_name',
              'organism_id']
    result = uniprot_utils.search_uniprot(query, fields)
    return json.dumps(result)


@app.route('/export', methods=['POST'])
def export_order():
    """Save export file, returning the url."""
    data = json.loads(request.data)['designs']
    file_id = str(uuid.uuid4()).replace('-', '_')
    dir_name = os.path.join(os.path.join(_STATIC_FOLDER, 'export'), file_id)

    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

    df = pd.DataFrame([[part['name'], part['seq']] for part in data],
                      columns=['Name', 'Sequence'])

    df.to_csv(os.path.join(dir_name, file_id + '.csv'), index=False)

    zip_file = os.path.join(dir_name + '.zip')

    with zipfile.ZipFile(zip_file, 'w') as zf:
        for dirpath, _, filenames in os.walk(dir_name):
            for filename in filenames:
                zf.write(os.path.join(dirpath, filename), filename)

    return json.dumps({'path': 'export/' + file_id + '.zip'})


@app.errorhandler(Exception)
def handle_error(error):
    """Handles errors."""
    app.logger.error('Exception: %s', (error))
    traceback.print_exc()

    response = jsonify({'message': traceback.format_exc()})
    response.status_code = 500
    return response


def main(argv):
    """main method."""
    if argv:
        app.run(host='0.0.0.0', threaded=True, port=int(argv[0]),
                use_reloader=False)
    else:
        app.run(host='0.0.0.0', threaded=True, use_reloader=False)


if __name__ == '__main__':
    main(sys.argv[1:])
