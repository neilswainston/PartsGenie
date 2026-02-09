"""
PartsGenie (c) GeneGenie Bioinformatics 2025

@author: neilswainston
"""
import json
import locale
import os
import time

from partsgenie.parts import PartsThread


class Listener:
    def event_fired(self, event):
        """Responds to event being fired."""
        self.status = event['update']['status']


def test_submit_simple():
    """Tests submit method with simple query."""
    _test_submit('simple_query.json')

def test_submit_complex():
    """Tests submit method with complex query."""
    _test_submit('complex_query.json')

def test_submit_promoter():
    """Tests submit method with promoter query."""
    _test_submit('promoter_query.json')

def test_submit_multiple():
    """Tests submit method with simple query."""
    _test_submit('multiple_query.json')

def _test_submit(filename):
    """Tests submit method."""
    directory = os.path.dirname(os.path.realpath(__file__))
    filename = os.path.join(directory, filename)

    with open(filename, encoding=locale.getpreferredencoding()) as fle:
        query = json.load(fle)

    # Do job in new thread, return result when completed:
    thread = PartsThread(query, idx=0, verbose=True)
    listener = Listener()
    thread.add_listener(listener)
    thread.start()

    while True:
        if listener.status in ['finished', 'unfinished', 'cancelled',
                                'error']:
            break

        time.sleep(1)

    assert listener.status in ['finished', 'unfinished']
