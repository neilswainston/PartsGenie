"""
PartsGenie (c) GeneGenie Bioinformatics 2025

@author: neilswainston
"""
# pylint: disable=attribute-defined-outside-init
import json
import locale
import os
import time

from partsgenie.parts import PartsThread


class TestPartsThread():
    """Test class for PartsThread."""

    def test_submit_simple(self):
        """Tests submit method with simple query."""
        self.__test_submit('simple_query.json')

    def test_submit_complex(self):
        """Tests submit method with complex query."""
        self.__test_submit('complex_query.json')

    def test_submit_promoter(self):
        """Tests submit method with promoter query."""
        self.__test_submit('promoter_query.json')

    def test_submit_multiple(self):
        """Tests submit method with simple query."""
        self.__test_submit('multiple_query.json')

    def event_fired(self, event):
        """Responds to event being fired."""
        self.__status = event['update']['status']

    def __test_submit(self, filename):
        """Tests submit method."""
        self.__status = None
        directory = os.path.dirname(os.path.realpath(__file__))
        filename = os.path.join(directory, filename)

        with open(filename, encoding=locale.getpreferredencoding()) as fle:
            query = json.load(fle)

        self.__test_submit_query(query)

    def __test_submit_query(self, query):
        """Tests submit method."""
        # Do job in new thread, return result when completed:
        thread = PartsThread(query, idx=0, verbose=True)
        thread.add_listener(self)
        thread.start()

        while True:
            if self.__status in ['finished', 'unfinished', 'cancelled',
                                 'error']:
                break

            time.sleep(1)

        assert self.__status in ['finished', 'unfinished']


def _get_filename(filename):
    """Get filename."""
    directory = os.path.dirname(os.path.realpath(__file__))
    return os.path.join(directory, filename)
