"""
PartsGenie (c) GeneGenie Bioinformatics 2025

@author: neilswainston
"""
# pylint: disable=too-few-public-methods
import json
import locale
import time

from partsgenie.parts import PartsThread


class Listener:
    """Default Listener class."""
    def __init__(self):
        self.status = None

    def event_fired(self, event):
        """Responds to event being fired."""
        self.status = event['update']['status']


def submit(query):
    """Submit."""

    if isinstance(query, str):
        with open(query, encoding=locale.getpreferredencoding()) as fle:
            query = json.load(fle)

    # Do job in new thread, return result when completed:
    thread = PartsThread(query, idx=0, verbose=True)
    listener = Listener()
    thread.add_listener(listener)
    thread.start()

    while True:
        if listener.status in ['finished', 'unfinished', 'cancelled', 'error']:
            break

        time.sleep(1)

    if listener.status == 'finished':
        return thread.solution

    raise ValueError(f'PartsGenie exited with status: {listener.status}')
