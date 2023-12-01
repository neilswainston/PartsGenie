'''
PartsGenie 2023

Additional utils from I2SysBio (CSIC/UV) and AI2 (UPV)

@author:  Pablo Carbonell
'''
import requests
import pandas as pd
from io import StringIO


def get_uniprot_sequence(uniprot_id):
    query = ("https://rest.uniprot.org/uniprotkb/search?query="+
             uniprot_id+
             "&format=tsv&fields=id,sequence")
    tab = pd.read_table(StringIO(requests.get(query).text))
    return tab.loc[0,'Sequence']

