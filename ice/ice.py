'''
PartsGenie (c) University of Liverpool 2020

All rights reserved.

@author:  neilswainston
'''
# pylint: disable=too-few-public-methods
# pylint: disable=too-many-arguments
#from liv_utils import dna_utils, ice_utils, net_utils
from liv_utils import dna_utils, net_utils
from thread_utils import AbstractThread


class IceThread(AbstractThread):
    '''Runs a save-to-ICE job.'''

    def __init__(self, query):
        AbstractThread.__init__(self, query)

        group_name = self._query['ice'].get('groups', None)
        self.__group_names = [group_name] if group_name else []

    def run(self):
        '''Saves results.'''
        keys = ['part', 'plasmid', 'strain']

        try:
            ice_client = ice_utils.get_ice_client(
                self._query['ice']['url'],
                self._query['ice']['username'],
                self._query['ice']['password'],
                group_names=self.__group_names)

            dna_writer = ice_utils.DNAWriter(ice_client)

            iteration = 0

            self._fire_designs_event(
                'running', iteration, 'Connecting to ICE...')

            url = self._query['ice']['url']
            self._query['ice']['url'] = url[:-1] if url[-1] == '/' else url

            for result in self._query['designs']:
                # Append links to results:
                links = {key: {'link': url + '/entry/' + str(entry_id),
                               'ice_id': str(entry_id)}
                         for key, entry_id in zip(
                             keys,
                             self.__write_design(result, ice_client,
                                                 dna_writer))
                         if entry_id}

                self._results.append(links)
                iteration += 1
                self._fire_designs_event('running', iteration, 'Saving...')

            if self._cancelled:
                self._fire_designs_event('cancelled', iteration,
                                         message='Job cancelled')
            else:
                self._fire_designs_event('finished', iteration,
                                         message='Job completed')
        except net_utils.NetworkError as err:
            self._fire_designs_event('error', iteration,
                                     message=err.get_text())
        finally:
            ice_client.close()

    def __write_design(self, result, ice_client, dna_writer):
        '''Write an individual design.'''
        ice_id, typ = dna_writer.submit(result)
        plasmid_id = ice_id if typ == 'PLASMID' else None
        strain_id = None

        if typ == 'PART' and self._query['ice'].get('plasmid', None):
            # Write plasmid.
            plasmid, _, _ = \
                write_ice_entry(ice_client,
                                ice_id,
                                self._query['ice']['plasmid'],
                                'PLASMID',
                                True,
                                self.__group_names)
            plasmid_id = plasmid.get_ice_id()

        if plasmid_id and self._query['ice'].get('strain', None):
            # Write strain.
            strain, _, _ = \
                write_ice_entry(ice_client,
                                plasmid_id,
                                self._query['ice']['strain'],
                                'STRAIN',
                                False,
                                self.__group_names)
            strain_id = strain.get_ice_id()\

        return ice_id, plasmid_id, strain_id


def write_ice_entry(ice_client, ice_id1, ice_id2, typ, write_seq, group_names):
    '''Write a composite ICE entry (part in plasmid, or plasmid in strain).'''
    comp1 = ice_client.get_ice_entry(ice_id1)
    comp2 = ice_client.get_ice_entry(ice_id2)

    name = comp1.get_metadata()['name'] + \
        ' (' + comp2.get_metadata()['name'] + ')'

    product = ice_utils.ICEEntry(typ=typ)
    product.set_values({'name': name[:127], 'shortDescription': name})

    taxonomy = comp1.get_parameter('Taxonomy')

    if taxonomy:
        product.set_parameter('Taxonomy', taxonomy)

    ice_client.set_ice_entry(product)
    ice_client.add_link(product.get_ice_id(), comp1.get_ice_id())
    ice_client.add_link(product.get_ice_id(), comp2.get_ice_id())

    if write_seq:
        product.set_dna(dna_utils.concat(
            [comp1.get_dna(), comp2.get_dna()]))

    ice_client.set_ice_entry(product)

    if group_names:
        groups = ice_client.get_groups()

        for group_name in group_names:
            ice_client.add_permission(product.get_ice_id(),
                                      groups[group_name])

    return product, comp1, comp2
