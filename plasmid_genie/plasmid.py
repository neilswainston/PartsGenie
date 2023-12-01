'''
PartsGenie (c) University of Liverpool 2020

All rights reserved.

@author:  neilswainston
'''
# pylint: disable=wrong-import-order
from liv_utils import dna_utils, pairwise, seq_utils
#from liv_utils import dna_utils, ice_utils, pairwise, seq_utils
from thread_utils import AbstractThread


class PlasmidThread(AbstractThread):
    '''Runs a PlasmidGenie job.'''

    def __init__(self, query):
        AbstractThread.__init__(self, query)

    def run(self):
        '''Designs dominoes (bridging oligos) for LCR.'''
        iteration = 0

        if 'components' not in self._query:
            self.__get_components()

        self._fire_designs_event('running', iteration, 'Running...')

        for dsgn in self._query['designs']:
            orig_comps = [comp.copy() for comp in dsgn['components']]

            # Apply restriction site digestion:
            dsgn['components'] = \
                [_apply_restricts(dna, self._query['restr_enzs'])
                 for dna in dsgn['components']]

            # Generate plasmid DNA object:
            dna = dna_utils.concat(dsgn['components'])
            dna['name'] = dsgn['name']
            dna['typ'] = dna_utils.SO_PLASMID
            dna['parameters']['Design id'] = self._query['design_id']
            dna['children'].extend(orig_comps)

            # Generate domino sequences:
            cmps = dsgn['components'] + [dsgn['components'][0]] \
                if self._query['circular'] else dsgn['components']

            dna['children'].extend([self.__get_domino(pair)
                                    for pair in pairwise(cmps)])

            self._results.append(dna)

            iteration += 1
            self._fire_designs_event('running', iteration, 'Running...')

        if self._cancelled:
            self._fire_designs_event('cancelled', iteration,
                                     message='Job cancelled')
        else:
            self._fire_designs_event('finished', iteration,
                                     message='Job completed')

    def analyse_dominoes(self):
        '''Analyse sequences for similarity using BLAST.'''
        for design in self._query['designs']:
            ids_seqs = dict(zip(design['design'], design['seqs']))
            analysis = seq_utils.do_blast(ids_seqs, ids_seqs)

            for result in analysis:
                for alignment in result.alignments:
                    for hsp in alignment.hsps:
                        if result.query != alignment.hit_def:
                            return hsp

        return None

    def __get_components(self):
        '''Gets DNA components from ICE.'''
        iteration = 0

        self._fire_designs_event('running', iteration,
                                 'Extracting sequences from ICE...')

        for design in self._query['designs']:
            design['components'] = \
                [self.__get_component(ice_id)
                 for ice_id in design['design']
                 if ice_id]

            iteration += 1
            self._fire_designs_event('running', iteration,
                                     'Extracting sequences from ICE...')

    def __get_component(self, ice_id):
        '''Gets a DNA component from ICE.'''
        try:
            ice_client = ice_utils.get_ice_client(
                self._query['ice']['url'],
                self._query['ice']['username'],
                self._query['ice']['password'],
                group_names=self._query['ice'].get('groups', None))

            ice_entry = ice_client.get_ice_entry(ice_id)
            dna = ice_entry.get_dna()
            dna['desc'] = ice_id
            return dna
        finally:
            ice_client.close()

    def __get_domino(self, pair):
        '''Gets a domino from a pair of DNA objects.'''
        dna = dna_utils.concat([self.__get_domino_branch(pair[0], False),
                                self.__get_domino_branch(pair[1])])
        dna['parameters']['Type'] = 'DOMINO'
        return dna

    def __get_domino_branch(self, comp, forward=True):
        '''Gets domino branch from DNA object.'''
        target_melt_temp = self._query['melt_temp']
        reag_concs = self._query.get('reagent_concs', None)

        seq, melt_temp = seq_utils.get_seq_by_melt_temp(comp['seq'],
                                                        target_melt_temp,
                                                        forward,
                                                        reag_concs)

        # Flip name and description for Dominoes:
        dna = dna_utils.DNA(name=comp['desc'],
                            desc=comp['name'],
                            seq=seq,
                            typ=dna_utils.SO_ASS_COMP,
                            forward=True)

        dna['features'].append(dna.copy())
        dna['parameters']['Tm'] = float('{0:.3g}'.format(melt_temp))

        return dna


def _apply_restricts(dna, restr_enz):
    '''Apply restruction enzyme.'''
    if not restr_enz:
        return dna

    restrict_dnas = dna_utils.apply_restricts(dna, restr_enz)

    # This is a bit fudgy...
    # Essentially, return the longest fragment remaining after digestion.
    # Assumes prefix and suffix are short sequences that are cleaved off.
    restrict_dnas.sort(key=lambda x: len(x['seq']), reverse=True)
    return restrict_dnas[0]
