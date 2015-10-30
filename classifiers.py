"""
Define classifiers which return (multiple) groups related to a gene record
"""

import sys, json, os.path

from filters import *
from utils import BaseClassifier, extract_gene_name, load_gene_annotations, annotate_seq_records


class GeneNameClassifier(BaseClassifier):
    """ Divide genes into known and unknown ones by checking if they have a name
    """
    data_file = 'dicty_primary_cds'

    def __init__(self):
        super().__init__()

        self.retroelements = set(['_'.join(gene.split('_')[:-1]) for gene in json.load(open(os.path.join(BaseClassifier.RESULTS_DIR, 'dicty_rte_list.json'), 'r'))])
        self.keywords = ['translation', 'transcription', 'stress response', 'cell cycle', 'rnai', 'cell signaling', 'splicing', 'cytokinesis', 'dna recombination', 'replication', 'transport', 'endocytosis']

    def get_groupname(self, record):
        """ Choose "best" annotation out of list of possible ones
        """
        gnames = []
        gnames.append('all')

        egn = extract_gene_name(record)
        if egn.endswith('_RTE'): gnames.append('rte')

        annos = ' | '.join(record.annotations['manual'])

        for kw in self.keywords:
            if kw.lower() in annos.lower():
                gnames.append(kw)

        if len(gnames) == 1: gnames.append('other')

        return gnames

    @staticmethod
    def preprocess(genes):
        fname = os.path.join(BaseClassifier.RESULTS_DIR, 'annotated_genes.json')

        if not os.path.isfile(fname):
            print('[ERROR] no annotations available')
            sys.exit(1)

        return annotate_seq_records(genes, json.load(open(fname, 'r')))

class RTEClassifier(BaseClassifier):
    """ Check if gene is part of a retroelement
    """
    data_file = 'DD_ComplexRepeats.fa'

    def __init__(self):
        super().__init__()

        self.rtes = ['_'.join(gene.split('_')[:-1]) for gene in json.load(open(os.path.join(BaseClassifier.RESULTS_DIR, 'dicty_rte_list.json'), 'r'))]

    def get_groupname(self, record):
        return [record.id.split()[0]]
