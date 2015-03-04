import json, os.path

from filters import *
from utils import BaseClassifier, extract_gene_name, load_gene_annotations, annotate_seq_records


class GeneNameClassifier(BaseClassifier):
    """ Divide genes into known and unknown ones by checking if they have a name
    """
    data_file = 'dicty_primary_cds'

    def __init__(self):
        super().__init__()

        self.result = ['known_genes', 'rte']
        self.rules = [
            {
                'condition': lambda record: extract_gene_name(record).endswith('_RTE'),
                'datafield': 'rte'
            },
            {
                'condition': lambda record: extract_gene_name(record).startswith('DDB_G'),
                'datafield': 'unknown_genes'
            },
            {
                'condition': lambda record: True,
                'datafield': 'known_genes'
            }
        ]

        self.retroelements = set(['_'.join(gene.split('_')[:-1]) for gene in json.load(open('results/dicty_rte_list.json', 'r'))])
        self.keywords = ['translation', 'transcription', 'stress response', 'cell cycle', 'rnai', 'cell signaling', 'splicing', 'cytokinesis', 'dna recombination']

    def get_groupname(self, record):
        """ Choose "best" annotation out of list of possible ones
        """
        egn = extract_gene_name(record)
        if egn.endswith('_RTE'): return 'rte'

        annos = ' | '.join(record.annotations['manual'])

        for kw in self.keywords:
            if kw.lower() in annos.lower():
                return kw

        return 'other' #record.annotations['manual'][0]

    @staticmethod
    def preprocess(genes):
        fname = 'results/annotated_genes.json'

        if not os.path.isfile(fname):
            load_gene_annotations(genes, fname)

        return annotate_seq_records(genes, json.load(open(fname, 'r')))

class RTEClassifier(BaseClassifier):
    """ Check if gene is part of a retroelement
    """
    data_file = 'DD_ComplexRepeats.fa'

    def __init__(self):
        super().__init__()

        self.rtes = ['_'.join(gene.split('_')[:-1]) for gene in json.load(open('results/dicty_rte_list.json', 'r'))]
        self.result = ['rte']
        self.rules = [
            {
                'condition': lambda record: True,
                'datafield': 'rte'
            }
        ]

    def get_groupname(self, record):
        return record.id.split()[0]
