import json, os.path

from filters import *
from utils import BaseClassifier, extract_gene_name, load_gene_annotations, annotate_seq_records


class GeneNameClassifier(BaseClassifier):
    """ Divide genes into known and unknown ones by checking if they have a name
    """
    requires_annotations = True

    def __init__(self):
        super().__init__()

        self.result = ['known_genes']
        self.rules = [
            {
                'condition': lambda record: extract_gene_name(record).startswith('DDB_G'),
                'datafield': 'unknown_genes'
            },
            {
                'condition': lambda record: True,
                'datafield': 'known_genes'
            }
        ]

    def get_groupname(self, record):
        """ Choose "best" annotation out of list of possible ones
        """
        keywords = json.load(open('keywords.json'))
        annos = ' | '.join(record.annotations['manual'])

        for kw in keywords:
            if kw.lower() in annos.lower():
                return kw

        return record.annotations['manual'][0]

    @staticmethod
    def preprocess(genes):
        fname = 'results/annotated_genes.json'

        if not os.path.isfile(fname):
            load_gene_annotations(genes, fname)

        return annotate_seq_records(genes, json.load(open(fname, 'r')))

class RTEClassifier(BaseClassifier):
    """ Check if gene is part of a retroelement
    """

    def __init__(self):
        super().__init__()
        
        self.rtes = []
        self.products = {}
        for rte in json.load(open('results/dicty_rte_list.json', 'r')):
            gene, product = rte.split(':')
            gene = '_'.join(gene.split('_')[:-1])

            self.rtes.append(gene)
            self.products[gene] = product

        def rte_detector(record):
            parts = record.id.split('|')

            if len(parts) == 2:
                return parts[1] in self.rtes

            return False

        self.result = ['rte']
        self.rules = [
            {
                'condition': rte_detector,
                'datafield': 'rte'
            }
        ]
        self.skip_filter = [FunctionalGroupFilter]

    def get_groupname(self, record):
        gene = record.id.split('|')[1]
        prod = self.products[gene]

        prod = prod.split()[0]

        return prod
