from utils import BaseClassifier, extract_gene_name


class GeneNameClassifier(BaseClassifier):
    def __init__(self):
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
