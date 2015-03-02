from utils import BaseClassifier, extract_gene_name


class GeneNameClassifier(BaseClassifier):
    requires_annotations = True

    def __init__(self):
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
