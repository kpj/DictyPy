from utils import BaseClassifier


class GeneNameClassifier(BaseClassifier):
    def __init__(self):
        self.rules = [
            {
                'condition': lambda record: self.extract_gene_name(record).startswith('DDB_G'),
                'datafield': 'unknown_genes'
            },
            {
                'condition': lambda record: True,
                'datafield': 'known_genes'
            }
        ]

    def extract_gene_name(self, record):
        """ Extract gene name from FastA entry description
        """
        return record.description.split('|')[3].split()[1] # regex anyone?
