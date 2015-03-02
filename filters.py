import json

from Bio.SeqUtils import CodonUsage

from utils import BaseFilter


class LengthFilter(BaseFilter):
    """ Skip all sequences which are shorter than specified threshold
    """

    def __init__(self):
        self.threshold = 100

    def apply(self, record):
        return len(record.seq) > self.threshold

class MultipleStopCodonFilter(BaseFilter):
    """ Skip all sequences which contain more than one STOP codon
    """
    skip = True

    def __init__(self):
        self.stop_codons = CodonUsage.SynonymousCodons['STOP']

    def apply(self, record):
        cods = [record.seq[i:i+3] for i in range(0, len(record.seq), 3)]
        return sum([cods.count(cod) for cod in self.stop_codons]) <= 1

class FunctionalGroupFilter(BaseFilter):
    """ Skip all sequences which are not annotated with one of the given keywords
    """
    post_annotation = True

    def __init__(self):
        self.keywords = set(json.load(open('keywords.json')))

    def apply(self, record):
        return record.annotations['group_name'] in self.keywords
