import json

from Bio.SeqUtils import CodonUsage

from sequence_analyzer import DNAAnalyzer
from utils import BaseFilter


class LengthFilter(BaseFilter):
    """ Skip all sequences which are shorter than specified threshold
    """

    def __init__(self):
        self.threshold = 250

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

class LysineAbundanceFilter(BaseFilter):
    def __init__(self):
        self.dnaa = DNAAnalyzer(strict=False)

    def apply(self, record):
        res = self.dnaa._count_codons(str(record.seq))
        lysin = res['AAA'] + res['AAG']
        norm = lysin * 1000 / len(record.seq)

        return norm > 76.6
