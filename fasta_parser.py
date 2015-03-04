import sys
import os.path
import collections

from Bio import SeqIO

from utils import parse_filters


class FastaParser(object):
    DATA_DIR = 'data'

    def __init__(self, fname):
        self.records = SeqIO.parse(os.path.join(FastaParser.DATA_DIR, fname), format='fasta')
        self.filters = parse_filters()

        self.data = None
        self.classifier = None

    def get_result(self):
        """ Return result(s) as specified by classifier
        """
        if self.classifier is None:
            raise RuntimeError('No classifier specified, did you forget to call parse(..) first?')
        if self.data is None:
            raise RuntimeError('No data generated, maybe the parse(..) function screwed up?')

        res = []
        for key in self.classifier.result:
            res.extend(self.data[key])

        return res

    def parse(self, Classifier, verbose=True):
        """ Parse FastA file and group sequences into known and unknown ones.
            A sequence is known, if its name is given, i.e. doesn't start with 'DDB_G'
        """
        self.classifier = Classifier()
        self.data = collections.defaultdict(list)
        skipped = []

        for r in self.records:
            skip = False
            for f in self.filters:
                if not f.skip and not f().apply(r):
                    skip = True
            if skip:
                skipped.append(r)
                continue

            for rule in self.classifier.rules:
                if rule['condition'](r):
                    self.data[rule['datafield']].append(r)
                    break

        self.data = dict(self.data)

        if verbose:
            for k, v in self.data.items():
                print('%i -> %s' % (len(v), k))
            print('%i skipped' % len(skipped))
