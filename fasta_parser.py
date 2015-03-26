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

    def parse(self, verbose=True):
        """ Parse FastA file and group sequences into known and unknown ones.
            A sequence is known, if its name is given, i.e. doesn't start with 'DDB_G'
        """
        filter_stats = collections.defaultdict(int)
        data = []
        rec_num = 0
        for r in self.records:
            rec_num += 1

            skip = False
            for f in self.filters:
                if not f.skip and not f().apply(r):
                    filter_stats[f.__name__] += 1
                    skip = True
            if skip: continue

            data.append(r)

        print('%d entries parsed' % rec_num)
        if len(filter_stats) > 0: print('Pre-Annotation filters:')
        for k, v in filter_stats.items(): print(' ', k, '->', v)

        print('%d entries remaining' % len(data))
        return data
