"""
Parse FastA files and apply fitting filters
"""

import sys
import os.path
import collections

from Bio import SeqIO

from utils import parse_filters


class FastaParser(object):
    DATA_DIR = 'data'

    def __init__(self, fname):
        """ If `fname` is not found, look for it in `DATA_DIR` before failing
        """
        if os.path.isfile(fname):
            full_path = fname
        else:
            full_path = os.path.join(FastaParser.DATA_DIR, fname)

        self.records = SeqIO.parse(full_path, format='fasta')
        self.filters = parse_filters()

    def parse(self, verbose=True):
        """ Parse FastA file and group sequences into known and unknown ones.
            A sequence is known, if its name is given, i.e. doesn't start with 'DDB_G'
        """
        filter_stats = collections.defaultdict(int)
        data = []
        cod_lens = []
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
            cod_lens.append(len(r.seq)/3)

        print('%d entries parsed' % rec_num)
        if len(filter_stats) > 0:
            print('Pre-Annotation filters:')
            for k, v in filter_stats.items():
                print(' ', k, '->', v)
        else:
            print(' ', 'No filters applied')

        print('%d entries remaining' % len(data))
        print(' > min/avg/max gene length:', '%d/%d/%d' % \
            (min(cod_lens), sum(cod_lens)/len(cod_lens), max(cod_lens)))

        return data
