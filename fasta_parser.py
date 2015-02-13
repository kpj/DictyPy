import sys
import os.path
import collections

from Bio import SeqIO


def parse_filters():
    """ Return list of filters given in `filter.py`
    """
    import filters
    classes = [getattr(filters, x) for x in dir(filters) if isinstance(getattr(filters, x), type) and x != 'BaseFilter']
    return classes

class FastaParser(object):
    DATA_DIR = 'data'

    def __init__(self, fname):
        self.records = SeqIO.parse(os.path.join(FastaParser.DATA_DIR, fname), format='fasta')
        self.filters = parse_filters()

    def parse(self, classifier, verbose=True):
        """ Parse FastA file and group sequences into known and unknown ones.
            A sequence is known, if its name is given, i.e. doesn't start with 'DDB_G'
        """
        data = collections.defaultdict(list)
        skipped = []

        for r in self.records:
            skip = False
            for f in self.filters:
                if not f.skip and not f().apply(r):
                    skip = True
            if skip:
                skipped.append(r)
                continue

            for rule in classifier.rules:
                if rule['condition'](r):
                    data[rule['datafield']].append(r)
                    break

        data = dict(data)

        if verbose:
            for k, v in data.items():
                print('%i -> %s' % (len(v), k))
            print('%i skipped' % len(skipped))

        return data
