import os.path

from io import StringIO

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

from utils import frint


class Blaster(object):
    """ Blast 'em up
    """

    CACHE_DIR = 'cache'
    E_VALUE_THRESHOLD = 0.04

    def __init__(self, genes):
        self.genes = genes

    def crawl(self):
        for r in self.genes:
            self.blast(r)

    def blast(self, record):
        """ Use BLAST to find similar sequences in order to gain additional information
        """

        frint('Blasting %s' % record.id, end='')
        fname = os.path.join(Blaster.CACHE_DIR, 'blast_%s.dat' % record.id)

        if os.path.isfile(fname):
            frint(' (cached)...', end=' ')
            with open(fname, 'r') as fd:
                txt = StringIO(fd.read())
        else:
            frint('...', end=' ')
            result_handle = NCBIWWW.qblast('blastn', 'nt', record.format('fasta'))
            blast_results = result_handle.read()

            with open(fname, 'w') as fd:
                fd.write(blast_results)

            txt = StringIO(blast_results)

        results = NCBIXML.parse(txt)
        frint('Done')

        for r in results:
            for alignment in r.alignments:
                for hsp in alignment.hsps:
                    e_value = hsp.expect
                    if e_value > Blaster.E_VALUE_THRESHOLD:
                        continue

                    print(alignment.title)
