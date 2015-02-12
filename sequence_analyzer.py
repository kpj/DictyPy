import collections

from Bio.SeqUtils import CodonUsage

import utils


class DNAAnalyzer(object):
    """ The functions 'get_codon_usage' and '_count_codons' are adapted from http://biopython.org/DIST/docs/api/Bio.SeqUtils.CodonUsage-pysrc.html
    """
    def get_avg_codon_usage(self, genes):
        """ Average codon usage over multiple genes
        """

        #
        # Note: this doesn't give reasonable results for now
        #

        usages = []
        for gene in genes:
            try:
                codu = self.get_codon_usage(str(gene.seq))
                usages.append(codu)
            except TypeError:
                print('invalid gene %s' % gene.id)

        avg_codu = collections.defaultdict(int)
        for i, usa in enumerate(usages):
            for codon, usage in usa.items():
                if usage is None:
                    avg_codu[codon] = None if avg_codu[codon] == 0 else avg_codu[codon]
                else:
                    avg_codu[codon] = utils.next_cma(usage, i, 0 if avg_codu[codon] is None else avg_codu[codon])

        return dict(avg_codu)

    def get_codon_usage(self, seq):
        """ Calculate codon usage in given sequence
        """
        codon_usage = {}

        seq = seq.upper()
        codon_count = self._count_codons(seq)

        for aa in CodonUsage.SynonymousCodons:
            total = 0.0
            rcsu = []
            codons = CodonUsage.SynonymousCodons[aa]

            for codon in codons:
                total += codon_count[codon]

            for codon in codons:
                if total == 0:
                    codon_usage[codon] = None
                else:
                    codon_usage[codon] = codon_count[codon] / total

        return codon_usage

    def _count_codons(self, seq):
        """ Count codon occurences in given sequence
        """
        codon_count = CodonUsage.CodonsDict.copy()

        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            if codon in codon_count:
                codon_count[codon] += 1
            else:
                raise TypeError('illegal codon %s' % codon)

        return codon_count
