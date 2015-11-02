"""
Compute single/average/cumulative codon usage for sequences
"""

import collections

from Bio.SeqUtils import CodonUsage

import utils


class DNAAnalyzer(object):
    """ The functions 'get_codon_usage' and '_count_codons' are adapted from http://biopython.org/DIST/docs/api/Bio.SeqUtils.CodonUsage-pysrc.html
    """
    def __init__(self, strict=True):
        self.strict = strict

    def get_avg_codon_usage(self, genes):
        """ Average codon usage over multiple genes
        """
        usages = []
        for gene in genes:
            try:
                codu = self.get_codon_usage(str(gene.seq))
                usages.append(codu)
            except TypeError as e:
                print('invalid gene %s (%s)' % (gene.id, str(e)))

        avg_codu = collections.defaultdict(int)
        cma_len_diff = collections.defaultdict(int)
        for i, usa in enumerate(usages):
            for codon, usage in usa.items():
                if usage is None:
                    cma_len_diff[codon] += 1
                    avg_codu[codon] = None if avg_codu[codon] == 0 else avg_codu[codon]
                else:
                    if avg_codu[codon] is None:
                        avg_codu[codon] = usage
                    else:
                        avg_codu[codon] = utils.next_cma(usage, i-cma_len_diff[codon], avg_codu[codon])

        return dict(avg_codu)

    def get_cum_codon_usage(self, genes):
        """ Get cumulative codon usage over multiple genes
        """
        cum_codu = collections.defaultdict(list)
        for gene in genes:
            try:
                codu = self.get_codon_usage(str(gene.seq))

                for codon, usage in codu.items():
                    if not usage is None:
                        cum_codu[codon].append(usage)
            except TypeError as e:
                print('invalid gene %s (%s)' % (gene.id, str(e)))

        return dict(cum_codu)

    def get_gene_codon_usages(self, genes):
        """ Get codon usage for individual genes
        """
        gene_codu = {}
        for gene in genes:
            try:
                codu = self.get_codon_usage(str(gene.seq))
                gene_codu[gene] = codu
            except TypeError as e:
                print('invalid gene %s (%s)' % (gene.id, str(e)))

        return gene_codu

    def get_codon_freqs(self, genes):
        """ Compute codon frequencies in given genes
        """
        # TODO: usa cma

        counts = []
        for gene in genes:
            seq = str(gene.seq)
            codon_count = self._count_codons(seq)
            tot_cod_num = len(seq) / 3
            codon_count.update((cod, cou / tot_cod_num) for cod, cou in codon_count.items())
            counts.append(codon_count)

        avg_counts = collections.defaultdict(int)
        for c in counts:
            for k, v in c.items():
                avg_counts[k] += v
        avg_counts.update((cod, s / len(counts)) for cod, s in avg_counts.items())

        return avg_counts

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
            codon = seq[i:i+3]
            if codon in codon_count:
                codon_count[codon] += 1
            else:
                if self.strict:
                    raise TypeError('illegal codon %s' % codon)

        return codon_count
