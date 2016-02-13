"""
Figure out if rare codons aggregate at certain gene positions
"""

import sys
import itertools
import collections

from Bio.SeqUtils import CodonUsage
import matplotlib.pylab as plt

sys.path.insert(0, '.')

from fasta_parser import FastaParser
from sequence_analyzer import DNAAnalyzer
from utils import find_all_positions


def get_rare_codons(codu):
    """ Extract rarest codon for each amino acid
    """
    res = {}
    for aa, codons in sorted(CodonUsage.SynonymousCodons.items()):
        min_codon = min(codons, key=lambda c: codu[c] if not codu[c] is None else float('inf'))
        res[aa] = min_codon
    return res

def get_codon_positions(codons, genes):
    """ Get codon positions in given genes
    """
    res = collections.defaultdict(list)
    for g in genes:
        for c in codons:
            pos = find_all_positions(str(g.seq), c, force_orf=True)
            res[c].extend(pos)
    return res

def plot_positions(positions):
    """ Plot given positions as histogram
    """
    pos = list(itertools.chain(*positions))

    n, bins, patches = plt.hist(
        pos, 100,
        facecolor='khaki')

    plt.title('Rare codon position overview')
    plt.xlabel('relative nucleotide position in gene')
    plt.ylabel('count')

    plt.savefig('rarest_codon_positions.pdf')
    #plt.show()

def main():
    """ Generate overview
    """
    farser = FastaParser(sys.argv[1])
    genes = farser.parse()

    dnana = DNAAnalyzer(strict=False)
    codu = dnana.get_avg_codon_usage(genes)

    rarest = get_rare_codons(codu)
    pos = get_codon_positions(rarest.values(), genes)

    plot_positions(pos.values())


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: %s <fasta file>' % sys.argv[0])
        sys.exit(1)

    main()
