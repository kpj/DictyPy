"""
Figure out if rare codons aggregate at certain gene positions
"""

import sys
import itertools
import collections

from Bio.SeqUtils import CodonUsage

import numpy as np
import matplotlib.pylab as plt

sys.path.insert(0, '.')

from fasta_parser import FastaParser
from sequence_analyzer import DNAAnalyzer
from utils import find_all_positions
from gene_expression_analysis import extract_expression_levels, group_expression_levels


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
            res[c].extend([p/len(g.seq) for p in pos])
    return res

def plot_positions(positions, label):
    """ Plot given positions as histogram
    """
    pos = list(itertools.chain(*positions))

    n, bins, patches = plt.hist(
        pos, np.arange(0, 1.0001, 0.01),
        facecolor='khaki')

    plt.title('Rare codon position overview (%s)' % label)
    plt.xlabel('relative position in gene')
    plt.ylabel('count')

    plt.xlim((0, 1))

    plt.savefig('rarest_codon_positions.pdf')
    #plt.show()

def show_codon_share(pos_dict, resolution=2):
    """ Show which triplet appears how often at which position
    """
    # transform data
    tmp = collections.defaultdict(list)
    for trip, pos in pos_dict.items():
        for p in pos:
            tmp[round(p, resolution)].append(trip)

    # save data
    with open('rare_codon_positions.txt', 'w') as fd:
        for pos, codons in sorted(tmp.items()):
            count = collections.Counter(codons)
            mc = sorted(count.most_common())

            fd.write('%.2f\n' % pos)
            for cod, am in mc:
                fd.write(' %s: %d\n' % (cod, am))
            fd.write('\n')

def get_codu(genes, group):
    """ Compute codon usage for all genes or only for certain expression group if file is given
    """
    exprs = extract_expression_levels(sys.argv[2]) if len(sys.argv) == 3 else None
    groups = {'all': genes} if exprs is None else group_expression_levels(genes, exprs)
    select = 'all' if exprs is None else group

    dnana = DNAAnalyzer(strict=False)
    codu = dnana.get_avg_codon_usage(groups[select])

    return codu, select

def main():
    """ Generate overview
    """
    farser = FastaParser(sys.argv[1])
    genes = farser.parse()

    codu, label = get_codu(genes, 'weak')
    rarest = get_rare_codons(codu)
    pos = get_codon_positions(rarest.values(), genes)

    show_codon_share(pos)
    plot_positions(pos.values(), label)


if __name__ == '__main__':
    if len(sys.argv) != 2 and len(sys.argv) != 3:
        print('Usage: %s <fasta file> [expression file]' % sys.argv[0])
        sys.exit(1)

    main()
