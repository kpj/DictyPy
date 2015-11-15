"""
Print table with simple codon usage overview
"""

import sys

from Bio.SeqUtils import CodonUsage

sys.path.insert(0, '.')

from fasta_parser import FastaParser
from sequence_analyzer import DNAAnalyzer
from gene_expression_analysis import extract_expression_levels, group_expression_levels


def output_data(group_codu, out_stream=sys.stdout):
    """ Output nicely formatted data
    """
    def frmt(codu, head):
        """ Format codon usage float according to header
        """
        if codu is None:
            return '-' * len(head)
        else:
            flt = '{:.{width}f}'.format(codu, width=len(head)-2)
            return flt

    labels = sorted(group_codu.keys())
    space = '   '

    for aa, codons in sorted(CodonUsage.SynonymousCodons.items()):
        out_stream.write(space.join([aa] + labels) + '\n')
        for c in codons:
            data = space.join([c] + [frmt(group_codu[l][c], l) for l in labels])
            out_stream.write('%s\n' % data)
        out_stream.write('\n')

def get_codu_per_group(groups):
    """ Compute codon usage per group
    """
    group_codu = {}
    dnana = DNAAnalyzer(strict=False)

    for label, genes in groups.items():
        group_codu[label] = dnana.get_avg_codon_usage(genes)

    return group_codu

def main():
    """ Generate overview
    """
    farser = FastaParser(sys.argv[1])
    genes = farser.parse()

    exprs = extract_expression_levels(sys.argv[2]) if len(sys.argv) == 3 else None
    if exprs is None:
        groups = {'all': genes}
    else:
        groups = group_expression_levels(genes, exprs)

    group_codu = get_codu_per_group(groups)
    with open('results/codon_usage_table.txt', 'w') as fd:
        output_data(group_codu, fd)


if __name__ == '__main__':
    if not len(sys.argv) in [2, 3]:
        print('Usage: %s <fasta file> [expr file]' % sys.argv[0])
        sys.exit(1)

    main()
