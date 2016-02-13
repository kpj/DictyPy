"""
Print amino acids and respective codon usage
"""

import sys

from Bio.SeqUtils import CodonUsage

sys.path.insert(0, '.')

from fasta_parser import FastaParser
from sequence_analyzer import DNAAnalyzer


def output_data(codu, out_stream=sys.stdout):
    """ Output nicely formatted data
    """
    for aa, codons in sorted(CodonUsage.SynonymousCodons.items()):
        out_stream.write(aa + '\n')
        for c in codons:
            out_stream.write('%s   %s\n' % (c, str(codu[c])))
        out_stream.write('\n')

def main():
    """ Generate overview
    """
    farser = FastaParser(sys.argv[1])
    genes = farser.parse()

    dnana = DNAAnalyzer(strict=False)
    codu = dnana.get_avg_codon_usage(genes)

    with open('results/plain_codon_usage_table.txt', 'w') as fd:
        output_data(codu, fd)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: %s <fasta file>' % sys.argv[0])
        sys.exit(1)

    main()
