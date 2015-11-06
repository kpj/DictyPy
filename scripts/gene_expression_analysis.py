"""
Relate gene expression levels to codon usage patterns
"""

import sys, csv
import pprint
import collections

sys.path.insert(0, '.')

from fasta_parser import FastaParser
from sequence_analyzer import DNAAnalyzer


def extract_expression_levels(fname):
    """ Extract normalized gene expression levels from CSV
    """
    def get_norm_expr(row):
        """ Extract normalized expression value from given row
        """
        cur = [int(e) for e in row[1:4]]
        return sum(cur) / len(cur)

    # read raw data
    with open(fname, 'r') as fd:
        creader = csv.reader(fd, delimiter=';')
        data = {}

        creader.__next__() # skip header
        for row in creader:
            gene = row[0]
            expr = get_norm_expr(row)

            data[gene] = expr

    # normalize levels
    max_level = max(data.values())
    data.update((gene, count / max_level) for gene, count in data.items())

    return data

def group_expression_levels(genes, exprs):
    """ Group genes according to their expression levels
    """
    # [(label, expr_threshold)]
    group_spec = [
        ('lowest', 0.1),
        ('weak', 0.3),
        ('middle', 0.6),
        ('strong', 0.9),
        ('highest', 1)
    ]

    groups = collections.defaultdict(list)
    parsed_num = len(genes)
    for gene in genes:
        try:
            name = gene.id.split('|')[1]
            expr = exprs[name]
        except (IndexError, KeyError):
            parsed_num -= 1
            continue

        for label, thres in group_spec:
            if expr < thres:
                groups[label].append(gene)
                break

    print(parsed_num, '/', len(genes), 'genes parsed')
    for k, v in groups.items(): print(' ', k, '->', len(v))

    return dict(groups)

def generate_codon_usage_summary(groups, out_fname):
    """ Save codon usage tables per group
    """
    dnana = DNAAnalyzer(strict=False)
    with open(out_fname, 'w') as fd:
        for label, genes in groups.items():
            codu = dnana.get_avg_codon_usage(genes)

            fd.write(label + '\n')
            pprint.pprint(codu, fd)
            fd.write('\n')

def main():
    """ Read and extract data
    """
    farser = FastaParser(sys.argv[1])
    genes = farser.parse()

    exprs = extract_expression_levels('data/gene_expression.csv')
    groups = group_expression_levels(genes, exprs)
    generate_codon_usage_summary(groups, 'results/grouped_gene_expressions.txt')


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: %s <fasta file>' % sys.argv[0])
        sys.exit(1)

    main()
