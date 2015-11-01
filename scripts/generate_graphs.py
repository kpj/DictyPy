"""
Generate some histograms related to codon usage and repetetive codon sequence lengths
"""

import sys, json, csv
import subprocess, re

import numpy as np

sys.path.insert(0, '.') # evil hack for importing scripts from parent directory

from fasta_parser import FastaParser
from sequence_analyzer import DNAAnalyzer
from utils import extract_gene_name


def handle_codon_usage(genes):
    """ Generate codon usage histograms
    """
    def extract(marker, data):
        """ Extract marked codon usage from each gene
        """
        out = []
        for gene, codu in data.items():
            if not codu[marker] is None:
                out.append(codu[marker])
        return out

    print('Computing codon statistics')
    dnana = DNAAnalyzer(strict=False)
    data = dnana.get_gene_codon_usages(genes)

    plot_data = []
    bin_width = 0.01
    for marker in ['AAA', 'GAA', 'CAA']:
        cur = {}
        usage = extract(marker, data)

        cur['marker'] = marker

        counts, edges = np.histogram(usage, bins=np.arange(0, 1+bin_width, bin_width))
        cur['counts'] = counts.tolist()
        cur['edges'] = edges.tolist()[1:]

        plot_data.append(cur)
    json.dump(plot_data, open('results/gene_codon_usages.json', 'w'))

    print('Plotting')
    subprocess.check_call(['Rscript', 'plotting/histogram_maker.R'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def store_low_CAA_genes(genes):
    """ Create list of genes where CAA usage < 0.9
    """
    # compute codon usage
    print('Computing codon statistics')
    dnana = DNAAnalyzer(strict=False)
    data = dnana.get_gene_codon_usages(genes)

    def compute_norm(gene, *args):
        """ Compute normalized occurrence frequency of aa
        """
        res = dnana._count_codons(str(gene.seq))
        aa = sum([res[blub] for blub in args])
        norm = aa * 1000 / len(gene.seq)
        return norm

    # filter for genes
    low_CAA_genes = []
    for gene, codu in data.items():
        if not codu['CAA'] is None and codu['CAA'] < 0.9:
            lys_freq = compute_norm(gene, 'AAA', 'AAG')
            glu_freq = compute_norm(gene, 'GAA', 'GAG')
            gln_freq = compute_norm(gene, 'CAA', 'CAG')

            low_CAA_genes.append((gene.id, extract_gene_name(gene), lys_freq, codu['AAA'], glu_freq, codu['GAA'], gln_freq, codu['CAA']))

    # store results
    with open('results/low_CAA_genes.csv', 'w') as fd:
        wrtr = csv.writer(fd)
        wrtr.writerow(['ID', 'name', 'LYS freq', 'CU: AAA', 'GLU freq', 'CU: GAA', 'GLN freq', 'CU: CAA'])

        for entry in low_CAA_genes:
            wrtr.writerow(entry)

def find_longest_stretch(genes):
    """ Find longest stretches
    """
    def get_longest_stretch(gene, codon):
        """ Find longtest stretch in given gene
        """
        pat = re.compile(r'((?:' + codon + ')+)')
        stretches = pat.findall(str(gene.seq))
        longest = max(stretches, key=len) if len(stretches) > 0 else ''
        return longest

    data = {'aaa_len': [], 'caa_len': []}
    for gene in genes:
        aaa_stretch = get_longest_stretch(gene, r'AAA')
        caa_stretch = get_longest_stretch(gene, r'CAA')

        data['aaa_len'].append(len(aaa_stretch))
        data['caa_len'].append(len(caa_stretch))

    with open('results/longest_stretches.json', 'w') as fd:
        json.dump(data, fd)


def main():
    """ Read and extract data
    """
    farser = FastaParser(sys.argv[1])
    genes = farser.parse()

    handle_codon_usage(genes)
    #store_low_CAA_genes(genes)
    #find_longest_stretch(genes)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: %s <fasta file>' % sys.argv[0])
        sys.exit(1)

    main()
