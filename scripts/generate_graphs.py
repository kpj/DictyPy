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


def do_binning(data, bin_width, bin_max=1):
    """ Bin data
    """
    counts, edges = np.histogram(data, bins=np.arange(0, bin_max+bin_width, bin_width))
    return counts.tolist(), edges.tolist()[1:]

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
        cur['counts'], cur['edges'] = do_binning(usage, bin_width)

        plot_data.append(cur)
    json.dump(plot_data, open('results/gene_codon_usages.json', 'w'))

    print('Plotting')
    subprocess.check_call(['Rscript', 'plotting/codon_usage_histogram_maker.R'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

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
        all_codon_num = dnana._count_codons(str(gene.seq))
        aa_num = sum([all_codon_num[codon] for codon in args])
        norm = aa_num * 1000 / len(gene.seq)
        return norm

    avg_codon_freqs = dnana.get_codon_freqs(genes)
    print(
        '  LYS freq: %f\n' % (avg_codon_freqs['AAA'] + avg_codon_freqs['AAG']) +
        '  GLU freq: %f\n' % (avg_codon_freqs['GAA'] + avg_codon_freqs['GAG']) +
        '  GLN freq: %f' % (avg_codon_freqs['CAA'] + avg_codon_freqs['CAG'])
    )

    # filter for genes
    low_CAA_genes = []
    for gene, codu in data.items():
        if not codu['CAA'] is None and codu['CAA'] < 0.9:
            lys_freq = (compute_norm(gene, 'AAA', 'AAG') / 1000) / (avg_codon_freqs['AAA'] + avg_codon_freqs['AAG'])
            glu_freq = (compute_norm(gene, 'GAA', 'GAG') / 1000) / (avg_codon_freqs['GAA'] + avg_codon_freqs['GAG'])
            gln_freq = (compute_norm(gene, 'CAA', 'CAG') / 1000) / (avg_codon_freqs['CAA'] + avg_codon_freqs['CAG'])

            low_CAA_genes.append((gene.id, extract_gene_name(gene), lys_freq, codu['AAA'], glu_freq, codu['GAA'], gln_freq, codu['CAA']))

    # store results
    with open('results/low_CAA_genes.csv', 'w') as fd:
        wrtr = csv.writer(fd)
        wrtr.writerow(['ID', 'name', 'LYS rel freq', 'CU: AAA', 'GLU rel freq', 'CU: GAA', 'GLN rel freq', 'CU: CAA'])

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

    bin_width = 1
    data = []
    for codon in ['AAA', 'CAA']:
        stretch_lens = []
        for gene in genes:
            stretch = get_longest_stretch(gene, codon)
            stretch_lens.append(len(stretch))

        counts, edges = do_binning(stretch_lens, bin_width, bin_max=max(stretch_lens))

        data.append({
            'codon': codon,
            'counts': counts,
            'edges': edges
        })

    with open('results/longest_stretches.json', 'w') as fd:
        json.dump(data, fd)

    print('Plotting')
    subprocess.check_call(['Rscript', 'plotting/stretch_histogram.R'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def main():
    """ Read and extract data
    """
    farser = FastaParser(sys.argv[1])
    genes = farser.parse()

    #handle_codon_usage(genes)
    store_low_CAA_genes(genes)
    #find_longest_stretch(genes)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: %s <fasta file>' % sys.argv[0])
        sys.exit(1)

    main()
