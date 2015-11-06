"""
Generate some histograms related to codon usage and repetetive codon sequence lengths
"""

import sys, json, csv
import subprocess

import numpy as np
import regex as re

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
    """ Generate 2D-Histogram of stretch length and relative position in gene
    """
    def get_stretches(gene, codon):
        """ Find stretches in ORF of given gene and relative position
        """
        pat = re.compile(r'((?:' + codon + ')+)')
        stretches = pat.finditer(str(gene.seq), overlapped=True)

        strtchs = []
        pstns = []
        for stretch in stretches:
            if stretch.start() % 3 == 0:
                strtchs.append(stretch.group())
                pstns.append(stretch.start() / len(gene.seq))

        return strtchs, pstns

    data = []
    for codon in ['AAA', 'CAA', 'AAT']:
        stretch_lens = []
        stretch_pos = []

        for gene in genes:
            stretches, rel_pos = get_stretches(gene, codon)
            stretch_lens.extend([len(stretch) / 3. for stretch in stretches])
            stretch_pos.extend(rel_pos)

        # make 2D-Histogram
        xedges = np.arange(0, max(stretch_lens)+1, 1)
        yedges = np.arange(0, 1+0.01, 0.01)

        counts, xedges, yedges = np.histogram2d(stretch_lens, stretch_pos, bins=(xedges, yedges))
        xedges, yedges = xedges[1:], yedges[1:]

        coords = []
        for i, xe in enumerate(xedges):
            for j, ye in enumerate(yedges):
                coords.append({
                    'x': xe,
                    'y': ye,
                    'z': counts[i, j]
                })

        data.append({
            'codon': codon,
            'data': coords
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
    #store_low_CAA_genes(genes)
    find_longest_stretch(genes)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: %s <fasta file>' % sys.argv[0])
        sys.exit(1)

    main()
