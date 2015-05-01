import collections, csv, json, os.path, subprocess
from pprint import pprint

import numpy as np

from fasta_parser import FastaParser
from sequence_analyzer import DNAAnalyzer
from classifiers import GeneNameClassifier, RTEClassifier
from gene_grouper import GeneGrouper
from utils import parse_filters


def group_genes(Classifier, genes, fname_out):
    """ Group genes given in filename and save results elsewhere
    """
    gegro = GeneGrouper(Classifier)
    genes = Classifier.preprocess(genes)
    groups = gegro.group(genes)

    foo = []
    filter_stats = collections.defaultdict(int)
    dnana = DNAAnalyzer(strict=False)
    for group_name, group_genes in groups.items():
        # apply post-annotation filters
        filters = parse_filters(post_annotation=True)
        genes = []
        for gene in group_genes:
            skip = False
            for f in filters:
                if not f.skip and not f().apply(gene):
                    filter_stats[f.__name__] += 1
                    skip = True
            if skip: continue

            genes.append(gene)

        if len(genes) == 0: continue

        # compute codon usage
        cum_codu = dnana.get_cum_codon_usage(genes)
        foo.append({
            'group': group_name,
            'cumulative_codon_usage': cum_codu
        })

    if len(filter_stats) > 0: print('Post-Annotation filters:')
    for k, v in filter_stats.items(): print(' ', k, '->', v)

    json.dump(foo, open(os.path.join(Classifier.RESULTS_DIR, fname_out), 'w'))
    #pprint(foo)

def plot_grouped_genes():
    print('Plotting...')
    subprocess.check_call(['Rscript', 'plotter.R'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def save_to_csv():
    aa_dict = {
        'V': ['GTT', 'GTA', 'GTC', 'GTG'],
        'A': ['GCT', 'GCC', 'GCG', 'GCA'],
        'C': ['TGT', 'TGC'],
        'Y': ['TAT', 'TAC'],
        'W': ['TGG'],
        'N': ['AAT', 'AAC'],
        'T': ['ACT', 'ACC', 'ACG', 'ACA'],
        'Q': ['CAG', 'CAA'],
        'L': ['CTC', 'CTA', 'CTT', 'CTG', 'TTA', 'TTG'],
        'S': ['AGC', 'TCT', 'TCC', 'TCG', 'AGT', 'TCA'],
        'E': ['GAA', 'GAG'],
        'R': ['CGC', 'CGA', 'CGT', 'AGG', 'CGG', 'AGA'],
        'H': ['CAC', 'CAT'],
        'F': ['TTC', 'TTT'],
        'D': ['GAT', 'GAC'],
        'M': ['ATG'],
        'K': ['AAA', 'AAG'],
        '*': ['TAG', 'TAA', 'TGA'],
        'P': ['CCA', 'CCC', 'CCT', 'CCG'],
        'I': ['ATT', 'ATC', 'ATA'],
        'G': ['GGT', 'GGG', 'GGA', 'GGC'],
    }

    with open('results/grouped_genes.json', 'r') as fd:
        content = json.load(fd)

    with open('out.csv', 'w') as fd:
        writer = csv.writer(fd)
        writer.writerow(['group', 'amino_acid', 'codon', 'codon_usage'])

        for entry in content:
            group = entry['group']
            cuco = entry['cumulative_codon_usage']

            for aa, cods in aa_dict.items():
                for c in cods:
                    codu = np.mean(cuco[c])

                    writer.writerow([group, aa, c, codu])

def apply_procedure(Classifier):
    farser = FastaParser(Classifier.data_file)
    genes = farser.parse()

    group_genes(Classifier, genes, 'grouped_genes.json')

    #plot_grouped_genes()
    save_to_csv()

def main():
    apply_procedure(GeneNameClassifier)
    #apply_procedure(RTEClassifier)


if __name__ == '__main__':
    main()
