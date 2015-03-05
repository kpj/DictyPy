import collections, json, os.path, subprocess
from pprint import pprint

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
    subprocess.check_call(['Rscript', 'plotter.R'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def apply_procedure(Classifier):
    farser = FastaParser(Classifier.data_file)
    genes = farser.parse()

    group_genes(Classifier, genes, 'grouped_genes.json')
    plot_grouped_genes()

def main():
    apply_procedure(GeneNameClassifier)
    #apply_procedure(RTEClassifier)


if __name__ == '__main__':
    main()
