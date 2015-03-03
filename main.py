import json, subprocess
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
    dnana = DNAAnalyzer()
    for group_name, group_genes in groups.items():
        # apply post-annotation filters
        filters = parse_filters(post_annotation=True)
        genes = []
        for gene in group_genes:
            skip = False
            for f in filters:
                if not f.skip and not f().apply(gene):
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
    json.dump(foo, open(fname_out, 'w'))
    #pprint(foo)

def plot_grouped_genes():
    subprocess.check_call(['Rscript', 'plotter.R'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def apply_procedure(Classifier):
    farser = FastaParser('dicty_primary_cds')
    farser.parse(Classifier)
    genes = farser.get_result()

    group_genes(Classifier, genes, 'results/grouped_genes.json')
    plot_grouped_genes()

def main():
    #apply_procedure(GeneNameClassifier)
    apply_procedure(RTEClassifier)


if __name__ == '__main__':
    main()
