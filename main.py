import json, os.path, subprocess
from pprint import pprint

from fasta_parser import FastaParser
from sequence_analyzer import DNAAnalyzer
from classifiers import GeneNameClassifier
from gene_annotator import GeneAnnotator
from gene_grouper import AnnotationGrouper
from utils import parse_filters


def annotate_genes(genes, fname):
    """ Save annotated genes to given filename
    """
    foo = {}
    errors = 0
    ganno = GeneAnnotator()
    for gene in genes:
        try:
            anno = ganno.get_direct_annotation(gene)
            foo[gene.id] = anno
        except:
            print('Error:', gene, '\n')
            errors += 1
    json.dump(foo, open(fname, 'w'))
    print(errors, 'errors')

def group_genes(GeneGrouper, genes, fname_in, fname_out):
    """ Group genes given in filename and save results elsewhere
    """
    gegro = GeneGrouper()
    rec_list = gegro.transform(genes, json.load(open(fname_in, 'r')))
    groups = gegro.group(rec_list)

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

def apply_procedure(Classifier, Grouper):
    farser = FastaParser('dicty_primary_cds')
    farser.parse(Classifier)
    genes = farser.get_result()

    if Classifier.requires_annotations:
        if not os.path.isfile('results/annotated_genes.json'):
            annotate_genes(genes, 'results/annotated_genes.json')

    group_genes(Grouper, genes, 'results/annotated_genes.json', 'results/grouped_genes.json')
    plot_grouped_genes()

def main():
    apply_procedure(GeneNameClassifier, AnnotationGrouper)


if __name__ == '__main__':
    main()
