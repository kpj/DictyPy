import json
from pprint import pprint

from fasta_parser import FastaParser
from blaster import Blaster
from sequence_analyzer import DNAAnalyzer
from classifiers import GeneNameClassifier
from gene_annotator import GeneAnnotator
from gene_grouper import GeneGrouper


def annotate_genes(genes, fname):
    """ Save annotated genes to given filename
    """
    foo = {}
    errors = 0
    ganno = GeneAnnotator()
    for gene in genes:
        try:
            anno = ganno.get_direct_annotation(gene)
        except:
            print('Error:', gene, '\n')
            errors += 1
        foo[gene.id] = anno
    json.dump(foo, open(fname, 'w'))
    print(errors, 'errors')

def group_genes(genes, fname_in, fname_out):
    """ Group genes given in filename and save results elsewhere
    """
    gego = GeneGrouper()
    transf = gego.transform(genes, json.load(open(fname_in, 'r')))
    groups = gego.group(transf)

    foo = []
    dnana = DNAAnalyzer()
    for group_name, group_genes in groups.items():
        avg_codu = dnana.get_avg_codon_usage(group_genes)
        foo.append({
            'group': group_name,
            'average_codon_usage': avg_codu
        })
    json.dump(foo, open(fname_out, 'w'))
    pprint(foo)

def main():
    farser = FastaParser('dicty_primary_cds')
    data = farser.parse(GeneNameClassifier())
    genes, ugenes = data['known_genes'], data['unknown_genes']

    #blaster = Blaster(ugenes)
    #blaster.crawl()

    #annotate_genes(genes, 'results/annotated_genes.json')
    group_genes(genes, 'results/annotated_genes.json', 'results/grouped_genes.json')


if __name__ == '__main__':
    main()
