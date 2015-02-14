import json

from fasta_parser import FastaParser
from blaster import Blaster
from sequence_analyzer import DNAAnalyzer
from classifiers import GeneNameClassifier
from gene_annotator import GeneAnnotator


def main():
    farser = FastaParser('dicty_primary_cds')
    data = farser.parse(GeneNameClassifier())
    genes, ugenes = data['known_genes'], data['unknown_genes']

    #blaster = Blaster(ugenes)
    #blaster.crawl()

    #dnana = DNAAnalyzer()
    #avg_codu = dnana.get_avg_codon_usage(genes)
    #print(avg_codu)

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
    json.dump(foo, open('annotated_genes.json', 'w'))
    print(errors, 'errors')

if __name__ == '__main__':
    main()
