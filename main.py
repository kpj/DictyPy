from fasta_parser import FastaParser
from blaster import Blaster
from sequence_analyzer import DNAAnalyzer
from classifiers import GeneNameClassifier


def main():
    farser = FastaParser('dicty_primary_cds')
    data = farser.parse(GeneNameClassifier())
    genes, ugenes = data['known_genes'], data['unknown_genes']

    #blaster = Blaster(ugenes)
    #blaster.crawl()

    dnana = DNAAnalyzer()
    avg_codu = dnana.get_avg_codon_usage(genes)
    print(avg_codu)

if __name__ == '__main__':
    main()
