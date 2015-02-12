from fasta_parser import FastaParser
from blaster import Blaster
from sequence_analyzer import DNAAnalyzer


def main():
    farser = FastaParser('dicty_primary_cds')
    genes, ugenes = farser.parse()

    #blaster = Blaster(ugenes)
    #blaster.crawl()

    dnana = DNAAnalyzer()
    avg_codu = dnana.get_avg_codon_usage(genes)
    print(avg_codu)

if __name__ == '__main__':
    main()
