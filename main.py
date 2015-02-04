from utils import frint
from fasta_parser import FastaParser
from blaster import Blaster


def main():
    farser = FastaParser('dicty_primary_cds')
    genes, ugenes = farser.parse()

    blaster = Blaster(ugenes)
    blaster.crawl()

if __name__ == '__main__':
    main()
