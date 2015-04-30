import sys
import difflib

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from progressbar import ProgressBar

from fasta_parser import FastaParser


THRESHOLD = 0.9

if len(sys.argv) != 4:
    print('Usage: %s <input1.fa> <input2.fa> <output.fa>' % sys.argv[0])
    sys.exit(1)

def intersect(gene_vec1, gene_vec2):
    """ Find "approximate" intersection of two gene sets
    """
    inter = []

    pbar = ProgressBar(maxval=len(gene_vec1)*len(gene_vec2))
    pbar.start()
    counter = 0
    for g1 in gene_vec1:
        for g2 in gene_vec2:
            rat = difflib.SequenceMatcher(
                a=str(g1.seq).lower(),
                b=str(g2.seq).lower()
            ).ratio()

            if rat > THRESHOLD:
                inter.append((g1, g2)) # add both as they might be slightly different

            pbar.update(counter)
            counter += 1
    pbar.finish()

    return inter

def writeFastA(struct):
    with open(sys.argv[3], 'w') as fd:
        for g1, g2 in struct:
            g1.annotations['origin'] = sys.argv[1]
            g2.annotations['origin'] = sys.argv[2]

            SeqIO.write(g1, fd, 'fasta')
            SeqIO.write(g2, fd, 'fasta')

def main():
    genes1 = FastaParser(sys.argv[1]).parse()
    genes2 = FastaParser(sys.argv[2]).parse()

    res = intersect(genes1, genes2)
    writeFastA(res)


if __name__ == '__main__':
    main()
