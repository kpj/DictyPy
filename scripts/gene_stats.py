import re, json, operator

from fasta_parser import FastaParser
from sequence_analyzer import DNAAnalyzer
import utils


def find_special_AAA_freqs(genes):
    #id_filter = ['DDB0305421|DDB_G0276433', 'DDB0347990|DDB_G0289359', 'DDB0347948|DDB_G0270662', 'DDB0349097|DDB_G0279651', 'DDB0306784|DDB_G0293038', 'DDB0218505|DDB_G0283527', 'DDB0348150|DDB_G0285779', 'DDB0347690|DDB_G0286087'] # AAA=0
    #id_filter = ['DDB0230164|DDB_G0293360', 'DDB0186263|DDB_G0284929', 'DDB0232396|DDB_G0282423', 'DDB0238636|DDB_G0269008', 'DDB0234236|DDB_G0289721', 'DDB0229439|DDB_G0270122'] # AAA=1
    id_filter = ['DDB0348668|DDB_G0276223', 'DDB0307442|DDB_G0269954', 'DDB0307413|DDB_G0269350', 'DDB0308362|DDB_G0269090', 'DDB0216219|DDB_G0269132'] # long AAA=1

    def get_record(gene_id):
        for gene in genes:
            if gene.id == gene_id:
                return gene
        return None

    dnaa = DNAAnalyzer()
    for gid in id_filter:
        rec = get_record(gid)

        if not rec is None:
            print(rec.id)
            print(' ', 'gene length:', len(rec.seq))
            coco = dnaa._count_codons(str(rec.seq))
            print(' ', 'AAA:', coco['AAA'])
            print(' ', 'AAG:', coco['AAG'])
            print()

def find_longest_A_stretch(genes):
    threshold = 8
    regex = re.compile("(A+A)")

    data = []
    A_lengths = []
    for gene in genes:
        stretches = regex.findall(str(gene.seq))
        A_lengths.extend([len(s) for s in stretches])

        longest_stretch = max(stretches)
        if len(longest_stretch) > threshold:
            data.append({
                'gene_id': gene.id,
                'gene_name': utils.extract_gene_name(gene),
                'gene_len': len(gene.seq),
                'stretch': longest_stretch,
                'stretch_len': len(longest_stretch)
            })

    data = list(reversed(sorted(data, key=operator.itemgetter('stretch'))))

    with open('longest_A_stretches.json', 'w') as fd:
        json.dump(data, fd)

    with open('all_A_stretch_lengths.json', 'w') as fd:
        json.dump(A_lengths, fd)

def main():
    farser = FastaParser('dicty_primary_cds')
    genes = farser.parse()

    #find_special_AAA_freqs(genes)
    find_longest_A_stretch(genes)


if __name__ == '__main__':
    main()
