import sys
import os.path

from io import StringIO

from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML


def frint(s, **kwargs):
    print(s, **kwargs)
    sys.stdout.flush()

def parse_filters():
    """ Return list of filters given in `filter.py`
    """
    import filters
    classes = [getattr(filters, x) for x in dir(filters) if isinstance(getattr(filters, x), type) and x != 'BaseFilter']
    return classes

def parse_fasta_file(fname, verbose=True):
    """ Parse FastA file and group sequences into known and unknown ones.
        A sequence is known, if its name is given, i.e. doesn't start with 'DDB_G'
    """
    def extract_gene_name(desc):
        """ Extract gene name from FastA entry description
        """
        return desc.split('|')[3].split()[1] # regex anyone?

    filters = parse_filters()

    known_genes = []
    skipped_genes = []
    unknown_genes = []

    records = SeqIO.parse(os.path.join('data', fname), format='fasta')
    for r in records:
        name = extract_gene_name(r.description)

        skip = False
        for f in filters:
            if not f().apply(r):
                skip = True
        if skip:
            skipped_genes.append(r)
            continue

        if name.startswith('DDB_G'):
            unknown_genes.append(r)
        else:
            known_genes.append(r)

    if verbose:
        frint('%i known genes' % len(known_genes))
        frint('%i unknown genes' % len(unknown_genes))
        frint('%i skipped genes' % len(skipped_genes))

    return known_genes, unknown_genes

def check_unknown_genes(unknown_genes):
    """ Use BLAST to find similar sequences in order to gain additional information
    """
    E_VALUE_THRESHOLD = 0.04

    def blaster(record):
        frint('Blasting %s' % record.id, end='')
        fname = 'blast_%s.dat' % record.id

        if os.path.isfile(fname):
            frint(' (cached)...', end=' ')
            with open(fname, 'r') as fd:
                txt = StringIO(fd.read())
        else:
            frint('...', end=' ')
            result_handle = NCBIWWW.qblast('blastn', 'nt', record.format('fasta'))
            blast_results = result_handle.read()

            with open(fname, 'w') as fd:
                fd.write(blast_results)

            txt = StringIO(blast_results)

        results = NCBIXML.parse(txt)
        frint('Done')

        for r in results:
            for alignment in r.alignments:
                for hsp in alignment.hsps:
                    e_value = hsp.expect
                    if e_value > E_VALUE_THRESHOLD:
                        continue

                    print(alignment.title)


    for r in unknown_genes:
        blaster(r)


def main():
    genes, ugenes = parse_fasta_file('dicty_primary_cds')
    check_unknown_genes(ugenes)

if __name__ == '__main__':
    main()
