import os.path

from Bio import SeqIO


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
        print('%i known genes' % len(known_genes))
        print('%i unknown genes' % len(unknown_genes))
        print('%i skipped genes' % len(skipped_genes))

    return known_genes, unknown_genes


def main():
    genes, ugenes = parse_fasta_file('dicty_primary_protein')

if __name__ == '__main__':
    main()
