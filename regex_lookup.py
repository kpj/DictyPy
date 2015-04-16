import csv, re

from progressbar import ProgressBar

from fasta_parser import FastaParser


def get_gene_name(desc):
    return re.search(r'gene: ([\w\-\./]+) ?', desc).group(1)

def get_position(desc):
    match = re.search(r'(chromosome:.*)', desc)
    return match.group(1) if match else 'no position available'

def lookUp(patterns, data_file):
    """ Try to match given pattern in all fasta entries specified in given file
    """
    def match(record):
        ret = []
        seq = str(record.seq)

        for pat in patterns:
            match = re.search(pat, seq)
            if not match: continue
            ret.append((record.description, ', '.join(match.groups()), match.span()))

        return ret

    # parse genome data
    farser = FastaParser(data_file)
    genes = farser.parse()

    print('Matching "%s"' % '", "'.join(patterns))

    pbar = ProgressBar(maxval=len(genes))
    pbar.start()

    # generate result
    res = []
    for i, record in enumerate(genes):
        foo = match(record)
        if foo: res.extend(foo)

        pbar.update(i)
    pbar.finish()

    # sort result in natural order
    def natural_keys(text):
        def atoi(text):
            return int(text) if text.isdigit() else text
        return [atoi(c) for c in re.split('(\d+)', text)]
    res = sorted(res, key=lambda e: natural_keys(get_position(e[0])))

    # save result
    with open('results/regex_lookup.csv', 'w') as fd:
        cwriter = csv.writer(fd)

        cwriter.writerow(['gene name', 'position in genome', 'position in match', 'matched sequence'])
        for desc, seq, span in res:
            name = get_gene_name(desc)
            pos = get_position(desc)

            cwriter.writerow([name, pos, span, seq])

if __name__ == '__main__':
    regexprs = [
        '({N}{{{z}}}AGTC{N}{{{y_min},{y_max}}}GTAGT{N}{{{x}}})',
        '({N}{{{x}}}TGATG{N}{{{y_min},{y_max}}}CTGA{N}{{{z}}})'
    ]
    lookUp(
        [s.format(x=5, y_min=30, y_max=120, z=5, N='[AGTCN]') for s in regexprs],
        'dicty_primary_genomic'
    )
