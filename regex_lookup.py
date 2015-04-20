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

    print('Matching \n > "%s"' % '"\n > "'.join(patterns))

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
        '(TGATG[AT]{N}{{{x_min},{x_max}}}(?:(?:TGA)|(?:CGA)|(?:TGT)|(?:TTA)|(?:TGC)){N}{{{y_min},{y_max}}}T{N}{{{z_min},{z_max}}}CTGA)',
        '(TCAG{N}{{{z_min},{z_max}}}A{N}{{{y_min},{y_max}}}(?:(?:TAA)|(?:ACA)|(?:TCG)|(?:TCA)|(?:GCA)){N}{{{x_min},{x_max}}}[AT]CATCA)'
    ]
    lookUp(
        [s.format(x_min=10, x_max=50, y_min=2, y_max=40, z_min=12, z_max=45, N='[AGTCN]') for s in regexprs],
        'dicty_primary_genomic'
    )
