import sys

from selenium import webdriver
from bs4 import BeautifulSoup
from Bio import SeqIO

import utils


class GeneAnnotator(object):
    def __init__(self):
        self.browser = None

    def __del__(self):
        if not self.browser is None:
            self.browser.close()

    def get_direct_annotation(self, record):
        """ Extract direct annotation information from amigo query result
        """
        soup = self._query_amigo(utils.extract_gene_name(record))

        table = soup.find('table', attrs={'class': 'bbop-js-search-pane-results-table'})
        table_body = table.find('tbody')
        rows = table_body.find_all('tr')

        # only check first row
        annotations = []
        for ele in rows[0].find_all('a'):
            if ele['href'].startswith('http://amigo.geneontology.org/amigo/term/GO:'):
                annotations.append(ele.text.lower())

        return annotations

    def _query_amigo(self, gene_name):
        """ Query amigo database over webinterface which requires javascript
        """
        if self.browser is None:
            self.browser = webdriver.Firefox()

        _base_url = 'http://amigo.geneontology.org/amigo/search/bioentity?q=%s'

        self.browser.get(_base_url % gene_name)
        content = self.browser.page_source
        soup = BeautifulSoup(content)

        return soup


def load_gene_annotations(genes, fname):
    """ Load gene annotations from the web
    """
    foo = {}
    errors = 0
    ganno = GeneAnnotator()
    for gene in genes:
        try:
            anno = ganno.get_direct_annotation(gene)
            foo[gene.id] = anno
        except:
            print('Error:', gene, '\n')
            errors += 1
    json.dump(foo, open(fname, 'w'))
    print(errors, 'errors')

def main():
    if len(sys.argv) != 3:
        print('Usage: %s <fasta file> <annotation output>' % sys.argv[0])
        sys.exit()

    genes = SeqIO.parse(sys.argv[1], format='fasta')

    load_gene_annotations(genes, sys.argv[2])

if __name__ == '__main__':
    main()
