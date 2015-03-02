import os, os.path
import urllib.request

from bs4 import BeautifulSoup


_base_url = 'http://dictybase.org/db/cgi-bin/search/results.pl?class=dicty::UI::Search::Gene&query=*_RTE&page=%d'
_pages = range(1, 6)
_save_dir = 'dicty_rte_data'

# load search results from dictybase
if not os.path.isdir(_save_dir):
    if not os.path.isdir(_save_dir): os.mkdir(_save_dir)

    for page in _pages:
        print('Loading page %d' % page)
        dat = urllib.request.urlopen(_base_url % page).read()
        with open(os.path.join(_save_dir, 'dicty_rte_%d.html' % page), 'w') as fd:
            fd.write(str(dat))

# parse search results
rte = []
for page in _pages:
    with open(os.path.join(_save_dir, 'dicty_rte_%d.html' % page), 'r') as fd:
        soup = BeautifulSoup(fd.read())
        for row in soup.find_all('tr', {'class': 'row2'}):
            name = row.find('a').getText()
            if name.startswith('DDB_G') and name.endswith('_RTE'):
                rte.append(name)

# save results
with open('dicty_rte_list.txt', 'w') as fd:
    for g in rte:
        fd.write(g + '\n')
