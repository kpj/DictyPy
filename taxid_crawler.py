import urllib.parse, urllib.request

from bs4 import BeautifulSoup


_base_url = 'http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=%d'

taxids = []
for page in [35325]: #439488
    dat = urllib.request.urlopen(_base_url % page).read()
    soup = BeautifulSoup(dat)

    for virus in soup.find_all('a', {'title': 'species'}):
        href = virus['href']
        res = urllib.parse.urlparse(href)
        query = urllib.parse.parse_qs(res.query)

        taxid = query['id'][0]
        taxids.append(taxid)

with open('taxids.txt', 'w') as fd:
    for tid in taxids:
        fd.write('%s\n' % tid)
