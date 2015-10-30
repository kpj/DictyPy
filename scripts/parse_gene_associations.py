"""
Convert gene_association.dictyBase into dict
README: ftp://ftp.geneontology.org/pub/go/gene-associations/readme/dictyBase.README
"""

import json
import collections


result = collections.defaultdict(list)
fname = 'data/gene_association.dictyBase'

with open(fname, 'r') as fd:
    content = fd.read()

for line in content.split('\n'):
    if len(line) == 0 or line.startswith('!'):
        continue

    parts = line.split('\t')
    DB_Object_ID = parts[1]
    DB_Object_Symbol = parts[2]
    Aspect = parts[8]
    DB_Object_Name = parts[9]
    Synonym = parts[10]

    if len(DB_Object_Name) > 0:
        #print(DB_Object_ID, DB_Object_Symbol, Aspect, DB_Object_Name, Synonym)
        if not DB_Object_Name in result[DB_Object_ID]:
            print(DB_Object_ID, DB_Object_Name)
            result[DB_Object_ID].append(DB_Object_Name)

print(len(result))
json.dump(result, open('results/annotated_genes.json', 'w'))
