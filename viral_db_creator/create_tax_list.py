import csv
import collections


specs = []
spec_to_tax = {}

# get species names of wanted viruses
viral_counter = collections.defaultdict(int)
with open('species_group.csv') as fd:
    creader = csv.reader(fd)

    for row in creader:
        species = row[5]
        group = row[8]

        if group in ['ssRNA(+)', 'ssRNA(-)', 'dsRNA']:
            specs.append((species, group))
            viral_counter[group] += 1

print('Group results:')
for k, v in viral_counter.items(): print(' >', k, v)

# create species -> taxonomy id mapping
with open('species_tax.dmp') as fd:
    for line in fd.read().split('\n'):
        parts = line.split('|')
        if len(parts) <= 1: continue

        species = parts[1].strip()
        tax_id = parts[0].strip()

        spec_to_tax[species] = tax_id

# create tax id list
no_mapping_count = collections.defaultdict(int)
with open('taxids.txt', 'w') as fd:
    for spc, grp in specs:
        try:
            fd.write(spec_to_tax[spc] + '\n')
        except KeyError:
            no_mapping_count[grp] += 1

print('No mapping found for:')
for k, v in no_mapping_count.items(): print(' >', k, v)
