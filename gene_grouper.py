import collections, json


class GeneGrouper(object):
    def __init__(self, Classifier):
        self.classifier = Classifier()
        for f in self.classifier.skip_filter: f.skip = True

    def group(self, record_list):
        """ Group genes according to their annotation
            'gene_dict' is a dict from gene id to a SeqRecord
            This function returns a dict with keys corresponding to the annotation and values being genes as biopythons SeqRecords
        """
        groups = collections.defaultdict(list)
        stats = collections.defaultdict(int)

        for record in record_list:
            gnames = self.classifier.get_groupname(record)
            record.annotations['group_name'] = gnames

            for name in gnames:
                groups[name].append(record)
                stats[name] += 1

        for k, v in stats.items(): print(k, '->', v)

        return dict(groups)
