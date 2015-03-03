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

        for record in record_list:
            gname = self.classifier.get_groupname(record)

            record.annotations['group_name'] = gname
            groups[gname].append(record)

        return dict(groups)
