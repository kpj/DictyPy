import collections


class GeneGrouper(object):
    def choose_annotation(self, record):
        """ Choose "best" annotation out of list of possible ones
        """
        return record.annotations[0]

    def group(self, record_list):
        """ Group genes according to their annotation
            'gene_dict' is a dict from gene id to a SeqRecord
            This function returns a dict with keys corresponding to the annotation and values being genes as biopythons SeqRecords
        """
        groups = collections.defaultdict(list)

        for record in record_list:
            groups[self.choose_annotation(record)].append(record)

        return dict(groups)

    def transform(self, genes, gene_dict):
        """ Tranform 'gene_dict' to list of (annotated) SeqRecords
        """
        def get_record(gene_id):
            for gene in genes:
                if gene.id == gene_id:
                    return gene
            return None

        record_list = []
        for gene_id, annotation in gene_dict.items():
            rec = get_record(gene_id)
            rec.annotations = annotation

            record_list.append(rec)
        return record_list
