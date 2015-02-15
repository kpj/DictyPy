import collections


class GeneGrouper(object):
    def group(self, gene_dict):
        """ Group genes according to their annotation
            'gene_dict' is a dict from gene id to a SeqRecord
            This function returns a dict with keys corresponding to the annotation and values being genes as biopythons SeqRecords
        """
        groups = collections.defaultdict(list)

        for gene_id, record in gene_dict.items():
            groups[record.annotations].append(record)

        return dict(groups)

    def transform(self, genes, gene_dict):
        """ Tranforms 'gene_dict' such that it actually contains SeqRecords
        """
        def get_record(gene_id):
            for gene in genes:
                if gene.id == gene_id:
                    return gene
            return None

        foo = {}
        for gene_id, annotation in gene_dict.items():
            rec = get_record(gene_id)
            rec.annotations = annotation

            foo[gene_id] = rec
        return foo
