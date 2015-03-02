import collections, json


class AnnotationGrouper(object):
    def choose_annotation(self, record):
        """ Choose "best" annotation out of list of possible ones
        """
        keywords = json.load(open('keywords.json'))
        annos = ' | '.join(record.annotations['manual'])

        for kw in keywords:
            if kw.lower() in annos.lower():
                return kw

        return record.annotations['manual'][0]

    def group(self, record_list):
        """ Group genes according to their annotation
            'gene_dict' is a dict from gene id to a SeqRecord
            This function returns a dict with keys corresponding to the annotation and values being genes as biopythons SeqRecords
        """
        groups = collections.defaultdict(list)

        for record in record_list:
            gname = self.choose_annotation(record)

            record.annotations['group_name'] = gname
            groups[gname].append(record)

        return dict(groups)

    def transform(self, genes, gene_dict):
        """ Tranform 'gene_dict' to list of (annotated) SeqRecords
        """
        def get_record(gene_id):
            for gene in genes:
                if gene.id == gene_id:
                    return gene

            raise RuntimeError('%s not found' % gene_id)

        record_list = []
        for gene_id, annotation in gene_dict.items():
            rec = get_record(gene_id)
            rec.annotations['manual'] = annotation

            record_list.append(rec)
        return record_list
