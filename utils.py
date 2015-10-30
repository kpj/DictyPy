import utils


class BaseFilter(object):
    """ Base class of every filter.
        Filters are the first measure to reduce the amount of data processed by this program and are executed on every entry present in the given fasta file. They should thus not be too resource intensive.
    """
    skip = False
    post_annotation = False # apply this filter after annotation step

    def apply(self):
        """ Return True if record should be kept, False if record should be discarted
        """
        raise NotImplementedError('Filter apply method not implemented')

class BaseClassifier(object):
    """ Base class of every classifier.
        Classifiers are used to group genes into their respective group. This is their logical group they will considered to be in in later processing steps.
    """
    RESULTS_DIR = 'results'
    data_file = ''

    def __init__(self):
        """ Init some variables with their default values
        """
        self.skip_filter = []

    def get_groupname(self, record):
        """ Generate groupname for given record
        """
        raise NotImplementedError('Groupname detection was not implemented')

    @staticmethod
    def preprocess(genes):
        """ Preprocess genes in some way
        """
        return genes

def next_cma(new_value, list_len, old_cma):
    """ Calculate next cumulative moving average
        'list_len' is the length of the currently being averaged list before adding the new value
    """
    return (new_value + list_len * old_cma) / (list_len + 1)

def extract_gene_name(record):
    """ Extract gene name from FastA entry description
    """
    return record.description.split('|')[3].split()[1] # regex anyone?

def parse_filters(post_annotation=False):
    """ Return list of filters given in `filter.py`
    """
    import filters
    classes = [getattr(filters, x) for x in dir(filters) if isinstance(getattr(filters, x), type) and x != 'BaseFilter' and getattr(filters, x).__bases__[0] == utils.BaseFilter and getattr(filters, x).post_annotation == post_annotation]
    return classes

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

def annotate_seq_records(genes, annotations):
    """ Annotate SeqRecord entries
    """
    def get_record(gene_id):
        for gene in genes:
            if gene.id == gene_id:
                return gene
        return None

    record_list = []
    processed_names = []
    for gene_id, annotation in annotations.items():
        rec = get_record(gene_id)

        if not rec is None:
            rec.annotations['manual'] = annotation

            record_list.append(rec)
            processed_names.append(rec.id)

    print(len(genes) - len(processed_names), '/', len(genes), 'annotated')

    # don't neglect entries which have no annotation
    for gene in genes:
        if not gene.id in processed_names:
            gene.annotations['manual'] = []
            record_list.append(gene)

    return record_list
