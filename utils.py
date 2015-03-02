class BaseFilter(object):
    """ Base class of every filter
    """
    skip = False
    post_annotation = False # apply this filter after annotation step

    def apply(self):
        """ Return True if record should be kept, False if record should be discarted
        """
        raise NotImplementedError('Filter apply method not implemented')

class BaseClassifier(object):
    """ Base class of every classifier
    """
    requires_annotations = False

    def __init__(self):
        """ self.rules is a list of dicts, whereeach entry has the two keys
            'condition' and 'datafield'.
            If 'condition' yields True, the respective gene is stored in 'datafield'.
            'condition' is a lambda function which takes single argument, a record.
            E.g. "lambda record: True" would work as a default case (only use as last element in list).

            self.result is the field which should be processed eventually
        """
        self.result = []
        self.rules = []

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
    classes = [getattr(filters, x) for x in dir(filters) if isinstance(getattr(filters, x), type) and x != 'BaseFilter' and getattr(filters, x).post_annotation == post_annotation]
    return classes
