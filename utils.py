class BaseFilter(object):
    """ Base class of every filter
    """
    skip = False

    def apply(self):
        """ Return True if record should be kept, False if record should be discarted
        """
        raise NotImplementedError('Filter apply method not implemented')

class BaseClassifier(object):
    """ Base class of every classifier
    """
    def __init__(self):
        # self.rules is a list of dicts, whereeach entry has the two keys
        # 'condition' and 'datafield'.
        # If 'condition' yields True, the respective gene is stored in 'datafield'.
        # 'condition' is a lambda function which takes single argument, a record.
        # E.g. "lambda record: True" would work as a default case (only use as last element in list).
        self.rules = []

def next_cma(new_value, list_len, old_cma):
    """ Calculate next cumulative moving average
        'list_len' is the length of the currently being averaged list before adding the new value
    """
    return (new_value + list_len * old_cma) / (list_len + 1)
