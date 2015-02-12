class BaseFilter(object):
    """ Base class of every filter
    """
    skip = False

    def __init__(self):
        pass

    def apply(self):
        """ Return True if record should be kept, False if record should be discarted
        """
        raise Exception('Filter apply method not implemented')

def next_cma(new_value, list_len, old_cma):
    """ Calculate next cumulative moving average
    """
    return (new_value + list_len * old_cma) / (list_len + 1)
