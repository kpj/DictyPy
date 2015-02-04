import sys


def frint(s, **kwargs):
    """ Print with auto-flush
    """

    print(s, **kwargs)
    sys.stdout.flush()

class BaseFilter(object):
    """ Base class of every filter
    """

    def __init__(self):
        pass

    def apply(self):
        """ Return True is record should be kept, False if record should be discarted
        """
        raise Exception('Filter apply method not implemented')
