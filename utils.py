class BaseFilter(object):
    def __init__(self):
        pass

    def apply(self):
        """ Return True is record should be kept, False if record should be discarted
        """
        raise Exception('Filter apply method not implemented')
