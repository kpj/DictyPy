from utils import BaseFilter


class LengthFilter(BaseFilter):
    """ Skip all sequences which are shorter than specified threshold
    """

    def __init__(self):
        self.threshold = 100

    def apply(self, record):
        return len(record.seq) > self.threshold
