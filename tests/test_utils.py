from unittest import TestCase

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import utils


class TestBaseClasses(TestCase):
    def test_basefilter(self):
        with self.assertRaises(NotImplementedError):
            utils.BaseFilter().apply()

    def test_baseclassifier(self):
        self.assertEqual(len(utils.BaseClassifier().rules), 0)

class TestMethods(TestCase):
    def test_cma(self):
        self.assertEqual(utils.next_cma(4, 1, 2), 3)
        self.assertEqual(utils.next_cma(25, 4, 2.5), 7)

    def test_gene_name_extraction(self):
        record = SeqRecord(Seq('AGTC'), description='DDB0191165|DDB_G0267380 |DNA coding sequence|gene: argE on chromosome: 1 position 414980 to 416538')

        self.assertEqual(utils.extract_gene_name(record), 'argE')
