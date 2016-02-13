from unittest import TestCase

import numpy as np

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import utils


class TestBaseClasses(TestCase):
    def test_basefilter(self):
        with self.assertRaises(NotImplementedError):
            utils.BaseFilter().apply()

class TestMethods(TestCase):
    def test_cma(self):
        self.assertEqual(utils.next_cma(4, 1, 2), 3)
        self.assertEqual(utils.next_cma(25, 4, 2.5), 7)

    def test_gene_name_extraction(self):
        record = SeqRecord(Seq('AGTC'), description='DDB0191165|DDB_G0267380 |DNA coding sequence|gene: argE on chromosome: 1 position 414980 to 416538')

        self.assertEqual(utils.extract_gene_name(record), 'argE')

    def test_2d_binning(self):
        xdat = [3,3,3,3, 4,4,4]
        ydat = [1,1,1,1, 2,2,2]
        coords = utils.do_2d_binning(
            xdat, ydat,
            0.5, 0.5,
            5, 5
        )

        self.assertEqual(len(coords), 100)
        self.assertIn({'x': 3.5, 'y': 1.5, 'z': 4.0}, coords)
        self.assertIn({'x': 4.5, 'y': 2.5, 'z': 3.0}, coords)
        for i, y in enumerate(np.arange(0.5, 5.5, 0.5)):
            self.assertEqual(coords[i], {'x': 0.5, 'y': y, 'z': 0.0})

    def test_find_all_positions(self):
        pos = utils.find_all_positions('AAAA', 'AA')
        self.assertEqual(pos, [0, 2])

        pos = utils.find_all_positions('AGAAAGAAAAAA', 'AAA', force_orf=True)
        self.assertEqual(pos, [6, 9])
