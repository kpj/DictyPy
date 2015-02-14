from unittest import TestCase

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
