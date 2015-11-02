from unittest import TestCase

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from sequence_analyzer import DNAAnalyzer


class TestCodonUsage(TestCase):
    def setUp(self):
        self.dnana = DNAAnalyzer()
        self.seq = 'AAAAAGAAA'
        self.genes = [SeqRecord(Seq('AAAAAAAAG')), SeqRecord(Seq('AAAAAGAAA')), SeqRecord(Seq('TTTTTCTTT')), SeqRecord(Seq('TTCTTTTTC'))]

    def test_codon_counter(self):
        count = self.dnana._count_codons(self.seq)

        self.assertEqual(count['AAA'], 2)
        self.assertEqual(count['AAG'], 1)
        self.assertEqual(count['AAT'], 0)
        self.assertEqual(count['AAC'], 0)

    def test_codon_usage(self):
        codu = self.dnana.get_codon_usage(self.seq)

        self.assertEqual(round(codu['AAA'], 3), round(0.6666, 3))
        self.assertEqual(round(codu['AAG'], 3), round(0.3333, 3))
        self.assertEqual(codu['AAT'], None)
        self.assertEqual(codu['AAC'], None)

    def test_average_codon_usage(self):
        avg_codu = self.dnana.get_avg_codon_usage(self.genes)

        self.assertEqual(round(avg_codu['AAA'], 3), round(0.6666, 3))
        self.assertEqual(round(avg_codu['AAG'], 3), round(0.3333, 3))
        self.assertEqual(avg_codu['AAT'], None)
        self.assertEqual(avg_codu['AAC'], None)
        self.assertEqual(round(avg_codu['TTT'], 3), round(0.5, 3))
        self.assertEqual(round(avg_codu['TTC'], 3), round(0.5, 3))

    def test_codon_frequencies(self):
        avg_cod_freqs = self.dnana.get_codon_freqs(self.genes)

        self.assertEqual(round(avg_cod_freqs['AAA'], 3), 0.333)
        self.assertEqual(round(avg_cod_freqs['AAG'], 3), 0.167)
        self.assertEqual(round(avg_cod_freqs['TTT'], 3), 0.25)
        self.assertEqual(round(avg_cod_freqs['TTC'], 3), 0.25)
