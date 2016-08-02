#!/usr/bin/env python

import unittest
import os.path
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from cigar import Cigar

class TestCigar(unittest.TestCase):
    def setUp(self):
        self.test_cigar_str = "5S5=3I15=4D10=3X5S"
        self.test_cigar_arr = [('S', 5), ('=', 5), ('I', 3), ('=', 15),
                          ('D', 4), ('=', 10), ('X', 3), ('S', 5)]
        self.c = Cigar(self.test_cigar_str)

    def test_constructor(self):
        c = self.c
        self.assertEqual(c.cigar_str, self.test_cigar_str)
        self.assertEqual(c.cigar_arr, self.test_cigar_arr)

    def test_equality(self):
        another_c = Cigar(self.test_cigar_str)
        self.assertEqual(self.c, another_c)

    def test_base_count(self):
        c = self.c
        self.assertEqual(c.num_total_bases(),         50)
        self.assertEqual(c.num_match(),               30)
        self.assertEqual(c.num_mismatch(),             3)
        self.assertEqual(c.num_insertion(),            3)
        self.assertEqual(c.num_deletion(),             4)
        self.assertEqual(c.num_softclip(),        (5, 5))
        self.assertEqual(c.num_query_bases(),         46)
        self.assertEqual(c.num_query_bases_in_aln(),  36)
        self.assertEqual(c.num_target_bases_in_aln(), 37)

    def test_reverse(self):
        reversed_cigar_str = "5S3X10=4D15=3I5=5S"
        self.assertEqual(self.c.reverse().cigar_str, reversed_cigar_str)

    def test_strip_softclip(self):
        c = self.c
        stripped_c = Cigar("5=3I15=4D10=3X")
        self.assertEqual(c.strip_softclip(), stripped_c)

    def test_substring(self):
        c = self.c
        self.assertEqual(c.slice_left(10).cigar_str,  "5S5=")
        self.assertEqual(c.slice_right(10).cigar_str, "2=3X5S")

    def test_alignment_stats(self):
        c = self.c
        self.assertEqual(c.aln_bgn(), 6)
        self.assertEqual(c.aln_end(), 45)
        self.assertEqual(c.aln_len(), 40)
        self.assertEqual(c.aln_identity(), 0.75)

    def test_to_legacy_match(self):
        c = self.c.to_legacy_match()
        self.assertEqual(c.cigar_str, "5S5M3I15M4D13M5S")
        c_str = Cigar.convert_to_legacy_match(self.test_cigar_str)
        self.assertEqual(c_str, "5S5M3I15M4D13M5S")

    def test_arr_to_str(self):
        c_str = Cigar.arr_to_str(self.test_cigar_arr)
        self.assertEqual(c_str, self.test_cigar_str)

    def test_create_from_array(self):
        c = Cigar.create_from_array(self.test_cigar_arr)
        self.assertEqual(c.cigar_str, self.test_cigar_str)
        self.assertEqual(c.cigar_arr, self.test_cigar_arr)


if __name__ == '__main__':
    unittest.main()
