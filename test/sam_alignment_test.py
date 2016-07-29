#!/usr/bin/env python

import unittest
import os.path
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from sam_alignment import SamAlignment

class TestSamAlignment(unittest.TestCase):
    def setUp(self):
        with open(os.path.join(os.path.dirname(__file__), 'test.sam'), 'r') as f:
            line = f.next()
            self.sa = SamAlignment(line)

    def test_mandatory_fields(self):
        m = self.sa.mandatory_fields
        self.assertEqual(m['qname'], 'm160124_190915_42172_c100907682550000001823206404291616_s1_p0/1630/0_21782')
        self.assertEqual(m['flag'],  '16')
        self.assertEqual(m['rname'], '7')
        self.assertEqual(m['pos'],   '11044')
        self.assertEqual(m['mapq'],  '254')
        self.assertEqual(m['cigar'], '5S5=3I15=4D10=3X5S')
        self.assertEqual(m['rnext'], '*')
        self.assertEqual(m['pnext'], '0')
        self.assertEqual(m['tlen'],  '0')
        self.assertEqual(m['seq'],   'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC')
        self.assertEqual(m['qual'],  '*')

    def test_optional_fields(self):
        o = self.sa.optional_fields
        self.assertEqual(o['np'], '1')
        self.assertEqual(o['qs'], '0')
        self.assertEqual(o['qe'], '21782')
        self.assertEqual(o['rq'], '0.841891')

    def test_alignment_attributes(self):
        sa = self.sa
        self.assertEqual(sa.q_name, 'm160124_190915_42172_c100907682550000001823206404291616_s1_p0/1630/0_21782')
        self.assertEqual(sa.q_bgn,  6)
        self.assertEqual(sa.q_end,  41)
        self.assertEqual(sa.q_len,  46)
        self.assertEqual(sa.t_name, '7')
        self.assertEqual(sa.t_bgn,  11044)
        self.assertEqual(sa.t_end,  11080)
        self.assertEqual(sa.aln_len, 40)
        self.assertEqual(sa.aln_identity, 0.75)


if __name__ == '__main__':
    unittest.main()
