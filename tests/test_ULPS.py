"""
Created on 25/05/2013

@author: thom
"""

from ULPS import Float_t

import unittest


class TestULPS(unittest.TestCase):

    def testFloat_t(self):
        self.assertTrue(Float_t(-1).negative())
        self.assertFalse(Float_t(1).negative())

    def testULPS(self):
        self.assertTrue(Float_t.almost_equal(1, 1))
        self.assertFalse(Float_t.almost_equal(1, 2))
        self.assertTrue(Float_t.almost_equal(1E-11, 1E-11))
        self.assertTrue(Float_t.almost_equal(1E10, 1E10))
        self.assertFalse(Float_t.almost_equal(1000, 1001))
        self.assertFalse(Float_t.almost_equal(67329.234, 67329.242))
        self.assertFalse(Float_t.almost_equal(1234567, 1234568))
        self.assertFalse(Float_t.almost_equal(1234567.12, 1234567.13))
        self.assertFalse(Float_t.almost_equal(123456789.12, 123456789.13))
        self.assertTrue(Float_t.almost_equal(123456789.12, 123456789.12))
        self.assertTrue(Float_t.almost_equal(-123456789.12, -123456789.12))
        self.assertTrue(Float_t.almost_equal(-123959.99999999942, -123959.99999999952))

        self.assertFalse(Float_t.almost_equal(381.4348779113786, 381.43487791113927, max_diff=1E-10))  # previous max_diff => false
        self.assertTrue(Float_t.almost_equal(381.4348779113786, 381.43487791113927))  # current max_diff for real momentum calculation result

        self.assertTrue(Float_t.almost_equal(240000.0, 240000.000201, maxUlpsDiff=1E7))  # current for real KE calculations

    def Struct(self):
        a = Float_t(1E-11).i
        b = Float_t(2E-11).i
        print("diff({0},{1}) = {2}".format(a, b, abs(a - b)))

if __name__ == "__main__":
    unittest.main()
