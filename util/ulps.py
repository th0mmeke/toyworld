"""
Created on 25/05/2013

@author: thom
"""

import struct
import logging


class Ulps(object):

    def __init__(self, f):

        rep = struct.pack('@d', f)  # use Python native sizes instead of 32-bit C standard float
        self.i = struct.unpack('@Q', rep)[0]

    def negative(self):
        return (self.i >> 63) != 0

    @classmethod
    def almost_equal(cls, a, b, max_diff=1E-7, maxUlpsDiff=3000):
        """http://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
        ...also used in Google Test
        """
        a = float(a)
        b = float(b)

        # Check if the numbers are really close -- needed when comparing numbers near zero
        if abs(a - b) <= max_diff:
            return True

        uA = Ulps(a)
        uB = Ulps(b)
        # Different signs means they do not match
        if uA.negative() != uB.negative():
            logging.info("NEGATIVE: {}={},{}={}".format(a, uA.i, b, uB.i))
            return False

        # Find the difference in ULPs

        ulpsDiff = abs(uA.i - uB.i)
        if ulpsDiff < maxUlpsDiff:
            return True

        logging.info("ULPSDIFF: {}={},{}={}".format(a, uA.i, b, uB.i))
        return False
