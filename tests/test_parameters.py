"""
Created on 22 Aug 2013

@author: thom
"""
import unittest
import config
import os
import random

import xml.etree.cElementTree as ElementTree

from parameters import Parameters


class Test(unittest.TestCase):

    def setUp(self):

        self.xml = """<xml>
    <One b='False'>1</One>
    <Two a='a'>2</Two>
    <Energy>
        <Start>100</Start>
        <End>200</End>
    </Energy>
    <Energy2>100</Energy2>
</xml>"""

    def test_get_param(self):
        defaults = {'One': 2, 'Three': 2, 'Missing': 2}
        p = Parameters(ElementTree.fromstring(self.xml), defaults)
        self.assertEqual(1, int(p.get('One')))
        self.assertEqual(2, int(p.get('Missing')))
        self.assertIsNone(p.get('Missing', default=False))  # fail if no default value provided
        self.assertIsNone(p.get('Fred'))  # fail if no default value provided

        p = Parameters(ElementTree.fromstring(self.xml))
        self.assertEqual(p.get("Energy").tag, "Energy")

    def test_get_attrib(self):
        p = Parameters(ElementTree.fromstring(self.xml).find('Two'))
        self.assertEqual('a', p.get_attrib('a'))

        p = Parameters(ElementTree.fromstring(self.xml).find('One'))
        self.assertEqual(False, p.get_attrib('b'))

        self.assertIsNone(p.get_attrib('c'))

    def test_get_filename(self):
        p = Parameters(ElementTree.fromstring('<xml><Fred>/test</Fred></xml>'))
        self.assertEqual(os.path.normpath('/test'), p.get_filename("Fred"))
        p = Parameters(ElementTree.fromstring('<xml><Fred>test</Fred></xml>'))
        self.assertEqual(os.path.join(config.DataDir, 'test'), p.get_filename("Fred"))

    def test_common(self):
        p = Parameters(ElementTree.fromstring('<xml><Common><Fred>1</Fred></Common></xml>'))
        self.assertIsNone(p.get('Fred'))  # We assume that any common and factor blocks have already been expanded (e.g., by Experiment.factors_to_parameters)

    def test_seed(self):
        for i in range(10):
            random.seed(1)
            a = random.randint(0, 10000)
            random.seed(1)
            b = random.randint(0, 10000)
            self.assertEqual(a, b)

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.test_init']
    unittest.main()
