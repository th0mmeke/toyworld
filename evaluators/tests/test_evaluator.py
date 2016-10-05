"""
Created on 6/03/2014

@author: thom
"""
import unittest
import os

import xml.etree.cElementTree as ElementTree

from evaluators.evaluator import Evaluator
import config
from parameters import Parameters


class EvaluatorSubClass(Evaluator):

    def evaluate(self):
        pass

    def get_result_titles(self):
        pass


class TestEvaluator(unittest.TestCase):

    def testIsPartitioned(self):
        e = EvaluatorSubClass()
        self.assertFalse(e.is_partitioned())
        e = EvaluatorSubClass(True)
        self.assertTrue(e.is_partitioned())

    def test_checkComplete(self):
        self.assertTrue(Evaluator.check_complete(os.path.join(config.TestDir, 'experiment-test-01.data')))
        self.assertFalse(Evaluator.check_complete(os.path.join(config.TestDir, 'experiment-incomplete-01.data')))

    def test_CheckPartitions(self):
        results_filename = os.path.join(config.TestDir, 'experiment-test-01.data')
        blocksize = 50

        # Test default of no specified partitions
        expected = Evaluator.load_results(results_filename)['reactions']
        actual = 0
        for block in Evaluator.incr_load_results(results_filename):
            actual += len(block['reactions'])
        self.assertEqual(len(expected), actual)

        partitions = [{'start': 50, 'end': blocksize + 50}]
        for block in Evaluator.incr_load_results(results_filename, partitions):
            for reaction in block['reactions']:
                self.assertTrue(reaction['iteration'] >= partitions[0]['start'] and reaction['iteration'] <= partitions[0]['end'])
                last_reaction = reaction['iteration']
        self.assertEqual(partitions[-1]['end'], last_reaction)

        partitions = [{'start': 50, 'end': 52}, {'start': 55, 'end': blocksize + 60}]
        for block in Evaluator.incr_load_results(results_filename, partitions):
            for reaction in block['reactions']:
                result = reaction['iteration'] >= partitions[0]['start'] and reaction['iteration'] <= partitions[0]['end']
                result = result or reaction['iteration'] >= partitions[1]['start'] and reaction['iteration'] <= partitions[1]['end']
                self.assertTrue(result)
                last_reaction = reaction['iteration']
        self.assertEqual(partitions[-1]['end'], last_reaction)


if __name__ == "__main__":
    unittest.main()
