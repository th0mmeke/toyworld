"""
Created on 6/03/2013

@author: thom
"""


import logging
import cPickle
import os
import copy

from abc import ABCMeta, abstractmethod


class Evaluator(object):

    __metaclass__ = ABCMeta

    def __init__(self, partition=False):
        self._partition = partition

    def is_partitioned(self):
        return self._partition

    @abstractmethod
    def evaluate(self, results_filename, **kwargs): pass

    @abstractmethod
    def get_result_titles(self): pass

    @classmethod
    def load_results(cls, filename):
        """Load the complete set of results from the file with the given filename.

        :rtype: Dictionary with an entry with key 'reactions' for all reaction entries; otherwise key:value pairs as found in results file."""

        results = {}
        for block in Evaluator.incr_load_results(filename):
            # print(block)
            for k, v in block.iteritems():
                if k in results.keys():
                    for reaction in block['reactions']:
                        results['reactions'].append(reaction)
                else:
                    results[k] = v

        return results

    @classmethod
    def check_complete(cls, filename):
        """Check if the data in a file is complete - meaning that the Experiment that generated the data finished without error."""

        logging.info("Checking if the repeat data at {} is complete (that is, finished without error)".format(filename))

        try:
            Evaluator.get_final_summary(filename)
        except:
            return False
        else:
            return True

    @classmethod
    def get_initial_parameters(cls, filename):
        for block in Evaluator.incr_load_results(filename):
            initial_parameters = block
            break
        return initial_parameters

    @classmethod
    def get_final_summary(cls, filename):
        """Get the summary information from an experiment run - iterations completed, final timestamp and so on.
        Everything that isn't reaction data...

        :rtype: Dictionary"""

        if not os.path.isfile(filename):
            raise ValueError

        summary_information = {}

        with open(filename, "rb") as f:
            while True:
                try:
                    r = cPickle.load(f)
                except:
                    if 'iterations_completed' in summary_information.keys():
                        return summary_information  # only return summary_information if complete
                    else:
                        raise ValueError
                else:
                    for key, value in r.iteritems():

                        # special treatment for reaction blocks as these are written out sequentially
                        if key != 'block':
                            summary_information[key] = value

    @classmethod
    def incr_load_results(cls, filename, selected_partitions=None):
        """Load the next block of results from a results file.

        Dictionary with an entry with key 'reactions' for all reaction entries; otherwise key:value pairs as found in results file."""

        logging.info("Loading results from {}, selected partitions = {}".format(filename, selected_partitions))

        if selected_partitions is None:
            partitions = [{'start': 0, 'end': -1}]
        else:
            partitions = copy.copy(selected_partitions)  # make a copy as we'll modify this list as we go

        with open(filename, "rb") as f:
            finished = False
            end_block = -1
            results = {}
            partition = partitions.pop(0)
            while not finished:
                try:
                    r = cPickle.load(f)
                except:
                    logging.info("Finished loading results")
                    finished = True
                else:
                    results = {'reactions': []}
                    for key, value in r.iteritems():
                        # special treatment for reaction blocks as these are written out sequentially
                        if key != 'block':
                            results[key] = value
                        else:
                            if value['start_block'] != end_block + 1:
                                raise ValueError("Inconsistency in reaction blocks - previous block ended at {}, this block starts at {}".format(end_block, value['start_block']))

                            for reaction in value['reactions']:
                                if reaction['iteration'] >= partition['start']:
                                    if partition['end'] < 0 or reaction['iteration'] <= partition['end']:
                                        results['reactions'].append(reaction)
                                    else:
                                        try:
                                            partition = partitions.pop(0)  # move onto the next partition
                                        except:
                                            finished = True
                            end_block = value['end_block']

                    yield results

    @classmethod
    def incr_load_states(cls, filename):
        """Load the next block of state information from a states file."""

        logging.info("State filename = {}".format(filename))
        with open(filename, "rb") as f:
            finished = False
            while not finished:
                try:
                    r = cPickle.load(f)
                except:
                    finished = True
                else:
                    yield r
