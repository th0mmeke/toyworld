"""
Created on 6/05/2013

@author: thom
"""

from evaluator import Evaluator
from reaction_network import ReactionNetwork

import networkx as nx

import logging
import collections
import itertools
import string
from abc import ABCMeta, abstractmethod


class EvaluatorCycles(Evaluator):

    __metaclass__ = ABCMeta

    def get_result_titles(self):
        return ["Number of cycles", "Length of longest cycle", "Count of most common cycle"]

    @abstractmethod
    def evaluate(self, results_filename, selected_partitions=None, min_length=3, min_count=2, cycles=None, **kwargs):

        if len(cycles) == 0:
            longest_cycle = 0
            common_cycle = 0
        else:
            if selected_partitions is not None:
                logging.info("Selected iterations: {}".format(string.join(["{}-{}".format(i['start'], i['end']) for i in selected_partitions], ",")))
            logging.info("Top 10 reaction cycles sorted by cycle count then cycle length (length >= {} and count >= {})".format(min_length, min_count))

            sorted_cycles = collections.OrderedDict(sorted(cycles.iteritems(), key=lambda t: len(t[0]) + t[1] * len(cycles), reverse=True))  # sorted by count
            common_cycles = collections.OrderedDict(itertools.islice(sorted_cycles.iteritems(), 10))  # top ten

            for cycle, count in common_cycles.iteritems():
                logging.info("{:<10} {}".format(count, cycle))

            logging.info("Top 10 reaction cycles sorted by cycle length then cycle count (length >= {} and count >= {})".format(min_length, min_count))
            max_cycles = max(cycles.itervalues())
            sorted_cycles = collections.OrderedDict(sorted(cycles.iteritems(), key=lambda t: len(t[0]) * max_cycles + t[1], reverse=True))  # sorted by length
            longest_cycles = collections.OrderedDict(itertools.islice(sorted_cycles.iteritems(), 10))  # top ten

            for cycle, count in longest_cycles.iteritems():
                logging.info("{:<10} {:<10} {}".format(count, len(cycle), cycle))

            longest_cycle = common_cycle = 0
            if len(longest_cycles) > 0:
                longest_cycle = len(longest_cycles.keys()[0])
            if len(common_cycles) > 0:
                common_cycle = common_cycles.values()[0]

        logging.info("{} unique cycles discovered, length of longest cycle = {}, count of most common cycle = {}".format(len(cycles), longest_cycle, common_cycle))
        return len(cycles), longest_cycle, common_cycle


class EvaluatorActualCycles(EvaluatorCycles):

    """Summarizes information about the actual cycles present in the experiment data. A cycle is one made up of molecules, rather than, say,
    molecule types. As an example, a cycle would not include an edge from 'H' to 'OH' to 'O' if the actual 'OH' molecules in the 'H'->'OH' and
    'OH'->'O' reactions aren't the same molecule.

    :rtype: number of cycles discovered, length of longest cycle, count of most common cycle"""

    def evaluate(self, results_filename, selected_partitions=None, min_length=3, min_count=2, **kwargs):
        logging.info("Beginning EvaluatorActualCycles")
        cycles = ReactionNetwork.discover_actual_reaction_cycles(
            results_filename, node_type="smiles", selected_partitions=selected_partitions, min_length=min_length, min_count=min_count)  # find all cycles
        return super(EvaluatorActualCycles, self).evaluate(results_filename, selected_partitions, min_length, min_count, cycles)


class EvaluatorPotentialCycles(EvaluatorCycles):

    """Summarizes information about the potential cycles present in the experiment data.

    :rtype: number of cycles discovered, length of longest cycle, count of most common cycle"""

    def evaluate(self, results_filename, selected_partitions=None, min_length=3, min_count=2, **kwargs):
        logging.info("Beginning EvaluatorPotentialCycles")
        cycles = ReactionNetwork.discover_potential_reaction_cycles(
            results_filename, node_type="smiles", selected_partitions=selected_partitions, min_length=min_length, min_count=min_count)  # find all cycles
        return super(EvaluatorPotentialCycles, self).evaluate(results_filename, selected_partitions, min_length, min_count, cycles)


class EvaluatorNXCycles(EvaluatorCycles):

    """Summarizes information about the potential cycles present in the experiment data.

    :rtype: number of cycles discovered, length of longest cycle, count of most common cycle"""

    def evaluate(self, results_filename, selected_partitions=None, min_length=3, min_count=2, **kwargs):
        reaction_network = ReactionNetwork._build_reaction_network(results_filename, selected_partitions=selected_partitions, node_type="smiles")
        logging.info("Running nx.simple_cycles()")
        cycles = nx.simple_cycles(reaction_network)
        logging.info("nx.simple_cycles() complete")
        return super(EvaluatorPotentialCycles, self).evaluate(results_filename, selected_partitions, min_length, min_count, cycles)


class EvaluatorActualCyclesDistribution(EvaluatorCycles):

    """Return a dictionary representing the distribution of cycle lengths in the reaction network.

    :rtype: dictionary of length cycle: count of cycles of this length"""

    def evaluate(self, results_filename, selected_partitions=None, min_length=3, min_count=2, **kwargs):

        cycles = ReactionNetwork.discover_actual_reaction_cycles(results_filename, selected_partitions=selected_partitions, min_length=min_length, min_count=min_count)  # find all cycles
        sorted_cycles = collections.OrderedDict(sorted(cycles.iteritems(), key=lambda t: len(t[0]) * len(cycles) + t[1], reverse=True))  # sorted by length
        d = {}
        for cycle, count in sorted_cycles.iteritems():
            try:
                d[len(cycle)] += count
            except:
                d[len(cycle)] = count

        return d
