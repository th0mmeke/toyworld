"""
Created on 6/05/2013

@author: thom
"""


import logging

from evaluator import Evaluator
from reaction_network import ReactionNetwork


class EvaluatorIterations(Evaluator):

    """Calculate some summary information about the number of reactions and collisions completed."""

    def get_result_titles(self):
        return ["Iterations", "Reactions", "Collisions"]

    def evaluate(self, results_filename, **kwargs):
        """:rtype: number of iterations completed, number of reactions, and number of non-reaction collisions"""

        iteration = collisions = 0
        for block in Evaluator.incr_load_results(results_filename):
            for reaction in block['reactions']:
                if ReactionNetwork._is_reaction(reaction):
                    iteration += 1
                else:
                    collisions += 1

        logging.info("The simulation ran for {} iterations ({} reactions and {} simple collisions)".format(iteration + collisions, iteration, collisions))

        return iteration + collisions, iteration, collisions
