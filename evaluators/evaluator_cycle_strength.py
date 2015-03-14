"""
Created on 6/05/2013

@author: thom
"""

import logging

from evaluators.evaluator import Evaluator
from evaluators.reaction_network import ReactionNetwork


class EvaluatorCycleStrength(Evaluator):

    """To compare cycles, we need average strength and variation for each observed cycle
    Hypothesis: that cycle strength (proxy for efficiency) increases over time if molecular evolution occurring
    """

    def get_result_titles(self):
        return ["Actual likelihood", "Expected likelihood"]

    def evaluate(self, results, min_length=3, min_count=2, **kwargs):
        """:rtype: list of actual cycle likelihood, list of expected cycle likelihood"""

        def calculate_probabilities(cycles, reaction_network, actual=True):
            cycle_probabilities = {}
            for cycle, count in cycles.iteritems():

                # logging.info("Analysing cycle {}".format(cycle))
                cumulative_probability = 1
                for mol_idx in range(len(cycle)):
                    mol = cycle[mol_idx]
                    fan_out = sum(attributes['weight'] for attributes in reaction_network[mol].itervalues())
                    # logging.info("{} = {:<.2f}% ({}/{})".format(mol, (count * 1.0 / fan_out) * 100.0, count, fan_out))
                    if actual:
                        cumulative_probability *= (count * 1.0 / fan_out)
                    else:
                        cumulative_probability *= 1.0 / len(reaction_network[mol])

                cycle_probabilities[cycle] = cumulative_probability

            return cycle_probabilities

        super(EvaluatorCycleStrength, self).evaluate(results, **kwargs)
        cycles = ReactionNetwork.discover_potential_reaction_cycles(results, min_length=min_length, min_count=min_count)  # find all cycles
        reaction_network = ReactionNetwork.build_smiles_reaction_network(results)

        actual = calculate_probabilities(cycles, reaction_network, actual=True)
        expected = calculate_probabilities(cycles, reaction_network, actual=False)

        strong_cycles = [(cycle, actual_probability, expected[cycle]) for cycle, actual_probability in actual.iteritems() if actual_probability > expected[cycle]]
        sorted_strong_cycles = sorted(strong_cycles, key=lambda t: (t[1] / t[2]), reverse=True)

        logging.info("Top 10 reaction cycles sorted by strength (actual probability/expected probability) (length >= {} and count >= {})".format(min_length, min_count))
        for cycle, actual_probability, expected_probability in sorted_strong_cycles[:10]:  # top 10
            logging.info("cycle count={:<10} len={:<10} strength={:<10.6f} actual={:<10.6f} expected={:<10.6f}".format(
                cycles[cycle], len(cycle), actual_probability / expected[cycle], actual_probability, expected[cycle]))

        return actual, expected
