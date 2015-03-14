"""
Created on 6/05/2013

@author: thom
"""


from evaluator import Evaluator
from molecular_population import MolecularPopulation


class EvaluatorPopulationSummary(Evaluator):

    def get_result_titles(self):
        return ["t", "Number of molecular species"]

    def evaluate(self, results_filename, **kwargs):
        """Calculates some interesting statistics for the experimental run.

        :rtype: Number of unique reactants, number of reactant types, maximum times a molecule was a reactant, length of longest molecule"""

        results = Evaluator.load_results(results_filename)
        population = MolecularPopulation(population=results['initial_population'], reactions=results['reactions'], size=100)

        result = []
        for t in population.get_times():
            p = population.get_slice_by_time([t])
            number_of_items = len(set(item for item in p.get_items() if p.get_quantity(item) > 0))
            result.append([t, number_of_items])

        # print(result)
        return result
