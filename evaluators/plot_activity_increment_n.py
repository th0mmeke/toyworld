"""
Created on 6/05/2013

@author: thom
"""

import matplotlib.colors as colors

from atoms.molecular_population import MolecularPopulation
from evaluator import Evaluator
from plot import Plot


class PlotActivityIncrementN(Plot):

    def draw_figure(self, f1, results_filename, **kwargs):

        results = Evaluator.load_results(results_filename)

        # Reduce number of iterations in population to something that makes sense to plot...

        population = MolecularPopulation(population=results['initial_population'], reactions=results['reactions'], size=250)

        activity_increment = []
        iterations = []
        expected_proportion = {}

        for t in population.get_times():
            iterations.append(t)

            population_slice = population.get_slice_by_time([t])
            population_size = population_slice.get_population_size()
            actual_proportion = {}
            current_activity_increment = []

            for molecular_type in population_slice.get_items():
                actual_proportion[molecular_type] = population_slice.get_quantity(molecular_type) * 1.0 / population_size
                actual = 0
                if molecular_type in expected_proportion.keys():
                    if actual_proportion[molecular_type] > expected_proportion[molecular_type]:
                        actual = 1.0 * population_size * (actual_proportion[molecular_type] - expected_proportion[molecular_type]) ** 2
                current_activity_increment.append(actual)

            expected_proportion.update(actual_proportion)
            activity_increment.append(current_activity_increment)

        ax = f1.add_subplot(1, 1, 1)  # one row, one column, first plot
        ax.set_title('Evolutionary Activity by Time')

        ax.set_xlabel('Time')
        ax.set_ylabel('$a^N$')

        ax.plot(iterations, activity_increment, color=colors.cnames['slategray'])

        ax.grid()
