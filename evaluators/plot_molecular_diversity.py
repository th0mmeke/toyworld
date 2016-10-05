"""
Created on 6/05/2013

@author: thom
"""

import matplotlib.colors as colors

from atoms.molecular_population import MolecularPopulation
from evaluator import Evaluator
from plot import Plot


class PlotMolecularDiversity(Plot):

    def draw_figure(self, f1, results_filename, **kwargs):

        results = Evaluator.load_results(results_filename)
        population = MolecularPopulation(population=results['initial_population'], reactions=results['reactions'])  # , size=500)

        diversity = []
        for t in population.get_times():
            slice = population.get_slice_by_time([t])
            quantities = [slice.get_quantity(item) for item in slice.get_items() if slice.get_quantity(item) > 0]
            iteration_diversity = (1.0 * len(quantities)) / (1.0 * sum(quantities))
            diversity.append(iteration_diversity)

        ax = f1.add_subplot(1, 1, 1)  # one row, one column, first plot
        ax.set_title('Diversity of Molecular Types')

        ax.set_xlabel('Time')
        ax.set_ylabel('Diversity')
        ax.set_xlim(left=0, right=population.get_times()[-1])

        ax.plot(diversity, color=colors.cnames['slategray'])
        ax.grid()
