"""
Created on 6/05/2013

@author: thom
"""

from plot import Plot
from evaluator import Evaluator
from molecular_population import MolecularPopulation

import matplotlib.colors as colors

import logging


class PlotMoleculeTypes(Plot):

    def draw_figure(self, f1, results_filename, **kwargs):

        results = Evaluator.load_results(results_filename)

        # Reduce number of iterations in population to something that makes sense to plot...sample to give a maximum of 1000 iterations to plot
        population = MolecularPopulation(population=results['initial_population'], reactions=results['reactions'], size=1000)

        # Check that results from reactions match our provided final numbers
        # for i in population.get_items():
        #    assert population.get_quantity(i) == results['final_population'].get_quantity(i)

        molecular_types_count = []
        for t in population.get_times():
            slice = population.get_slice_by_time([t])
            molecular_types = set(item for item in slice.get_items() if slice.get_quantity(item) > 0)
            molecular_types_count.append(len(molecular_types))

        ke = results['initial_kinetic_energy'] / (results['initial_population'].get_population_size() * 1.0)
        logging.info("Initial average KE={}".format(ke))

        ax = f1.add_subplot(1, 1, 1)  # one row, one column, first plot
        ax.set_title('Count of Molecular Types by time (initial average KE={:.2f})'.format(ke))

        ax.set_xlabel('Time')
        ax.set_ylabel('Quantity')

        ax.plot(population.get_times(), molecular_types_count, color=colors.cnames['slategray'])

        ax.grid()
