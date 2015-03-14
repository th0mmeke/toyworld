"""
Created on 6/05/2013

@author: thom
"""

from plot import Plot
from evaluator import Evaluator

import matplotlib.colors as colors
import logging


class PlotNewMoleculeTypes(Plot):

    def draw_figure(self, f1, results_filename, **kwargs):

        iterations = [0]
        molecular_types_difference = [0]
        cumulative_types = set()

        count = 0
        iteration = 0
        for block in Evaluator.incr_load_results(results_filename):
            for reaction in block['reactions']:
                for product in reaction['products']:
                    if product['smiles'] not in cumulative_types:
                        count += 1
                    cumulative_types.add(product['smiles'])

                iteration += 1

                if iteration % 50 == 0:
                    molecular_types_difference.append(count)
                    iterations.append(iteration)

        # Append final values
        if iteration % 50 != 0:
            molecular_types_difference.append(count)
            iterations.append(iteration)

        initial_parameters = Evaluator.get_initial_parameters(results_filename)

        ke = initial_parameters['initial_kinetic_energy'] / initial_parameters['initial_population'].get_population_size() * 1.0
        logging.info("Initial average KE={}".format(ke))

        ax = f1.add_subplot(1, 1, 1)  # one row, one column, first plot
        ax.set_title('New Molecular Types by time (initial average KE={:.2f})'.format(ke))

        ax.set_xlabel('Iteration')
        ax.set_ylabel('Quantity')

        ax.set_xlim(left=0, right=iterations[-1])

        ax.plot(iterations, molecular_types_difference, color=colors.cnames['slategray'])
        ax.grid()
