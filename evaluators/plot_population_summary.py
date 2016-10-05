"""
Created on 6/05/2013

@author: thom
"""

import logging
import random

import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt

from atoms.molecular_population import MolecularPopulation
from evaluator import Evaluator
from plot import Plot


class PlotPopulationSummary(Plot):

    def draw_figure(self, f1, results_filename, **kwargs):

        results = Evaluator.load_results(results_filename)
        num_labels = 6
        faint_but_visible = 5

        colormap = plt.get_cmap('Greys')
        cNorm = colors.Normalize(vmin=0, vmax=num_labels + 10)
        scalar_map = cm.ScalarMappable(norm=cNorm, cmap=colormap)

        # Reduce number of iterations in population to something that makes sense to plot...sample to give a maximum of 1000 iterations to plot
        population = MolecularPopulation(population=results['initial_population'], reactions=results['reactions'], size=1000)

        max_of_items = tuple()
        for item in population.get_items():
            values = population.get_slice_by_items([item])._population
            max_idx = values.argmax(axis=0)[0]
            max_of_items = max_of_items + ((values[max_idx][0], max_idx, item),)

        max_x = population.get_times()[-1]
        colours = {}
        point_maximums = {}
        count = 0
        for max_value, max_idx, item in sorted(max_of_items, reverse=True):
            if count < num_labels:
                colours[item] = (item.__hash__() % num_labels) + 10  # in the darker portion of the range
                if max_idx < max_x * 0.1:
                    # special handling for supplied values and others that start high and decrease - labels otherwise all overlap...
                    t = random.choice(population.get_times())
                    point_maximums[item] = (t, population.get_slice_by_time([t]).get_quantity(item))
                else:
                    point_maximums[item] = (max_idx, max_value)
            else:
                colours[item] = faint_but_visible
            count += 1

        ke = results['initial_kinetic_energy'] / results['initial_population'].get_population_size() * 1.0
        logging.info("Initial population size={}".format(results['initial_population'].get_population_size()))
        logging.info("Initial average KE={}".format(ke))

        ax = f1.add_subplot(1, 1, 1)  # one row, one column, first plot
        ax.set_title('Item quantity by time (initial average KE={:.2f})'.format(ke))

        ax.set_xlabel('Time')
        ax.set_ylabel('Quantity')

        for item in population.get_items():
            colorVal = scalar_map.to_rgba(colours[item])
            r = population.get_slice_by_items([item])._population[0:len(population.get_times()), 0:len(population.get_items())]
            line, = ax.plot(population.get_times(), r, color=colorVal)
            if colours[item] != faint_but_visible:
                line.set_label(item)
                x, y = point_maximums[item]
                h_alignment = 'center'
                if x > max_x * 0.95:
                    x = max_x - 5
                    h_alignment = 'right'
                if x < max_x * 0.1:
                    x = 5
                    h_alignment = 'left'
                ax.annotate(item, xy=(x, y + 1), color='black', backgroundcolor='white', horizontalalignment=h_alignment, size='small')
            else:
                line.set_linestyle('--')

        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc='best', fontsize='small', ncol=2)
        ax.grid()
