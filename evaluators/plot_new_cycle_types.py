"""
Created on 6/05/2013

@author: thom
"""

from plot import Plot
from evaluator import Evaluator
from reaction_network import ReactionNetwork

import matplotlib.colors as colors
from matplotlib.ticker import MaxNLocator


class PlotNewCycleTypes(Plot):

    def draw_figure(self, f1, results_filename, **kwargs):

        final_summary = Evaluator.get_final_summary(results_filename)
        partition_length = 250
        number_of_partitions = 50
        partition_generator = range(final_summary['iterations_completed'] / number_of_partitions,
                                    final_summary['iterations_completed'] + 1, final_summary['iterations_completed'] / number_of_partitions)
        experiment_partitions = [[{'start': i - partition_length, 'end': i}] for i in partition_generator if i > partition_length]

        new_cycle_types_count = []
        cycle_types_count = []
        cumulative_cycles = set()

        for partition in experiment_partitions:
            cycles = ReactionNetwork.discover_actual_reaction_cycles(results_filename, node_type="smiles",
                                                                     selected_partitions=partition)
            cycle_types = set(cycles.iterkeys())
            new_cycle_types_count.append(len(cycle_types.difference(cumulative_cycles)))
            cumulative_cycles = cumulative_cycles.union(cycle_types)
            cycle_types_count.append(len(cumulative_cycles))

        ax = f1.add_subplot(1, 1, 1)  # one row, one column, first plot
        ax.set_title('New Cycle Types by time')

        ax.set_xlim(left=0, right=partition_generator[-1])

        ax.set_xlabel('Time')
        ax.set_ylim(bottom=0, top=max(new_cycle_types_count) + 1)
        ya = ax.get_yaxis()

        ya.set_major_locator(MaxNLocator(integer=True))
        ax.set_ylabel('Number of new cycles')
        ax.bar([p[0]['start'] for p in experiment_partitions], new_cycle_types_count, edgecolor=colors.cnames['lightgray'], color=colors.cnames['lightgray'], align='center', linewidth=10)

        ax2 = ax.twinx()

        ax2.set_ylim(bottom=0, top=max(cycle_types_count) + 1)
        ax2.set_ylabel('Cumulative number of new cycles')
        ax2.set_xlim(left=0, right=partition_generator[-1])
        ax2.plot([p[0]['start'] for p in experiment_partitions], cycle_types_count, linestyle='--', color=colors.cnames['slategray'])
