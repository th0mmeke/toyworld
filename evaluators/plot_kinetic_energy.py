"""
Created on 6/05/2013

@author: thom
"""

from plot import Plot
from evaluator import Evaluator
import matplotlib.colors as colors
import logging


class PlotKineticEnergy(Plot):

    def draw_figure(self, f1, results_filename, **kwargs):

        states_filename = kwargs['states_filename']
        iterations = []
        average_ke = []

        for block in Evaluator.incr_load_states(states_filename):

            # Each block consists of {'t': time, 'state': state }
            # Each state is {'locations': {molecule id:position}, 'molecule_states':[mol.get_state()]}

            iterations.append(block['iteration'])
            state = block['state']

            total_ke = 0
            for molecule_state in state['molecule_states']:
                total_ke += molecule_state['ke']
            average_ke.append(total_ke / len(state['molecule_states']))

        final_summary = Evaluator.get_final_summary(results_filename)
        iterations_completed = final_summary['iterations_completed']

        if iterations[-1] < iterations_completed:  # final state entry may not correspond to the end of the experiment, so add in final state if required
            logging.info("State information for {} out of {} iterations".format(iterations[-1], iterations_completed))
            final_ke = final_summary['final_kinetic_energy'] / final_summary['number_molecules']
            average_ke.append(final_ke)
            iterations.append(iterations_completed)

        ax = f1.add_subplot(1, 1, 1)  # one row, one column, first plot
        ax.set_title('Kinetic Energy over Time')

        ax.set_xlabel('Iteration')
        ax.set_ylabel('Kinetic Energy')
        ax.set_ylim(bottom=0, top=max(average_ke) + 10)

        ax.plot(iterations, average_ke, color=colors.cnames['slategray'])
        ax.grid()
