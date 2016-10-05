"""
Created on 6/05/2013

@author: thom
"""

import logging

import numpy as np
import pysal

from evaluator import Evaluator
from reactors.spatial_reactor import SpatialReactor
from abc import ABCMeta


class EvaluatorSpatialCorrelation(Evaluator):

    __metaclass__ = ABCMeta

    def get_result_titles(self):
        return ["Gamma"]

    def evaluate(self, results_filename, reaction_network=None, **kwargs):
        """Calculate Gamma index for spatial autocorrelation for each time point in the reaction file

        Algorithm:

        for each time that a reaction occurred
            calculate Gamma for the molecules
        """

        states_filename = kwargs['states_filename']
        number_of_cells = 2 * SpatialReactor.reaction_vessel_size
        final_t = Evaluator.get_final_summary(results_filename)['t']
        delta_t = final_t / 5.0
        logging.info("Final t = {}, delta_t = {}".format(final_t, delta_t))
        next_t = 0

        weights = pysal.lat2W(number_of_cells, number_of_cells, rook=False)

        results_gamma = []

        for block in Evaluator.incr_load_states(states_filename):

            # Each block consists of {'t': time, 'state': state }
            # Each state is {'locations': {molecule id:position}, 'molecule_states':[mol.get_state()]}

            t = block['t']
            state = block['state']

            if t >= next_t:
                locations = state['locations']
                cell_locations = self._get_spatial_matrix(locations, number_of_cells)

                logging.info("Starting Gamma evaluation for t={}".format(t))
                gamma = pysal.Gamma(cell_locations, weights)
                results = {"t": t, "g": gamma.g, "g_z": gamma.g_z, "p_sim_g": gamma.p_sim_g, "min_g": gamma.min_g, "mean_g": gamma.mean_g, "max_g": gamma.max_g, "n": len(locations)}
                logging.info("Gamma for {}={}".format(t, results))
                # logging.info("Starting Geary evaluation for t={}".format(last_reaction['t']))
                # geary = pysal.Geary(cell_locations, weights).z_norm
                # logging.info("t={}, gamma={} geary={}".format(last_reaction['t'],gamma,geary))
                results_gamma.append(results)
                # results_geary.append(geary)

                next_t += delta_t

        return results_gamma

    def _get_spatial_matrix(self, locations, number_of_cells):
        grid = np.zeros((number_of_cells * number_of_cells))
        for location in locations.values():
            if abs(location[0]) > SpatialReactor.reaction_vessel_size or abs(location[1]) > SpatialReactor.reaction_vessel_size:
                print(location)
            else:
                cell = self._to_grid_coordinates(*location)
                grid[cell[1] * number_of_cells + cell[0]] = 1

        return grid

    def _to_grid_coordinates(self, x, y):
        """map [-reaction_vessel_size,reaction_vessel_size] to [0,2*reaction_vessel_size]"""
        x = int(x + SpatialReactor.reaction_vessel_size)
        y = int(y + SpatialReactor.reaction_vessel_size)
        return x, y
