"""
Created on 6/05/2013

@author: thom
"""

import random
import copy

import numpy as np
import pysal
import matplotlib.pyplot as plt
from scipy.cluster import vq

from reactors.spatial_reactor import SpatialReactor


class TestSpatialCorrelation(object):

    @classmethod
    def save_grid(cls, locations, filename):

        number_of_cells = 100
        cell_locations = cls._get_spatial_matrix(locations, number_of_cells)
        frame2D = cell_locations.reshape((number_of_cells, number_of_cells))
        plt.imshow(frame2D, cmap=plt.cm.hot)
        filename = filename + ".eps"
        print("Drawing frame to path {}".format(filename))
        plt.savefig(filename, format='eps')

    @classmethod
    def evaluate_grid(cls, locations):
        number_of_cells = 100
        weights = pysal.lat2W(number_of_cells, number_of_cells, rook=False)
        cell_locations = cls._get_spatial_matrix(locations, number_of_cells)
        # print("Starting Gamma evaluation")
        gamma = pysal.Gamma(cell_locations, weights)
        results = {"g": gamma.g, "g_z": gamma.g_z, "p_sim_g": gamma.p_sim_g, "min_g": gamma.min_g, "mean_g": gamma.mean_g, "max_g": gamma.max_g}
        print("G gamma={}".format(results))

    @classmethod
    def evaluate_distance(cls, locations):
        np_locations = np.array(locations)
        weights = pysal.threshold_continuousW_from_array(np_locations, 0.75, alpha=-2.0)
#         for i in range(0,len(locations)):
#             print(weights.weights[i])
        attributes = np.ones((np_locations.shape[0], 1))
        print(attributes)
        gamma = pysal.Gamma(attributes, weights)
        results = {"g": gamma.g, "g_z": gamma.g_z, "p_sim_g": gamma.p_sim_g, "min_g": gamma.min_g, "mean_g": gamma.mean_g, "max_g": gamma.max_g}
        print("D gamma={}".format(results))
        geary = pysal.Geary(attributes, weights)
        print("Geary={} EC={}".format(geary.C, geary.EC))

    @classmethod
    def _get_spatial_matrix(cls, locations, number_of_cells):
        cell_size = 2.0 / number_of_cells
        grid = np.zeros((number_of_cells * number_of_cells))
        for location in locations:
            cell = SpatialReactor._get_grid_cell(location, cell_size)
            i = cell[1] * number_of_cells + cell[0]
            try:
                grid[i] = 1
            except:
                print(cell)
        return grid


if __name__ == '__main__':

    n = 800
    locations_clustered = []
    for i in range(0, n / 4):
        centre = [random.uniform(-1, 1), random.uniform(-1, 1)]
        locations_clustered.append(centre)
        for j in range(0, 3):  # three clustered points around center
            new = copy.deepcopy(centre)
            for dim in range(0, len(centre)):
                new[dim] += random.uniform(-0.01, 0.01)
                if new[dim] < -1.0:
                    new[dim] = -1.0
                elif new[dim] > 1.0:
                    new[dim] = 1.0
            locations_clustered.append(new)
    print("Clustered distribution of {} molecules".format(n))

    np_locations = np.array(locations_clustered)
    whitened = vq.whiten(np_locations)
    TestSpatialCorrelation.evaluate_grid(locations_clustered)

    n = 800
    print("Random distribution of {} molecules".format(n))
    locations_random = [[random.uniform(-1, 1), random.uniform(-1, 1)] for i in range(0, n)]
    TestSpatialCorrelation.evaluate_grid(locations_random)
    TestSpatialCorrelation.evaluate_distance(locations_random)

    n = 400
    print("Random distribution of {} molecules".format(n))
    locations_random = [[random.uniform(-1, 1), random.uniform(-1, 1), random.uniform(-1, 1)] for i in range(0, n)]
    TestSpatialCorrelation.evaluate_grid(locations_random)
    TestSpatialCorrelation.evaluate_distance(locations_random)
