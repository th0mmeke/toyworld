"""
Created on 6/05/2013

@author: thom
"""

import logging
import os
import argparse

import config
from evaluators.evaluator import Evaluator
from experiment import Experiment
import xml.etree.cElementTree as ElementTree
from reactor_model.spatial_reaction_vessel import SpatialReactionVessel


class Locations(object):

    _dimension = 2

    @classmethod
    def _is_reaction(cls, reaction):
        reactant_smiles = [x['smiles'] for x in reaction['reactants']]
        product_smiles = [x['smiles'] for x in reaction['products']]
        return sorted(reactant_smiles) != sorted(product_smiles)

    def get_frame(self, frame, number_of_cells, id_to_smiles, id_to_index):
        """Output the molecule locations in a specific frame, one molecule per line. The format is: UID, SMILES, x, y
        (x,y start from 1.)"""

        cell_size = 2.0 / number_of_cells
        print("{}".format(len(frame)))
        for global_id, location in frame.iteritems():
            index = id_to_index[global_id]  # 0...number of unique molecules-1
            smiles = id_to_smiles[global_id]
            x, y = SpatialReactionVessel._get_grid_cell(location, cell_size)  # map to grid - each value 0,1...
            print("{},{},{},{}".format(index, smiles, x + 1, y + 1))

    def build_reaction_locations(self, results_filename, number_of_cells, start_t, delta_t):
        """Output the molecule locations at a series of frames (a frame group) taken from the original datastream.

        :param results_filename: Filename of the output experimental data
        :param start_t: Start time for the frame group
        :param delta_t: Size of frame group - end time = start_t + delta_t

        The output format is:
        <number of frames, start time, end time>
        <number of locations in this frame>
        <output from get_frame() for this frame>
        <number of locations in next frame>
        <output from get_frame() for this next frame>...
        """

        frame_group = {}
        id_to_smiles = {}

        finished = False

        for block in Evaluator.incr_load_results(results_filename):
            if 'xml_parameters' in block.keys():
                e = Experiment(None, ElementTree.fromstring(block['xml_parameters']))
                if e._dimension != 2:
                    raise ValueError("Must be 2-dimensional data")

            if 'smiles_map' in block.keys():
                id_to_smiles = {global_id: smiles for global_id, smiles in block['smiles_map'].iteritems()}  # short-cut

            for reaction in block['reactions']:
                for reaction_item in reaction['reactants'] + reaction['products']:
                    assert reaction_item['id'] not in id_to_smiles.keys() or id_to_smiles[reaction_item['id']] == reaction_item['smiles']
                    if reaction_item['id'] not in id_to_smiles.keys():
                        id_to_smiles[reaction_item['id']] = reaction_item['smiles']
                try:
                    t = reaction['t']
                except:
                    raise ValueError("Results data must include reaction time-stamps - commit 2b1f5a8 or later is required")

                if t >= start_t and t < start_t + delta_t:
                    frame_group[t] = reaction['locations']
                    assert len(reaction['locations'].keys()) == len(set(reaction['locations'].keys()))
                else:
                    finished = True

            if finished:
                # write out locations...
                if len(frame_group) > 0:
                    id_to_index = {global_id: index for global_id, index in zip(id_to_smiles.iterkeys(), range(len(id_to_smiles)))}
                    print("{},{},{},{},{}".format(len(frame_group), number_of_cells, len(id_to_index) - 1, start_t, start_t + delta_t))

                    logging.info("Saving {} slices from {} to {}.".format(len(frame_group), start_t, start_t + delta_t))

                    for frame in frame_group.itervalues():
                        self.get_frame(frame, number_of_cells, id_to_smiles, id_to_index)
                break

        logging.info("Built reaction locations")

    def get_all_times(self, results_filename):
        """Show all unique timestamps of items in datafile."""

        for block in Evaluator.incr_load_results(results_filename):
            t = -1  # a time known not to be in the results
            for reaction in block['reactions']:
                if reaction['t'] != t:
                    print("{}".format(reaction['t']))
                    t = reaction['t']

if __name__ == "__main__":

    kwargs = {'format': '%(asctime)s %(levelname)s:%(message)s', 'datefmt': '%m/%d/%Y %I:%M:%S %p', 'level': getattr(logging, 'FATAL', None)}
    logging.basicConfig(**kwargs)

    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--start_t', default=0, type=float, help='Start point of the frame group, in seconds')
    parser.add_argument('-d', '--delta_t', default=0.01, type=float, help='Size of the frame group in seconds')
    parser.add_argument('--end_time', action='store_true', help='Show timestamp of final item in data file')
    parser.add_argument('--all_times', action='store_true', help='Show all timestamps of items in data file')
    parser.add_argument('--num_cells', type=int, default=50, help='Number of cells on each dimension of each frame')
    parser.add_argument('-l', '--log_level', default='INFO', help="Set the logging level")
    parser.add_argument('filename', help="Filename for source data")
    args = parser.parse_args()

    results_filename = os.path.join(config.DataDir, args.filename)

    if args.end_time:
        results = Evaluator.load_results(results_filename)
        print("{}".format(results['t']))
    else:
        locations = Locations()
        if args.all_times:
            locations.get_all_times(results_filename)
        else:
            locations.build_reaction_locations(results_filename, args.num_cells, args.start_t, args.delta_t)
