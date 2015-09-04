"""
Created on 6/05/2013

@author: thom
"""

import logging

import pygame
import os

from evaluator import Evaluator
from charged_molecule import ChargedMolecule
from reactor_model.spatial_reaction_vessel import SpatialReactionVessel

from abc import ABCMeta


class EvaluatorForceIntegral(Evaluator):

    __metaclass__ = ABCMeta

    def __init__(self, partition=False):
        super(EvaluatorVisualiser, self).__init__(partition)

    def _get_force(self, locations, molecule_state):


            smiles = state['smiles']
            id = state['id']
            position = locations[id]

            if 'orientation' in state.keys():
                orientation = state['orientation']
            else:
                orientation = None
            mol = ChargedMolecule(smiles, orientation=orientation)

            for shape in mol.body.shapes:
                x = position.x + shape.offset.x * 2.0
                y = position.y + shape.offset.y * 2.0

                # map [-reaction_vessel_size/2,reaction_vessel_size/2] to [0,self._screen.get_height()]"""
                screen_x = int((x + SpatialReactionVessel.reaction_vessel_size) * EvaluatorVisualiser.ratio_vessel_screen)
                screen_y = int((y + SpatialReactionVessel.reaction_vessel_size) * EvaluatorVisualiser.ratio_vessel_screen)

                pygame.draw.circle(self._screen, pygame.color.THECOLORS["black"], (screen_x, screen_y), int(shape.radius), 0)

    def get_result_titles(self):
        return ["Force Integral"]

    def evaluate(self, results_filename, reaction_network=None, **kwargs):
        '''
        Sum forces acting on each molecule. Calculate force vector at each interval, and keep running total for each molecule.
        '''

        states_filename = kwargs['states_filename']
        dirname = os.path.dirname(states_filename)

        try:
            final_t = Evaluator.get_final_summary(results_filename)['t']
            delta_t = final_t / 50
            logging.info("Final t = {}, delta_t = {}".format(final_t, delta_t))
        except:
            delta_t = 1.0

        next_t = 0
        forces = {}  # dictionary of molecule['id'] : cumulative force or force integral

        for block in Evaluator.incr_load_states(states_filename):

            # Each block consists of {'t': time, 'state': state }
            # Each state is {'locations': {molecule id:position}, 'molecule_states':[mol.get_state()]}

            t = block['t']

            state = block['state']

            if t >= next_t:
                locations = state['locations']
                molecule_states = state['molecule_states']

                self._screen.fill(pygame.color.THECOLORS["white"])
                for molecule_state in molecule_states:

            body.force_vector = []
            nearest_points = self._space.nearest_point_query(mol.get_position(), ChargeReactionVessel.distance_cutoff)  # world co-ordinates
            # e.g., 'distance': 36.06281003619938, 'shape': <pymunk.Circle object at 0x107d7a610>, 'point': Vec2d(-27.2849762796, 322.918629605)

            for shape in body.shapes:
                # Find all shapes closer than distance_cutoff
                shape_position = body.local_to_world(shape.offset)  # world co-ordinates

                for near_point in nearest_points:
                    other_shape = near_point['shape']
                        near_shape_position = near_point['point']
                        magnitude = other_shape.charge * shape.charge * self._dipole_force_constant

                        if magnitude > ChargeReactionVessel.magnitude_cutoff:
                            minimum_distance = self._bodies[near_point['shape'].body].get_size() + mol.get_size()
                            force_vector = ChargeReactionVessel._get_force_vector(shape_position, near_shape_position, magnitude, minimum_distance)

                            body.apply_impulse(force_vector, shape_position)  # apply force at an offset (world-coordinates)
                            body.force_vector.append((shape_position, force_vector))

                    forces[state['id']] += self._get_force(state)

                next_t += delta_t

        print(forces)

        return {forces}
