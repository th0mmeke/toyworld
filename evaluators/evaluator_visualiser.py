"""
Created on 6/05/2013

@author: thom
"""

import logging

import pygame
import os

from evaluator import Evaluator
from charged_molecule import ChargedMolecule
from reactors.spatial_reactor import SpatialReactor

from abc import ABCMeta


class EvaluatorVisualiser(Evaluator):

    __metaclass__ = ABCMeta

    ratio_vessel_screen = 0.75

    def __init__(self, partition=False):
        super(EvaluatorVisualiser, self).__init__(partition)

    def _show_state(self, locations, molecule_states):

        for state in molecule_states:
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
                screen_x = int((x + SpatialReactor.reaction_vessel_size) * EvaluatorVisualiser.ratio_vessel_screen)
                screen_y = int((y + SpatialReactor.reaction_vessel_size) * EvaluatorVisualiser.ratio_vessel_screen)

                pygame.draw.circle(self._screen, pygame.color.THECOLORS["black"], (screen_x, screen_y), int(shape.radius), 0)

    def get_result_titles(self):
        return ["Visualiser"]

    def evaluate(self, results_filename, reaction_network=None, **kwargs):
        '''Visualise molecules in reaction vessel'''

        states_filename = kwargs['states_filename']
        dirname = os.path.dirname(states_filename)

        try:
            final_t = Evaluator.get_final_summary(results_filename)['t']
            delta_t = final_t / 50
            logging.info("Final t = {}, delta_t = {}".format(final_t, delta_t))
        except:
            delta_t = 1.0

        pygame.init()
        screen_size = int(SpatialReactor.reaction_vessel_size * EvaluatorVisualiser.ratio_vessel_screen)
        self._screen = pygame.display.set_mode((screen_size * 2, screen_size * 2))

        next_t = 0

        for block in Evaluator.incr_load_states(states_filename):

            # Each block consists of {'t': time, 'state': state }
            # Each state is {'locations': {molecule id:position}, 'molecule_states':[mol.get_state()]}

            t = block['t']

            state = block['state']

            if t >= next_t:
                locations = state['locations']
                molecule_states = state['molecule_states']

                self._screen.fill(pygame.color.THECOLORS["white"])
                self._show_state(locations, molecule_states)
                image_filename = os.path.join(dirname, "{}-{}.jpg".format(os.path.splitext(results_filename)[0], t))
                logging.info("Writing image {}".format(image_filename))
                pygame.image.save(self._screen, image_filename)
                pygame.display.flip()

                next_t += delta_t

        pygame.display.quit()

        return {}
