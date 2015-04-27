"""
Created on 6/05/2013

@author: thom
"""

import logging

import numpy as np
import pygame

from evaluator import Evaluator
from reactor_model.spatial_reaction_vessel import SpatialReactionVessel
from charged_molecule import ChargedMolecule
from reactor_model.spatial_reaction_vessel import SpatialReactionVessel

from abc import ABCMeta


class EvaluatorVisualiser(Evaluator):

    __metaclass__ = ABCMeta

    ratio_vessel_screen = 0.75

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
                screen_x = int((x + SpatialReactionVessel.reaction_vessel_size) * EvaluatorVisualiser.ratio_vessel_screen)
                screen_y = int((y + SpatialReactionVessel.reaction_vessel_size) * EvaluatorVisualiser.ratio_vessel_screen)

                pygame.draw.circle(self._screen, pygame.color.THECOLORS["black"], (screen_x, screen_y), int(shape.radius), 0)

    def get_result_titles(self):
        return ["Visualiser"]

    def evaluate(self, results_filename, reaction_network=None, **kwargs):
        '''Visualise molecules in reaction vessel'''

        states_filename = kwargs['states_filename']
        final_t = Evaluator.get_final_summary(results_filename)['t']
        delta_t = final_t / 50
        logging.info("Final t = {}, delta_t = {}".format(final_t, delta_t))
        next_t = 0

        pygame.init()
        screen_size = int(SpatialReactionVessel.reaction_vessel_size * EvaluatorVisualiser.ratio_vessel_screen)
        self._screen = pygame.display.set_mode((screen_size * 2, screen_size * 2))


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
                image_filename = "save_{}.jpg".format(t)
                logging.info("Writing save image {}".format(image_filename))
                pygame.image.save(self._screen, image_filename)
                pygame.display.flip()

                next_t += delta_t

        pygame.display.quit()

        return {}
