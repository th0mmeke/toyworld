"""
Created on 22/03/2013

@author: thom
"""

import random
import logging
import importlib
import os

from rdkit.Chem import AllChem as Chem

import pymunk as pm
import pygame
from pygame.locals import *

from reactor_model.reaction_vessel import ReactionVessel
from ULPS import Float_t
from kinetics_2D import Kinetics2D
import config
from parameters import Parameters


class SpatialReactionVessel(ReactionVessel):

    """
    n-dimensional structure of fixed size (ranging [-reaction_vessel_size,reaction_vessel_size] in each dimension).
    Molecules bounce off virtual walls.
    With KineticMolecules of size 1, ratio of vessel size:molecule size = 10^6:1
    """

    reaction_vessel_size = 500  # -500->500
    ratio_vessel_screen = 0.75  # 0.75
    molecule_collision_type = 1
    wall_collision_type = 2

    def __init__(self, chemistry, population=None, parameters=Parameters(), product_selection_strategy="energy",
                 results_filename=os.devnull, states_filename=os.devnull):

        self.initial_average_ke = int(parameters.get('Energy'))
        self._visualize = parameters.get('Visualize')
        logging.info("Visualize is {}".format('ON' if self._visualize else 'OFF'))
        self._show_force_vectors = parameters.get('ShowForceVectors')
        self._show_orientation = parameters.get('ShowOrientation')
        if self._visualize:
            pygame.init()
            screen_size = int(SpatialReactionVessel.reaction_vessel_size * SpatialReactionVessel.ratio_vessel_screen)
            self._screen = pygame.display.set_mode((screen_size * 2, screen_size * 2))

        self._bodies = {}  # dictionary body:mol
        self._nothing_happening = 0
        self._space = pm.Space()
        self._space.gravity = pm.Vec2d(0, 0)

        wall_thickness = 100  # nice and thick so that fast moving molecules don't tunnel straight through
        wall_end_point = SpatialReactionVessel.reaction_vessel_size + wall_thickness
        self._walls = [pm.Segment(self._space.static_body, (-wall_end_point, -wall_end_point), (-wall_end_point, wall_end_point), wall_thickness),
                       pm.Segment(self._space.static_body, (-wall_end_point, wall_end_point), (wall_end_point, wall_end_point), wall_thickness),
                       pm.Segment(self._space.static_body, (wall_end_point, wall_end_point), (wall_end_point, -wall_end_point), wall_thickness),
                       pm.Segment(self._space.static_body, (-wall_end_point, -wall_end_point), (wall_end_point, -wall_end_point), wall_thickness)
                       ]
        for wall in self._walls:  # can't set these in the constructor
            wall.collision_type = SpatialReactionVessel.wall_collision_type
            wall.elasticity = 0.9999
            wall.friction = 0
        self._space.add(self._walls)

        super(SpatialReactionVessel, self).__init__(chemistry, population, parameters=parameters,
                                                    product_selection_strategy=product_selection_strategy,
                                                    results_filename=results_filename, states_filename=states_filename)

        self._space.add_collision_handler(SpatialReactionVessel.molecule_collision_type,
                                          SpatialReactionVessel.molecule_collision_type,
                                          begin=SpatialReactionVessel._begin_molecule_collision_handler,
                                          pre_solve=None,
                                          post_solve=None,
                                          separate=SpatialReactionVessel._end_separation_handler,
                                          context3=self)

        wall_energy_parameter = parameters.get('VesselWallEnergy')
        if wall_energy_parameter is not None:
            if wall_energy_parameter:  # could be False or True or a string...
                # could be True or a string...
                if type(wall_energy_parameter) == bool:  # True
                    wall_energy = self.initial_average_ke
                else:
                    wall_energy = float(wall_energy_parameter)

                if wall_energy > 0:  # cannot be zero
                    logging.info("Vessel Wall Energy = {} (molecules bouncing off the container will be set to this KE)".format(wall_energy))
                    self._space.add_collision_handler(SpatialReactionVessel.molecule_collision_type,
                                                      SpatialReactionVessel.wall_collision_type,
                                                      begin=SpatialReactionVessel._begin_wall_collision_handler,
                                                      pre_solve=None, post_solve=None, separate=None, context4=self,
                                                      wall_energy=wall_energy)

        if population is not None:

            mols = [getattr(importlib.import_module(self._molecule_module), self._molecule_class)(smiles, kinetic_energy=self.initial_average_ke) for smiles in population.get_population()]
            self.add_molecules(mols)

        self._energy_input = self._energy_output = 0  # reset after adding molecules

        logging.info("Spatial Reaction Vessel with initial KE = {}".format(self.initial_average_ke))
        self._write_initial(parameters, population)

    def add_molecules(self, molecules, locations=None):

        logging.debug(
            "Adding {} molecules of type {}.{} to a vessel with {} pre-existing molecules".format(len(molecules),
                                                                                                  self._molecule_module,
                                                                                                  self._molecule_class,
                                                                                                  len(self._space.bodies)))

        super(SpatialReactionVessel, self).add_molecules(molecules)
        if locations is None:  # generate some locations
            locations = [[random.uniform(-SpatialReactionVessel.reaction_vessel_size, SpatialReactionVessel.reaction_vessel_size) for i in range(2)] for mol in molecules]
        # Add them into the main space using the newly separated locations
        for mol, (x, y) in zip(molecules, locations):
            self._bodies[mol.body] = mol
            mol.set_position(x, y)
            for shape in mol.body.shapes:
                shape.collision_type = SpatialReactionVessel.molecule_collision_type
                self._space.add(shape)
            mol.color = 0, 0, 0
            self._space.add(mol.body)
            mol.tag = 1

    def remove_molecules(self, molecules):
        super(SpatialReactionVessel, self).remove_molecules(molecules)

        for molecule in molecules:
            for shape in molecule.body.shapes:
                self._space.remove(shape)
                del shape
            self._space.remove(molecule.body)
            del self._bodies[molecule.body]
            # objgraph.show_backrefs(molecule)
            del molecule

    def get_molecules(self):
        return self._bodies.itervalues()

    def get_state(self):
        state = {'locations': {mol.global_id: body.position for body, mol in self._bodies.iteritems()}}
        state['molecule_states'] = [mol.get_state() for mol in self.get_molecules()]
        return state

    def step(self):

        if self._visualize:
            for event in pygame.event.get():  # must call event.get() to prevent window freezing...
                if event.type == QUIT:  # window close
                    quit()

        self._t += self._delta_t

        self._space.step(self._delta_t)  # automatically trigger collision handler if required
        self._apply_energy_model(self._energy_object, self._delta_t)

        if self._visualize:
            self._screen.fill(pygame.color.THECOLORS["white"])
            self._show_state()
            image_filename = "save_{}.jpg".format(self._t)
            logging.info("Writing save image {}".format(image_filename))
            pygame.image.save(self._screen, image_filename)
            pygame.display.flip()

        self._write_data(self._reactions)
        ke = self.get_total_ke()
        if ke / (1.0 * len(self._bodies)) < 10.0:  # we're not goin' to make it...
            self._nothing_happening += 1
        else:  # reset the count
            self._nothing_happening = 0

        if self.iteration >= self.end_iteration or self._nothing_happening > 10:
            self._write_final(self.iteration)
            if self._visualize:
                pygame.display.quit()
            del self._space
            if self._nothing_happening > 10:
                logging.info("We'vvvve runnnnn outtttt offfff ennnnnnergyyyyy...")
                raise ValueError

    def _show_state(self):
        for body, mol in self._bodies.iteritems():

            for shape in mol.body.shapes:
                x = body.position.x + shape.offset.x * 2.0
                y = body.position.y + shape.offset.y * 2.0

                # map [-reaction_vessel_size/2,reaction_vessel_size/2] to [0,self._screen.get_height()]"""
                screen_x = int((x + SpatialReactionVessel.reaction_vessel_size) * SpatialReactionVessel.ratio_vessel_screen)
                screen_y = int((y + SpatialReactionVessel.reaction_vessel_size) * SpatialReactionVessel.ratio_vessel_screen)

                pygame.draw.circle(self._screen, mol.color, (screen_x, screen_y), int(shape.radius), 0)

    def _show_ke(self):
        ke = "{:.2f}".format(self.get_total_ke() / self.get_number_molecules())

        font = pygame.font.Font(None, 36)
        text = font.render(ke, 1, (10, 10, 10))
        textpos = text.get_rect()
        textpos.centerx = self._screen.get_rect().centerx
        self._screen.blit(text, textpos)

    def discover_reaction(self, reactant_mols):
        """
        Return a list of all feasible reaction options for these molecules, where feasible means without exceeding the available energy for the reaction, and
        preserving overall energy (potential + kinetic + internal) and conserving momentum. If no reaction is feasible, then just 'bounce' the reactant molecules -
        create a collision where the reactants bounce off each other.

        initial_PE + initial_KE + initial_IE = final_PE + final_KE + final_IE, and given:

        initial_PE + rxn.get_energy_delta() = final_PE, then:

        initial_PE + initial_KE + initial_IE = initial_PE + rxn.get_energy_delta() + final_KE + final_IE, which implies that:

        final_KE + final_IE = initial_KE + initial_IE - rxn.get_energy_delta()

        Note that bonds have NEGATIVE energy, so adding bonds REDUCES PE. Natural state is towards reduced PE - atoms want to be together...
        Reducing PE means increased KE - creating a bond RELEASES energy; breaking a bond ABSORBS energy

        :rtype: List of Reaction - each option is a Reaction which, when fired, results in a list of Molecule products
        """

        initial_pe = sum([mol.get_potential_energy(self.chemistry) for mol in reactant_mols])
        initial_ke = sum([mol.get_kinetic_energy() for mol in reactant_mols])
        initial_ie = sum([mol.get_internal_energy() for mol in reactant_mols])

        # Energy available for a reaction = internal energy + energy of collision - energy of centre of mass
        available_energy_for_reaction = initial_ie + initial_ke - Kinetics2D.get_CM_energy(reactant_mols)
        logging.debug("Available energy for reaction = IE {} + KE {} - CM {} = {}".format(initial_ie, initial_ke,
                                                                                          Kinetics2D.get_CM_energy(
                                                                                              reactant_mols),
                                                                                          available_energy_for_reaction))

        reaction_options = self.reaction_model.get_reaction_options(reactant_mols)

        rxn = self._select_reaction(reaction_options, available_energy_for_reaction)
        if rxn is not None:

            logging.debug("RXN energy delta = {}".format(rxn.get_energy_delta()))  # finalPE = initial_PE + energy_delta

            # Adjust velocities and internal energies of product molecules to maintain overall system energy and momentum

            product_mols = rxn.fire()  # these product_mols will initially have zero internal and zero kinetic energy
            out_v, final_ie = Kinetics2D.inelastic_collision(reactant_mols, product_mols,
                                                             rxn.get_energy_delta())  # delta = increase in PE, therefore decrease in KE
            total_mass_product_mols = sum([mol.get_mass() for mol in product_mols])
            for mol, v in zip(product_mols, out_v):
                mol.set_velocity(*v)
                mol.set_internal_energy(
                    final_ie * mol.get_mass() / total_mass_product_mols)  # split total IE across product mols

            # ################
            # Post-conditions
            final_pe = sum([x.get_potential_energy(self.chemistry) for x in product_mols])
            final_ke = sum([mol.get_kinetic_energy() for mol in product_mols])
            final_ie = sum([mol.get_internal_energy() for mol in product_mols])
            assert Float_t.almost_equal(final_pe, initial_pe + rxn.get_energy_delta(), max_diff=config.EnergyTolerance)
            assert Float_t.almost_equal(sum([mol.get_mass() for mol in product_mols]),
                                        sum([mol.get_mass() for mol in reactant_mols]))  # conservation of mass
            assert Float_t.almost_equal(initial_pe + initial_ke + initial_ie, final_pe + final_ke + final_ie,
                                        max_diff=config.EnergyTolerance)  # conservation of energy

            # logging.info("Reactant KEs = {}".format([mol.get_kinetic_energy() for mol in reactant_mols]))
            logging.debug("In_PE + In_KE + In_IE = {} + {} + {} = {}".format(initial_pe, initial_ke, initial_ie,
                                                                             initial_pe + initial_ke + initial_ie))
            logging.debug("Out_PE + Out_KE + Out_IE = {} + {} + {} = {}".format(final_pe, final_ke, final_ie,
                                                                                final_pe + final_ke + final_ie))

        return rxn

    @classmethod
    def _end_separation_handler(cls, space, arbiter, context3,  *args, **kwargs):
        mols = []
        for shape in arbiter.shapes:
            if shape.body in context3._bodies.keys():
                mols.append(context3._bodies[shape.body])
        for mol in mols:
            mol.tag = 0
        return False

    @classmethod
    def _begin_molecule_collision_handler(cls, space, arbiter, context3, *args, **kwargs):

        if context3.iteration >= context3.end_iteration:
            return False

        try:
            reactant_mols = [context3._bodies[shape.body] for shape in arbiter.shapes]
        except:
            return False  # one or more reactants were in collision with another molecule previously in this timestep

        tags = set([mol.tag for mol in reactant_mols])
        if len(tags) > 1 or 0 not in tags:
            return False

        rxn = context3.discover_reaction(reactant_mols)

        if rxn is None:
            return True  # let the standard collision handler take over - bounce these molecules

        context3.iteration += 1
        product_mols = rxn.fire()

        for mol in product_mols:
            mol.tag = 1

        collision_location = list(sum([shape.body.position for shape in arbiter.shapes]) / len(arbiter.shapes))
        context3.remove_molecules(reactant_mols)  # remove molecules from reactor and from space
        for mol in product_mols:
            mol_location = collision_location + mol.get_velocity() * context3._delta_t * 0.01
            context3.add_molecules([mol], locations=[mol_location])  # add it in

        logging.info("t={}, iteration={} (average ke={:.2f}): Reaction between {} giving {}".format(context3._t, context3.iteration,
                                                                                                    context3.get_total_ke() / context3.get_number_molecules(),
                                                                                                    [str(mol) for mol in reactant_mols],
                                                                                                    [str(mol) for mol in product_mols]))
        reactants = [{'id': mol.global_id, 'smiles': Chem.MolToSmiles(mol), 'ke': mol.get_kinetic_energy()} for mol in reactant_mols]
        products = [{'id': mol.global_id, 'smiles': Chem.MolToSmiles(mol), 'ke': mol.get_kinetic_energy()} for mol in product_mols]
        reaction = {'iteration': context3.iteration, 't': context3._t, 'reaction_site': collision_location, 'reactants': reactants, 'products': products}

        context3._reactions.append(reaction)

        return False

    @classmethod
    def _begin_wall_collision_handler(cls, space, arbiter, context4, *args, **kwargs):
        for shape in arbiter.shapes:
            if shape.body in context4._bodies.keys():
                mol = context4._bodies[shape.body]
                initial_ke = mol.get_kinetic_energy()
                mol.set_kinetic_energy(kwargs['wall_energy'])
                delta = initial_ke - kwargs['wall_energy']

                if delta > 0:
                    context4._energy_output += delta
                else:
                    context4._energy_input += delta

        return True  # now let the standard collision handler take over to bounce this molecule off the wall
