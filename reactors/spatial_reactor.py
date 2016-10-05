"""
Created on 22/03/2013

@author: thom
"""

import importlib
import logging
import random

import pymunk as pm
from rdkit.Chem import AllChem as Chem

import config
from kinetics_2D import Kinetics2D
from parameters import Parameters
from reactors.reactor import Reactor
from util.ulps import Ulps


class SpatialReactor(Reactor):

    """
    n-dimensional structure of fixed size (ranging [-reaction_vessel_size,reaction_vessel_size] in each dimension).
    Molecules bounce off virtual walls.
    With KineticMolecules of size 1, ratio of vessel size:molecule size = 10^6:1
    """

    reaction_vessel_size = 500  # -500->500
    REACTANT = 1
    PRODUCT = 2
    WALL = 3
    bodies = {}  # dictionary body:mol
    space = pm.Space()
    space.gravity = pm.Vec2d(0, 0)

    def __init__(self, chemistry, population=None, parameters=Parameters(), product_selection_strategy="energy"):

        self.initial_average_ke = int(parameters.get('Energy'))

        wall_thickness = 100  # nice and thick so that fast moving molecules don't tunnel straight through
        wall_end_point = SpatialReactor.reaction_vessel_size + wall_thickness
        self._walls = [pm.Segment(self.space.static_body, (-wall_end_point, -wall_end_point), (-wall_end_point, wall_end_point), wall_thickness),
                       pm.Segment(self.space.static_body, (-wall_end_point, wall_end_point), (wall_end_point, wall_end_point), wall_thickness),
                       pm.Segment(self.space.static_body, (wall_end_point, wall_end_point), (wall_end_point, -wall_end_point), wall_thickness),
                       pm.Segment(self.space.static_body, (-wall_end_point, -wall_end_point), (wall_end_point, -wall_end_point), wall_thickness)
                       ]
        for wall in self._walls:  # can't set these in the Segment constructor
            wall.collision_type = SpatialReactor.WALL
            wall.elasticity = 0.9999
            wall.friction = 0

        super(SpatialReactor, self).__init__(chemistry, population, parameters=parameters,
                                             product_selection_strategy=product_selection_strategy)

        self._end_iteration = int(parameters.get('Iterations'))
        h = self.space.add_collision_handler(SpatialReactor.REACTANT,
                                             SpatialReactor.REACTANT)

        h.begin = SpatialReactor._begin_handler
        h.pre_solve = None
        h.post_solve = None
        h.separate = SpatialReactor._end_handler

        if population is not None:
            mols = [getattr(importlib.import_module(self._molecule_module), self._molecule_class)(smiles, kinetic_energy=self.initial_average_ke) for smiles in population.get_population()]
            self.add_molecules(mols)

        self._energy_input = self._energy_output = 0  # reset after adding molecules

        logging.info("Spatial Reaction Vessel with initial KE = {}".format(self.initial_average_ke))

    def add_molecules(self, molecules, locations=None):

        logging.debug(
            "Adding {} molecules of type {}.{} to a vessel with {} pre-existing molecules".format(len(molecules),
                                                                                                  self._molecule_module,
                                                                                                  self._molecule_class,
                                                                                                  len(self.space.bodies)))

        super(SpatialReactor, self).add_molecules(molecules)
        if locations is None:  # generate some locations
            locations = [[random.uniform(-SpatialReactor.reaction_vessel_size, SpatialReactor.reaction_vessel_size) for i in range(2)] for mol in molecules]
        # Add them into the main space using the newly separated locations
        for mol, (x, y) in zip(molecules, locations):
            self.bodies[mol.body] = mol
            mol.set_position(x, y)
            for shape in mol.body.shapes:
                shape.collision_type = SpatialReactor.REACTANT  # Mark molecule as potential reactant
                self.space.add(shape)
            self.space.add(mol.body)

    @classmethod
    def remove_molecules(cls, molecules):
        for molecule in molecules:
            for shape in molecule.body.shapes:
                cls.space.remove(shape)
                del shape
            cls.space.remove(molecule.body)
            del cls.bodies[molecule.body]
            del molecule

    def get_molecules(self):
        return self.bodies.itervalues()

    def get_state(self):
        '''Returns a snapshot of the current state for the reaction vessel
        'locations': dictionary of molecule_id:position,
        'molecule_states': list of molecule states (from mol.get_state())
        '''

        state = {'locations': {mol.global_id: body.position for body, mol in self.bodies.iteritems()}}
        state['molecule_states'] = [mol.get_state() for mol in self.get_molecules()]
        return state

    def step(self):

        self.t += self._delta_t
        self.space.step(self._delta_t)  # automatically trigger collision handler if required
        self._apply_energy_model(self._energy_object, self._delta_t)

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
                                                                                          Kinetics2D.get_CM_energy(reactant_mols),
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
            assert Ulps.almost_equal(final_pe, initial_pe + rxn.get_energy_delta(), max_diff=config.EnergyTolerance)
            assert Ulps.almost_equal(sum([mol.get_mass() for mol in product_mols]),
                                     sum([mol.get_mass() for mol in reactant_mols]))  # conservation of mass
            assert Ulps.almost_equal(initial_pe + initial_ke + initial_ie, final_pe + final_ke + final_ie,
                                     max_diff=config.EnergyTolerance)  # conservation of energy

            # logging.info("Reactant KEs = {}".format([mol.get_kinetic_energy() for mol in reactant_mols]))
            logging.debug("In_PE + In_KE + In_IE = {} + {} + {} = {}".format(initial_pe, initial_ke, initial_ie,
                                                                             initial_pe + initial_ke + initial_ie))
            logging.debug("Out_PE + Out_KE + Out_IE = {} + {} + {} = {}".format(final_pe, final_ke, final_ie,
                                                                                final_pe + final_ke + final_ie))

        return rxn

    @classmethod
    def _end_handler(cls, arbiter, space, data):
        '''Called when two molecules separate. Mark them as potential reactants.
        :param arbiter:
        :param space:
        :param data:
        :return:
        '''

        for shape in arbiter.shapes:
            shape.collision_type = SpatialReactor.REACTANT

        return False

    @classmethod
    def _begin_handler(cls, arbiter, space, data):

        try:
            reactant_mols = [cls.bodies[shape.body] for shape in arbiter.shapes]
        except:
            return False  # one or more reactants were in collision with another molecule previously in this timestep

        rxn = data.discover_reaction(reactant_mols)

        if rxn is None:
            return True  # let the standard collision handler take over - bounce these molecules

        data['iteration'] += 1
        product_mols = rxn.fire()

        collision_location = list(sum([shape.body.position for shape in arbiter.shapes]) / len(arbiter.shapes))
        cls.remove_molecules(reactant_mols)  # remove molecules from reactor and from space
        for mol in product_mols:
            mol_location = collision_location + mol.get_velocity() * data._delta_t * 0.01
            data.add_molecules([mol], locations=[mol_location])  # add it in
            for shape in mol.body.shapes:
                shape.collision_type = SpatialReactor.PRODUCT  # as a product, won't react again until _end_handler triggered

        logging.info("iteration={}: Reaction between {} giving {}".format(data['iteration'],
                                                                          [str(mol) for mol in reactant_mols],
                                                                          [str(mol) for mol in product_mols]))
        reactants = [{'id': mol.global_id, 'smiles': Chem.MolToSmiles(mol), 'ke': mol.get_kinetic_energy()} for mol in reactant_mols]
        products = [{'id': mol.global_id, 'smiles': Chem.MolToSmiles(mol), 'ke': mol.get_kinetic_energy()} for mol in product_mols]
        reaction = {'iteration': data.iteration, 'reaction_site': collision_location, 'reactants': reactants, 'products': products}

        data['reactions'].append(reaction)

        return False