"""
Created on 22/03/2013

@author: thom
"""

import copy
import logging
import random

from rdkit.Chem import AllChem as Chem

from atoms.molecule import Molecule
from parameters import Parameters
from reactions.reaction import Reaction
from reactor import Reactor
from util.ulps import Ulps


class AspatialReactor(Reactor):

    """
    Zero-dimensional reaction vessel modelling a well-mixed container.
    """

    def __init__(self, chemistry, population, parameters=Parameters(), product_selection_strategy="energy"):

        self.initial_average_ke = int(parameters.get('Energy'))
        super(AspatialReactor, self).__init__(chemistry, population, parameters=parameters, product_selection_strategy=product_selection_strategy)

        logging.info("Aspatial Reaction Vessel with initial KE = {}".format(self.initial_average_ke))

        self._molecules = []
        if population is not None:
            self.add_molecules([Molecule(smiles, kinetic_energy=self.initial_average_ke) for smiles in population.get_population()])
            self._energy_input = 0

    def add_molecules(self, molecules):
        super(AspatialReactor, self).add_molecules(molecules)
        self._molecules.extend(molecules)

    def remove_molecules(self, molecules):
        for mol in molecules:
            self._molecules.remove(mol)

    def get_molecules(self):
        return list(self._molecules)

    def step(self):
        """Will save no less than end_iteration, if possible, but may return additional reactions - all of those that
        took place at the time of the final iteration."""

        logging.debug("ke={},pe={},ie={}, input={}, output={}".format(self.get_total_ke(), self.get_total_pe(), self.get_total_ie(), self.get_total_energy_input(), self.get_total_energy_output()))

        self.t += self._delta_t
        self.iteration += 1

        reactant_mols = random.sample(self._molecules, 2)
        rxn = self.discover_reaction(reactant_mols)
        if rxn is None:
            return True  # let the standard collision handler take over - bounce these molecules

        product_mols = rxn.fire()

        self.add_molecules(product_mols)
        self.remove_molecules(reactant_mols)  # remove molecules from reactor

        logging.info("{}: Reaction between {} giving {}".format(self.iteration, [str(mol) for mol in reactant_mols], [str(mol) for mol in product_mols]))
        reactants = [{'id': mol.global_id, 'smiles': Chem.MolToSmiles(mol), 'ke': mol.get_kinetic_energy()} for mol in reactant_mols]
        products = [{'id': mol.global_id, 'smiles': Chem.MolToSmiles(mol), 'ke': mol.get_kinetic_energy()} for mol in product_mols]
        reaction = {'iteration': self.iteration, 't': self.t, 'reactants': reactants, 'products': products}

        self._reactions.append(reaction)

    def get_state(self):
        return {'molecule_states': [mol.get_state() for mol in self.get_molecules()]}

    def adjust_ke(self, absolute, relative):
        for mol in self._molecules:
            energy_output = mol.get_kinetic_energy() * relative
            mol.set_kinetic_energy(mol.get_kinetic_energy() + absolute - energy_output)
            self._energy_input += absolute
            self._energy_output += energy_output

    def get_speeds(self):
        return []

    def discover_reaction(self, reactant_mols):
        """
        Return a list of all feasible reaction options for these molecules, where feasible means without exceeding the available energy for the reaction, and
        preserving overall energy (potential + kinetic) and conserving momentum.

        :rtype: List of Reaction - each option is a Reaction which, when fired, results in a List of Molecule products
        """

        initial_PE = sum([mol.get_potential_energy(self.chemistry) for mol in reactant_mols])
        initial_KE = sum([mol.get_kinetic_energy() for mol in reactant_mols])
        initial_IE = sum([mol.get_internal_energy() for mol in reactant_mols])

        available_energy_for_reaction = initial_IE + initial_KE  # Energy available for a reaction = internal energy + energy of collision - energy of centre of mass

        reaction_options = self.reaction_model.get_reaction_options(reactant_mols)
        rxn = self._select_reaction(reaction_options, available_energy_for_reaction)

        if rxn is not None:
            # distribute the KE across the products
            product_mol = reactant_mols[0].combine_molecules(rxn.fire())  # create one combined KineticMolecule from reaction product Molecules
            product_mol.set_kinetic_energy(initial_KE - rxn.get_energy_delta())  # assign correct post_reaction KE
            product_mols = product_mol.split_molecule()  # split to assign proportional KE to products
            rxn = Reaction(product_mols, energy_delta=rxn.get_energy_delta())
        else:  # in case all the reaction options for these reactants need more energy than we have...
            product_mols = [copy.deepcopy(mol) for mol in reactant_mols]
            rxn = Reaction(product_mols, energy_delta=0, is_reaction=False)

        delta_pe = initial_PE - sum([x.get_potential_energy(self.chemistry) for x in product_mols])
        assert Ulps.almost_equal(delta_pe, -rxn.get_energy_delta())

        return rxn
