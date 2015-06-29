"""
Created on 22/03/2013

@author: thom
"""

import random
import logging
import importlib
import os
from abc import ABCMeta, abstractmethod
import cPickle

from rdkit.Chem import AllChem as Chem

from parameters import Parameters


class ReactionVessel(object):

    """
    Abstract class for Reaction Vessels.
    Methods for quantities and concentrations of molecules, and for tracking energy inputs and outputs.
    """

    __metaclass__ = ABCMeta

    def __init__(self, chemistry, population=None, parameters=Parameters(), product_selection_strategy="energy"):

        self._product_selection_strategy = product_selection_strategy.lower()
        self.chemistry = chemistry
        self._energy_input = 0
        self._energy_output = 0

        reactions_module, reactions_class = parameters.get('Reactions').rsplit('.', 1)
        self.reaction_model = getattr(importlib.import_module(reactions_module), reactions_class)(chemistry)

        self._molecule_module, self._molecule_class = parameters.get('Molecule').rsplit(".", 1)
        energy_module, energy_class = parameters.get('EnergyModel').rsplit('.', 1)
        self._energy_object = getattr(importlib.import_module(energy_module), energy_class)(parameters)

        reactions_module, reactions_class = parameters.get('Reactions').rsplit('.', 1)
        self.reaction_model = getattr(importlib.import_module(reactions_module), reactions_class)(chemistry)

        self.iteration = self.t = 0
        self._delta_t = float(parameters.get('DeltaT'))

        self._reactions = []

    def get_total_energy_input(self):
        return self._energy_input

    def get_total_energy_output(self):
        return self._energy_output

    def get_number_molecules(self):
        return len(list(self.get_molecules()))

    def get_total_ke(self):
        return sum([mol.get_kinetic_energy() for mol in self.get_molecules()]) * 1.0

    def get_total_pe(self):
        return sum([mol.get_potential_energy(self.chemistry) for mol in self.get_molecules()]) * 1.0

    def get_total_ie(self):
        return sum([mol.get_internal_energy() for mol in self.get_molecules()]) * 1.0

    def get_speeds(self):
        return [mol.get_speed() for mol in self.get_molecules()]

    def get_mol_to_SMILES_map(self):
        return {mol.global_id: Chem.MolToSmiles(mol) for mol in self.get_molecules()}

    @abstractmethod
    def get_molecules(self):
        pass

    @abstractmethod
    def step(self):
        pass

    @abstractmethod
    def add_molecules(self, molecules):
        pass

    @abstractmethod
    def remove_molecules(self, molecules):
        pass

    def write_initial(self, population, parameters, f_data, f_states):
        '''Record the initial starting conditions.
        :param population: the initial population. Can't assume that the reaction vessel will have this as a class variable
        :param parameters: likewise the experiment parameters for this run with this reaction vessel
        :param f_data: the file handle of the output file
        :return:
        '''
        cPickle.dump({'xml_parameters': parameters.to_xml(),
                      'initial_population': population,
                      'initial_kinetic_energy': self.get_total_ke(),
                      'initial_potential_energy': self.get_total_pe(),
                      'smiles_map': self.get_mol_to_SMILES_map()}, f_data, 2)

        logging.info("Recording initial state of molecules in reaction vessel")
        cPickle.dump({'t': 0, 'iteration': 0, 'state': self.get_state()}, f_states)

    def write_reactions(self, f_data, last_write):
        '''Dump the reaction list to the data file

        Each block in the data file has the following structure:
        'block':{'start_block': int, 'end_block': int, 'reactions': list of reactions}
        '''

        cPickle.dump({'block': {'start_block': last_write, 'end_block': self.iteration, 'reactions': self._reactions}}, f_data, 2)
        f_data.flush()  # force python file write
        os.fsync(f_data.fileno())  # force OS file write
        del self._reactions
        self._reactions = []

    def write_state(self, f_states):
        '''Dump a state snapshot to states file

        Each entry in the states file has this structure:
        't': real, 'iteration': int, 'state': snapshot of state from self.get_state()'''

        cPickle.dump({'t': self.t, 'iteration': self.iteration, 'state': self.get_state()}, f_states)
        f_states.flush()
        os.fsync(f_states.fileno())

    def write_final(self, f_data):
        """If reached end successfully, write out final summary"""
        logging.info("Writing final summary to data file, ending at iteration {} and time {}".format(self.iteration, self.t))
        end_state = {'final_kinetic_energy': self.get_total_ke(),
                     'final_potential_energy': self.get_total_pe(),
                     'final_internal_energy': self.get_total_ie(),
                     'number_molecules': self.get_number_molecules(),
                     'energy_input': self.get_total_energy_input(),
                     'energy_output': self.get_total_energy_output(),
                     'iterations_completed': self.iteration,
                     't': self.t}
        cPickle.dump(end_state, f_data, 2)
        f_data.flush()  # force python file write
        os.fsync(f_data.fileno())  # force OS file write

    def _apply_energy_model(self, energy_model, t):
        """Adjust the KE of each molecule in the reaction vessel according to the energy model and radiation rate.
        The energy model's get_absolute(t) method provides an energy input; the get_relative(t) proportional to a molecule's current KE.
        """
        for mol in self.get_molecules():
            ie = mol.get_internal_energy()
            ke = mol.get_kinetic_energy()

            absolute_input = energy_model.get_absolute_input(ie, t)
            radiation_output = energy_model.get_relative_output(ke, t)

            mol.set_internal_energy(absolute_input + ie)
            mol.set_kinetic_energy(ke - radiation_output)

            self._energy_input += absolute_input
            self._energy_output += radiation_output

    @abstractmethod
    def discover_reaction(self, reactants, chemistry, reaction_model):
        """Return a Reaction chosen from the set of all reactions possible between the provided reactants.

        :param reactants: reactants for the discovered reaction
        :type reactants: Molecule (containing one or more molecules)
        :param reaction_model: model to discover possible interactions between reactants
        :type reaction_model: Reactions
        :rtype: Reaction or None if no reaction possible
        """
        pass

    def _select_reaction(self, reaction_options, reaction_energy):
        """Select a Reaction from a list of alternative Reactions, weighted towards least energy required"""

        def sample(population, k=1):
            """Population is {value:probability}"""

            total = 0
            for value, p in population.iteritems():
                total += p

            if total != 0:
                bound = random.random() * total

                cumulative_sum = 0
                for value, p in population.iteritems():
                    cumulative_sum += p
                    if cumulative_sum > bound:
                        return value

            return None

        logging.debug(
            "Selecting reaction from a population of {} items, with reaction energy of {}".format(len(reaction_options),
                                                                                                  reaction_energy))
        if len(reaction_options) == 1 and reaction_options[0].get_energy_delta() == 0:  # a non-reactive collision
            return None

        weighted_options = {}
        for option in reaction_options:
            if option.get_energy_delta() < 0:
                weighted_options[option] = -option.get_energy_delta()
            else:
                e = reaction_energy - option.get_energy_delta()
                if e > 0:
                    weighted_options[option] = e
        logging.debug("The weighted options are: {}".format(weighted_options))

        # TODO: do the uniform selection before the weighting for speed..
        if self._product_selection_strategy == "energy":
            return sample(weighted_options)
        else:
            try:
                return random.choice(weighted_options.keys())
            except:
                logging.info("No reaction possible - reaction energy = {}".format(reaction_energy))
                return None  # if no valid options
