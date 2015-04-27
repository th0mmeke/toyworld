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

    def __init__(self, chemistry, population=None, parameters=Parameters(), product_selection_strategy="energy", results_filename=os.devnull, states_filename=os.devnull):

        self._product_selection_strategy = product_selection_strategy.lower()
        self.chemistry = chemistry
        self._energy_input = 0
        self._energy_output = 0

        reactions_module, reactions_class = parameters.get('Reactions').rsplit('.', 1)
        self.reaction_model = getattr(importlib.import_module(reactions_module), reactions_class)(chemistry)

        self._state_record_rate = float(parameters.get('StateRecordRate'))
        self._molecule_module, self._molecule_class = parameters.get('Molecule').rsplit(".", 1)
        energy_module, energy_class = parameters.get('EnergyModel').rsplit('.', 1)
        self._energy_object = getattr(importlib.import_module(energy_module), energy_class)(parameters)

        reactions_module, reactions_class = parameters.get('Reactions').rsplit('.', 1)
        self.reaction_model = getattr(importlib.import_module(reactions_module), reactions_class)(chemistry)

        self.iteration = self._previous_write_iteration = 0
        self._iteration_blocksize = int(parameters.get('IterationBlocksize'))
        self._next_write_iteration = self._iteration_blocksize

        self._t = 0
        self._delta_t = float(parameters.get('DeltaT'))
        self._next_frame = self._state_record_rate  # Experiment has already written out the initial state, so we want to start off after that

        self._reactions = []
        self._f_data = open(results_filename, "wb")
        self._f_states = open(states_filename, "wb")

    def __del__(self):
        self._f_states.close()
        self._f_data.close()

    def set_end_iteration(self, end_iteration):
        self.end_iteration = end_iteration

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

    def _write_initial(self, parameters, population):
        cPickle.dump({'xml_parameters': parameters.to_xml(),
                      'initial_population': population,
                      'initial_kinetic_energy': self.get_total_ke(),
                      'initial_potential_energy': self.get_total_pe(),
                      'smiles_map': self.get_mol_to_SMILES_map()}, self._f_data, 2)

        logging.info("Recording initial state of molecules in reaction vessel")
        cPickle.dump({'t': 0, 'iteration': 0, 'state': self.get_state()}, self._f_states)

    def _write_data(self, reactions):
        '''Dump the current block of information - reaction list to the data file, and state snapshot to states file

        Each block in the data file has the following structure:
        'block':{'start_block': int, 'end_block': int, 'reactions': list of reactions}

        Each entry in the states file has this structure:
        't': real, 'iteration': int, 'state': snapshot of state from self.get_state()'''

        if self.iteration >= self._next_write_iteration or self.iteration >= self.end_iteration:
            self._next_write_iteration = self.iteration + self._iteration_blocksize
            # Dump data to file
            cPickle.dump({'block': {'start_block': self._previous_write_iteration, 'end_block': self.iteration, 'reactions': self._reactions}}, self._f_data, 2)
            cPickle.dump({'t': self._t, 'iteration': self.iteration, 'state': self.get_state()}, self._f_states)
            self._previous_write_iteration = self.iteration + 1
            self._f_data.flush()  # force python file write
            os.fsync(self._f_data.fileno())  # force OS file write
            self._f_states.flush()
            os.fsync(self._f_states.fileno())
            del self._reactions
            self._reactions = []

    def _write_final(self, end_iteration):
        """If reached end successfully, write out final summary"""
        logging.info("Writing final summary to data file, ending at iteration {}".format(end_iteration))
        end_state = {'final_kinetic_energy': self.get_total_ke(),
                     'final_potential_energy': self.get_total_pe(),
                     'final_internal_energy': self.get_total_ie(),
                     'number_molecules': self.get_number_molecules(),
                     'energy_input': self.get_total_energy_input(),
                     'energy_output': self.get_total_energy_output(),
                     'iterations_completed': end_iteration,
                     't': self._t}
        cPickle.dump(end_state, self._f_data, 2)
        self._f_data.flush()  # force python file write
        os.fsync(self._f_data.fileno())  # force OS file write

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
