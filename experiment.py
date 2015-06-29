"""
Created on 13 Aug 2013

@author: thom
"""

import logging
import random
import collections
import os
import xml.etree.cElementTree as ElementTree

from molecular_population import MolecularPopulation
from parameters import Parameters
from evaluators.evaluator import Evaluator
from reactor_model.reaction_vessel_factory import ReactionVesselFactory


class Experiment(object):
    """Run one experiment with same method, population and reactions, possibly with multiple repeats"""

    def __init__(self, dirname, experiment_parameters, common=None, factor_definitions=None):
        """ Current parameters used include:

        |  EnergyModel - addition and removal of energy during the simulation
        |  Reactions - model to determine the products of two interacting molecules
        |  Molecule - type of molecule
        |  MolecularRadius - any molecules passing within this distance will interact
        |  ReactionVesselDimension - the size of the reaction vessel (really only the ratio with the MoleculeInteractionRadius is of interest)

        :param dirname: Base directory for the results file
        :param experiment_parameters: Parameters for the simulation run, for example the energy model and reaction model to be used.
        :type experiment_parameters: ElementTree Element
        :param common: Common factors block
        :param factor_definitions: Factor definitions
        """

        self._dirname = dirname
        if factor_definitions is not None:
            experiment_parameters, self._factors, self._factor_titles = self.factors_to_parameters(experiment_parameters, common, factor_definitions)

        self._parameters = Parameters(experiment_parameters)
        self.name = self._parameters.get_attrib('name')
        self.end_iteration = int(self._parameters.get('Iterations'))
        self.end_time = int(self._parameters.get('Time'))
        self.repeats = int(self._parameters.get_attrib('repeats'))


    def factors_to_parameters(self, parameters, common, factor_definitions):
        """Convert a factor-based design to fully-specified parameters"""
        factors = {}  # save for later...
        factor_titles = {}  # save for later...
        factor = parameters.find('Factor')
        while factor is not None:
            parameters.remove(factor)
            key = factor.get('key')
            value = factor.get('value')
            if key in factor_definitions.keys() and value in factor_definitions[key].keys():
                factors[key] = value
                factor_titles[key] = factor_definitions[key]['Title']
                new_factor = factor_definitions[key][value]
                if len(new_factor) > 0:
                    try:
                        parameters.append(new_factor)
                    except:
                        parameters.extend(new_factor)

            factor = parameters.find('Factor')

        if common is not None:
            parameters.extend(common)
        logging.debug("Converted parameters = \n{}".format(ElementTree.tostring(parameters)))

        return parameters, factors, factor_titles

    def get_factors(self):
        """Return the factors for this experiment in standard order (ascending sort on factor key)"""
        return collections.OrderedDict(sorted(self._factors.items(), key=lambda t: t[0])).values()

    def get_factor_titles(self):
        return collections.OrderedDict(sorted(self._factor_titles.items(), key=lambda t: t[0])).values()

    def run(self):
        """Run this experiment"""

        experiment_random_seed = self._parameters.get_attrib('seed')
        if experiment_random_seed is not None:
            experiment_random_seed = float(experiment_random_seed)

        population_filename = self._parameters.get_filename("PopulationFilename")

        attempt_recovery = self._parameters.get_attrib('recover')

        with open(population_filename, "rU") as f:
            xml_population = f.read()
        original_population = MolecularPopulation(xml=xml_population)

        logging.info("Running experiment, repeated {} times".format(self.repeats))
        for repeat in range(self.repeats):
            if attempt_recovery and Evaluator.check_complete(self.get_results_filename(repeat)):
                logging.info("Skipping repeat {} as it seems we've already done it...".format(repeat + 1))
            else:
                self._run_repeat(repeat, original_population, experiment_random_seed)

    def _run_repeat(self, repeat, original_population, experiment_random_seed):

        if experiment_random_seed is None:
            logging.info("Using system random number seed")
        else:
            logging.info("Using random seed of {}".format(experiment_random_seed + repeat * 1.0))
            random.seed(experiment_random_seed + repeat * 1.0)

        results_filename = self.get_results_filename(repeat)
        states_filename = self.get_states_filename(repeat)
        iteration_blocksize = next_write = int(self._parameters.get('IterationBlocksize'))
        state_record_rate = next_state_write = float(self._parameters.get('StateRecordRate'))

        reaction_vessel = ReactionVesselFactory.new(population=original_population, parameters=self._parameters)

        logging.info("Running repeat #{} of {} iterations and maximum time of {}, and writing results to {} and states to {}".format(repeat + 1, self.end_iteration, self.end_time, results_filename, states_filename))

        initial_ke = reaction_vessel.get_total_ke()
        initial_pe = reaction_vessel.get_total_pe()
        initial_ie = reaction_vessel.get_total_ie()

        f_data = open(results_filename, "wb")
        f_states = open(states_filename, "wb")

        reaction_vessel.write_initial(original_population, self._parameters, f_data, f_states)

        while reaction_vessel.iteration < self.end_iteration and reaction_vessel.t < self.end_time:
            reaction_vessel.step()
            if reaction_vessel.iteration >= next_write:
                reaction_vessel.write_reactions(f_data, reaction_vessel.iteration - iteration_blocksize)
                logging.info("Writing reaction block at iteration {}-{}".format(reaction_vessel.iteration - iteration_blocksize, reaction_vessel.iteration))
                next_write += iteration_blocksize
            if reaction_vessel.t >= next_state_write:
                reaction_vessel.write_state(f_states)
                next_state_write = reaction_vessel.t + state_record_rate
                logging.info("Writing state at time {} (next write at {})".format(reaction_vessel.t, next_state_write))

        reaction_vessel.write_reactions(f_data, reaction_vessel.iteration - iteration_blocksize)
        reaction_vessel.write_final(f_data)

        final_energy = reaction_vessel.get_total_ke() + \
                       reaction_vessel.get_total_pe() + \
                       reaction_vessel.get_total_ie() - \
                       reaction_vessel.get_total_energy_input() + \
                       reaction_vessel.get_total_energy_output()

        logging.info('Finishing iteration = {}'.format(reaction_vessel.iteration))
        logging.info('{:>10}{:>15}{:>15}'.format('', 'initial', 'final'))
        logging.info('{:>10}{:>15.1f}{:>15.1f}'.format('PE', initial_pe, reaction_vessel.get_total_pe()))
        logging.info('{:>10}{:>15.1f}{:>15.1f}'.format('KE', initial_ke, reaction_vessel.get_total_ke()))
        logging.info('{:>10}{:>15.1f}{:>15.1f}'.format('IE', initial_ie, reaction_vessel.get_total_ie()))
        logging.info('{:>10}{:>15}{:>15.1f}'.format('+E', '', -reaction_vessel.get_total_energy_input()))
        logging.info('{:>10}{:>15}{:>15.1f}'.format('-E', '', reaction_vessel.get_total_energy_output()))
        logging.info('{:>10}{:>15.1f}{:>15.1f}'.format('Total', initial_ke + initial_pe + initial_ie, final_energy))

        del reaction_vessel

        logging.info("Finished repeat #{}".format(repeat + 1))

    def get_results_filename(self, suffix=0):
        return "{}-{:0>2}.data".format(os.path.join(self._dirname, self.name), suffix + 1)

    def get_states_filename(self, suffix=0):
        return "{}-{:0>2}.states".format(os.path.join(self._dirname, self.name), suffix + 1)
