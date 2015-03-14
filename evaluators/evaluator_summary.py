"""
Created on 6/05/2013

@author: thom
"""


import logging
import string

from rdkit.Chem import AllChem as Chem

from evaluator import Evaluator
from molecular_population import MolecularPopulation


class EvaluatorSummary(Evaluator):

    def get_result_titles(self):
        return ["Number of unique reactants", "Number of reactant types", "Maximum times a molecule was a reactant", "Length of longest molecule"]

    def evaluate(self, results_filename, **kwargs):
        """Calculates some interesting statistics for the experimental run.

        :rtype: Number of unique reactants, number of reactant types, maximum times a molecule was a reactant, length of longest molecule"""

        results = Evaluator.load_results(results_filename)
        population = MolecularPopulation(population=results['initial_population'], reactions=results['reactions'], size=100)

        initial_population = population.get_slice_by_time([0])
        initial_average_ke = results['initial_kinetic_energy'] / (initial_population.get_population_size() * 1.0)

        final_average_ke = results['final_kinetic_energy'] / (population.get_population_size() * 1.0)

        changing = len(population.get_changing_items())
        changing_percent = (changing * 1.0 / len(population.get_items())) * 100.0
        supplied_items = set(item for item in initial_population.get_items() if initial_population.get_quantity(item) > 0)
        final_items = set(item for item in population.get_items() if population.get_quantity(item) > 0)

        # TODO: change this to use Molecule construction rather than MolFromSmiles
        supplied_atom_count = sum([Chem.AddHs(Chem.MolFromSmiles(i)).GetNumAtoms() * initial_population.get_quantity(i) for i in supplied_items])
        final_atom_count = sum([Chem.AddHs(Chem.MolFromSmiles(i)).GetNumAtoms() * population.get_quantity(i) for i in final_items])

        iteration = collisions = 0
        reactant_ids = {}
        reactant_smiles = {}
        keys = set()

        for reaction in results['reactions']:

            keys = keys.union(set(reaction.keys()))

            # if ReactionNetwork._is_reaction(reaction):

            for reactant in reaction['reactants']:
                try:
                    reactant_ids[reactant['id']].append(iteration)
                except:
                    reactant_ids[reactant['id']] = [iteration]

                try:
                    reactant_smiles[reactant['smiles']].append(iteration)
                except:
                    reactant_smiles[reactant['smiles']] = [iteration]

            iteration += 1

        assert iteration + collisions == len(results['reactions'])

        logging.info("We began with {} types of items; there were {} active types at the end out of {} (both supplied and discovered) overall".format(
            len(supplied_items), len(final_items), len(population.get_items())))
        logging.info("The initial population size was {}; at the end, {}".format(initial_population.get_population_size(), population.get_population_size()))
        logging.info("The initial average ke was {}; at the end, {}".format(initial_average_ke, final_average_ke))
        logging.info("Approximately {:.2f}% of the item types ({}) changed quantity during the simulation, while {:.2f}% didn't change at all".format(
            changing_percent, changing, 100 - changing_percent))
        logging.info("The simulation ran for t={} and {} iterations ({} reactions and {} simple collisions)".format(results['t'], len(results['reactions']), iteration, collisions))
        logging.info("Supplied atoms = {}, final atoms = {}".format(supplied_atom_count, final_atom_count))
        if supplied_atom_count != final_atom_count:
            logging.info("Food items detected")
        logging.info("There were {} unique reactants and {} reactant types".format(len(reactant_ids), len(reactant_smiles)))
        sorted_items = sorted(final_items, key=lambda t: len(t), reverse=True)
        logging.info("The longest molecule type was {}".format(sorted_items[0]))

        reactant_smiles_occurences = [len(occurences) for id, occurences in reactant_smiles.iteritems()]
        logging.info("{} molecule types were a reactant in more than one reaction".format(len([x for x in reactant_smiles_occurences if x > 1])))
        logging.info("The minimum and maximum number of times a molecule type served as a reactant was {} and {}".format(min(reactant_smiles_occurences), max(reactant_smiles_occurences)))

        reactant_id_occurences = [len(occurences) for id, occurences in reactant_ids.iteritems()]
        logging.info("{} molecules were a reactant in more than one reaction (should be zero)".format(len([x for x in reactant_id_occurences if x > 1])))

        final_population_sorted = sorted([(item, population.get_quantity(item)) for item in final_items], key=lambda t: t[1], reverse=True)
        final_population_string = ["{} x {}".format(item, quantity) for item, quantity in final_population_sorted]
        logging.info("Final population = {}".format(string.join(final_population_string, ", ")))

        return len(reactant_ids), len(reactant_smiles), max(reactant_smiles_occurences), len(sorted_items[0])
