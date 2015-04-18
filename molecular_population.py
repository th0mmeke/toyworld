"""
Created on 30/01/2013

@author: thom
"""

import random
import logging

import numpy as np
from rdkit.Chem import AllChem as Chem
import xml.etree.cElementTree as ElementTree

from population import Population
from molecule import Molecule  # only used to get canonical SMILES representation for a molecule


class MolecularPopulation(Population):

    """Population subclass where items must be SMILES strings; method parameters expect either SMILES or Rdkit Mol values"""

    def __init__(self, xml=None, population=None, reactions=None, size=0):
        """As for Population, except we expect all population items to be in SMILES. We automatically standardise these -
        all hydrogens are made explicit, and the item is converted to canonical SMILES.

        :param size: attempt to resize to include only this many time-stamps"""

        if xml is not None:
            super(MolecularPopulation, self).__init__()
            self._load_population_from_xml(xml)
        else:
            super(MolecularPopulation, self).__init__(population=population)

        if reactions is not None:
            self._sample_reactions(reactions, size)

    def _load_population_from_xml(self, xml):

        logging.info("Loading population from XML")

        self.set_t(0)
        for xml in ElementTree.fromstring(xml).findall('setvalue'):
            item = xml.attrib['item'].encode("utf-8")
            # make items canonical via Molecule! Important! Otherwise N(=O)[O], say, in XML will be a different element from [O]N=O
            self.set_quantity(item=Chem.MolToSmiles(Molecule(item)), quantity=int(xml.attrib['quantity']))

        if len(self._index) == 0 or len(self._t) == 0:
            raise ValueError('No population found in XML')
        else:
            logging.info("Population loaded - {} unique items, {} total population size".format(len(self._index), self.get_population_size()))

    def _apply_reactions(self, reactions):

        count = 1

        for reaction in reactions:

            self.set_t(reaction['t'])  # use count rather than 'iteration' as collisions don't increment 'iteration' but do appear in reactions list
            for reactant in reaction['reactants']:
                self.set_quantity(reactant['smiles'], self.get_quantity(reactant['smiles']) - 1)
            for product in reaction['products']:
                self.set_quantity(product['smiles'], self.get_quantity(product['smiles']) + 1)

            if count % 1000 == 0:
                logging.debug("Added reactions up to {}".format(count))

            count += 1

    def _sample_reactions(self, const_reactions, size):
        """Apply a set of reactions in form of a list of reactions, where each
        reaction = {'iteration':iteration, 'reactants':[reactants], 'products':[products]}

        Quantities of each reactant and each product are assumed to be 1

        :param size: the maximum number of iterations that will be found in the final population.
        We conduct a sample at regular intervals. If None then no sampling is done.
        """

        reactions = const_reactions[:]
        logging.debug("Calculating unique elements...")
        elements = set(self.get_items())  # COPY initial items
        for reaction in reactions:
            elements.update([product['smiles'] for product in reaction['products']])
        number_of_unique_elements = len(elements)

        logging.debug("Resizing population to {} elements".format(number_of_unique_elements))
        incremental_shape = (0, max(0, number_of_unique_elements - self._population.shape[1]))
        self._expand_population(incremental_shape)

        logging.debug("Applying reactions...this may take some time...")

        count_t = len(set([reaction['t'] for reaction in reactions]))

        if size == 0 or count_t < size:
            self._apply_reactions(reactions)

        else:

            final_t = reactions[-1]['t']
            block_size = final_t / (size * 1.0)

            sample_population = MolecularPopulation(population=self.get_last_slice())  # start with an initial population
            self._population = self._population[0:len(self._t), :]  # because using vstack to add on samples start with just the initial population
            logging.debug("Shape is now {}, and sample block_size is {}".format(self._population.shape, block_size))

            reactions.reverse()
            block_end = block_size

            while len(reactions) > 0:

                reaction = reactions.pop()

                sample_population.set_t(reaction['t'])
                for reactant in reaction['reactants']:
                    sample_population.set_quantity(reactant['smiles'], sample_population.get_quantity(reactant['smiles']) - 1)
                for product in reaction['products']:
                    sample_population.set_quantity(product['smiles'], sample_population.get_quantity(product['smiles']) + 1)

                # print(reaction['t'], block_end)
                if reaction['t'] >= block_end or len(reactions) == 0:

                    block_end += block_size

                    # transfer final state of this block to the real population
                    # at this point we can be sure that the final_population and the full population (in self)
                    # are of the same 'y' dimension (elements) and that any elements in self will be in the same
                    # position in final_population (as it is an incremental copy) therefore we can do the transfer by
                    # 1) updating set_t of self, 2) appending the nparray from final_population to the nparray
                    # of self and 3) replacing the index of self by the updated index of final_population

                    final_population = sample_population.get_last_slice()

                    self._t.append(final_population.get_times()[-1])  # should only be one time in the final_population...
                    self._population = np.vstack((self._population, final_population._population))
                    self._index = final_population._index

                    # logging.info("Full population now of size {}".format(self._population.shape))

                    # reset the sample population, ready for next block of reactions
                    del sample_population
                    sample_population = MolecularPopulation(population=final_population.get_last_slice())

            logging.info("Reactions applied")

    def _get_molecules_matching_pattern(self, pattern):
        """Find all items in the population that contain the provided molecular pattern.

        :param pattern: A Mol descriptor of the pattern we are to attempt to match
        :type pattern: Mol
        :rtype:list of SMILES items in population matching the pattern"""

        matches = [x for x in self._index if Chem.MolFromSmiles(x).HasSubstructMatch(pattern)]
        return matches

    def count_molecules_matching_pattern(self, pattern):
        """Count the number of molecules which contain the given fragment.

        :param smarts_pattern: A Mol descriptor of the pattern we're interested in
        :type smarts_pattern: Mol
        :rtype: int
        """

        quantities = [self.get_quantity(mol) for mol in self._get_molecules_matching_pattern(pattern)]
        return sum(quantities)

    def choose_molecule(self, pattern):
        """Select one reactant molecule (type Mol) with probability proportional to a molecule's share in the population.

        :param pattern: pattern of interest
        :type pattern: Mol
        :rtype: Mol
        """

        logging.debug("Choosing molecule matching {}".format(Chem.MolToSmiles(pattern)))
        candidates = self._get_molecules_matching_pattern(pattern)
        upper_bound = 0
        for mol in candidates:
            upper_bound += self.get_quantity(mol)
        bound = random.random() * upper_bound
        cumulative_sum = 0
        for mol in candidates:
            cumulative_sum += self.get_quantity(mol)
            if cumulative_sum > bound:
                logging.debug("Molecule chosen")
                return Chem.MolFromSmiles(mol)
