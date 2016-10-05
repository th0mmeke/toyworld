"""
Created on 27/04/2013

@author: thom
"""
import math
import unittest

from atoms.molecule import Molecule
from rdkit.rdBase import DisableLog, EnableLog

from atoms.kinetic_molecule import KineticMolecule
from atoms.molecular_population import MolecularPopulation
from reactions.reaction import Reaction
from reactions.semi_realistic_chemistry import SemiRealisticChemistry
from reactors.aspatial_reactor import AspatialReactor
from reactors.spatial_reactor import SpatialReactor


class TestReactor(unittest.TestCase):

    def setUp(self):
        DisableLog('rdApp.error')
        self._chemistry = SemiRealisticChemistry()
        self._population = MolecularPopulation()
        self._population.set_quantity("O", 2)
        self._population.set_quantity("N", 1)
        self._vessel = [SpatialReactor(self._chemistry, self._population), AspatialReactor(self._chemistry, self._population)]

    def tearDown(self):
        EnableLog('rdApp.error')

    def testGetNumberMolecules(self):
        for rv in self._vessel:
            self.assertEqual(self._population.get_population_size(), rv.get_number_molecules())

    def testAddMolecules(self):
        for rv in self._vessel:
            self.assertEqual(0, rv._energy_input)

            mols = [KineticMolecule("C", kinetic_energy=10)]
            rv.add_molecules(mols)
            self.assertEqual(self._population.get_population_size() + len(mols), rv.get_number_molecules())

    def testGetReactants(self):
        pass

    def testAdjustKE(self):
        pass

    def testGetTotalPE(self):
        for rv in self._vessel:
            mols = [Molecule(smiles) for smiles in self._population.get_population()]
            pe = sum([mol.get_potential_energy(self._chemistry) for mol in mols])
            self.assertAlmostEqual(rv.get_total_pe(), pe)

    def testIntegration(self):
        pass

    def Weight(self):
        reaction_energy = 100  # more likely to form bond than break it at low energies
        self.assertTrue(self._vessel.weight(100, reaction_energy) < self._vessel.weight(-100, reaction_energy))
        reaction_energy = 700  # more likely to break bond than form it at high energies
        self.assertTrue(self._vessel.weight(100, reaction_energy) > self._vessel.weight(-100, reaction_energy))

    def RegressionOnlyPossibleReactions(self):
        self.assertEqual(0, self._vessel.weight(100, 100))  # not possible
        self.assertLessEqual(0, self._vessel.weight(-10, 100))  # possible
        self.assertEqual(0, self._vessel.weight(200, 100))  # not possible

    def test_SelectReaction(self):

        def weight(delta, reaction_energy):
            if delta < 0:
                return(-delta)
            else:
                return max(0, reaction_energy - delta)
        for rv in self._vessel:
            reaction_options = [Reaction('BREAK', 100), Reaction('FORM', -100)]
            count = 0
            for i in range(1000):
                rxn = rv._select_reaction(reaction_options, 700)
                count += 'BREAK' == rxn.fire()
            ratio = 1000.0 * weight(100, 700) / (weight(-100, 700) + weight(100, 700))
            self.assertAlmostEqual(math.floor(ratio) / 1000, count / 1000.0, 1)

            count = 0
            for i in range(1000):
                rxn = rv._select_reaction(reaction_options, 100)
                count += 'FORM' == rxn.fire()
            ratio = 1000.0 * weight(-100, 100) / (weight(-100, 100) + weight(100, 100))
            self.assertAlmostEqual(math.floor(ratio) / 1000, count / 1000.0, 1)


if __name__ == "__main__":
    unittest.main()
