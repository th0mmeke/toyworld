"""
Created on 27/04/2013

@author: thom
"""
import unittest

from rdkit.Chem import AllChem as Chem
from rdkit.rdBase import DisableLog, EnableLog

from molecule import Molecule
from kinetic_molecule import KineticMolecule
from chemistry_model.chemistry_factory import ChemistryFactory
from chemistry_model.emergent_reactions import EmergentReactions


class TestEmergentReactions(unittest.TestCase):

    def setUp(self):
        DisableLog('rdApp.error')
        self._chemistry = ChemistryFactory.new()
        self._reactions_model = EmergentReactions(self._chemistry)

    def tearDown(self):
        EnableLog('rdApp.error')

    def testDuplicateOptions(self):
        self.assertEqual(8, len(self._reactions_model.get_reaction_options([KineticMolecule('O=C=O', kinetic_energy=300), KineticMolecule('C', kinetic_energy=200)])))

    def testAtomsMaintained(self):
        mol = Molecule('O', kinetic_energy=0)
        for rxn in self._reactions_model.get_reaction_options([mol]):
            products = rxn.fire()
            self.assertEqual(mol.GetNumAtoms(onlyExplicit=False), sum(m.GetNumAtoms(onlyExplicit=False) for m in products))
            for m in products:
                self.assertEqual(Chem.AddHs(m).GetNumAtoms(), m.GetNumAtoms())

    def testGetReactionOptions(self):
        # Use KE=300 to give adequate energy for all bond options to be available

        # mol = Molecule(Chem.MolFromSmiles('[CH2-2].[CH2-2]'))
        # options should include formation of single C-C bond, and double C=C bond
        options = self._reactions_model.get_reaction_options([Molecule('O=C=O', kinetic_energy=300)])
        self.assertEqual(4, len(options))  # four options: break and drop to single bond from O=C bonds

        # This requires three mini-steps - first, recognition that the ion has free unbonded electrons, second that can use those to form
        # bond between O and H, and third, that O and H are in different components of the combined molecule
        # This only applies if two molecules - if only [OH-] without [H] or other molecule then don't have this option
        options = self._reactions_model.get_reaction_options([Molecule('[H+].[OH-]', kinetic_energy=300)])
        self.assertEqual(2, len(options))  # bond from H+ to O- (tricky - needs ion manipulation), break bond between O and H-; no bond possible between H+ and H-!

        options = self._reactions_model.get_reaction_options([Molecule('O', kinetic_energy=300)])
        self.assertEqual(2, len(options))  # two options, both breaks of H bonds

        mol = Molecule('[OH-]', kinetic_energy=300)
        self.assertEqual(2, mol.GetNumAtoms(onlyExplicit=False))
        self.assertEqual(1, mol.GetNumBonds(onlyHeavy=False))
        options = self._reactions_model.get_reaction_options([mol])
        self.assertEqual(1, len(options))  # break H bond

        self.assertEqual(3, len(self._reactions_model.get_reaction_options([Molecule('[C].[C]', kinetic_energy=300)])))  # 3 types of bond formation - single, double, triple
        self.assertEqual(3, len(self._reactions_model.get_reaction_options([Molecule('[C].[C]', kinetic_energy=0)])))  # 3 types of bond formation - single, double, triple

        self.assertEqual(6, len(self._reactions_model.get_reaction_options([Molecule('C=C', kinetic_energy=300)])))  # five complete breaks, and one drop from double to single
        self.assertEqual(2, len(self._reactions_model.get_reaction_options([Molecule('[O].[O]', kinetic_energy=300)])))  # oxygen ions...pretty rare in nature - single and double bonds
        self.assertEqual(1, len(self._reactions_model.get_reaction_options([Molecule('[H].[O]', kinetic_energy=300)])))  # oxygen ion and proton...pretty rare in nature

    def testEnergyChanges(self):
        mol = Molecule('[H+].[OH-]', kinetic_energy=300)
        initial_PE = mol.get_potential_energy(self._chemistry)
        options = self._reactions_model.get_reaction_options([mol])
        # self.assertEqual(2, len(options))  # bond from H+ to O- (tricky - needs ion manipulation), break bond between O and H-; no bond possible between H+ and H-!
        for rxn in options:
            final_PE = sum([mol.get_potential_energy(self._chemistry) for mol in rxn.fire()])
            self.assertEqual(initial_PE + rxn.get_energy_delta(), final_PE)

    def testMultipleProductMoleculesFromReaction(self):
        # Expect reations to return a list of separate products molecules, not one combined product
        mol = Molecule('[OH-]', kinetic_energy=300)
        options = self._reactions_model.get_reaction_options([mol])
        self.assertEqual(1, len(options))  # break H bond
        self.assertEqual(2, len(options[0].fire()))

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testInit']
    unittest.main()
