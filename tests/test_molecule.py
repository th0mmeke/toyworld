"""
Created on 27/04/2013

@author: thom
"""
import unittest
import string
import copy

from rdkit.Chem import AllChem as Chem
from rdkit.rdBase import DisableLog, EnableLog

from molecule import Molecule
from chemistry_model.chemistry_factory import ChemistryFactory


class Test(unittest.TestCase):

    def setUp(self):
        DisableLog('rdApp.error')

    def tearDown(self):
        EnableLog('rdApp.error')

    def testInit(self):
        mol = Molecule("O=C=O", internal_energy=30)
        self.assertEqual(3, mol.GetNumAtoms())
        self.assertEqual(30, mol.get_internal_energy())

        mol = Molecule(Chem.MolFromSmiles("O=C=O"), internal_energy=30)
        self.assertEqual(3, mol.GetNumAtoms())
        self.assertEqual(30, mol.get_internal_energy())
        self.assertEqual(44.009, mol.get_mass())

        with self.assertRaises(ValueError):
            Molecule(Chem.MolFromSmiles("O"), -10)

        mol = Molecule('[O][H-]')
        self.assertEqual('[H][O]', Chem.MolToSmiles(mol))  # RDKit converts to this...
        mol = Molecule('[H+].[OH-]')
        self.assertEqual('[H+].[H][O-]', Chem.MolToSmiles(mol))  # RDKit converts to this...

    def test_global_id(self):
        mol1 = Molecule('O')
        mol2 = Molecule('O')
        self.assertNotEqual(mol1.global_id, mol2.global_id)
        del mol2
        mol2 = Molecule('O')
        self.assertNotEqual(mol1.global_id, mol2.global_id)

    def test_deepcopy(self):
        mol = Molecule("O=C=O", kinetic_energy=30)
        mol2 = copy.deepcopy(mol)
        self.assertEqual(mol.get_internal_energy(), mol2.get_internal_energy())
        self.assertEqual(Chem.MolToSmiles(mol), Chem.MolToSmiles(mol2))
        self.assertEqual(mol.GetNumAtoms(), mol2.GetNumAtoms())

        mol = Molecule("[H][O-]", kinetic_energy=30)
        mol2 = copy.deepcopy(mol)
        self.assertEqual(mol.get_internal_energy(), mol2.get_internal_energy())
        self.assertEqual(Chem.MolToSmiles(mol), Chem.MolToSmiles(mol2))
        self.assertEqual(mol.GetNumAtoms(), mol2.GetNumAtoms())

        mol = Molecule("O=C=O", kinetic_energy=30, canonize=False)
        mol2 = copy.deepcopy(mol)
        self.assertEqual(mol.get_internal_energy(), mol2.get_internal_energy())
        self.assertEqual(Chem.MolToSmiles(mol), Chem.MolToSmiles(mol2))
        self.assertEqual(mol.GetNumAtoms(), mol2.GetNumAtoms())

    def testPotentialEnergy(self):
        chem = ChemistryFactory.new()
        self.assertEqual(0, Molecule('[O].[O]').get_potential_energy(chem))
        self.assertEqual(chem.get_bond_energy(Chem.Atom('O'), Chem.Atom('O'), end_bond_type=2), Molecule('[O]=[O]').get_potential_energy(chem))
        self.assertTrue(Molecule('[O]=[O]').get_potential_energy(chem) < Molecule('[O].[O]').get_potential_energy(chem))  # bonds have NEGATIVE energy, so adding bonds REDUCES PE
        O1H_bond = chem.get_bond_energy(Chem.Atom('O'), Chem.Atom('H'), end_bond_type=1)
        O2O_bond = chem.get_bond_energy(Chem.Atom('O'), Chem.Atom('O'), end_bond_type=2)
        self.assertEqual(O2O_bond, Molecule('[O]=[O]').get_potential_energy(chem))
        self.assertEqual(2 * O1H_bond, Molecule('[H]O[H]').get_potential_energy(chem))
        self.assertEqual(O1H_bond, Molecule('[H+].[OH-]').get_potential_energy(chem))

    def testGetStronglyConnectedComponents(self):
        mol = Molecule("O=C=O.O")
        self.assertEqual(6, mol.GetNumAtoms())
        self.assertEqual([[0, 1, 2], [3, 4, 5]], mol._get_strongly_connected_components())
        mol = Molecule("[H]")
        self.assertEqual([[0]], mol._get_strongly_connected_components())

    def testSameComponents(self):
        mol = Molecule("O=C=O.[H]O[H]")
        self.assertFalse(mol.same_component(0, 5))
        self.assertTrue(mol.same_component(0, 1))
        mol = Molecule("CO.N", components=[set(range(6)), set(range(6, 10))])
        self.assertTrue(mol.same_component(0, 2))
        self.assertFalse(mol.same_component(0, 8))

    def testAssignFormalCharge(self):
        # [H]O[H] -> [OH-]+[H+] Hydroxyl plus proton
        mol = Molecule('[H].[OH]')
        mol._assign_formal_charge()
        self.assertEqual('[H+].[H][O-]', Chem.MolToSmiles(mol))

    def testSplitAndCombineMolecule(self):
        split_smiles = ('[H]O[H]', '[H]O[H]', 'O=C=O')
        mols = [Molecule(smiles, kinetic_energy=10) for smiles in split_smiles]
        combined_mols = mols[0].combine_molecules(mols)
        target = Molecule(string.join(split_smiles, "."), 10)
        self.assertEqual(Chem.MolToSmiles(target), Chem.MolToSmiles(combined_mols))
        self.assertEqual(split_smiles, tuple(Chem.MolToSmiles(mol) for mol in combined_mols.split_molecule()))

        split_smiles = ('[H+]', '[OH-]')
        mols = [Molecule(smiles, kinetic_energy=10) for smiles in split_smiles]
        combined_mols = mols[0].combine_molecules(mols)
        target = '[H+].[H][O-]'
        self.assertEqual(target, Chem.MolToSmiles(combined_mols))

        split_smiles = ['[H+].[OH-]']
        mols = [Molecule(smiles, kinetic_energy=10) for smiles in split_smiles]
        combined_mols = mols[0].combine_molecules(mols)
        target = '[H+].[H][O-]'
        self.assertEqual(target, Chem.MolToSmiles(combined_mols))

    def testRegressionCombineMoleculeComponents(self):
        split_smiles = ('[H+]', '[OH-]')
        mols = [Molecule(smiles, kinetic_energy=10) for smiles in split_smiles]
        combined_mols = mols[0].combine_molecules(mols)
        self.assertFalse(combined_mols.same_component(0, 2))

        mol = Molecule('[H+].[OH-]', kinetic_energy=0)
        self.assertFalse(mol.same_component(0, 2))

        mol2 = mol.combine_molecules([mol])  # must correctly transfer across *existing* components from mol
        self.assertFalse(mol2.same_component(0, 2))

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testInit']
    unittest.main()
