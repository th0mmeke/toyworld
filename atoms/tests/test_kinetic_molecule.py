"""
Created on 27/04/2013

@author: thom
"""
import copy
import math
import string
import unittest

from rdkit.Chem import AllChem as Chem
from rdkit.rdBase import DisableLog, EnableLog

from atoms.kinetic_molecule import KineticMolecule
from util.ulps import Ulps


class TestKineticMolecule(unittest.TestCase):

    def setUp(self):
        DisableLog('rdApp.error')

    def tearDown(self):
        EnableLog('rdApp.error')

    def testInit(self):
        mol = KineticMolecule("O=C=O")
        self.assertEqual(3, mol.GetNumAtoms())

        mol = KineticMolecule(Chem.MolFromSmiles("O=C=O"))
        self.assertEqual(3, mol.GetNumAtoms())
        self.assertEqual(44.009, mol.get_mass())

        self.assertEqual(18.015, KineticMolecule(Chem.MolFromSmiles("[H]O[H]")).get_mass())
        self.assertEqual(18.015, KineticMolecule(Chem.MolFromSmiles("O")).get_mass())

        mol = KineticMolecule("O", kinetic_energy=30)
        s = mol.get_speed()
        v = mol.get_velocity()
        e = mol.get_kinetic_energy()
        self.assertTrue(Ulps.almost_equal(30, e))
        self.assertAlmostEqual(s, math.sqrt(v[0] * v[0] + v[1] * v[1]))
        self.assertAlmostEqual(e, 0.5 * mol.get_mass() * s * s)

        mol = KineticMolecule("O")
        mol.set_kinetic_energy(0)

    def testGetSet(self):

        mol = KineticMolecule('O')
        v = 3, 4
        mol.set_velocity(*v)
        speed = math.sqrt(3 * 3 + 4 * 4)
        self.assertEqual(speed, mol.get_speed())

    def test_global_id(self):
        mol1 = KineticMolecule('O')
        mol2 = KineticMolecule('O')
        self.assertNotEqual(mol1.global_id, mol2.global_id)
        del mol2
        mol2 = KineticMolecule('O')
        self.assertNotEqual(mol1.global_id, mol2.global_id)

    def test_deepcopy(self):
        mol = KineticMolecule("O=C=O", 30)
        mol.set_velocity(5, 6)
        mol2 = copy.deepcopy(mol)
        self.assertEqual(mol.get_kinetic_energy(), mol2.get_kinetic_energy())
        self.assertEqual(mol.get_velocity(), mol2.get_velocity())
        self.assertEqual(mol.get_speed(), mol2.get_speed())
        self.assertEqual(Chem.MolToSmiles(mol), Chem.MolToSmiles(mol2))
        self.assertEqual(mol.GetNumAtoms(), mol2.GetNumAtoms())

        mol = KineticMolecule("[H][O-]", 30)
        mol.set_velocity(8, 6)
        mol2 = copy.deepcopy(mol)
        self.assertEqual(mol.get_kinetic_energy(), mol2.get_kinetic_energy())
        self.assertEqual(Chem.MolToSmiles(mol), Chem.MolToSmiles(mol2))
        self.assertEqual(mol.GetNumAtoms(), mol2.GetNumAtoms())

        mol = KineticMolecule("O=C=O", 30, canonize=False)
        mol.set_velocity(8, 6)
        mol2 = copy.deepcopy(mol)
        self.assertEqual(mol.get_kinetic_energy(), mol2.get_kinetic_energy())
        self.assertEqual(Chem.MolToSmiles(mol), Chem.MolToSmiles(mol2))
        self.assertEqual(mol.GetNumAtoms(), mol2.GetNumAtoms())
        self.assertNotEqual(mol.global_id, mol2.global_id)  # we use deepcopy in get_interaction_options to bounce molecules when no reaction is possible

    def testSplitAndCombineMolecule(self):
        split_smiles = ('[H]O[H]', '[H]O[H]', 'O=C=O')
        mols = [KineticMolecule(smiles, 10) for smiles in split_smiles]

        initial_ke = 0
        for mol in mols:
            mol.set_velocity(8, 6)
            initial_ke += mol.get_kinetic_energy()
        combined_mols = mols[0].combine_molecules(mols)

        target = KineticMolecule(string.join(split_smiles, "."), 10)
        self.assertEqual(Chem.MolToSmiles(target), Chem.MolToSmiles(combined_mols))
        split_mols = combined_mols.split_molecule()
        self.assertEqual(split_smiles, tuple(Chem.MolToSmiles(mol) for mol in split_mols))

        # Now check energy is conserved through both combining and splitting
        self.assertAlmostEqual(initial_ke, combined_mols.get_kinetic_energy())
        split_ke = 0
        for mol in split_mols:
            split_ke += mol.get_kinetic_energy()
        self.assertAlmostEqual(initial_ke, split_ke)

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testInit']
    unittest.main()
