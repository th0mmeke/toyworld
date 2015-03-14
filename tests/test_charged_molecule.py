"""
Created on 27/04/2013

@author: thom
"""
import unittest
import math

from rdkit.Chem import AllChem as Chem
from rdkit.rdBase import DisableLog, EnableLog
import pymunk as pm

from charged_molecule import ChargedMolecule
from ULPS import Float_t


class TestChargedMolecule(unittest.TestCase):

    def setUp(self):
        DisableLog('rdApp.error')

    def tearDown(self):
        EnableLog('rdApp.error')

    def testInit(self):
        mol = ChargedMolecule("O=C=O")
        self.assertEqual(3, mol.GetNumAtoms())

        mol = ChargedMolecule(Chem.MolFromSmiles("O=C=O"))
        self.assertEqual(3, mol.GetNumAtoms())
        self.assertEqual(44.009, mol.get_mass())

        self.assertEqual(18.015, ChargedMolecule(Chem.MolFromSmiles("[H]O[H]")).get_mass())
        self.assertEqual(18.015, ChargedMolecule(Chem.MolFromSmiles("O")).get_mass())

        mol = ChargedMolecule("O", kinetic_energy=30)
        s = mol.get_speed()
        v = mol.get_velocity()
        e = mol.get_kinetic_energy()
        self.assertTrue(Float_t.almost_equal(30, e))
        self.assertAlmostEqual(s, math.sqrt(v[0] * v[0] + v[1] * v[1]))
        self.assertAlmostEqual(e, 0.5 * mol.get_mass() * s * s)

        mol = ChargedMolecule("O")
        mol.set_kinetic_energy(0)

    def testGetSet(self):

        mol = ChargedMolecule('O')
        v = 3, 4
        mol.set_velocity(*v)
        speed = math.sqrt(3 * 3 + 4 * 4)
        self.assertEqual(speed, mol.get_speed())

    def testKE(self):
        mol = ChargedMolecule('O')
        mol.set_kinetic_energy(100)
        self.assertAlmostEqual(100, mol.get_kinetic_energy())
        self.assertAlmostEqual(100, 0.5 * mol.get_mass() * mol.get_speed() * mol.get_speed())
        v = pm.Vec2d(3, 4)
        mol.set_velocity(*v)
        ke = 0.5 * mol.get_mass() * v.get_length() * v.get_length()
        self.assertAlmostEqual(ke, mol.get_kinetic_energy())

    def test_get_clusters(self):
        mol = ChargedMolecule("C1CCCC1")
        mol.get_clusters()

        mol = ChargedMolecule("O")
        mol.get_clusters()

        mol = ChargedMolecule("[O-]")
        mol.get_clusters()

    def test_locations(self):
        mol = ChargedMolecule("C1CCCC1")

        for cluster_ids in mol.get_clusters():
            cluster_centre = mol._get_cluster_centre(cluster_ids)

        mol1 = ChargedMolecule("CC(=O)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@@H]4[C@@]3(CCC(=O)C4)C)C")
        for cluster_ids in mol1.get_clusters():
            cluster_centre = mol1._get_cluster_centre(cluster_ids)

if __name__ == "__main__":
    unittest.main()
