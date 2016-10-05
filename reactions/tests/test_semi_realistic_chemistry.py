"""
Created on 27/04/2013

@author: thom
"""
import unittest
import xml.etree.cElementTree as ElementTree

from rdkit.Chem import AllChem as Chem

from atoms.molecule import Molecule
from reactions.chemistry_factory import ChemistryFactory
from parameters import Parameters


class TestSemiRealisticChemistry(unittest.TestCase):

    def testBondEnergy(self):
        # energy is energy REQUIRED => - means releases energy
        chem = ChemistryFactory.new()
        self.assertEqual(-38.4, chem.get_bond_energy(Chem.Atom('N'), Chem.Atom('N'), end_bond_type=1))  # create single bond
        self.assertEqual(-38.4, chem.get_bond_energy(Chem.Atom('N'), Chem.Atom('N'), start_bond_type=0, end_bond_type=1))  # create single bond
        self.assertEqual(38.4, chem.get_bond_energy(Chem.Atom('N'), Chem.Atom('N'), start_bond_type=1, end_bond_type=0))  # destroy single bond
        self.assertEqual(149 - 38.4, chem.get_bond_energy(Chem.Atom('N'), Chem.Atom('N'), start_bond_type=2, end_bond_type=1))  # from double to single
        self.assertEqual(38.4 - 149, chem.get_bond_energy(Chem.Atom('N'), Chem.Atom('N'), start_bond_type=1, end_bond_type=2))  # from single to double
        self.assertEqual(38.4, chem.get_bond_energy(Chem.Atom('N'), Chem.Atom('N'), start_bond_type=1))  # delete single bond

    def testBondEnergyWithCustomChemistry(self):
        chem = ChemistryFactory.new(parameters=ElementTree.fromstring('<xml></xml>'))  # non-null, to check non-null parameter pathway - should default to standard bond energies
        self.assertEqual(-38.4, chem.get_bond_energy(Chem.Atom('N'), Chem.Atom('N'), end_bond_type=1))  # create single bond
        bond_formation_xml = '<BondFormationEnergies><Single>10</Single><Double>20</Double><Triple>30</Triple></BondFormationEnergies>'
        bond_break_xml = '<BondBreakEnergies><Single>11</Single><Double>22</Double><Triple>33</Triple></BondBreakEnergies>'
        chem1 = ChemistryFactory.new(parameters=Parameters(ElementTree.fromstring('<xml>' + bond_formation_xml + '</xml>')))
        self.assertEqual(-10, chem1.get_bond_energy(Chem.Atom('N'), Chem.Atom('N'), end_bond_type=1))  # create single bond
        self.assertEqual(20, chem1.get_bond_energy(Chem.Atom('N'), Chem.Atom('N'), start_bond_type=2))  # destroy double bond
        chem2 = ChemistryFactory.new(parameters=Parameters(ElementTree.fromstring('<xml>' + bond_formation_xml + bond_break_xml + '</xml>')))
        self.assertEqual(-10, chem2.get_bond_energy(Chem.Atom('N'), Chem.Atom('N'), end_bond_type=1))  # create single bond
        self.assertEqual(22, chem2.get_bond_energy(Chem.Atom('N'), Chem.Atom('N'), start_bond_type=2))  # destroy double bond
        # Check safe handling of quadruple bonds
        self.assertEqual(33, chem2.get_bond_energy(Chem.Atom('N'), Chem.Atom('N'), start_bond_type=4))  # destroy double bond

    def testBondEnergyManyAtoms(self):
        chem = ChemistryFactory.new(parameters=Parameters(ElementTree.fromstring('<xml><atoms><Atom>N</Atom><Atom>S</Atom><Atom>Q</Atom></atoms></xml>')))
        actual = chem.get_bond_energy(Chem.Atom('S'), Chem.Atom('S'), start_bond_type=2)
        self.assertEqual(actual, chem.get_bond_energy(Chem.Atom('S'), Chem.Atom('S'), start_bond_type=2))  # same as for known bond energy
        with self.assertRaises(Exception):
            chem.get_bond_energy(Chem.Atom('O'), Chem.Atom('O'), start_bond_type=1)
        with self.assertRaises(Exception):
            chem.get_bond_energy(Chem.Atom('S'), Chem.Atom('O'), start_bond_type=1)

    def testGetBondPotential(self):
        chem = ChemistryFactory.new()
        mol = Molecule('[CH2-2].[CH2-2]')
        self.assertEqual(4, chem.get_bond_potential(mol.GetAtoms()[0]))

        mol = Molecule('O')  # H2O
        self.assertEqual(3, mol.GetNumAtoms())
        self.assertEqual(8, mol.GetAtoms()[0].GetAtomicNum())
        self.assertEqual(0, chem.get_bond_potential(mol.GetAtoms()[0]))

        mol = Molecule('O')  # H2O, implicit Hs
        self.assertEqual(3, mol.GetNumAtoms())  # implicit no longer...
        self.assertEqual(8, mol.GetAtoms()[0].GetAtomicNum())
        self.assertEqual(0, chem.get_bond_potential(mol.GetAtoms()[0]))

        mol = Molecule('[H]')
        self.assertEqual(1, mol.GetAtoms()[0].GetAtomicNum())
        self.assertEqual(1, chem.get_bond_potential(mol.GetAtoms()[0]))

        mol = Molecule('O=C=O')  # CO2 - bond potentials of zero all round (full octets)
        for atom in mol.GetAtoms():
            self.assertEqual(0, chem.get_bond_potential(atom))

        mol = Molecule('[OH-]')  # = [H][O-] Hydroxl anion
        self.assertEqual(8, mol.GetAtoms()[0].GetAtomicNum())
        self.assertEqual(1, mol.GetAtoms()[1].GetAtomicNum())
        self.assertEqual(2, chem.get_bond_potential(mol.GetAtoms()[0]))
        self.assertEqual(0, chem.get_bond_potential(mol.GetAtoms()[1]))

        mol = Molecule('[H].[OH-]')
        self.assertEqual(1, mol.GetAtoms()[0].GetAtomicNum())  # the H in [OH-]
        self.assertEqual(1, chem.get_bond_potential(mol.GetAtoms()[0]))
        self.assertEqual(1, mol.GetAtoms()[2].GetAtomicNum())  # the H in [OH-]
        self.assertEqual(0, chem.get_bond_potential(mol.GetAtoms()[2]))

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testInit']
    unittest.main()
