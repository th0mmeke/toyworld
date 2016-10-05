"""
Created on 27/04/2013

@author: thom
"""
import random
import unittest

import numpy as np

from atoms.kinetic_molecule import KineticMolecule
from reactions.chemistry_factory import ChemistryFactory
from reactions.emergent_reactions import EmergentReactions
from kinetics_2D import Kinetics2D


class TestKinetics2D(unittest.TestCase):

    def setUp(self):
        # logging.basicConfig(level=logging.DEBUG)
        self._kinetics = Kinetics2D

    def test_xyz_tofrom_radial(self):
        for i in range(0, 10):
            x = random.uniform(-10, 10)
            y = random.uniform(-10, 10)
            x_prime, y_prime = self._kinetics.radial_to_xy(*self._kinetics.xy_to_radial(x, y))
            self.assertAlmostEqual(x, x_prime)
            self.assertAlmostEqual(y, y_prime)
        x, y = 10, 0
        theta, r = self._kinetics.xy_to_radial(x, y)
        self.assertAlmostEqual(0, theta)
        self.assertAlmostEqual(10, r)

    def test_reaction_energy(self):
        mol0 = KineticMolecule("O=C=O", internal_energy=10)
        mol1 = KineticMolecule("O=C=O", internal_energy=5)
        v = random.uniform(-10, 10)
        mol0.set_velocity(v, 0)
        mol1.set_velocity(-v, 0)  # equal and opposite
        self.assertAlmostEqual(mol0.get_kinetic_energy() + mol0.get_kinetic_energy(), self._kinetics.get_CM_energy([mol0, mol0]))  # no collision energy and should ignore internal energy
        self.assertAlmostEqual(0, self._kinetics.get_CM_energy([mol0, mol1]))

    def test_get_CM_velocity(self):
        mol0 = KineticMolecule("O=C=O", internal_energy=10)
        mol1 = KineticMolecule("O=C=O", internal_energy=5)
        v = random.uniform(-10, 10)
        mol0.set_velocity(v, 0)
        mol1.set_velocity(-v, 0)  # equal and opposite
        theta, velocity = self._kinetics.xy_to_radial(*self._kinetics.get_CM_velocity([mol0, mol1]))
        self.assertAlmostEqual(0, velocity)  # velocity = 0 given rounding

    def test_inelastic_collision_angles(self):
        reactant_mols = [KineticMolecule('O=C=O', kinetic_energy=300), KineticMolecule('C', kinetic_energy=200)]
        product_mols = [KineticMolecule('[H]', kinetic_energy=0), KineticMolecule('[H][C]([H])[H]', kinetic_energy=0), KineticMolecule('O=C=O', kinetic_energy=0)]
        # ['[H]', '[H][C]([H])[H]', 'O=C=O']

        out_v, IE = self._kinetics.inelastic_collision(reactant_mols, product_mols, 0)

        for mol, v in zip(product_mols, out_v):
            mol.set_velocity(*v)

        # test that momentum conserved
        in_CM = self._kinetics.get_CM_velocity(reactant_mols)
        out_CM = self._kinetics.get_CM_velocity(product_mols)
        for in_, out_ in zip(in_CM, out_CM):
            self.assertAlmostEqual(in_, out_)

        # test that the out molecules do not all have the same velocity... possible to have correct CM and energy if all along CoM, but boring!
        out_v_radial = [self._kinetics.xy_to_radial(*mol.get_velocity()) for mol in product_mols]
        components = np.array(out_v_radial).transpose()
        self.assertNotEqual(0, np.count_nonzero([np.std(components[i]) for i in range(2)]))

    def test_inelastic_collision(self):
        reactant_mols = [KineticMolecule('O=C=O', kinetic_energy=300), KineticMolecule('C', kinetic_energy=250)]  # different energies so will collide
        options = EmergentReactions(ChemistryFactory.new()).get_reaction_options(reactant_mols)
        for rxn in options:
            self._kinetics.inelastic_collision(reactant_mols, rxn.fire(), 20)

        mol0 = KineticMolecule("O=C=O", kinetic_energy=0, internal_energy=10)
        mol1 = KineticMolecule("O", kinetic_energy=0, internal_energy=5)
        with self.assertRaises(Exception):
            self._kinetics.inelastic_collision([mol0, mol1], [mol0, mol1], 20)  # shouldn't fail even though initial KE = 0

    def test_inelastic_collision_stick(self):
        # Reaction is not possible - if momentum conserved, then energy cannot match without specific energy released in reaction
        reactant_mols = [KineticMolecule('[H+]', kinetic_energy=300), KineticMolecule('[OH-]', kinetic_energy=250)]
        product_mol = [KineticMolecule('[H]O[H]', kinetic_energy=0)]
        self._kinetics.inelastic_collision(reactant_mols, product_mol, 111)

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testInit']
    unittest.main()
