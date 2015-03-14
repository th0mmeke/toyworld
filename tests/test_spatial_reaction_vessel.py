"""
Created on 27/04/2013

@author: thom
"""
import unittest
import logging
import random

from chemistry_model.chemistry_factory import ChemistryFactory
from kinetic_molecule import KineticMolecule
from reactor_model.spatial_reaction_vessel import SpatialReactionVessel
from ULPS import Float_t


class TestSpatialReactionVessel(unittest.TestCase):

    def setUp(self):
        logging.basicConfig(level=logging.INFO)
        self._vessel = SpatialReactionVessel(chemistry=ChemistryFactory.new())

    def _init_vessel(self, num_molecules=2):
        """Guaranteed not to bounce...each molecule (other than first one) travels shortest path between two points"""

        mr = self._vessel._base_molecule_radius
        mol0 = KineticMolecule("O=C=O", internal_energy=10)
        mol0.set_velocity(0, 0, 0)
        p0 = [random.uniform(-1 + mr, 1 - mr), random.uniform(-1 + mr, 1 - mr), random.uniform(-1 + mr, 1 - mr)]
        mols = {mol0: p0}

        delta_t = 1E-2
        collisions = {}

        for i in range(num_molecules - 1):
            p1 = [random.uniform(-1 + mr, 1 - mr), random.uniform(-1 + mr, 1 - mr), random.uniform(-1 + mr, 1 - mr)]
            num_steps = random.randint(15, 100)
            delta_x = (p0[0] - p1[0]) / num_steps / delta_t
            delta_y = (p0[1] - p1[1]) / num_steps / delta_t
            delta_z = (p0[2] - p1[2]) / num_steps / delta_t

            for p, v, actual_p in zip(p1, (delta_x, delta_y, delta_z), p0):
                assert Float_t.almost_equal(p + v * num_steps * delta_t, actual_p)
            mol1 = KineticMolecule("O=C=O", internal_energy=10)
            mol1.set_velocity(delta_x, delta_y, delta_z)
            mols[mol1] = p1
            try:
                collisions[num_steps].extend(mol1)
            except:
                collisions[num_steps] = [mol1]

        self._vessel.initialize(chemistry=ChemistryFactory.new())
        self._vessel.add_molecules(mols)

        return delta_t, mol0, collisions

    def test_get_state(self):
        #  self.assertEqual({}, self._vessel.get_state())
        mols = [KineticMolecule("C", kinetic_energy=10), KineticMolecule("O", kinetic_energy=10)]
        self._vessel.add_molecules(mols)
        state = self._vessel.get_state()
        self.assertEqual(2, len(state))
        for mol in mols:
            self.assertIn(mol.global_id, state['locations'].keys())

    def test_add_molecules(self):
        mols = [KineticMolecule("C", kinetic_energy=10), KineticMolecule("O", kinetic_energy=10)]
        self._vessel.add_molecules(mols)
        self.assertEqual(len(mols), self._vessel.get_number_molecules())
        for mol in mols:
            self.assertTrue(mol in self._vessel.get_molecules())

        self._vessel.add_molecules([KineticMolecule("C", kinetic_energy=0)])  # should be able to add a stationary molecule

    def test_discover_reaction(self):
        # Test handled by post-condition assertions in discover_reaction
        for i in range(10):  # stochastic test - multiple reaction options
            mol0 = KineticMolecule("O=C=O", kinetic_energy=100, internal_energy=10)
            mol1 = KineticMolecule("O", kinetic_energy=100, internal_energy=5)
            self._vessel.discover_reaction([mol0, mol1])

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testInit']
    unittest.main()
