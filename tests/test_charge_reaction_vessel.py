"""
Created on 27/04/2013

@author: thom
"""
import unittest

import pymunk as pm
import numpy as np

from chemistry_model.chemistry_factory import ChemistryFactory
from charged_molecule import ChargedMolecule
from reactor_model.charge_reaction_vessel import ChargeReactionVessel


class TestChargeReactionVessel(unittest.TestCase):

    def setUp(self):
        # logging.basicConfig(level=logging.INFO)
        self._vessel = ChargeReactionVessel(chemistry=ChemistryFactory.new())

    def test_get_charge_forces(self):
        self.assertEqual(pm.Vec2d(0, 0), self._vessel._get_force_vector(pm.Vec2d(0, 0), pm.Vec2d(1, 0), 0))
        force1 = self._vessel._get_force_vector(pm.Vec2d(0, 0), pm.Vec2d(2, 0), -1)  # - = attraction
        self.assertTrue(force1[0] > 0)  # force of p2 on p1 => movement to right
        self.assertAlmostEqual(0, force1[1])

        force2 = self._vessel._get_force_vector(pm.Vec2d(0, 0), pm.Vec2d(1, 0), -1)  # half the distance...
        self.assertAlmostEqual(force1[0] * 4.0, force2[0])
        self.assertAlmostEqual(0, force2[1])

        p2 = pm.Vec2d(3, 4)
        force3 = self._vessel._get_force_vector(pm.Vec2d(0, 0), p2, -1)  # + = repulsion
        d_squared = p2.get_length_sqrd()

        f = pm.Vec2d(p2.normalized()[0] / d_squared, p2.normalized()[1] / d_squared)
        self.assertAlmostEqual(f[0], force3[0])
        self.assertAlmostEqual(f[1], force3[1])

    def test_position(self):
        mol1 = ChargedMolecule("[H]")
        mol1.set_orientation(0)
        for shape in mol1.body.shapes:
            shape_position = mol1.body.local_to_world(shape.offset)  # wc
            self.assertEqual(pm.Vec2d(0, 0), shape_position)
        mol1.set_orientation(1)
        for shape in mol1.body.shapes:
            shape_position = mol1.body.local_to_world(shape.offset)  # wc
            self.assertEqual(pm.Vec2d(0, 0), shape_position)

        body = pm.Body()
        body.angle = 0
        body.position = pm.Vec2d(10, 1)
        shape = pm.Circle(body, 1, (1, 0))
        shape_position = body.local_to_world(shape.offset)  # wc
        self.assertEqual(pm.Vec2d(11, 1), shape_position)
        body.angle = np.pi
        shape_position = body.local_to_world(shape.offset)  # wc
        self.assertAlmostEqual(9.0, shape_position[0])
        self.assertAlmostEqual(1.0, shape_position[1])


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testInit']
    unittest.main()
