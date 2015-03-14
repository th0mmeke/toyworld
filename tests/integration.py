"""
Created on 10 Sep 2013

@author: thom
"""
import unittest
import os
import logging

import xml.etree.cElementTree as ElementTree

from doe import FactorialDesign
from runner import Runner
import config


class TestIntegration(unittest.TestCase):

    def setUp(self):
        logging.basicConfig(level=logging.FATAL)

    def test_DOE(self):
        directory = os.path.join(config.DataDir, 'test')
        experiment = 'test_integration'
        seed = 123123
        iterations = 5
        population = os.path.join(config.TestDir, 'test_population.xml')

        kwargs = {'format': '%(asctime)s %(levelname)s:%(message)s', 'datefmt': '%m/%d/%Y %I:%M:%S %p', 'level': getattr(logging, 'DEBUG', None)}
        logging.basicConfig(**kwargs)

        if not os.path.lexists(directory):
            os.makedirs(directory)

        factors = []
        factors.append({'key': 'dimensionality',
                        'Title': 'Dimensionality',
                        'Low': "<Molecule>molecule.Molecule</Molecule><Vessel>reactor_model.aspatial_reaction_vessel.AspatialReactionVessel</Vessel><ReactionVesselDimension>0</ReactionVesselDimension>",
                        'High': "<Molecule>kinetic_molecule.KineticMolecule</Molecule><Vessel>reactor_model.spatial_reaction_vessel.SpatialReactionVessel</Vessel><ReactionVesselDimension>3</ReactionVesselDimension>"})
        factors.append({'key': 'energy', 'Title': 'Energy', 'Low': "<Energy>100</Energy>", 'High': "<Energy>300</Energy>"})
        factors.append({'key': 'foodset', 'Title': 'FoodSet', 'Low': "<FoodSet>False</FoodSet>", 'High': "<FoodSet>True</FoodSet>"})
        factors.append({'key': 'bonds', 'Title': 'Bond Energies', 'Low': "", 'High': "<BondFormationEnergies><Single>50</Single><Double>100</Double><Triple>200</Triple></BondFormationEnergies>"})

        experiment_design = FactorialDesign.design(directory, experiment, seed, iterations, population, factors, repeats=1, recover=False)
        runner = Runner(ElementTree.fromstring(experiment_design), directory)
        runner.run()

        self.assertTrue(True)

    def test_DOE_2D_3D(self):
        directory = os.path.join(config.DataDir, 'test')
        experiment = 'test_integration'
        seed = 123123
        iterations = 5
        population = os.path.join(config.TestDir, 'test_population.xml')

        kwargs = {'format': '%(asctime)s %(levelname)s:%(message)s', 'datefmt': '%m/%d/%Y %I:%M:%S %p', 'level': getattr(logging, 'DEBUG', None)}
        logging.basicConfig(**kwargs)

        if not os.path.lexists(directory):
            os.makedirs(directory)

        factors = []
        factors.append({'key': 'dimensionality',
                        'Title': 'Dimensionality',
                        'Low': "<Molecule>kinetic_molecule.KineticMolecule</Molecule><Vessel>reactor_model.spatial_reaction_vessel.SpatialReactionVessel</Vessel><ReactionVesselDimension>2</ReactionVesselDimension>",
                        'High': "<Molecule>kinetic_molecule.KineticMolecule</Molecule><Vessel>reactor_model.spatial_reaction_vessel.SpatialReactionVessel</Vessel><ReactionVesselDimension>3</ReactionVesselDimension>"})
        factors.append({'key': 'energy', 'Title': 'Energy', 'Low': "<Energy>100</Energy>", 'High': "<Energy>300</Energy>"})
        factors.append({'key': 'foodset', 'Title': 'FoodSet', 'Low': "<FoodSet>False</FoodSet>", 'High': "<FoodSet>True</FoodSet>"})
        factors.append({'key': 'bonds', 'Title': 'Bond Energies', 'Low': "", 'High': "<BondFormationEnergies><Single>50</Single><Double>100</Double><Triple>200</Triple></BondFormationEnergies>"})

        experiment_design = FactorialDesign.design(directory, experiment, seed, iterations, population, factors, repeats=1, recover=False)
        runner = Runner(ElementTree.fromstring(experiment_design), directory)
        runner.run()

        self.assertTrue(True)

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
