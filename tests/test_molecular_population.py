"""
Created on 10/01/2013

@author: thom
"""
import unittest
from molecular_population import MolecularPopulation
from molecule import Molecule
from rdkit.Chem import AllChem as Chem


class MolecularPopulationTest(unittest.TestCase):

    def setUp(self):
        self._xml = """<?xml version="1.0" ?>
<population>
    <setvalue item="CC(=O)O" quantity="100" />
    <setvalue item="CN" quantity="50" food='true'/>
    <setvalue item="O=C(O)CCOC(=O)O" quantity="25" />
    <setvalue item="C=CO" quantity="100" food='True'/>
    <setvalue item="C=CC(=C)N" quantity="50" food='false'/>
</population>
        """

    def tearDown(self):
        del self._xml

    def test_copy_constructor(self):
        p = MolecularPopulation(xml=self._xml)
        q = MolecularPopulation(population=p)
        self.assertEqual(len(p.get_items()), len(q.get_items()))
        self.assertEqual(len(p.get_times()), len(q.get_times()))

    def test_population_count_molecules(self):
        p = MolecularPopulation(xml=self._xml)
        self.assertEqual(100, p.count_molecules_matching_pattern(Chem.MolFromSmiles("CN")))
        self.assertEqual(25, p.count_molecules_matching_pattern(Chem.MolFromSmiles('O=C(O)CCOC(=O)O')))
        self.assertEqual(225, p.count_molecules_matching_pattern(Chem.MolFromSmiles('O')))

    def test_population_count_molecules_R1(self):  # regresssion - handle case where population doesn't contain the requested molecule
        p = MolecularPopulation(xml=self._xml)
        self.assertEqual(0, p.count_molecules_matching_pattern(Chem.MolFromSmiles('Br')))

    def test_apply_reactions(self):
        xml = """<?xml version="1.0" ?>
<population>
    <setvalue item="[H][H]" quantity="100" />
    <setvalue item="O=O" quantity="100" />
    <setvalue item="O" quantity="200" /> <!-- H2O -->
    <setvalue item="[O-][N+](=O)[N+]([O-])=O" quantity="100" />
    <setvalue item="N(=O)[O]" quantity="100" />
    <setvalue item="O=C=O" quantity="200" />
</population>
        """

        p1 = MolecularPopulation(xml=xml)
        items = p1.get_items()  # KineticMolecules
        self.assertEqual(6, len(items))

        reactions = []
        reactions.append({'iteration': 1000, 't': 0, 'reactants': [{'smiles': items[0]}, {'smiles': items[1]}], 'products': [{'smiles': items[2]}]})
        reactions.append({'iteration': 1001, 't': 1, 'reactants': [{'smiles': items[0]}, {'smiles': items[2]}], 'products': [{'smiles': items[3]}]})
        reactions.append({'iteration': 1002, 't': 1, 'reactants': [{'smiles': Chem.MolToSmiles(Molecule(items[4]))}, {'smiles': items[1]}], 'products': [{'smiles': items[2]}]})
        reactions.append({'iteration': 1003, 't': 2, 'reactants': [{'smiles': items[3]}, {'smiles': items[1]}], 'products': [{'smiles': items[2]}]})
        reactions.append({'iteration': 1005, 't': 3, 'reactants': [{'smiles': items[4]}, {'smiles': items[1]}], 'products': [{'smiles': items[2]}]})

        p1 = MolecularPopulation(xml=xml, reactions=reactions)
        self.assertEqual([0, 1, 2, 3], p1.get_times())
        self.assertEqual(6, len(p1.get_items()))  # Test for BI76 - if we get 7 we're wrong...

        p2 = MolecularPopulation(xml=xml, reactions=reactions, size=2)
        self.assertEqual(6, len(p2.get_items()))
        self.assertEqual([0, 2, 3], p2.get_times())
        self.assertEqual(3, p2.get_times()[-1])

        for item in p1.get_items():
            self.assertEqual(p1.get_quantity(item), p2.get_quantity(item))
            self.assertEqual(p1.get_quantity(item), p2.get_quantity(Chem.MolToSmiles(Molecule(item))))  # Test for BI76

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
