"""
Created on 3/05/2013

@author: thom
"""
import unittest

from rdkit.Chem import AllChem as Chem

from reactions.reaction import Reaction


class Test(unittest.TestCase):

    def testInit(self):
        result = Reaction(Chem.MolFromSmiles('O=C=O'), 40)
        self.assertEqual(40, result.get_energy_delta())

if __name__ == "__main__":
    unittest.main()
