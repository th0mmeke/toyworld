"""
Created on 6/05/2013

@author: thom
"""

from evaluator import Evaluator
from molecule import Molecule
import logging

from rdkit.Chem import Draw
from rdkit.Chem import AllChem as Chem


class DrawMolecules(Evaluator):

    """Draw to a file the longest molecule found in the final population.
    """

    @classmethod
    def evaluate(self, results, save_filebase, min_length=10, **kwargs):

        max_atoms = 0
        for smiles in results['final_population'].get_items():
            mol = Molecule(smiles)
            if mol.GetNumAtoms() > max_atoms:
                max_mol = mol
                max_atoms = mol.GetNumAtoms()
        logging.info("Longest molecule found is {}".format(Chem.MolToSmiles(max_mol)))
        Chem.Compute2DCoords(max_mol)
        Draw.MolToFile(max_mol, save_filebase + "-molecule.png", size=(2000, 2000))
