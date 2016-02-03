"""
Created on 6/05/2013

@author: thom
"""

from evaluator import Evaluator
from molecule import Molecule
from molecular_population import MolecularPopulation
import logging
import os

from rdkit.Chem import Draw
from rdkit.Chem import AllChem as Chem


class DrawMolecules(Evaluator):

    """Draw to a file the longest molecule found in the final population.
    """

    def get_result_titles(self):
        return []

    @classmethod
    def evaluate(self, results_filename, **kwargs):

        results = Evaluator.load_results(results_filename)
        population = MolecularPopulation(population=results['initial_population'], reactions=results['reactions'], size=100)
        final_items = set(item for item in population.get_items() if population.get_quantity(item) > 0)

        max_atoms = 0
        for smiles in final_items:
            mol = Molecule(smiles)
            if mol.GetNumAtoms() > max_atoms:
                max_mol = mol
                max_atoms = mol.GetNumAtoms()
        logging.info("Longest molecule found is {}".format(Chem.MolToSmiles(max_mol)))
        Chem.Compute2DCoords(max_mol)
        Draw.MolToFile(max_mol, os.path.splitext(results_filename)[0] + "-molecule.png", size=(2000, 2000))
