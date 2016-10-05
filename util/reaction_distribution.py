"""
Created on 13 Aug 2013

@author: thom
"""

import logging
import os

from chemistry_factory import ChemistryFactory
from emergent_reactions import EmergentReactions
from physics_factory import PhysicsFactory
from rdkit.Chem import AllChem as Chem

import config
from atoms.molecule import Molecule

if __name__ == '__main__':

    logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s', datefmt='%m/%d/%Y %I:%M:%S %p',
                        level=logging.INFO)
    samples = 1000

    import matplotlib.pyplot as plt
    import collections

    test_data = [("[O]", "[C]=O", "[O]N=O"), ("[O-][N+](=O)[N+]([O-])=O", "O=C=O"), ("N(=O)[O]", "O=C=O"),
                 ("[H][H]", "N(=O)[O]")]
    f1 = plt.figure()
    if len(test_data) == 1:
        grid_columns = 1
    else:
        grid_columns = 2

    grid_rows = int(round(len(test_data) / (grid_columns * 1.0)))
    plt.gray()
    f1, ax = plt.subplots(nrows=grid_rows, sharex=True, sharey=True, ncols=grid_columns,
                          figsize=(10 * grid_columns, 7.5 * grid_rows), dpi=300, squeeze=False)
    plot_index = 0
    for reactant_smiles in test_data:
        axis = ax[plot_index / 2][plot_index % 2]
        axis.set_xlabel('KE')
        axis.set_ylabel('Product proportion')
        axis.set_title('Product proportion by KE for reactants {}'.format(reactant_smiles))
        plot_index += 1
        count = {}

        for _ke in range(0, 161, 5):
            reactant_mols = [Molecule(smiles, _ke) for smiles in reactant_smiles]
            # logging.info("{} at {}".format(reactant_smiles,_ke))

            for i in range(samples):

                reactions_object = EmergentReactions(chemistry=ChemistryFactory.new(), physics=PhysicsFactory.new(3))
                rxn = reactions_object.discover_reaction(Molecule.combine_molecules(reactant_mols))
                if rxn is not None:
                    product_mols = rxn.fire()
                    products = tuple(Chem.MolToSmiles(mol) for mol in product_mols)
                    if products not in count.keys():
                        count[products] = {}
                    try:
                        count[products][_ke] += 1
                    except:
                        count[products][_ke] = 1

            for product in count.keys():
                data = count[product]
                data = {k: v * 1.0 / samples for k, v in data.iteritems()}
                sorted_data = collections.OrderedDict(sorted(data.iteritems(), key=lambda t: t[0]))

                line, = axis.plot(sorted_data.keys(), sorted_data.values())

    file_name = "{}.eps".format(os.path.join(config.DataDir, 'emergent_reactions_plot'))
    logging.info("Saving figures to {}".format(file_name))
    f1.savefig(file_name, format='eps')
