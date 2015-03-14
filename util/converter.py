import config
import os

from rdkit.Chem import AllChem as Chem

if __name__ == "__main__":

    """Population of molecules (in SMILES format: http://www.daylight.com/dayhtml/doc/theory/theory.smiles.html)"""
    """perl -e 'print sort { length $a<=>length $b || $a cmp $b } <>' mols.smi > sorted-mols.smi"""

    smiles_database_filename = os.path.join(config.DataDir, 'desc_sorted.gdb13.rand1M.smi')
    population_filename = os.path.join(config.DataDir, 'converted_population.xml')

    in_file = open(smiles_database_filename, "rU")
    out_file = open(population_filename, "w")

    out_file.write('<?xml version="1.0" ?>\n')
    out_file.write("<population>\n")

    count = 0
    for smiles in in_file:
        item = Chem.MolToSmiles(Chem.MolFromSmiles(smiles.rstrip()))
        entry = '<setvalue item="{}" quantity="{}" />\n'.format(item, "1")
        out_file.write(entry)
        count += 1
        print ("{:<10}{}".format(count, entry.rstrip()))
        if count >= 10:
            break

    out_file.write("</population>\n")
    in_file.close()
    out_file.close()
