"""
@author: thom
"""

import logging
import time

from rdkit.Chem import AllChem as Chem
import networkx as nx


class Molecule(Chem.Mol):

    """A base representation of a Molecule with potential and kinetic energy.

    * Potential energy is the energy required to form all bonds in the Molecule (therefore always negative as bond formation releases energy).
    * Kinetic energy is, as usual, equal to 1/2 * mass * velocity ^ 2."""

    def __init__(self, source, internal_energy=0, kinetic_energy=None, canonize=True, components=None, **kwargs):
        """
        :param internal_energy: initial internal energy for this Molecule
        :type internal_energy: float
        :param kinetic_energy: initial kinetic energy for this Molecule
        :type kinetic_energy: float
        :param canonize: make implicit Hs explicit? Default True, but when copying Molecules we don't want to chance that these changes might be introduced
        :type canonize: bool
        :param components: for molecules that consist of multiple disjoint components, a mapping of atoms to component
        :type components: list of set of indexes of atoms within the molecule, where each set represents a different component
        """

        # just has to be unique over lifetime of simulation - id() only guarantees unique over lifetime of object, and time.clock() includes integer component which can overlap
        self.global_id = "{}.{}".format(id(self), time.clock())

        # Make all H explicit for later processing

        if not isinstance(source, Chem.Mol):
            source = Chem.MolFromSmiles(source)

        if canonize:
            source = Chem.AddHs(source)

        Chem.Mol.__init__(self, source.ToBinary())

        if components is None and Chem.MolToSmiles(self).find(".") == -1:  # if simple molecule with only one component, set components manually
            components = [set(range(source.GetNumAtoms()))]

        self._components = components
        self._mass = sum([atom.GetMass() for atom in self.GetAtoms()])
        self.set_internal_energy(internal_energy)
        if kinetic_energy is not None:
            self.set_kinetic_energy(kinetic_energy)

    def get_state(self):
        state = {'ke': self.get_kinetic_energy(), 'ie': self.get_internal_energy()}
        return state

    def get_mass(self):
        return self._mass

    def get_potential_energy(self, chemistry):
        """Return the energy required to form all of the molecule's bonds (therefore a negative quantity as bond formation releases energy)
        WARNING: assumes formation energy = energy of breaking (symmetrical)

        :rtype: float"""

        return sum([chemistry.get_bond_energy(bond.GetBeginAtom(), bond.GetEndAtom(), end_bond_type=int(bond.GetBondType())) for bond in self.GetBonds()])

    def get_internal_energy(self):
        return self._internal_energy

    def set_internal_energy(self, value):
        if value < 0:
            raise ValueError
        self._internal_energy = value

    def set_kinetic_energy(self, value):
        if value < 0:
            raise ValueError
        self._kinetic_energy = value

    def get_kinetic_energy(self):
        return self._kinetic_energy

    def split_molecule(self):
        """Allocate the initial energy proportional to the square of the mass of each resulting molecule.

        A rather simplified calculation as we can't easily work out the transfer of energy as the reaction is changing
        the interacting molecules

        :rtype: list of Molecule or subclass"""

        split_mols = [Molecule(smiles) for smiles in Chem.MolToSmiles(self).split(".")]
        total_mass_squared = sum([mol.get_mass() ** 2 for mol in split_mols])

        for mol in split_mols:
            mol.set_internal_energy(self._internal_energy * (mol.get_mass() ** 2 / total_mass_squared))
            mol.set_kinetic_energy(self.get_kinetic_energy() * (mol.get_mass() ** 2 / total_mass_squared))

        return split_mols

    def combine_molecules(self, mols):
        """Combine a number of molecules into one. Bookkeeping rather than Chemistry - does not connect molecules with bonds,
        just groups them in RDKit. The kinetic energy of the combined molecule is assumed to be preserved - that is, we assume
        a head-on collision.

        :rtype: Molecule or subclass"""

        combined_mol = mols[0]
        combined_IE = mols[0].get_internal_energy()
        combined_KE = mols[0].get_kinetic_energy()
        num_atoms = mols[0].GetNumAtoms()
        components = mols[0]._components
        for mol in mols[1:]:
            combined_mol = Chem.CombineMols(combined_mol, mol)
            combined_IE += mol.get_internal_energy()
            combined_KE += mol.get_kinetic_energy()
            components.append(set(range(num_atoms, num_atoms + mol.GetNumAtoms())))
            num_atoms = num_atoms + mol.GetNumAtoms()

        return Molecule(combined_mol, internal_energy=combined_IE, kinetic_energy=combined_KE, components=components)

    def same_component(self, idx1, idx2):
        # Components is a list of sets of atom indexes in molecule - same set = same component
        if self._components is None:  # Catch-all - components should have been set in init() and multiple components initialized in previous call to combine_molecules
            self._components = self._get_strongly_connected_components()

        for component in self._components:
            if idx1 in component and idx2 in component:
                return True
        return False

    def get_total_formal_charge(self):
        return sum([i.GetFormalCharge() for i in self.GetAtoms()])

    def _assign_formal_charge(self):
        """OpenEye Charge Model - http://www.eyesopen.com/docs/toolkits/current/html/OEChem_TK-python/valence.html#subsection-valence-openeye-charge
        The OpenEye formal charge model assigns formal charges to elements based upon their total valence.
        In OEChem, this functionality is invoked by the OEAssignFormalCharges function.
        If the formal charge on an atom is non-zero, it is left unchanged.

        Hydrogen
        If the valence isn't one, the formal charge is +1.
        Carbon
        If the valence is three, the formal charge is +1 if the atom has a polar neighbor, i.e. N, O or S, and formal charge -1 otherwise.
        Nitrogen
        If the valence is two, the formal charge is -1, and if the valence is four the formal charge is +1.
        Oxygen
        If the valence is one, the formal charge is -1, and if the valence is three the formal charge is +1."""

        for i in self.GetAtoms():

            valence = i.GetDegree()
            formal_charge = 0

            if i.GetAtomicNum() == 1:  # H
                if valence != 1:
                    formal_charge = 1

            if i.GetAtomicNum() == 6:  # C
                if valence == 3:
                    formal_charge = -1
                    for neighbor in i.GetNeighbors():
                        if neighbor.GetAtomicNum() == 7 or neighbor.GetAtomicNum() == 8:
                            formal_charge = 1

            if i.GetAtomicNum() == 7:  # N
                if valence == 2:
                    formal_charge = -1
                if valence == 3:
                    formal_charge = 1

            if i.GetAtomicNum() == 8:  # O
                if valence == 1:
                    formal_charge = -1
                elif valence == 3:
                    formal_charge = 1

            i.SetFormalCharge(formal_charge)

    def _get_strongly_connected_components(self):
        logging.info("Call to _get_strongly_connected_components")
        g = nx.Graph()
        for bond in self.GetBonds():
            g.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        connected_components = list(nx.connected_components(g))
        # Add single atoms as independent components
        for idx in range(self.GetNumAtoms()):
            if len(self.GetAtomWithIdx(idx).GetBonds()) == 0:
                connected_components.append([idx])

        return connected_components

    def __deepcopy__(self, memo):
        # don't mess with current structure - just leave it exactly as it is
        return Molecule(self, internal_energy=self.get_internal_energy(), kinetic_energy=self.get_kinetic_energy(), canonize=False)

    def __str__(self):
        return Chem.MolToSmiles(self)
