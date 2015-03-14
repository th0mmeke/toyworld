"""
Created on 8/01/2013

@author: thom
"""


import logging
import copy

from rdkit.Chem import AllChem as Chem

from chemistry_model.reactions import Reactions
from chemistry_model.reaction import Reaction
from ULPS import Float_t
from molecule import Molecule  # only used in assertion to check conservation of potential energy


class EmergentReactions(Reactions):

    """Emergent reactions are those that are created on-the-fly based on properties of the interacting molecules.

    Covalent bonds mean shared electrons - single bond means one from each atom, two electrons shared, each atom appears to gain an electron
    A higher bond energy (or a higher bond order or shorter bond length) means that a bond is less likely to break apart. In other words,
    it is more stable than a molecule with a lower bond energy. With Lewis Structures then, the structure with the higher bond energy is  more likely to occur.
    Breaking a bond takes energy; making a bond releases energy.
    Bond order (single, double, triple) inversely related to bond length
    Double bond takes/releases less than 2x single bond energy
    Enthalpy = sum of bond energy changes. Positive = endothermic (reaction requires energy)
    The Octet Rule requires all atoms in a molecule to have 8 valence electrons--either by sharing, losing or gaining electrons--to become stable.
    For Covalent bonds, atoms tend to share their electrons with each other to satisfy the Octet Rule. It requires 8 electrons because that is the
    amount of electrons needed to fill a s- and p- orbital (electron configuration); also known as a noble gas configuration. Each atom wants to become as
    stable as the noble gases that have their outer valence shell filled because noble gases have a test of 0. Although it is important to remember
    the "magic number", 8, note that there are many Octet rule exceptions.

    - http://chemwiki.ucdavis.edu/Theoretical_Chemistry/Chemical_Bonding/General_Principles/Bond_Energies

    http://en.wikipedia.org/wiki/Formal_charge for formal test formula:
        FC = V - N - B/2, where:
        V is the number of valence electrons of the atom in isolation (atom in ground state)
        N is the number of non-bonding valence electrons on this atom in the molecule; and
        B is the total number of electrons shared in covalent bonds with other atoms in the molecule."""

    def __init__(self, chemistry):
        self._chemistry = chemistry

    def get_reaction_options(self, reactants):
        """Discover all possible options for reactions between the separate components of this Molecule, based on the
        Chemistry provided. The options are based on possible bond changes; each option is returned in a Reaction object as a list of product Molecules along
        with the associated change in potential energy.

        That is, self.get_potential_energy() + rxn.get_energy_delta() = sum(product molecule potential energy)

        If any Hs have been added earlier this will preserve them.

        Note: this is surprisingly difficult in RDKit. RDKit's transformation method - EditableMol - doesn't play well in subclasses as it GetMol() method
        returns an object of class Mol rather than the original type. So using EditableMol results in a different type from the original.
        This means that we either 1) avoid using EditableMol, 2) restore the type, or 3) base everything on Mol rather than subclasses
        3) is impractical - KineticMolecule and Molecule are natural extensions of Mol
        1) is impossible - we need to break bonds, which requires EditableMol (SetBondType() can't remove a bond as Chem.BondType doesn't have a 'nothing' option)
        2) is difficult. We need transformed copies of our original reactants (which are subclasses of Mol), but Mol cannot be copied (no copy method, doesn't
        play with copy.copy or copy.deepcopy as doesn't support Pickle)

        :param reactants: List of reactants
        :type reactants: List of Molecule
        :rtype: list of Reaction"""

        def _adjust_formal_charge(begin_atom, end_atom, bond_order):
            """Bring the functional charges on the atoms in a bond to as close to zero as possible..."""

            begin_atom_fc = begin_atom.GetFormalCharge()
            end_atom_fc = end_atom.GetFormalCharge()
            # logging.debug("Adjusting test on {} ({}) and {} ({})".format(begin_atom.GetSymbol(), begin_atom_fc, end_atom.GetSymbol(), end_atom_fc))

            if (cmp(begin_atom_fc, 0) != cmp(end_atom_fc, 0)):  # must be opposite sign for this to work
                adjustment = min(abs(begin_atom_fc), abs(end_atom_fc), bond_order)
                begin_atom.SetFormalCharge(begin_atom_fc - adjustment * cmp(begin_atom_fc, 0))
                end_atom.SetFormalCharge(end_atom_fc - adjustment * cmp(end_atom_fc, 0))

        def _change_bond(mol, begin_atom_idx, end_atom_idx, new_bond_order):
            """Change type of bond - throw exception if cannot legitimately be changed.

            :rtype: same type as caller"""

            # logging.debug("Attempting to change bond between {} and {} to {}".format(begin_atom_idx, end_atom_idx, new_bond_order))

            if new_bond_order > 3:
                raise ValueError  # to meet RDKit restriction from organic chemistry that maximum likely bond is triple-bond

            bond = mol.GetBondBetweenAtoms(begin_atom_idx, end_atom_idx)
            e_mol = Chem.EditableMol(mol)
            if bond is not None:  # remove any existing bond
                e_mol.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
            if new_bond_order != 0:  # add in any new bond
                e_mol.AddBond(begin_atom_idx, end_atom_idx, Chem.BondType.values[new_bond_order])
            new_mol = e_mol.GetMol()

            # Adjust formal charges if using bond based on charge_order...goal is to make fc of both zero
            _adjust_formal_charge(new_mol.GetAtomWithIdx(begin_atom_idx), new_mol.GetAtomWithIdx(end_atom_idx), new_bond_order)

            try:
                new_mol = type(mol)(new_mol, kinetic_energy=mol.get_kinetic_energy())
            except:
                new_mol = type(mol)(new_mol)

            Chem.SanitizeMol(new_mol, catchErrors=False)  # possibly not required
            assert mol.GetNumAtoms(onlyExplicit=False) == new_mol.GetNumAtoms(onlyExplicit=False)  # reaction must preserve matter
            return new_mol  # return same type as caller, self defined in parent method

        reactant_mol = reactants[0].combine_molecules(reactants)
        reaction_options = []

        # Reaction options:
        # 1. breaking of existing bonds
        # 2. addition of bonds
        # 3. switch one bond type for another
        # In each of these cases we update the formal charges to reflect the changes, and test the resulting product for sanity

        # Bond breaking and substitution options...
        for bond in reactant_mol.GetBonds():
            begin_atom_idx = bond.GetBeginAtomIdx()
            end_atom_idx = bond.GetEndAtomIdx()
            old_bond_order = int(bond.GetBondType())

            for new_bond_order in range(0, old_bond_order):  # all bond options between none (break bond) and one-less-than-current
                reactants_copy = copy.deepcopy(reactant_mol)  # operate on copy of self, so options are not cumulative
                try:
                    new_mol = _change_bond(reactants_copy, begin_atom_idx, end_atom_idx, new_bond_order)

                except:
                    pass  # just ignore if this option isn't possible
                else:  # no exception, so managed to execute _change_bond()...
                    bond_energy = self._chemistry.get_bond_energy(new_mol.GetAtomWithIdx(begin_atom_idx), new_mol.GetAtomWithIdx(
                        end_atom_idx), start_bond_type=old_bond_order, end_bond_type=new_bond_order)  # change of energy from bond of GetBondType() to bond_order
                    assert Float_t.almost_equal(reactant_mol.get_potential_energy(
                        self._chemistry) + bond_energy, new_mol.get_potential_energy(self._chemistry))
                    logging.debug("Reaction option: {} to {} - change bond {}-{} from {} to {} (e={})".format(Chem.MolToSmiles(reactant_mol),
                                                                                                              Chem.MolToSmiles(new_mol), begin_atom_idx, end_atom_idx, old_bond_order, new_bond_order, bond_energy))
                    reaction_options.append(Reaction(new_mol.split_molecule(), bond_energy))

        # Bond addition options...
        bond_potential = []
        for atom in reactant_mol.GetAtoms():
            bond_potential.append(self._chemistry.get_bond_potential(atom))
            #logging.debug('Bond potential for {}={}'.format(atom.GetSymbol(), self._chemistry.get_bond_potential(atom)))

        # Check for possible bonds between all atoms with the potential for one or more additional bonds...
        for begin_atom_idx in range(len(bond_potential)):
            for end_atom_idx in range(begin_atom_idx + 1, len(bond_potential)):
                #logging.debug('{} and {} are in the same component = {}'.format(begin_atom_idx, end_atom_idx, reactant_mol.same_component(begin_atom_idx, end_atom_idx)))
                if not reactant_mol.same_component(begin_atom_idx, end_atom_idx):

                    # Add options for bonds of various order between begin_atom and end_atom...
                    begin_atom_valence = bond_potential[begin_atom_idx]
                    end_atom_valence = bond_potential[end_atom_idx]
                    max_bond_order = min(begin_atom_valence, end_atom_valence)
                    # to meet RDKit restriction from organic self._chemistry that maximum likely bond is triple-bond
                    max_bond_order = min(max_bond_order, 3)

                    # logging.debug("Checking potential bond of order {} between {} and {}".format(max_bond_order, reactant_mol.GetAtomWithIdx(begin_atom_idx).GetSymbol(), reactant_mol.GetAtomWithIdx(end_atom_idx).GetSymbol()))

                    if max_bond_order > 0:

                        for bond_order in range(1, max_bond_order + 1):  # all bond options up to and including max_bond_order
                            reactants_copy = copy.deepcopy(reactant_mol)  # operate on fresh copy so options don't accumulate
                            try:
                                new_mol = _change_bond(reactants_copy, begin_atom_idx, end_atom_idx, bond_order)
                            except:
                                pass  # just ignore invalid options
                            else:
                                bond_energy = self._chemistry.get_bond_energy(new_mol.GetAtomWithIdx(begin_atom_idx), new_mol.GetAtomWithIdx(
                                    end_atom_idx), end_bond_type=bond_order)  # bond creation of order bond_order
                                logging.debug("Reaction option: {} to {} - addition of {} bond from {} to {} (e={})".format(
                                    Chem.MolToSmiles(reactant_mol), Chem.MolToSmiles(new_mol), bond_order, begin_atom_idx, end_atom_idx, bond_energy))
                                assert Float_t.almost_equal(reactant_mol.get_potential_energy(
                                    self._chemistry) + bond_energy, Molecule(new_mol).get_potential_energy(self._chemistry))
                                reaction_options.append(Reaction(new_mol.split_molecule(), bond_energy))

        logging.debug("{} reaction options found".format(len(reaction_options)))
        return reaction_options
