"""
Created on 14 Aug 2013

@author: thom
"""

from rdkit.Chem import AllChem as Chem


class SemiRealisticChemistry(object):

    """A simple Chemistry based on real-world reactions."""

    def __init__(self, parameters=None):
        """:param parameters: Parameters object"""

        # Single : 77.7 = 777.1/10 = 104.2 + 83 + 38.4 + 35 + 99 + 93 + 111 + 73 + 85.5 + 55
        # Double : 148.2 = 889/6 = 185 + 146 + 149 + 119 + 147 + 143
        # Triple : 224.3 = 897/4 = 258 + 200 + 226 + 213
        default_data = {
            'H1H': 104.2,
            'C1C': 83,
            'N1N': 38.4,
            'O1O': 35,
            'H1C': 99, 'C1H': 99,
            'H1N': 93, 'N1H': 93,
            'H1O': 111, 'O1H': 111,
            'C1N': 73, 'N1C': 73,
            'C1O': 85.5, 'O1C': 85.5,
            'N1O': 55, 'O1N': 55,
            'C2O': 185, 'O2C': 185,  # rough average of range
            'C2C': 146,
            'N2N': 149,
            'O2O': 119,
            'C2N': 147, 'N2C': 147,
            'N2O': 143, 'O2N': 143,
            'C3O': 258, 'O3C': 258,
            'C3C': 200,
            'N3N': 226,
            'C3N': 213, 'N3C': 213,
            'C4C': 200  # theoretically possible from valences, but in nature forms a C2C bond instead
        }
        count = {}
        default_bond_energies = {}
        for bond, energy in default_data.iteritems():
            key = int(bond[1])
            try:
                count[key] += 1
                default_bond_energies[key] += energy
            except:
                count[key] = 1
                default_bond_energies[key] = energy
        for i in (1, 2, 3):
            default_bond_energies[i] = default_bond_energies[i] / count[i]

        self._atoms = ['C', 'N', 'O', 'H']
        self._bond_formation_energies = {}
        self._bond_break_energies = {}
        formation_energies = None
        break_energies = None

        if parameters is not None:  # Parameters object
            atoms = parameters.get('atoms')
            if atoms is not None:
                self._atoms = []
                for atom in atoms.findall('Atom'):
                    self._atoms.append(atom.text)

            formation_energies = parameters.get('BondFormationEnergies')
            break_energies = parameters.get('BondBreakEnergies')

        for atom_1 in self._atoms:
            for atom_2 in self._atoms:
                for bond_type, xml_key in {1: 'Single', 2: 'Double', 3: 'Triple'}.iteritems():
                    key = "{}{}{}".format(atom_1, bond_type, atom_2)
                    if formation_energies is None:
                        if key in default_data.keys():
                            self._bond_formation_energies[key] = default_data[key]
                        else:
                            self._bond_formation_energies[key] = default_bond_energies[bond_type]
                    else:
                        self._bond_formation_energies[key] = float(formation_energies.find(xml_key).text)
                    if break_energies is None:
                        self._bond_break_energies[key] = self._bond_formation_energies[key]
                    else:
                        self._bond_break_energies[key] = float(break_energies.find(xml_key).text)

    def get_bond_potential(self, atom):
        """Requires Explicit Hs!

        Simple method based on standard Lewis dot-structures e.g., http://library.thinkquest.org/C006669/data/Chem/bonding/lewis.html

        Bond calculation:
        FC = V - N - B (where single bond = 1, double = 2, etc) to make N = 8
        """
        if atom.GetAtomicNum() == 1:
            if len(atom.GetBonds()) == 0:  # if not already bound...
                return 1
            else:
                return 0
        else:
            bonded_electrons = 0
            for bond in atom.GetBonds():
                bonded_electrons += bond.GetBondType()  # relies on Chem.BondType mapping to int...

            valence_electrons = Chem.GetPeriodicTable().GetNOuterElecs(atom.GetAtomicNum())
            # logging.info("Bond potential: 8 - ({} + {} + {})".format(valence_electrons, bonded_electrons, atom.GetFormalCharge()))
            return 8 - (valence_electrons + bonded_electrons + atom.GetFormalCharge())

    def get_bond_energy(self, atom_1, atom_2, end_bond_type=0, start_bond_type=0):
        """Returns the energy REQUIRED to make the bond change from start_bond_type (or existing type if not provided) to end_bond_type.
        Creation of a bond requires -e; breaking the bond +e
        Energies taken from http://www.cem.msu.edu/~reusch/OrgPage/bndenrgy.htm - Average Bond Dissociation Enthalpies in kcal per mole

        :param atom_1: One end of the bond
        :type atom_1: Chem.Atom
        :param atom_2: Other end of the bond
        :type atom_2: Chem.Atom
        :param bond_type: Type of the bond, corresponding to index into Chem.BondType.values
        :type bond_type: int
        :rtype: int
        """

        # Energy to release current bond state
        if start_bond_type == 0:
            start_energy = 0
        else:
            start_energy = self._bond_break_energies[atom_1.GetSymbol() + str(min(3, start_bond_type)) + atom_2.GetSymbol()]

        # Energy to create desired bond state
        if end_bond_type == 0:
            end_energy = 0
        else:
            end_energy = self._bond_formation_energies[atom_1.GetSymbol() + str(min(3, end_bond_type)) + atom_2.GetSymbol()]

        return start_energy - end_energy
