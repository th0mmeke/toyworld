"""
@author: thom
"""

import math

from rdkit.Chem import AllChem as Chem
import pymunk as pm

from kinetics_2D import Kinetics2D
from molecule import Molecule


class KineticMolecule(Molecule):

    """A base representation of a Molecule with potential, kinetic and internal energy.

    * Potential energy is the energy required to form all bonds in the Molecule (therefore always negative as bond formation releases energy).
    * Kinetic energy is, as usual, equal to 1/2 * mass * velocity ^ 2.
    * Internal energy is a molecular store topped up by radiation, depleted by entropy, and convertible between kinetic and potential energies during reactions.

    Reaction model - on a collision between two molecules:

    * Energy available for a reaction = internal energy + kinetic energy - energy of centre of mass
    * A number of reaction options are then generated from the colliding molecules
    * Post-reaction energy = Energy available for a reaction - energy of reaction (change of potential energy)
    * Post-reaction energy is then split between kinetic energy and internal energy - the assumption is that internal energy is excess energy
      stored in atoms, so can be completely exhausted in a reaction

    Radiation/Entropy model:

    * incoming radiation is assigned to internal energy, available later for reactions or kinetic energy on collisions
    * entropy depletes internal energy by a proportion (so never results in negative internal energy.)
    """

    def __init__(self, source, internal_energy=0, kinetic_energy=0, canonize=True, components=None, base_molecule_radius=1):
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

        self._base_molecule_radius = base_molecule_radius
        super(KineticMolecule, self).__init__(source, internal_energy=internal_energy, canonize=canonize, components=components)  # set mass
        inertia = pm.moment_for_circle(self.get_mass(), 0, self.get_size(), (0, 0))
        self.body = pm.Body(self.get_mass(), inertia)
        radius = self.get_size()
        shape = pm.Circle(self.body, radius, (0, 0))
        shape.elasticity = 0.999  # required for standard collision handler to do a 'perfect' bounce
        shape.friction = 0.0
        self._shapes = [shape]

        self.set_kinetic_energy(kinetic_energy)

    def get_state(self):
        state = super(KineticMolecule, self).get_state()
        state.update({'speed': self.get_speed(), 'size': self.get_size(), 'ke': self.get_kinetic_energy()})
        return state

    def set_position(self, x, y):
        self.body.position = x, y

    def get_position(self):
        return self.body.position

    def set_kinetic_energy(self, ke):
        new_speed = math.sqrt(2.0 * ke / self.get_mass())  # v = sqrt(2ke/m)
        try:
            self.body.velocity = self.body.velocity * new_speed / self.get_speed()
        except:
            self.body.velocity = Kinetics2D.radial_to_xyz(r=new_speed)  # random direction

    def get_kinetic_energy(self):
        # Three different ways that the kinetic energy of this molecule can be changed:
        # 1. Explicitly through set_kinetic_energy()
        # 2. By modifying the velocity of the molecule through set_velocity()
        # 3. Implicitly by pymunk integrating force impulses on the molecule
        # For these reasons, must be careful velocity has not changed if optimizing this calculation by caching a value
        speed = self.get_speed()
        return 0.5 * self.get_mass() * speed * speed

    def get_speed(self):
        return self.body.velocity.get_length()
        # return Kinetics2D.get_speed(*self.body.velocity)

    def get_size(self):
        return self._base_molecule_radius

    def get_velocity(self):
        """:return: vec2d"""
        return self.body.velocity

    def set_velocity(self, *args):
        """:param args: vec2d"""
        self.body.velocity = args

    def split_molecule(self):
        """Allocate the initial energy proportional to the square of the mass of each resulting molecule.

        A rather simplified calculation as we can't easily work out the transfer of energy as the reaction is changing
        the interacting molecules

        :rtype: list of Molecule or subclass"""

        split_mols = [type(self)(smiles) for smiles in Chem.MolToSmiles(self).split(".")]
        total_mass_squared = sum([mol.get_mass() ** 2 for mol in split_mols])

        for mol in split_mols:
            mol.set_internal_energy(self._internal_energy * (mol.get_mass() ** 2 / total_mass_squared))
            mol.set_velocity(*self.get_velocity())

        return split_mols

    def combine_molecules(self, mols):
        """Combine a number of molecules into one. Bookkeeping rather than Chemistry - does not connect molecules with bonds,
        just groups them in RDKit. Preserves momentum and KE.

        :rtype: Molecule or subclass"""

        combined_mol = mols[0]
        combined_IE = mols[0].get_internal_energy()
        num_atoms = mols[0].GetNumAtoms()
        components = mols[0]._components
        for mol in mols[1:]:
            combined_mol = Chem.CombineMols(combined_mol, mol)
            combined_IE += mol.get_internal_energy()
            components.append(set(range(num_atoms, num_atoms + mol.GetNumAtoms())))
            num_atoms = num_atoms + mol.GetNumAtoms()

        combined_mol = type(self)(combined_mol, internal_energy=combined_IE, components=components)  # same type as mols

        # Adjust velocity of combined molecule to preserve momentum and KE
        combined_momentum = []
        for dim in range(2):
            combined_momentum.append(sum(mol.get_mass() * mol.get_velocity()[dim] for mol in mols))

        combined_velocity = tuple(i / combined_mol.get_mass() for i in combined_momentum)
        combined_mol.set_velocity(*combined_velocity)
        return combined_mol

    def __deepcopy__(self, memo):
        # create copy, but don't set kinetic_energy as we need to set velocity explicitly anyway
        mol = KineticMolecule(self, internal_energy=self._internal_energy, canonize=False)  # don't mess with current structure - just leave it exactly as it is
        mol.set_velocity(self.body.velocity)
        return mol

    def __str__(self):
        return Chem.MolToSmiles(self)
