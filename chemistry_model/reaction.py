"""
Created on 8/01/2013

@author: thom
"""


class Reaction(object):

    """A reaction takes reactant molecules and transforms them into product molecules."""

    def __init__(self, product_molecules, energy_delta=0, is_reaction=True):
        """
        Each reaction is unique to particular implicit reactants, known at the time the Reaction is defined

        :param product_molecules: The resulting products of the reaction, one KineticMolecule per product
        :type product_molecules: list of KineticMolecule
        :param energy_delta: The difference in energy between products and reactants (-=exothermic)
        :type energy_delta: int
        :param reaction_energy: The energy available to power the reaction
        :type reaction_energy: int
        :param reaction: True if this represents a true reaction; False if just a "bounce" of the reactant molecules without chemical changes
        :type reaction: bool
        """

        self._product_molecules = product_molecules
        self._energy_delta = energy_delta
        self._is_reaction = is_reaction

    def fire(self):
        return self._product_molecules

    def get_energy_delta(self):
        return self._energy_delta

    def is_reaction(self):
        return self._is_reaction
