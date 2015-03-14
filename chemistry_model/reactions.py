"""
Created on 8/01/2013

@author: thom
"""

from abc import abstractmethod


class Reactions(object):

    """The set of all possible reactions between reactants."""

    @abstractmethod
    def get_reaction_options(self, reactants):
        """Return the list of all reactions possible between the provided reactants.

        :param reactants: list of reactants for the discovered reaction
        :type reactants: list of Molecule (containing one or more molecules)
        :rtype: list of Reaction, or None if no reaction possible
        """
