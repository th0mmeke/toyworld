"""
Created on 14 Aug 2013

@author: thom
"""

from semi_realistic_chemistry import SemiRealisticChemistry


class ChemistryFactory(object):

    """Factory to create new Chemistry objects."""

    @classmethod
    def new(cls, parameters=None):
        return SemiRealisticChemistry(parameters)
