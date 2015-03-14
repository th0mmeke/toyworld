"""
Created on 14 Aug 2013

@author: thom
"""

from chemistry_model.default_chemistry import DefaultChemistry


class ChemistryFactory(object):

    """Factory to create new Chemistry objects."""

    @classmethod
    def new(cls, parameters=None):
        return DefaultChemistry(parameters)
