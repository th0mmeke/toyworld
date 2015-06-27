"""
Created on 23 Aug 2013

@author: thom
"""
import importlib

from chemistry_model.chemistry_factory import ChemistryFactory


class ReactionVesselFactory(object):

    """Factory class to create Reaction Vessels."""

    @classmethod
    def new(cls, population, parameters):
        """
        :param population: the initial molecular population (will not be modified)
        :type population: MolecularPopulation
        :param parameters: additional simulation parameters
        :type parameters: Parameters
        """

        product_selection_strategy = parameters.get('ProductSelectionStrategy')
        chemistry = ChemistryFactory.new(parameters)

        vessel_module, vessel_class = parameters.get('Vessel').rsplit(".", 1)
        return getattr(importlib.import_module(vessel_module), vessel_class)(chemistry=chemistry,
                                                                             population=population,
                                                                             product_selection_strategy=product_selection_strategy,
                                                                             parameters=parameters)
