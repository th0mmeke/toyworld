"""
Created on 23 Aug 2013

@author: thom
"""
import importlib

from chemistry_model.chemistry_factory import ChemistryFactory


class ReactionVesselFactory(object):

    """Factory class to create Reaction Vessels."""

    @classmethod
    def new(cls, population, results_filename, states_filename, parameters):
        """
        :param parameters: additional simulation parameters
        :type parameters: Parameters
        :param population: the initial molecular population (will not be modified)
        :type population: MolecularPopulation
        :param initial_average_ke: starting energy level
        :type initial_average_ke: float
        """

        product_selection_strategy = parameters.get('ProductSelectionStrategy')
        chemistry = ChemistryFactory.new(parameters)

        vessel_module, vessel_class = parameters.get('Vessel').rsplit(".", 1)
        return getattr(importlib.import_module(vessel_module), vessel_class)(chemistry=chemistry,
                                                                             population=population,
                                                                             product_selection_strategy=product_selection_strategy,
                                                                             results_filename=results_filename,
                                                                             states_filename=states_filename,
                                                                             parameters=parameters)
