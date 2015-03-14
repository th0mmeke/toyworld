"""
Created on 21 Aug 2013

@author: thom
"""

import logging


class DefaultEnergyModel(object):

    """
    A simple energy model.
    """

    def __init__(self, parameters):
        self._relative_rate = float(parameters.get('RadiationRate'))
        self._absolute_rate = float(parameters.get('EnergyInput'))
        logging.info("Default Energy Model with radiation rate of {}/second and absolute rate of {}/second".format(self._relative_rate, self._absolute_rate))

    def get_absolute_input(self, initial_energy, t):
        """Relative rate is relative to initial value"""
        return self._absolute_rate * t

    def get_relative_output(self, initial_energy, t):
        return initial_energy * self._relative_rate * t
