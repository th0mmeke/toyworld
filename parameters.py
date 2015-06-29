"""
Created on 22 Aug 2013

@author: thom
"""
import config
import os

import xml.etree.cElementTree as ElementTree


class Parameters(object):

    """Simple wrapper to provide a smarter method to get experiment parameters than the raw XML API.
    """

    # Cannot extend ElementTree.Element as it's a C method - "cannot create 'builtin_function_or_method' instances", so rather than wrapping, just provide a static method

    def __init__(self, parameters_etree=None, defaults=None):
        """:param parameters_etree: Set of parameters as an ElementTree
        :param defaults: Dictionary of default values to use if any requested value is missing - {key:value}

        Energy: Initial KE
        EnergyModel:
        Iterations: Maximum number of model iterations to perform
        Time: Maximum model time
        """

        if defaults is not None:
            self._defaults = defaults
        else:
            self._defaults = {
                'Energy': 300,
                'EnergyModel': 'energy_model.default_energy_model.DefaultEnergyModel',
                'Iterations': 10000,
                'Time': 1000,
                'RadiationRate': 0.0,
                'EnergyInput': 0.0,
                'Reactions': 'chemistry_model.emergent_reactions.EmergentReactions',
                'ProductSelectionStrategy': 'energy',
                'Molecule': 'kinetic_molecule.KineticMolecule',
                'Vessel': 'reactor_model.spatial_reaction_vessel.SpatialReactionVessel',
                'DipoleForceConstant': 0.01,

                'repeats': 1,
                'DeltaT': 1.0,
                'IterationBlocksize': 50,
                'recover': True,
                'seed': None,
                'StateRecordRate': 10.0,
            }

        if parameters_etree is None:
            self._parameters = ElementTree.fromstring("<dummy></dummy>")
        else:
            self._parameters = parameters_etree

    def get(self, parameter, default=True):
        """
        Attempt to get a parameter value from the XML element; if not found, return a predefined default value, or None if none exists
        """

        result = self._parameters.find(parameter)
        if result is None:
            if default and parameter in self._defaults.keys():
                return self._defaults[parameter]
            else:
                return None

        if len(result) == 0:
            return Parameters.convert_if_boolean(result.text.rstrip())
        else:
            return result  # return the raw element

    def get_attrib(self, attrib, default=True):
        """
        Attempt to get a parameter value from the XML element; return a predefined default value if not found
        """

        try:
            result = self._parameters.get(attrib).rstrip()
        except:
            if default and attrib in self._defaults.keys():
                result = self._defaults[attrib]
            else:
                return None

        return Parameters.convert_if_boolean(result)

    def get_filename(self, filename):
        path = self._parameters.find(filename).text.rstrip()
        if not os.path.isabs(path):
            path = os.path.join(config.DataDir, path)
        return os.path.normpath(path)

    @classmethod
    def convert_if_boolean(cls, text):
        try:
            l = text.lower()
        except:
            return text
        else:
            if l == 'true':
                return True
            elif l == 'false':
                return False
            return text

    def to_xml(self):
        return ElementTree.tostring(self._parameters)
