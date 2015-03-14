"""
Created on 1/03/2013

@author: thom
"""

from experiment import Experiment
from parameters import Parameters

from rdkit.rdBase import DisableLog, EnableLog

import logging
import collections
import importlib


class Runner(object):

    """
    Run a series of experiments.
    """

    def __init__(self, parameters, dirname):
        self._parameters = parameters
        self._dirname = dirname

        self._factor_definitions = {}
        if self._parameters.find("Factors") is not None:
            # set up factors list
            factors_xml = self._parameters.find('Factors').findall('Factor')

            for factor_xml in factors_xml:
                factor = {}
                for f in factor_xml:
                    if f.tag == 'Title':
                        factor[f.tag] = f.text
                    else:
                        factor[f.tag] = list(f)
                self._factor_definitions[factor_xml.get('key')] = factor
        logging.info("Experiment factors: {}".format(self._factor_definitions))

    def get_factor_definitions(self):
        return self._factor_definitions

    def get_experiments(self):
        experiments = self._parameters.findall("Experiment")  # might have multiple experiment blocks
        return collections.deque([Experiment(self._dirname, xml, list(self._parameters.find("Common")), self._factor_definitions) for xml in experiments])

    def get_evaluations(self):
        evaluations = []

        for xml in self._parameters.find("Evaluation"):  # one evaluation block, but with multiple methods
            evaluator_module, evaluator_class = xml.text.rsplit(".", 1)
            evaluator_module = evaluator_module.rstrip()
            evaluator_class = evaluator_class.rstrip()
            partition_input = Parameters.convert_if_boolean(xml.get('partition'))  # set flag in evaluator if requested to partition the input into sample blocks
            evaluations.append(getattr(importlib.import_module(evaluator_module), evaluator_class)(partition=partition_input))
        return evaluations

    def run(self):
        """
        Expand DoE XML design into canonical experiment format
        DoE format has <Factors> block, defining the varying factors, and a <Common> block with static elements of the design
        Each experiment has a set of <factor key=<key> value=<value>/> elements which need to be expanded
        We also need to copy in the common elements
        """
        experiments = self.get_experiments()

        DisableLog('rdApp.error')  # disable rdKit error messages, in particular "Explicit valence..."
        logging.info("Running {} experiment{}".format(len(experiments), "s" if len(experiments) > 1 else ""))

        number_of_experiments = len(experiments)
        for i in range(number_of_experiments):
            experiment = experiments.popleft()
            logging.info("Starting experiment #{}".format(i + 1))
            experiment.run()
            del experiment

        logging.info("Finished running experiment{}".format("s" if len(experiments) > 1 else ""))
        EnableLog('rdApp.error')
