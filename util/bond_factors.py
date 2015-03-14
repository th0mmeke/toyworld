"""
Created on 11/03/2014

@author: thom
"""

import os
import config

import xml.etree.cElementTree as ElementTree
from runner import Runner

if __name__ == '__main__':
    dirname = os.path.join(config.DataDir, "experiments/bonds2")

    experiment_filename = os.path.join(dirname, "experiment_design.xml")
    parameters = ElementTree.parse(experiment_filename)
    runner = Runner(parameters, dirname)
    print("Experiment,Single,Double,Triple")
    for experiment in runner.get_experiments():
        formation_energies = experiment._parameters.find('BondFormationEnergies')
        values = {}
        for a in formation_energies:
            values[a.tag] = a.text
        print("{},{},{},{}".format(experiment._name, values['Single'], values['Double'], values["Triple"]))
