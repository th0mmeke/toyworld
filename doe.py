"""
Design a factorial experiment.

-e, --experiment  Name of experiment
-i, --iterations  Number of iterations for each repeat or run
-r, --repeats  Number of repeats for experiment
-s, --seed  Starting random number seed
-p, --population  Filename of population definition
-t, --preview  Only display design; do not write or overwrite

"""

import itertools
import os
import config
import argparse
import random
import xml.dom.minidom
import re


class FactorialDesign:

    @classmethod
    def design(cls, directory, experiment_name, seed, iterations, population, factors, repeats=1, recover=True):

        if recover:
            recover_attrib = 'True'
        else:
            recover_attrib = 'False'

        experiment_design = "<ExperimentDesign>".format(experiment_name)

        experiment_design += "<Common>"
        experiment_design += "<Iterations>{}</Iterations>".format(iterations)
        experiment_design += "<PopulationFilename>{}</PopulationFilename>".format(population)
        experiment_design += "<Reactions>chemistry_model.emergent_reactions.EmergentReactions</Reactions>"
        experiment_design += "</Common>"

        experiment_design += "<Factors>"
        design_factors = []
        for factor in factors:
            experiment_design += "<Factor key='{}'>".format(factor['key'])
            f = []
            for key, value in factor.iteritems():
                if key != 'key':
                    experiment_design += "<{0}>{1}</{0}>".format(key.capitalize(), value)
                    if key != 'Title':
                        f.append(value)
            experiment_design += "</Factor>"
            design_factors.append(f)

        experiment_design += "</Factors>"

        experiment_number = 1

        design = list(itertools.product(*design_factors))
        random.shuffle(design)

        for experiment_factors in design:
            experiment_design += "<Experiment name='{}-{}' repeats='{}' recover='{}' seed='{}'>".format(experiment_name, experiment_number, repeats, recover_attrib, seed)
            for factor in experiment_factors:
                experiment_design += factor
            experiment_design += "</Experiment>"

            experiment_number += 1

        experiment_design += "<Evaluation>"
        experiment_design += "<Method>evaluators.evaluator_summary.EvaluatorSummary</Method>"
        experiment_design += "<Method>evaluators.evaluator_cycles.EvaluatorActualCycles</Method>"
        experiment_design += "<Method>evaluators.evaluator_cycles.EvaluatorPotentialCycles</Method>"
        experiment_design += "<Method>evaluators.evaluator_iterations.EvaluatorIterations</Method>"
        experiment_design += "</Evaluation>"
        experiment_design += "</ExperimentDesign>"

        return experiment_design

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--experiment', default='experiment', help="Name of experiment")
    parser.add_argument('-i', '--iterations', default=10000)
    parser.add_argument('-r', '--repeats', default=1)
    parser.add_argument('-s', '--seed', default=None)
    parser.add_argument('-b', '--bonds', action='store_true')
    parser.add_argument('-p', '--population', default="mixed_population.xml")
    parser.add_argument('-t', '--preview', action='store_true')

    args = parser.parse_args()
    if args.seed == None:
        args.seed = random.randint(0, 1203123123123)

    directory = os.path.join(config.DataDir, args.experiment)

    if not args.bonds:
        factors = []
        factors.append({'key': 'dimensionality', 'Title': 'Dimensionality', 'Low': "<Molecule>molecule.Molecule</Molecule><Vessel>reactor_model.aspatial_reaction_vessel.AspatialReactionVessel</Vessel>",
                        'High': "<Molecule>kinetic_molecule.KineticMolecule</Molecule><Vessel>reactor_model.spatial_reaction_vessel.SpatialReactionVessel</Vessel>"})
        factors.append({'key': 'energy', 'Title': 'Energy', 'Low': "<Energy>100</Energy>", 'High': "<Energy>300</Energy>"})
        factors.append({'key': 'bonds', 'Title': 'Bond Energies', 'Low': "", 'High': "<BondFormationEnergies><Single>50</Single><Double>100</Double><Triple>200</Triple></BondFormationEnergies>"})
    else:
        factor = {'key': 'bonds', 'Title': 'Bond Energies'}
        for i in range(0, 32):
            factor['set' + str(i)] = "<BondFormationEnergies><Single>{}</Single><Double>{}</Double><Triple>{}</Triple></BondFormationEnergies>".format(random.randint(0, 300),
                                                                                                                                                       random.randint(0, 300), random.randint(0, 300))
        factors = [factor]

    experiment_design = FactorialDesign.design(directory, args.experiment, args.seed, args.iterations, args.population, factors, repeats=args.repeats, recover=True)

    if args.preview:
        # http://stackoverflow.com/questions/749796/pretty-printing-xml-in-python
        xml = xml.dom.minidom.parseString(experiment_design)
        uglyXml = xml.toprettyxml(indent='  ')

        text_re = re.compile('>\n\s+([^<>\s].*?)\n\s+</', re.DOTALL)
        prettyXml = text_re.sub('>\g<1></', uglyXml)
        print(prettyXml)

    else:

        if not os.path.lexists(directory):
            os.makedirs(directory)

        experiment_filename = 'experiment_design.xml'
        print("Writing experiment file {} to directory {}".format(experiment_filename, directory))
        with open(os.path.join(directory, experiment_filename), "w") as f:
            f.write(experiment_design)

        print("Now run 'python toyworld/main.py -e {1}/{0}'".format(experiment_filename, directory))
        print("And then 'python toyworld/evaluate.py -e {1}/{0}'".format(experiment_filename, directory))
