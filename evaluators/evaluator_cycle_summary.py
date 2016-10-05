"""
Created on 6/05/2013

@author: thom
"""

import logging

import networkx as nx

from atoms.kinetic_molecule import KineticMolecule
from evaluator import Evaluator
from reaction_network import ReactionNetwork


class EvaluatorCycleSummary(Evaluator):

    def get_result_titles(self):
        return ["Nodes", "Edges", "Cycles", "Components", "v-cycles", "v-other"]

    def evaluate(self, results_filename, selected_partitions=None, reaction_network=None, **kwargs):

        actual_cycles = ReactionNetwork.discover_actual_reaction_cycles(results_filename, selected_partitions, node_type='id', min_length=3, min_count=2)
        cycle_ids = set(self.flatten(actual_cycles))

        locations = {}
        velocities = {}
        # smiles = {}

        for block in Evaluator.incr_load_results(results_filename, selected_partitions):

            for reaction in block['reactions']:

                if 'reaction_site' in reaction.keys():
                    reaction_molecules = [reaction[y] for y in ['reactants', 'products']][0]
                    reaction_ids = set([x['id'] for x in reaction_molecules])

                    for mol_id in cycle_ids.intersection(reaction_ids):
                        try:
                            locations[mol_id].append([reaction['t'], reaction['reaction_site']])
                        except:
                            locations[mol_id] = [[reaction['t'], reaction['reaction_site']]]

                    for mol in reaction_molecules:
                        v = KineticMolecule(mol['smiles'], kinetic_energy=mol['ke']).get_speed()
                        mol_id = mol['id']
                        try:
                            velocities[mol_id].extend(v)
                        except:
                            velocities[mol_id] = [v]

        n = nx.DiGraph()
        for cycle in actual_cycles:
            for idx in range(len(cycle)):
                # print(cycle[idx],cycle[(idx + 1) % len(cycle)])
                n.add_edge(cycle[idx], cycle[(idx + 1) % len(cycle)])

        # nx.draw(n)
        # plt.savefig("{}-cycles.eps".format(output_filepath), format='eps')
        cycle_velocities = self.flatten([v for id, v in velocities.iteritems() if id in cycle_ids])
        non_cycle_velocities = self.flatten([v for id, v in velocities.iteritems() if id not in cycle_ids])
        try:
            cycle_average_v = sum(cycle_velocities) / len(cycle_velocities)
        except:
            cycle_average_v = 0
        try:
            other_average_v = sum(non_cycle_velocities) / len(non_cycle_velocities)
        except:
            other_average_v = 0

        cycles = nx.simple_cycles(n)
        components = nx.weakly_connected_components(n)
        logging.info("Average v for cycles = {}".format(cycle_average_v))
        logging.info("Average v for all other molecules = {}".format(other_average_v))
        logging.info("Nodes={}, edges={}".format(n.number_of_nodes(), n.number_of_edges()))
        logging.info("Cycles={}".format(list(cycles)))
        logging.info("Components={}".format(list(components)))

        return n.number_of_nodes(), n.number_of_edges(), list(cycles), list(components), cycle_average_v, other_average_v

    def flatten(self, l):
        return [y for x in l for y in x]
