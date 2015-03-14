"""
Created on 27/04/2013

@author: thom
"""
import unittest
import logging
import cPickle
import os

import networkx as nx

import config
from evaluators.reaction_network import ReactionNetwork


class TestEvaluatorNetwork(unittest.TestCase):

    filename = os.path.join(config.DataDir, 'test/test_evaluator_network.data')

    def writeToResultsFile(self, reactions):
        with open(TestEvaluatorNetwork.filename, "wb") as f:
            cPickle.dump({'block': {'start_block': 0, 'end_block': len(reactions) - 1, 'reactions': reactions}}, f, 2)
            f.flush()  # force python file write
            os.fsync(f.fileno())  # force OS file write

    def testStrongComponents(self):
        self.writeToResultsFile([
            {'iteration': 0, 'reactants': [{'id': 1, 'smiles': 'a'}, {'id': 2, 'smiles': 'b'}], 'products': [{'id': 3, 'smiles': 'c'}, {'id': 4, 'smiles': 'd'}]}
        ])
        reaction_network = ReactionNetwork.build_molecule_reaction_network(TestEvaluatorNetwork.filename, only_reactions=True)
        self.assertEqual(4, nx.number_strongly_connected_components(reaction_network))  # each node
        self.assertEqual(1, nx.number_weakly_connected_components(reaction_network))  # single component

        logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.INFO)

        self.writeToResultsFile([
            {'iteration': 0, 'reactants': [{'smiles': 'A'}], 'products': [{'smiles': 'B'}]},
            {'iteration': 1, 'reactants': [{'smiles': 'B'}], 'products': [{'smiles': 'A'}]},
        ])

        self.writeToResultsFile([
            {'iteration': 0, 'reactants': [{'smiles': 'A'}], 'products': [{'smiles': 'B'}, {'smiles': 'B'}]},
            {'iteration': 1, 'reactants': [{'smiles': 'B'}], 'products': [{'smiles': 'C'}]},
            {'iteration': 2, 'reactants': [{'smiles': 'C'}], 'products': [{'smiles': 'A'}]},
        ])  # simple cycle A->B->C->A
        reaction_network = ReactionNetwork.build_smiles_reaction_network(TestEvaluatorNetwork.filename, only_reactions=True)
        self.assertEqual(1, nx.number_strongly_connected_components(reaction_network))  # single cycle = single component
        self.assertEqual(1, nx.number_weakly_connected_components(reaction_network))  # single component

        self.writeToResultsFile([
            {'iteration': 0, 'reactants': [{'smiles': 'A'}], 'products': [{'smiles': 'B'}]},
            {'iteration': 1, 'reactants': [{'smiles': 'B'}], 'products': [{'smiles': 'A'}]},
            {'iteration': 2, 'reactants': [{'smiles': 'B'}], 'products': [{'smiles': 'A'}]},
            {'iteration': 3, 'reactants': [{'smiles': 'A'}], 'products': [{'smiles': 'B'}]},
        ])  # simple A->B->A
        reaction_network = ReactionNetwork.build_smiles_reaction_network(TestEvaluatorNetwork.filename, only_reactions=True)
        self.assertEqual(1, nx.number_strongly_connected_components(reaction_network))  # single cycle = single component
        self.assertEqual(1, nx.number_weakly_connected_components(reaction_network))  # single component

        self.writeToResultsFile([
            {'iteration': 0, 'reactants': [{'smiles': 'A'}], 'products': [{'smiles': 'B'}]},
            {'iteration': 1, 'reactants': [{'smiles': 'B'}], 'products': [{'smiles': 'C'}]},
            {'iteration': 2, 'reactants': [{'smiles': 'C'}], 'products': [{'smiles': 'A'}]},
        ])

        self.writeToResultsFile([
            {'iteration': 0, 'reactants': [{'smiles': 'A'}], 'products': [{'smiles': 'B'}]},
            {'iteration': 1, 'reactants': [{'smiles': 'B'}], 'products': [{'smiles': 'C'}]},
            {'iteration': 2, 'reactants': [{'smiles': 'C'}], 'products': [{'smiles': 'A'}]},
            {'iteration': 3, 'reactants': [{'smiles': 'B'}], 'products': [{'smiles': 'C'}]},
            {'iteration': 4, 'reactants': [{'smiles': 'C'}], 'products': [{'smiles': 'B'}]},
        ])
        reaction_network = ReactionNetwork.build_smiles_reaction_network(TestEvaluatorNetwork.filename, only_reactions=True)
        self.assertEqual(1, nx.number_strongly_connected_components(reaction_network))
        self.assertEqual(1, nx.number_weakly_connected_components(reaction_network))

        self.writeToResultsFile([
            {'iteration': 0, 'reactants': [{'smiles': 'A'}], 'products': [{'smiles': 'B'}]},
            {'iteration': 1, 'reactants': [{'smiles': 'B'}], 'products': [{'smiles': 'A'}]},
            {'iteration': 2, 'reactants': [{'smiles': 'A'}], 'products': [{'smiles': 'C'}]},
            {'iteration': 3, 'reactants': [{'smiles': 'C'}], 'products': [{'smiles': 'A'}]},
        ])
        reaction_network = ReactionNetwork.build_smiles_reaction_network(TestEvaluatorNetwork.filename, only_reactions=True)
        self.assertEqual(1, nx.number_strongly_connected_components(reaction_network))  # [['A', 'B', 'C']]
        self.assertEqual(1, nx.number_weakly_connected_components(reaction_network))  # single component

        self.writeToResultsFile([
            {'iteration': 0, 'reactants': [{'smiles': 'A'}], 'products': [{'smiles': 'B'}]},
            {'iteration': 1, 'reactants': [{'smiles': 'B'}], 'products': [{'smiles': 'A'}]},
            {'iteration': 2, 'reactants': [{'smiles': 'D'}], 'products': [{'smiles': 'C'}]},
            {'iteration': 3, 'reactants': [{'smiles': 'C'}], 'products': [{'smiles': 'D'}]},
        ])
        # self.assertEqual({('A', 'B'):1, ('D', 'C'):1}, ReactionNetwork.discover_potential_reaction_cycles(TestEvaluatorNetwork.filename, min_length=1, min_count=1))
        reaction_network = ReactionNetwork.build_smiles_reaction_network(TestEvaluatorNetwork.filename, only_reactions=True)
        self.assertEqual(2, nx.number_strongly_connected_components(reaction_network))  # single cycle = single component
        self.assertEqual(2, nx.number_weakly_connected_components(reaction_network))  # single component

        self.writeToResultsFile([
            {'iteration': 0, 'reactants': [{'smiles': 'A'}], 'products': [{'smiles': 'B'}]},
            {'iteration': 1, 'reactants': [{'smiles': 'B'}], 'products': [{'smiles': 'A'}]},
            {'iteration': 2, 'reactants': [{'smiles': 'D'}], 'products': [{'smiles': 'C'}]},
            {'iteration': 3, 'reactants': [{'smiles': 'C'}], 'products': [{'smiles': 'D'}]},
            {'iteration': 4, 'reactants': [{'smiles': 'A'}], 'products': [{'smiles': 'E'}]},
        ])
        # self.assertEqual({('A', 'B'):1, ('D', 'C'):1}, ReactionNetwork.discover_potential_reaction_cycles(TestEvaluatorNetwork.filename, min_length=1, min_count=1))
        reaction_network = ReactionNetwork.build_smiles_reaction_network(TestEvaluatorNetwork.filename, only_reactions=True)
#         nx.draw(reaction_network)
#         plt.savefig("network.png")
        self.assertEqual(3, nx.number_strongly_connected_components(reaction_network))  # [['A', 'B'], ['C', 'D'], ['E']]
        self.assertEqual(2, nx.number_weakly_connected_components(reaction_network))  # [['A', 'B', 'E'], ['C', 'D']]


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testSplitMolecule']
    unittest.main()
