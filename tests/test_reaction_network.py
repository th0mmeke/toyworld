"""
Created on 3/05/2013

@author: thom
"""
import unittest
import os
import cPickle

from evaluators.reaction_network import ReactionNetwork

import config


class TestReactionNetwork(unittest.TestCase):

    filename = os.path.join(config.DataDir, 'test/test_reaction_network.data')
    
    def tearDown(self):
        if os.path.exists(TestReactionNetwork.filename):
            os.remove(TestReactionNetwork.filename)

    def writeToResultsFile(self, reactions):
        with open(TestReactionNetwork.filename, "wb") as f:
            cPickle.dump({'block': {'start_block': 0, 'end_block': len(reactions) - 1, 'reactions': reactions}}, f, 2)
            f.flush()  # force python file write
            os.fsync(f.fileno())  # force OS file write

    def testIsReaction(self):
        self.assertFalse(ReactionNetwork._is_reaction({'reactants':[{'id':1,'smiles':'a'}],'products':[{'id':2,'smiles':'a'}]}))
        self.assertTrue(ReactionNetwork._is_reaction({'reactants':[{'id':1,'smiles':'b'}],'products':[{'id':2,'smiles':'a'}]}))
        self.assertFalse(ReactionNetwork._is_reaction({'reactants':[{'id':1,'smiles':'a'},{'id':2,'smiles':'b'}],'products':[{'id':3,'smiles':'b'},{'id':4,'smiles':'a'}]}))

    def testIsReactionRegression(self):
        self.writeToResultsFile([
                {'iteration':0, 'reactants':[{'smiles':'a'}, {'smiles':'b'}], 'products':[{'smiles':'c'}, {'smiles':'d'}]},
                {'iteration':1, 'reactants':[{'smiles':'d'}, {'smiles':'e'}], 'products':[{'smiles':'a'}]},
                {'iteration':2, 'reactants':[{'smiles':'f'}, {'smiles':'f'}], 'products':[{'smiles':'f'},{'smiles':'f'}]},
                {'iteration':3, 'reactants':[{'smiles':'f'}, {'smiles':'g'}], 'products':[{'smiles':'b'}]},
                ])
        reaction_network = ReactionNetwork.build_smiles_reaction_network(TestReactionNetwork.filename,only_reactions=True)
        full_network = ReactionNetwork.build_smiles_reaction_network(TestReactionNetwork.filename,only_reactions=False)
        self.assertNotEqual(reaction_network.edges(),full_network.edges()) # number of nodes
        self.assertTrue(0 < len(ReactionNetwork.build_smiles_reaction_network(TestReactionNetwork.filename,only_reactions=False)))

    def testSimple(self):
        # a->a->b->a => a->b
        self.writeToResultsFile([
                {'iteration':1, 'reactants':[{'id':2, 'smiles':'a'}], 'products':[{'id':3, 'smiles':'b'}]},
                {'iteration':2, 'reactants':[{'id':3, 'smiles':'b'}], 'products':[{'id':10, 'smiles':'a'}]}
                ])
        self.assertEqual({('a', 'b'): 1}, ReactionNetwork.discover_actual_reaction_cycles(TestReactionNetwork.filename, min_length=1, min_count=1))

        self.writeToResultsFile([
                {'iteration':0, 'reactants':[{'id':1, 'smiles':'A'}, {'id':2, 'smiles':'B'}], 'products':[{'id':3, 'smiles':'C'},{'id':4, 'smiles':'A'},{'id':5, 'smiles':'D'}]},
                {'iteration':1, 'reactants':[{'id':4, 'smiles':'A'}, {'id':6, 'smiles':'E'}], 'products':[{'id':7, 'smiles':'C'},{'id':8, 'smiles':'A'},{'id':9, 'smiles':'D'}]},
                {'iteration':2, 'reactants':[{'id':7, 'smiles':'C'}, {'id':15, 'smiles':'D'}], 'products':[{'id':10, 'smiles':'D'},{'id':11, 'smiles':'F'}]},
                {'iteration':3, 'reactants':[{'id':11, 'smiles':'F'}, {'id':12, 'smiles':'A'}], 'products':[{'id':13, 'smiles':'A'},{'id':14, 'smiles':'B'}]},
                ])
        # Actual is longer than shortest potential cycle, but actual is still a subset of potential cycles
        actual_cycles = ReactionNetwork.discover_actual_reaction_cycles(TestReactionNetwork.filename, min_length=2, min_count=1)
        self.assertEqual({('A', 'C', 'F'): 1, ('A', 'C', 'F', 'B'): 1},actual_cycles) # self-loops are length 2...

    def testShortestCycles(self):
        # a->a->b->a => a->b
        self.writeToResultsFile([
                {'iteration':0, 'reactants':[{'id':1, 'smiles':'a'}], 'products':[{'id':2, 'smiles':'a'}]},
                {'iteration':1, 'reactants':[{'id':2, 'smiles':'a'}], 'products':[{'id':3, 'smiles':'b'}]},
                {'iteration':2, 'reactants':[{'id':3, 'smiles':'b'}], 'products':[{'id':10, 'smiles':'a'}]}
                ])
        self.assertEqual({('a', 'b'): 1}, ReactionNetwork.discover_actual_reaction_cycles(TestReactionNetwork.filename, min_length=1, min_count=1))

    def testA(self):
        self.writeToResultsFile([
            {'iteration':0, 'reactants':[{'smiles':'A'}], 'products':[{'smiles':'B'}]},
            {'iteration':1, 'reactants':[{'smiles':'B'}], 'products':[{'smiles':'A'}]},
            {'iteration':2, 'reactants':[{'smiles':'B'}], 'products':[{'smiles':'A'}]},
            {'iteration':3, 'reactants':[{'smiles':'A'}], 'products':[{'smiles':'B'}]},
            ])
        self.assertEqual({('A', 'B'):2}, ReactionNetwork.discover_potential_reaction_cycles(TestReactionNetwork.filename, min_length=1, min_count=1))

    def testB(self):
        self.writeToResultsFile([
            {'iteration':0, 'reactants':[{'smiles':'A'}], 'products':[{'smiles':'B'}]},
            {'iteration':1, 'reactants':[{'smiles':'B'}], 'products':[{'smiles':'C'}]},
            {'iteration':2, 'reactants':[{'smiles':'C'}], 'products':[{'smiles':'A'}]},
            {'iteration':3, 'reactants':[{'smiles':'B'}], 'products':[{'smiles':'C'}]},
            {'iteration':4, 'reactants':[{'smiles':'C'}], 'products':[{'smiles':'B'}]},
            ])
        self.assertEqual({('A', 'B', 'C'):1, ('B', 'C'):1}, ReactionNetwork.discover_potential_reaction_cycles(TestReactionNetwork.filename, min_length=1, min_count=1))

    def testDetectCycles(self):
        # logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.INFO)
        self.writeToResultsFile([
                {'iteration':0, 'reactants':[{'id':1, 'smiles':'a'}, {'id':2, 'smiles':'b'}], 'products':[{'id':3, 'smiles':'c'}]},
                {'iteration':1, 'reactants':[{'id':4, 'smiles':'c'}, {'id':6, 'smiles':'e'}], 'products':[{'id':1, 'smiles':'a'}]},
                {'iteration':2, 'reactants':[{'id':3, 'smiles':'c'}, {'id':9, 'smiles':'g'}], 'products':[{'id':10, 'smiles':'d'}]},
                {'iteration':3, 'reactants':[{'id':10, 'smiles':'d'}, {'id':11, 'smiles':'h'}], 'products':[{'id':10, 'smiles':'a'}]},
                ])
        # Actual is longer than shortest potential cycle, but actual is still a subset of potential cycles
        actual_cycles = ReactionNetwork.discover_actual_reaction_cycles(TestReactionNetwork.filename, min_length=1, min_count=1)
        potential_cycles = ReactionNetwork.discover_potential_reaction_cycles(TestReactionNetwork.filename, min_length=1, min_count=1)
        self.assertEqual({('a', 'c', 'd'): 1}, actual_cycles)
        self.assertEqual({('a', 'c'): 1, ('a', 'c', 'd'): 1}, potential_cycles)
        self.assertTrue(set(actual_cycles.keys()).issubset(set(potential_cycles.keys())))

        self.writeToResultsFile([
                {'iteration':0, 'reactants':[{'id':1, 'smiles':'a'}, {'id':2, 'smiles':'b'}], 'products':[{'id':3, 'smiles':'c'}, {'id':4, 'smiles':'d'}]}
                ])
        self.assertEqual({}, ReactionNetwork.discover_actual_reaction_cycles(TestReactionNetwork.filename, min_length=1, min_count=1))
        self.writeToResultsFile([
                {'iteration':0, 'reactants':[{'id':1, 'smiles':'a'}, {'id':2, 'smiles':'b'}], 'products':[{'id':3, 'smiles':'c'}, {'id':4, 'smiles':'d'}]},
                {'iteration':1, 'reactants':[{'id':4, 'smiles':'d'}, {'id':6, 'smiles':'e'}], 'products':[{'id':7, 'smiles':'a'}]},
                {'iteration':2, 'reactants':[{'id':8, 'smiles':'f'}, {'id':9, 'smiles':'g'}], 'products':[{'id':10, 'smiles':'b'}]}
                ])
        self.assertEqual({('a', 'd'): 1}, ReactionNetwork.discover_actual_reaction_cycles(TestReactionNetwork.filename, min_length=1, min_count=1))
        self.assertEqual({('a', 'd'): 1}, ReactionNetwork.discover_potential_reaction_cycles(TestReactionNetwork.filename, min_length=1, min_count=1))
        self.writeToResultsFile([
                {'iteration':0, 'reactants':[{'id':1, 'smiles':'a'}, {'id':2, 'smiles':'b'}], 'products':[{'id':3, 'smiles':'c'}, {'id':4, 'smiles':'d'}]},
                {'iteration':1, 'reactants':[{'id':5, 'smiles':'d'}, {'id':6, 'smiles':'e'}], 'products':[{'id':7, 'smiles':'a'}]},  # different d, so not a cycle!
                {'iteration':2, 'reactants':[{'id':8, 'smiles':'f'}, {'id':9, 'smiles':'g'}], 'products':[{'id':10, 'smiles':'b'}]}
                ])
        self.assertEqual({}, ReactionNetwork.discover_actual_reaction_cycles(TestReactionNetwork.filename, min_length=1, min_count=1))

        self.writeToResultsFile([
            {'iteration':0, 'reactants':[{'smiles':'A'}], 'products':[{'smiles':'A'}]},
            ])
        self.assertEqual({}, ReactionNetwork.discover_potential_reaction_cycles(TestReactionNetwork.filename, min_length=1, min_count=1))

        self.writeToResultsFile([
            {'iteration':0, 'reactants':[{'smiles':'A'}], 'products':[{'smiles':'B'}]},
            {'iteration':1, 'reactants':[{'smiles':'B'}], 'products':[{'smiles':'A'}]},
            ])
        self.assertEqual({('A', 'B'):1}, ReactionNetwork.discover_potential_reaction_cycles(TestReactionNetwork.filename, min_length=1, min_count=1))

        self.writeToResultsFile([
            {'iteration':0, 'reactants':[{'smiles':'A'}], 'products':[{'smiles':'B'}, {'smiles':'B'}]},
            {'iteration':1, 'reactants':[{'smiles':'B'}], 'products':[{'smiles':'C'}]},
            {'iteration':2, 'reactants':[{'smiles':'C'}], 'products':[{'smiles':'A'}]},
            ])
        self.assertEqual({('A', 'B', 'C'):1}, ReactionNetwork.discover_potential_reaction_cycles(TestReactionNetwork.filename, min_length=1, min_count=1))

        self.writeToResultsFile([
            {'iteration':0, 'reactants':[{'smiles':'A'}], 'products':[{'smiles':'B'}]},
            {'iteration':1, 'reactants':[{'smiles':'B'}], 'products':[{'smiles':'A'}]},
            {'iteration':2, 'reactants':[{'smiles':'B'}], 'products':[{'smiles':'A'}]},
            {'iteration':3, 'reactants':[{'smiles':'A'}], 'products':[{'smiles':'B'}]},
            ])
        self.assertEqual({('A', 'B'):2}, ReactionNetwork.discover_potential_reaction_cycles(TestReactionNetwork.filename, min_length=1, min_count=1))

        self.writeToResultsFile([
            {'iteration':0, 'reactants':[{'smiles':'A'}], 'products':[{'smiles':'B'}]},
            {'iteration':1, 'reactants':[{'smiles':'B'}], 'products':[{'smiles':'C'}]},
            {'iteration':2, 'reactants':[{'smiles':'C'}], 'products':[{'smiles':'A'}]},
            ])
        self.assertEqual({('A', 'B', 'C'):1}, ReactionNetwork.discover_potential_reaction_cycles(TestReactionNetwork.filename, min_length=1, min_count=1))

        self.writeToResultsFile([
            {'iteration':0, 'reactants':[{'smiles':'A'}], 'products':[{'smiles':'B'}]},
            {'iteration':1, 'reactants':[{'smiles':'B'}], 'products':[{'smiles':'C'}]},
            {'iteration':2, 'reactants':[{'smiles':'C'}], 'products':[{'smiles':'A'}]},
            {'iteration':3, 'reactants':[{'smiles':'B'}], 'products':[{'smiles':'C'}]},
            {'iteration':4, 'reactants':[{'smiles':'C'}], 'products':[{'smiles':'B'}]},
            ])
        self.assertEqual({('A', 'B', 'C'):1, ('B', 'C'):1}, ReactionNetwork.discover_potential_reaction_cycles(TestReactionNetwork.filename, min_length=1, min_count=1))

        self.writeToResultsFile([
            {'iteration':0, 'reactants':[{'smiles':'A'}], 'products':[{'smiles':'B'}]},
            {'iteration':1, 'reactants':[{'smiles':'B'}], 'products':[{'smiles':'A'}]},
            {'iteration':2, 'reactants':[{'smiles':'A'}], 'products':[{'smiles':'C'}]},
            {'iteration':3, 'reactants':[{'smiles':'C'}], 'products':[{'smiles':'A'}]},
            ])
        self.assertEqual({('A', 'B'):1, ('A', 'C'):1}, ReactionNetwork.discover_potential_reaction_cycles(TestReactionNetwork.filename, min_length=1, min_count=1))

        self.writeToResultsFile([
            {'iteration':0, 'reactants':[{'smiles':'A'}], 'products':[{'smiles':'B'}]},
            {'iteration':1, 'reactants':[{'smiles':'B'}], 'products':[{'smiles':'B'}]},
            ])
        self.assertEqual({}, ReactionNetwork.discover_potential_reaction_cycles(TestReactionNetwork.filename, min_length=1, min_count=1))

        self.writeToResultsFile([
                {'iteration':0, 'reactants':[{'id':1, 'smiles':'a'}, {'id':2, 'smiles':'b'}], 'products':[{'id':3, 'smiles':'c'}, {'id':4, 'smiles':'d'}]},
                {'iteration':1, 'reactants':[{'id':4, 'smiles':'d'}, {'id':5, 'smiles':'e'}], 'products':[{'id':6, 'smiles':'f'}]},
                {'iteration':2, 'reactants':[{'id':6, 'smiles':'f'}, {'id':7, 'smiles':'g'}], 'products':[{'id':8, 'smiles':'a'}]},
                ])
        # reaction_network = ReactionNetwork.build_smiles_reaction_network(TestReactionNetwork.filename)
        self.assertEqual({('a', 'd', 'f'): 1}, ReactionNetwork.discover_potential_reaction_cycles(TestReactionNetwork.filename, min_length=1, min_count=1))

    def testIdNotSmiles(self):
        self.writeToResultsFile([
                {'iteration':0, 'reactants':[{'id':1, 'smiles':'A'}, {'id':2, 'smiles':'B'}], 'products':[{'id':3, 'smiles':'C'},{'id':4, 'smiles':'A'},{'id':5, 'smiles':'D'}]},
                {'iteration':1, 'reactants':[{'id':4, 'smiles':'A'}, {'id':6, 'smiles':'E'}], 'products':[{'id':7, 'smiles':'C'},{'id':8, 'smiles':'A'},{'id':9, 'smiles':'D'}]},
                {'iteration':2, 'reactants':[{'id':7, 'smiles':'C'}, {'id':15, 'smiles':'D'}], 'products':[{'id':10, 'smiles':'D'},{'id':11, 'smiles':'F'}]},
                {'iteration':3, 'reactants':[{'id':11, 'smiles':'F'}, {'id':12, 'smiles':'A'}], 'products':[{'id':13, 'smiles':'A'},{'id':14, 'smiles':'B'}]},
                ])
        actual_cycles = ReactionNetwork.discover_actual_reaction_cycles(TestReactionNetwork.filename, id_cycles=True, min_length=3, min_count=1)

    def testRegressionReactionCycles(self):
        # last link to complete a cycle appears two or more times as product in last reaction -> cycle added multiple times (once for each product)    
        self.writeToResultsFile([
            {'iteration':0, 'reactants':[{'smiles':'B'}], 'products':[{'smiles':'C'}]},
            {'iteration':1, 'reactants':[{'smiles':'C'}], 'products':[{'smiles':'A'}]},
            {'iteration':2, 'reactants':[{'smiles':'A'}], 'products':[{'smiles':'B'}, {'smiles':'B'}]},
            ])
        self.assertEqual({('A', 'B', 'C'):1}, ReactionNetwork.discover_potential_reaction_cycles(TestReactionNetwork.filename, min_count=1, min_length=1))

    def testRegressionReactionCycle2(self):
        # Test correct identification of multiple cycles completed by same reaction
        self.writeToResultsFile([
                {'iteration':0, 'reactants':[{'smiles':'A'}, {'smiles':'a'}], 'products':[{'smiles':'B'}]},
                {'iteration':1, 'reactants':[{'smiles':'B'}, {'smiles':'a'}], 'products':[{'smiles':'C'}]},
                {'iteration':2, 'reactants':[{'smiles':'A'}, {'smiles':'a'}], 'products':[{'smiles':'D'}]},
                {'iteration':3, 'reactants':[{'smiles':'D'}, {'smiles':'a'}], 'products':[{'smiles':'E'}]},
                {'iteration':4, 'reactants':[{'smiles':'E'}, {'smiles':'a'}], 'products':[{'smiles':'C'}]},
                {'iteration':5, 'reactants':[{'smiles':'C'}, {'smiles':'a'}], 'products':[{'smiles':'A'}]},
                ])
        potential_cycles = ReactionNetwork.discover_potential_reaction_cycles(TestReactionNetwork.filename, min_length=1, min_count=1)
        self.assertEqual({('A', 'B','C'): 1, ('A', 'D', 'E', 'C'): 1}, potential_cycles)

    def testRegressionReactionCycle3(self):
        # Live data scenario
        self.writeToResultsFile([
                {'iteration':0, 'reactants':[{'id':1, 'smiles':'A'}, {'id':2, 'smiles':'B'}], 'products':[{'id':3, 'smiles':'C'},{'id':4, 'smiles':'A'},{'id':5, 'smiles':'D'}]},
                {'iteration':1, 'reactants':[{'id':4, 'smiles':'A'}, {'id':6, 'smiles':'E'}], 'products':[{'id':7, 'smiles':'C'},{'id':8, 'smiles':'A'},{'id':9, 'smiles':'D'}]},
                {'iteration':2, 'reactants':[{'id':7, 'smiles':'C'}, {'id':15, 'smiles':'D'}], 'products':[{'id':10, 'smiles':'D'},{'id':11, 'smiles':'F'}]},
                {'iteration':3, 'reactants':[{'id':11, 'smiles':'F'}, {'id':12, 'smiles':'A'}], 'products':[{'id':13, 'smiles':'A'},{'id':14, 'smiles':'B'}]},
                ])

        actual_cycles = ReactionNetwork.discover_actual_reaction_cycles(TestReactionNetwork.filename, min_length=3, min_count=1)
        self.assertEqual({('A', 'C', 'F'): 1, ('A', 'C', 'F', 'B'): 1}, actual_cycles)

        potential_cycles = ReactionNetwork.discover_potential_reaction_cycles(TestReactionNetwork.filename, min_length=3, min_count=1)
        self.assertEqual({('A', 'B', 'D', 'F'): 1, ('A', 'C', 'F'): 1,
                         ('A', 'C', 'D', 'F'): 1,
                         ('A', 'B', 'C', 'F'): 1,
                         ('A', 'B', 'C', 'D', 'F'): 1,
                         ('B', 'D', 'F'): 1,
                         ('A', 'D', 'F'): 1,
                         ('A', 'C', 'F', 'B'): 1,
                         ('A', 'D', 'F', 'B'): 1,
                         ('A', 'C', 'D', 'F', 'B'): 1,
                         ('B', 'C', 'F'): 1,
                         ('B', 'C', 'D', 'F'): 1},
                         potential_cycles)

        actual_cycles_set = set(actual_cycles.keys())
        potential_cycles_set = set(potential_cycles.keys())
        self.assertEqual(set([]),actual_cycles_set-potential_cycles_set)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testInit']
    unittest.main()