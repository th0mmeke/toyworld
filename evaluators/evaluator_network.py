"""
Created on 6/05/2013

@author: thom
"""

from evaluator import Evaluator
from reaction_network import ReactionNetwork

import networkx as nx

import math
import collections
import string
from abc import ABCMeta, abstractmethod


class EvaluatorNetwork(Evaluator):

    __metaclass__ = ABCMeta

    def get_result_titles(self):
        return ["Nodes", "Edges", "Mean Degree", "Std.Dev of Degree.", "Min Degree", "Max Degree", "Number strongly connected components", "Number weakly connected components"]

    @abstractmethod
    def evaluate(self, results_filename, reaction_network=None, **kwargs):

        raw_degree = nx.degree(reaction_network)
        degree = raw_degree.values()
        non_zero_in_degree = [x for x in reaction_network.in_degree().values() if x > 0]
        average_in_degree = sum(non_zero_in_degree) * 1.0 / len(non_zero_in_degree)  # average number of REACTANTS
        non_zero_out_degree = [x for x in reaction_network.out_degree().values() if x > 0]
        average_out_degree = sum(non_zero_out_degree) * 1.0 / len(non_zero_out_degree)  # average number of PRODUCTS
        mean = sum(degree) * 1.0 / len(degree)
        variance = sum([x * x for x in degree]) / len(degree) - mean * mean
        sorted_degree = collections.OrderedDict(sorted(raw_degree.items(), key=lambda t: t[1]))
        min_degree = sorted_degree.popitem(last=False)
        max_degree = sorted_degree.popitem(last=True)
        nodes = reaction_network.number_of_nodes()
        edges = reaction_network.number_of_edges()

        print("Nodes = {}, Edges = {}".format(nodes, edges))
        print("In graph theory, the degree (or valency) of a vertex of a graph is the number of edges incident to the vertex, with loops counted twice. - Wikipedia")
        print("Average reactants = {}, average products = {}, average degree = {}, standard deviation = {}".format(average_in_degree, average_out_degree, mean, math.sqrt(variance)))
        print("Min degree of {} at node {}".format(min_degree[1], min_degree[0]))
        print("Max degree of {} at node {}".format(max_degree[1], max_degree[0]))
        products = [reaction_network.node[x[1]]['smiles'] for x in reaction_network.out_edges(max_degree[0])]
        print("Products of max degree node = {}".format(string.join(products, ", ")))

        number_strongly_connected_components = nx.number_strongly_connected_components(reaction_network)
        print(
            "In the mathematical theory of directed graphs, a graph is said to be strongly connected if every vertex is reachable \
            from every other vertex. The strongly connected components of an arbitrary directed graph form a partition into subgraphs that are themselves strongly connected. - Wikipedia")
        print("{} strongly connected components".format(number_strongly_connected_components))
        print(nx.strongly_connected_components(reaction_network))

        number_weakly_connected_components = nx.number_weakly_connected_components(reaction_network)
        print("A directed graph is called weakly connected if replacing all of its directed edges with undirected edges produces a connected \
            (undirected) graph. - Wikipedia")
        print(
            "A weakly connected component is a maximal subgraph of a directed graph such that for every pair of vertices  u, v in the subgraph,\
            there is an undirected path from u to v and a directed path from v to u - Wolfram MathWorld")
        print("{} weakly connected components".format(number_weakly_connected_components))

        undirected_reaction_network = nx.convert_node_labels_to_integers(nx.Graph(reaction_network))
        cliques = nx.find_cliques(undirected_reaction_network)
        print("A k-clique community is the union of all cliques of size k that can be reached through adjacent (sharing k-1 nodes) k-cliques.")
        for k in range(2, 10):
            print("{} {}-cliques".format(len(list(nx.k_clique_communities(undirected_reaction_network, k, cliques))), k))

        # nx.draw(undirected_reaction_network)#,pos=nx.spring_layout(undirected_reaction_network))
        # plt.savefig("network.png")

        return nodes, edges, mean, math.sqrt(variance), min_degree[1], max_degree[1], number_strongly_connected_components, number_weakly_connected_components


class EvaluatorActualNetwork(EvaluatorNetwork):

    """Show some standard graph metrics (connected components, degree and so on) for an actual reaction network."""

    def evaluate(self, results_filename, **kwargs):

        reaction_network = ReactionNetwork.build_molecule_reaction_network(results_filename, only_reactions=True)
        return super(EvaluatorActualNetwork, self).evaluate(results_filename, reaction_network)


class EvaluatorPotentialNetwork(EvaluatorNetwork):

    """Show some standard graph metrics (connected components, degree and so on) for a potential reaction network."""

    def evaluate(self, results_filename, **kwargs):

        reaction_network = ReactionNetwork.build_smiles_reaction_network(results_filename, only_reactions=True)
        return super(EvaluatorPotentialNetwork, self).evaluate(results_filename, reaction_network)
