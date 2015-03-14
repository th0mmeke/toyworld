"""
Created on 6/05/2013

@author: thom
"""

from collections import OrderedDict
import itertools
import logging
import networkx as nx

from evaluators.evaluator import Evaluator


class ReactionNetwork(object):

    """Various methods to build Reaction Networks from reaction data."""

    @classmethod
    def _make_canonical(cls, path):
        hashes = [hash(x) for x in path]
        start_idx = hashes.index(min(hashes))
        return tuple(path[(i + start_idx) % len(path)] for i in range(len(path)))

    @classmethod
    def _is_reaction(cls, reaction):
        reactant_smiles = [x['smiles'] for x in reaction['reactants']]
        product_smiles = [x['smiles'] for x in reaction['products']]
        return sorted(reactant_smiles) != sorted(product_smiles)  # A reaction if there is a change between reactants and products

    @classmethod
    def _build_reaction_network(cls, results_filename, node_type='smiles', selected_partitions=None, only_reactions=True):
        """Build complete reaction network. Result is represented as a networkx Directed Graph object with molecules as nodes, and edges representing a reaction
        from a reactant molecule to a product, with the weight being the number of times that that path has been taken.
        """

        reaction_network = nx.DiGraph()  # here we track all reactions paths, with weight = number of times we've encountered an edge in the reaction list

        logging.info("Building {} reaction network...".format(node_type.upper()))
        if selected_partitions is not None:
            logging.info("Using selected iterations: {}".format(selected_partitions))

        for block in Evaluator.incr_load_results(results_filename, selected_partitions):
            for reaction in block['reactions']:

                if not only_reactions or ReactionNetwork._is_reaction(reaction):
                    for reactant, product in itertools.product(reaction['reactants'], reaction['products']):
                        # if reactant[node_type] != product[node_type]: # no reason why we shouldn't track unchanging events
                        if reaction_network.has_edge(reactant[node_type], product[node_type]):
                            reaction_network[reactant[node_type]][product[node_type]]['weight'] += 1
                        else:
                            reaction_network.add_edge(reactant[node_type], product[node_type], weight=1)
                            reaction_network.node[product[node_type]]['smiles'] = product['smiles']  # add an attribute with smiles representation for easing later analysis

        logging.info("Built reaction network of {} nodes and {} edges".format(reaction_network.number_of_nodes(), reaction_network.number_of_edges()))

        return reaction_network

    @classmethod
    def build_smiles_reaction_network(cls, results_filename, selected_partitions=None, only_reactions=True):
        """Build complete reaction network. Result is represented as a networkx Directed Graph object with molecules as nodes, and edges representing a reaction
        from a reactant molecule to a product, with the weight being the number of times that that path has been taken.
        """

        return ReactionNetwork._build_reaction_network(results_filename, 'smiles', selected_partitions, only_reactions)

    @classmethod
    def build_molecule_reaction_network(cls, results_filename, selected_partitions=None, only_reactions=True):
        """Build complete reaction network. Result is represented as a networkx Directed Graph object with molecules as nodes, and edges representing a reaction
        from a reactant molecule to a product, with the weight being the number of times that that path has been taken.
        """

        return ReactionNetwork._build_reaction_network(results_filename, 'id', selected_partitions, only_reactions)

    @classmethod
    def discover_potential_reaction_cycles(cls, results_filename, selected_partitions=None, min_length=3, min_count=1, only_reactions=True, **kwargs):
        """Identify all potential reaction cycles from the reaction network graph. A potential cycle is a sequential set of reactions from a molecule
        returning back to any molecule of the same type. Duplicated reactants or products only count once; for example, if we have
        reactions A+B->C+C+D, then we only have one step A or B -> C, not two.

        The straightforward algorithm would be to use networkx.simple_cycles(), but graph sizes lead to recursion depth problems.
        So try reducing size by finding paths in strongly-connected components - cycles can only exist within a component
        See https://networkx.lanl.gov/trac/report/12?sort=reporter&asc=1&page=5 for a similar solution
        Unfortunately, our reaction networks may still contain one huge component...
        Final algorithm exploits an iterative reaction file by checking for the existence of a shortest-path as each reaction is examined in turn
        """

        reaction_network = nx.DiGraph()
        paths = set()

        logging.info("Evaluating reaction graph...")
        if selected_partitions is not None:
            logging.info("Using selected iterations: {}".format(selected_partitions))

        for block in Evaluator.incr_load_results(results_filename, selected_partitions):
            for reaction in block['reactions']:
                if not only_reactions or ReactionNetwork._is_reaction(reaction):

                    # remove duplicate products and reactants - otherwise get extra cycles, one for each duplicate
                    reactants = set([x['smiles'] for x in reaction['reactants']])
                    products = set([x['smiles'] for x in reaction['products']])

                    for reactant, product in itertools.product(reactants, products):

                        if not reaction_network.has_edge(reactant, product):
                            reaction_network.add_edge(reactant, product, weight=1)
                        else:
                            reaction_network[reactant][product]['weight'] += 1

                            # path = nx.shortest_path(reaction_network, product, reactant) # Not correct - may be more than one cycle completed, so need to check all paths!
                        simple_paths = nx.all_simple_paths(reaction_network, product, reactant)  # same results from shortest_paths and simple_paths, but shortest is 2x faster...

                        for path in simple_paths:
                            if len(path) >= min_length:
                                # add this discovered cycle to the list...
                                # first though, make canonical so we group duplicates. Do this by starting cycle at smallest (measured by hash) molecule
                                paths.add(cls._make_canonical(path))

        # Now count the number of times through each cycle - it is the minimum edge weight on the cycle
        cycles = {}
        for path in list(paths):
            count = min([reaction_network[path[idx]][path[(idx + 1) % len(path)]]['weight'] for idx in range(len(path))])
            cycles[path] = count

        compliant_cycles = {cycle: count for cycle, count in cycles.iteritems() if len(cycle) >= min_length and count >= min_count}
        logging.info("{}/{} potential reaction cycles after taking into account min_length and min_count".format(len(compliant_cycles), len(cycles)))
        return compliant_cycles

    @classmethod
    def discover_actual_reaction_cycles(cls, results_filename, selected_partitions=None, node_type="smiles", min_length=3, min_count=1, only_reactions=True, **kwargs):
        """Identify all reaction cycles from the reaction network graph. A cycle is a sequential set of reactions from a molecule
        returning back to a molecule of the same type, and where at least one product molecule (not just type) from one step is a reactant molecule of the
        next. This differs from a potential cycle where the requirement is simply that products and reactants are of the same *type*.

        Note that the final step in a cycle must be to return to a molecule of the same type, not the same molecule, as we started from. That earlier molecule
        won't still be around...
        """
        def find_shortest_paths(network, source, target):  # http://eddmann.com/posts/depth-first-search-and-breadth-first-search-in-python/
            stack = [(source, [source])]
            while stack:
                (vertex, path) = stack.pop()
                for next_node in set(network.predecessors(vertex)) - set(path):
                    if network.node[next_node]['smiles'] == target:
                        yield list(reversed(path + [next_node]))  # reverse order as path found in reverse direction
                    else:
                        stack.append((next_node, path + [next_node]))

        reaction_network = nx.DiGraph()
        cycles = {}
        id_cycles = {}
        active_mols = {}
        cycle_id = 0

        logging.info("Evaluating reaction graph for {}".format(node_type.upper()))
        if selected_partitions is not None:
            logging.info("Using selected iterations: {}".format(selected_partitions))

        for block in Evaluator.incr_load_results(results_filename, selected_partitions):

            for reaction in block['reactions']:

                # print(reaction['iteration'])

                if not only_reactions or ReactionNetwork._is_reaction(reaction):

                    for reactant, product in itertools.product(reaction['reactants'], reaction['products']):
                        # First, check if this combination extends an existing cycle
                        try:
                            candidate_cycles = active_mols[reactant['id']]
                        except:
                            pass
                        else:
                            # yes it does, and candidate_cycles contains pointers to all of them
                            for cycle_id in candidate_cycles:
                                expected_smiles, count = cycles[cycle_id].popitem(last=False)
                                active_mols[reactant['id']].remove(cycle_id)
                                if product['smiles'] == expected_smiles:

                                    # cycle continues - so update count for this mol, and advance to next position in the cycle
                                    cycles[cycle_id][expected_smiles] = count + 1
                                    try:
                                        active_mols[product['id']].append(cycle_id)
                                    except:
                                        active_mols[product['id']] = [cycle_id]
                                else:
                                    cycles[cycle_id][expected_smiles] = count  # just restore it if cycle broken to preserve counts

                        # Second, check if we have a new cycle
                        if reaction_network.has_node(reactant['id']):  # first iteration won't have reactants added yet

                            # check for all paths from reactant['id'] to any node in the network with same smiles as product (product['smiles'])
                            paths = find_shortest_paths(reaction_network, reactant['id'], product['smiles'])
                            for id_path in paths:
                                if len(id_path) >= min_length:
                                    new_cycle = OrderedDict()
                                    for x in id_path:
                                        new_cycle[reaction_network.node[x]['smiles']] = 1
                                    start_smiles, count = new_cycle.popitem(last=False)
                                    new_cycle[start_smiles] = count  # move the current location to last in list; first item is now the next one we'll expect
                                    # add the cycle to the list of active cycles
                                    cycle_id += 1
                                    cycles[cycle_id] = new_cycle  # cycles recorded by smiles
                                    id_cycles[cycle_id] = id_path  # cycles recorded by id
                                    # now add this cycle to the list of cycles of which this molecule is a participant
                                    try:
                                        active_mols[product['id']].append(cycle_id)
                                    except:
                                        active_mols[product['id']] = [cycle_id]

                        reaction_network.add_edge(reactant['id'], product['id'])

                        reaction_network.node[reactant['id']]['smiles'] = reactant['smiles']
                        reaction_network.node[product['id']]['smiles'] = product['smiles']

        # Now aggregate by canonical form
        compliant_cycles = {}
        for cycle_id in cycles.keys():
            if node_type == "id":
                cycle = id_cycles[cycle_id]
            else:
                cycle = cls._make_canonical(tuple(cycles[cycle_id].keys()))

            try:
                compliant_cycles[cycle] += 1
            except:
                compliant_cycles[cycle] = 1

        compliant_cycles = {cycle: count for cycle, count in compliant_cycles.iteritems() if len(cycle) >= min_length and count >= min_count}

        logging.info("{}/{} actual reaction cycles after taking into account min_length={} and min_count={}".format(len(compliant_cycles), len(cycles), min_length, min_count))
        return compliant_cycles
