"""
@author: thom
"""

import pymunk as pm

from rdkit.Chem import rdPartialCharges
from rdkit.Chem import AllChem as Chem

import random
import math

from kinetic_molecule import KineticMolecule


class ChargedMolecule(KineticMolecule):

    def __init__(self, source, internal_energy=0, kinetic_energy=0, canonize=True, components=None, orientation=None, base_molecule_radius=1):

        self._base_molecule_radius = base_molecule_radius

        self._shapes = []  # never used, just so that shapes aren't gc'd
        self._cluster_centre = {}
        self._cluster_size = {}

        super(KineticMolecule, self).__init__(source, internal_energy=internal_energy, canonize=canonize, components=components)  # set mass
        inertia = pm.moment_for_circle(self.get_mass(), 0, self.get_size(), (0, 0))
        self.body = pm.Body(self.get_mass(), inertia)
        
        if orientation is not None:
            self.set_orientation(orientation)
        else:
            self.set_orientation(random.uniform(0, 2.0 * math.pi))

        for cluster in self.get_clusters():
            center = pm.Vec2d(self._get_cluster_centre(cluster))
            radius = int(math.ceil(self._get_cluster_size(cluster)))
            center.rotate(self.get_orientation())
            shape = pm.Circle(self.body, radius, center)
            shape.elasticity = 0.999  # required for standard collision handler to do a 'perfect' bounce
            shape.friction = 0.0
            shape.charge = self._get_cluster_charge(cluster)
            self._shapes.append(shape)

        self.set_kinetic_energy(kinetic_energy)

    def get_orientation(self):
        return self.body.angle

    def set_orientation(self, orientation):
        self.body.angle = orientation

    def get_state(self):
        state = super(ChargedMolecule, self).get_state()
        state.update({'orientation': self.get_orientation()})
        return state

    def get_size(self):
        return self._get_cluster_size(range(self.GetNumAtoms()))  # size of a cluster containing all the atoms in the molecule...

    def split_molecule(self):
        split_mols = super(ChargedMolecule, self).split_molecule()
        for mol in split_mols:
            mol.set_orientation(mol.get_orientation())  # match orientation
        return split_mols

    def get_charge(self):

        try:
            self.GetAtomWithIdx(0).GetProp("_GasteigerCharge")
        except:
            rdPartialCharges.ComputeGasteigerCharges(self, throwOnParamFailure=True)

        try:
            return self._charge
        except:
            g = [float(a.GetProp("_GasteigerCharge")) for a in self.GetAtoms()]
            positive = len([a for a in g if a > 0.0])
            negative = len([a for a in g if a < 0.0])
            self._charge = (positive - negative) / (1.0 * len(g))

        return self._charge

    def get_clusters(self):
        """Divide the nodes of this molecule into a set of distinct clusters, where all the nodes in a cluster have a
        common Gasteiger polarity."""

        try:
            return self._clusters
        except:
            try:
                self.GetAtomWithIdx(0).GetProp("_GasteigerCharge")
            except:
                rdPartialCharges.ComputeGasteigerCharges(self, throwOnParamFailure=True)
            self._clusters = list()
            already_clustered = set()
            for idx in range(self.GetNumAtoms()):
                if idx not in already_clustered:
                    cluster = self._find_cluster(idx)
                    already_clustered = already_clustered.union(cluster)
                    self._clusters.append(cluster)

        return self._clusters

    def _get_cluster_charge(self, cluster):

        return sum([float(self.GetAtomWithIdx(a).GetProp("_GasteigerCharge")) for a in cluster])

    def _get_cluster_centre(self, cluster):

        cluster_idx = frozenset(cluster)

        try:
            return self._cluster_centre[cluster_idx]
        except:
            try:
                conformer = self.GetConformer(0)
            except:
                id = Chem.EmbedMolecule(self, clearConfs=True)
                assert id == 0
                Chem.UFFOptimizeMolecule(self, confId=0)
                conformer = self.GetConformer(0)

            # cluster centre is the average position of all atoms in the cluster

            self._cluster_centre[cluster_idx] = [0, 0]
            for idx in cluster:
                position = conformer.GetAtomPosition(idx)
                self._cluster_centre[cluster_idx][0] += position[0]
                self._cluster_centre[cluster_idx][1] += position[1]

            self._cluster_centre[cluster_idx][0] /= len(cluster)
            self._cluster_centre[cluster_idx][1] /= len(cluster)
            self._cluster_centre[cluster_idx][0] *= self._base_molecule_radius
            self._cluster_centre[cluster_idx][1] *= self._base_molecule_radius

        return self._cluster_centre[cluster_idx]

    def _get_cluster_size(self, cluster):
        """Model a cluster of molecules as a circle.
        :param cluster: list of nodes in this cluster"""

        cluster_idx = frozenset(cluster)

        try:
            return self._cluster_size[cluster_idx]
        except:
            if len(cluster) == 0:
                return 0
            elif len(cluster) == 1:
                return self._base_molecule_radius * 2.0  # diameter

            center = pm.Vec2d(self._get_cluster_centre(cluster))

            try:
                conformer = self.GetConformer(0)
            except:
                id = Chem.EmbedMolecule(self, clearConfs=True)
                assert id == 0
                Chem.UFFOptimizeMolecule(self, confId=0)
                conformer = self.GetConformer(0)

            max_distance = max([center.get_distance(pm.Vec2d(conformer.GetAtomPosition(idx))) for idx in cluster])

            self._cluster_size[cluster_idx] = self._base_molecule_radius * (max_distance + 2.0)  # extreme sides of molecules positioned distance*radius apart

        return self._cluster_size[cluster_idx]

    def _find_cluster(self, source):
        """Identify the cluster (a maximal set of adjacent nodes sharing a common Gasteiger polarity) of which the source node is a member."""

        cluster = set([source])
        stack = [source]
        cluster_charge = float(self.GetAtomWithIdx(source).GetProp("_GasteigerCharge")) < 0
        while stack:
            next_node = stack.pop()
            neighbors = set([a.GetIdx() for a in self.GetAtomWithIdx(next_node).GetNeighbors()]) - set(cluster)
            for next_node in neighbors:
                next_node_charge = float(self.GetAtomWithIdx(next_node).GetProp("_GasteigerCharge")) < 0
                if cluster_charge == next_node_charge:
                    cluster.add(next_node)
                    stack.append(next_node)
        return cluster

    def __deepcopy__(self, memo):
        # create copy, but don't set kinetic_energy as we need to set velocity explicitly anyway
        mol = ChargedMolecule(self, internal_energy=self._internal_energy, canonize=False)  # don't mess with current structure - just leave it exactly as it is
        mol.set_velocity(self.body.velocity)
        mol.set_orientation(self.get_orientation())
        return mol

    def __str__(self):
        return Chem.MolToSmiles(self)
