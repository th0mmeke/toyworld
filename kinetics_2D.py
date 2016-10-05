"""
Created on 22/03/2013

@author: thom
"""

import random
import math
import logging

from rdkit.Chem import AllChem as Chem

from ULPS import Float_t
import config


class Kinetics2D(object):

    @classmethod
    def get_ke(cls, m, x, y):
        return 0.5 * m * (x * x + y * y)

    @classmethod
    def get_speed(cls, x, y):
        return math.sqrt(x * x + y * y)

    @classmethod
    def radial_to_xy(cls, theta=None, r=None):
        """Always returns a 2-D x,y"""
        if theta is None:
            theta = random.uniform(0, 2.0 * math.pi)
        if r is None:
            r = random.uniform(0, 1)
        y = math.sin(theta) * r
        x = math.cos(theta) * r
        return x, y

    @classmethod
    def xy_to_radial(cls, x, y):
        """Always returns a 2-D theta,r"""
        r = math.hypot(x, y)
        theta = math.atan2(y, x)
        return theta, r

    @classmethod
    def get_distance(cls, l1, l2):
        return math.sqrt(sum([(_l1 - _l2) * (_l1 - _l2) for _l1, _l2 in zip(l1, l2)]))

    @classmethod
    def get_CM_energy(cls, mols):
        """Return KE of Centre of Mass: _ke = 1/2mv^2, where mv for the centre of mass = sum (mi * vi) for all particles i

        :param mols: list of Molecule"""

        total_mass = sum([mol.get_mass() for mol in mols])
        return cls.get_ke(total_mass, *cls.get_CM_velocity(mols))

    @classmethod
    def get_CM_velocity(cls, mols):
        """Return the momentum (mdx,mdy) of the centre of mass for these particles"""

        cm_momentum = [0, 0]
        total_mass = sum([mol.get_mass() for mol in mols])
        for mol in mols:
            cm_momentum += mol.get_velocity() * mol.get_mass()
        CM_velocity = cm_momentum / total_mass
        logging.debug("CM velocity = {}".format(CM_velocity))
        return CM_velocity

        # for mol in mols:
        #     cm_momentum[0] += mol.get_mass() * mol.get_velocity()[0]
        #     cm_momentum[1] += mol.get_mass() * mol.get_velocity()[1]
        # return [mv / total_mass for mv in cm_momentum]

    @classmethod
    def inelastic_collision(cls, reactant_mols, product_mols, energy_delta):
        """Determine velocities of product molecules following a collision of reactant molecules, for between one and three product molecules.

        Model as a collision, followed by an explosion, meaning that the total momentum of the system is conserved - if two particles, each has equal and opposite momentum in CoM frame
        Assume an impulse, or force, splitting the particles apart, acting equally on each particle
        Then impulse J = mv2-mv1 and so momentum change will be the same for all particles
        Implies that for two particles, equal and opposite mv in CoM frame, and for three particles, mv arranged in equilateral triangle

        Post-conditions:
        1. Sum in_mass = Sum out_mass - although #in_molecules ne #out_molecules
        2. Vector speed and direction of CoM remains constant
        3. in_KE + in_PE + in_IE = Sum out_KE + out_PE + out_IE or in_KE - delta_KE = out_KE

        :param reactant_mols: reactants - must have total KE > 0
        :type reactant_mols: list of Molecule
        :param product_mols: products of reaction - must be 1, 2 or 3 products only
        :type product_mols: list of Molecule
        :param energy_delta: final KE = initial KE - energy_delta
        """

        def total_mv(mv):
            totals = [0, 0]
            for mv_ in mv:
                for dim in range(len(totals)):
                    totals[dim] += mv_[dim]
            return totals

        if len(product_mols) < 1 or len(product_mols) > 3:
            raise ValueError()
        logging.debug("reactant_mols = {}, product_mols = {}".format([Chem.MolToSmiles(mol) for mol in reactant_mols], [Chem.MolToSmiles(mol) for mol in product_mols]))

        in_v = [mol.get_velocity() for mol in reactant_mols]
        in_mass = [mol.get_mass() for mol in reactant_mols]
        in_mv = [[m * v_ for v_ in v] for m, v in zip(in_mass, in_v)]
        in_ke = sum([mol.get_kinetic_energy() for mol in reactant_mols])
        in_ie = sum([mol.get_internal_energy() for mol in reactant_mols])

        # Velocity of centre of mass after collision
        # Momentums add to zero in the CoM frame

        out_mass = [mol.get_mass() for mol in product_mols]
        cm_in_v = cls.get_CM_velocity(reactant_mols)
        cm_in_radial_v = cls.xy_to_radial(*cm_in_v)

        # Bound energy_of_collision to above zero (rounding errors for small values)
        # consistent sense with that in discover_reaction - final_PE = initial_PE + energy_delta => final_KE = initial_KE - energy_delta
        energy_of_collision = max(0, in_ke + in_ie - energy_delta - cls.get_CM_energy(reactant_mols))
        if energy_of_collision <= 0:
            raise ValueError

        out_v_in_CoM_frame = []

        if len(out_mass) == 1:
            # One out particle is stationary in out_CoM frame
            IE = energy_of_collision  # inelastic collision -> loss of KE -> must go to IE
            out_v_in_CoM_frame.append([0, 0])

        elif len(out_mass) == 2:
            ke_in_CM_frame = random.uniform(0, energy_of_collision)
            IE = energy_of_collision - ke_in_CM_frame
            mv = math.sqrt((2.0 * ke_in_CM_frame * out_mass[0] * out_mass[1]) / (out_mass[0] + out_mass[1]))
            out_v_in_CoM_frame.append(cls.radial_to_xy(cm_in_radial_v[0] + math.pi * 0.5, mv))
            out_v_in_CoM_frame.append(cls.radial_to_xy(cm_in_radial_v[0] + math.pi * 1.5, mv))

        elif len(out_mass) == 3:
            # Sum of vector momentums = 0, and in centre of momentum frame arranged as equilateral triangle, side mv
            # Must then convert to velocities by dividing by particle mass, which means no longer equilateral...but unimportant, as only needed equilateral to initially arrange
            ke_in_CM_frame = random.uniform(0, energy_of_collision)  # The energy of the collision - over and above the energy of the centre of mass, which is invariant
            IE = energy_of_collision - ke_in_CM_frame
            mv = math.sqrt((2.0 * ke_in_CM_frame * out_mass[0] * out_mass[1] * out_mass[2]) / (out_mass[0] * out_mass[1] + out_mass[1] * out_mass[2] + out_mass[0] * out_mass[2]))
            out_v_in_CoM_frame.append(cls.radial_to_xy(cm_in_radial_v[0] + math.pi / 3.0, mv))
            out_v_in_CoM_frame.append(cls.radial_to_xy(cm_in_radial_v[0] - math.pi / 3.0, mv))
            out_v_in_CoM_frame.append(cls.radial_to_xy(cm_in_radial_v[0] + math.pi, mv))

        # Now convert from momentums to velocities by scaling by 1/mass
        out_v_in_CoM_frame = [[mv_component / mass for mv_component in particle_mv] for particle_mv, mass in zip(out_v_in_CoM_frame, out_mass)]
        # Finally convert back from CoM frame to lab frame
        out_v = [[v_ + cm_v_ for v_, cm_v_ in zip(v, cm_in_v)] for v in out_v_in_CoM_frame]

        #########################
        # Confirm post-conditions
        # 1. Mass
        assert Float_t.almost_equal(sum(in_mass), sum(out_mass))

        # 2. Momentum

        out_mv = [[m * v_ for v_ in v] for m, v in zip(out_mass, out_v)]
        in_mv_total = total_mv(in_mv)
        out_mv_total = total_mv(out_mv)
        logging.debug("IN MV = {}, OUT MV = {}".format(in_mv_total, out_mv_total))
        for in_, out_ in zip(in_mv_total, out_mv_total):
            assert Float_t.almost_equal(in_, out_)

        # 3. Energy
        out_ke = sum([cls.get_ke(m, *v) for m, v in zip(out_mass, out_v)])
        logging.debug("IN_KE + IN_IE = {}+{} = {}, OUT_KE + DELTA + IE = {} + {} + {} = {}".format(in_ke, in_ie, in_ke + in_ie, out_ke, energy_delta, IE, out_ke + energy_delta + IE))

        assert Float_t.almost_equal(in_ke + in_ie, out_ke + energy_delta + IE, max_diff=config.EnergyTolerance)

        return out_v, IE
