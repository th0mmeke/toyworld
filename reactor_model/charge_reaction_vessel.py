"""
Created on 22/03/2013

@author: thom
"""

import logging
import os

from parameters import Parameters
from reactor_model.spatial_reaction_vessel import SpatialReactionVessel
from ULPS import Float_t


class ChargeReactionVessel(SpatialReactionVessel):

    """
    n-dimensional structure of fixed size.  Molecules bounce off virtual walls.
    """

    distance_cutoff = 50.0
    magnitude_cutoff = 1E-5

    def __init__(self, chemistry, population=None, parameters=Parameters(), product_selection_strategy="energy"):

        super(ChargeReactionVessel, self).__init__(chemistry, population, parameters=parameters, product_selection_strategy=product_selection_strategy)
        self._dipole_force_constant = float(parameters.get('DipoleForceConstant'))
        logging.info("ChargeReactionVessel with dipole force of {} and delta_t of {}".format(self._dipole_force_constant, self._delta_t))

    def __del__(self):
        super(ChargeReactionVessel, self).__del__()

    def step(self):

        for body, mol in self._bodies.iteritems():
            body.force_vector = []
            nearest_points = self._space.nearest_point_query(mol.get_position(), ChargeReactionVessel.distance_cutoff)  # world co-ordinates
            # e.g., 'distance': 36.06281003619938, 'shape': <pymunk.Circle object at 0x107d7a610>, 'point': Vec2d(-27.2849762796, 322.918629605)

            for shape in body.shapes:
                # Find all shapes closer than distance_cutoff
                shape_position = body.local_to_world(shape.offset)  # world co-ordinates

                for near_point in nearest_points:
                    other_shape = near_point['shape']
                    if other_shape.body != body and other_shape.body is not self._space.static_body:  # not the molecule and not a wall/edge of reaction vessel
                        near_shape_position = near_point['point']
                        magnitude = other_shape.charge * shape.charge * self._dipole_force_constant

                        if magnitude > ChargeReactionVessel.magnitude_cutoff:
                            minimum_distance = self._bodies[near_point['shape'].body].get_size() + mol.get_size()
                            force_vector = ChargeReactionVessel._get_force_vector(shape_position, near_shape_position, magnitude, minimum_distance)

                            body.apply_impulse(force_vector, shape_position)  # apply force at an offset (world-coordinates)
                            body.force_vector.append((shape_position, force_vector))

        super(ChargeReactionVessel, self).step()

    @classmethod
    def _get_force_vector(cls, p1, p2, magnitude, minimum_distance=0):
        """Return force vector on p1 as a result of force with given magnitude from p2 - +magnitude = repulsion - according to Coulomb's Law"""

        distance_squared = max(minimum_distance * minimum_distance, p1.get_dist_sqrd(p2))
        direction = (p2 - p1).normalized()  # normalize as just for direction - length determined by distance_squared

        force_vector = -magnitude * direction / distance_squared

        return force_vector
