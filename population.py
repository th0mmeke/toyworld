"""
Created on 30/01/2013

@author: thom
"""

import copy
import logging
import collections

import numpy as np


class Population(object):

    """A representation of a generic population of items, where each unique item has an associated quantity.

    Optionally, copy a population object from the provided Population.

    Population is stored as a Numpy Array - _population[t,i] where t=iteration time, and i = element index (t goes downwards, elements across)

    :param population: Copy constructor - the existing population to copy from
    """

    def __init__(self, population=None):

        if population is not None:
            self._population = np.array(population._population, copy=True)
            self._index = copy.deepcopy(population._index)
            self._t = copy.deepcopy(population._t)

        else:
            self._population = np.zeros((100, 10), dtype=np.dtype(int))
            self._index = {}  # dictionary of item:row in _population
            self._t = []

    def __str__(self):
        result = {}
        for item in self.get_items():
            result[item] = self.get_quantity(item)
        return str(result)

    def set_t(self, t):
        """Add on new time-slice (a copy of the last slice)
        :param t: timestamp of new time-slice (must be greater than last timestamp)"""

        if len(self._t) > 0 and t <= self._t[-1]:
            return
        try:
            self._population[len(self._t), :] = self._population[len(self._t) - 1, :]
        except:
            # logging.info("Expanding population in set_t")
            self._expand_population((1000, 0))
            self._population[len(self._t), :] = self._population[len(self._t) - 1, :]
        self._t.append(t)

    def get_items(self):
        return self._index.keys()

    def get_times(self):
        """Return the full list of timestamps for this population"""
        return self._t

    def get_population_size(self):
        """Return total size of population at most recent time

        :rtype: int
        """

        return self._population[len(self._t) - 1, :].sum()

    def get_slice_by_items(self, items):
        """
        Slice the population by item, returning the quantities at every time for each item

        :param items: Items to slice by
        :type items: List of string
        :rtype: Population
        """

        if len(self._index) == 0:
            raise ValueError('No slices available')

        slice = type(self)()
        slice_index = {item: self._index[item] for item in items}
        # order the slice index by row to ensure relationship between item and position is consistent
        ordered_slice_index = collections.OrderedDict(sorted(slice_index.iteritems(), key=lambda t: t[1]))

        slice._population = np.copy(self._population[:, ordered_slice_index.values()])  # index by the row positions of the items in the slice - row order preserved
        # rebase the slice index to match the sliced population
        row = 0
        slice._index = {}
        for item, old_row in ordered_slice_index.iteritems():
            slice._index[item] = row
            row += 1
        slice._t = copy.copy(self._t)

        return slice
        # return Population(population_tuple=(self._population[:, slice_index], [self._index[i] for i in slice_index], self._t))

    def get_slice_by_time(self, times):
        """Slice the population by time, returning the quantities of all items at the times given.

        :param times: Times to slice by
        :type times: List of float - times are the recorded times, not the ordinal index of the time
        :rtype: Population"""

        if len(self._t) == 0:
            raise ValueError('No slices available')

        slice = type(self)()
        slice._t = times
        array_indices = [self._t.index(i) for i in times]
        slice._population = np.copy(self._population[array_indices, :])
        slice._index = copy.copy(self._index)
        return slice

    def get_last_slice(self):
        logging.debug("Getting last slice (at {})".format(self._t[-1]))
        return self.get_slice_by_time([self._t[-1]])

    def get_changing_items(self):
        """Return a set of all items that have changed quantity."""

        index = np.all(np.equal(self._population[0:len(self._t), :], self._population[0, :]), 0)
        return set([i for i in self._index if not index[self._index[i]]])

    def get_unchanging_items(self):
        """Return a set of all items that have NOT changed quantity."""

        index = np.all(np.equal(self._population[0:len(self._t), :], self._population[0, :]), 0)
        return set([i for i in self._index if index[self._index[i]]])

    def get_quantity(self, item):
        """Just tell us the quantity of this item...if it's not in the population, return 0 (rather than None).

        :param item: Item in question
        :type item: string
        """

        try:
            return self._population[len(self._t) - 1, self._index[item]]
        except:
            return 0

    def set_quantity(self, item, quantity):
        """Update the quantity of an item, while handling the special case of setting the quantity for a new item.

        :param item: Item to update
        :type item: string
        :param quantity: New quantity
        :type quantity: int
        """

#         if quantity < 0:
#             raise Exception

        if item not in self._index:
            self._index[item] = len(self._index)  # set the position of this new item equal to the end row in the population matrix

        try:
            self._population[len(self._t) - 1, self._index[item]] = quantity
        except:
            # Check if we need to upsize our population matrix...if we do, up it by a good whack as upsizing is expensive
            logging.debug("set_quantity(): expanding population size by (0,100)")
            self._expand_population((0, 100))
            self._population[len(self._t) - 1, self._index[item]] = quantity

    def get_population(self):
        """Return a list of all of the items in the population. An item of quantity n will appear n times in the return list.

        :rtype: list of string"""

        result = []
        for item in self.get_items():
            result.extend([item for i in range(self.get_quantity(item))])
        return result

    def _expand_population(self, shape_increment):
        """Increase the size of the population.

        :param shape_increment: number of additional iterations, number of additional population items
        :type shape_increment: tuple (int,int)"""

        if len(shape_increment) != 2:
            raise ValueError
        logging.debug("_expand_population(): expanding population from {} by {}".format(self._population.shape, shape_increment))
        self._population = np.vstack((self._population, np.zeros((shape_increment[0], self._population.shape[1]), dtype=np.dtype(int))))
        self._population = np.hstack((self._population, np.zeros((self._population.shape[0], shape_increment[1]), dtype=np.dtype(int))))
