"""
Created on 10/01/2013

@author: thom
"""
import unittest

from atoms.population import Population


class PopulationTest(unittest.TestCase):

    def setUp(self):
        self.p = Population()
        self.p.set_t(0)
        self.p.set_quantity("CC(=O)O", 100)
        self.p.set_quantity("CN", 50)
        self.p.set_quantity("O=C(O)CCOC(=O)O", 25)
        self.p.set_quantity("C=CO", 100)
        self.p.set_quantity("C=CC(=C)N", 50)

        self.q = Population()
        self.q.set_t(0)
        self.q.set_quantity('AA', 10)
        self.q.set_quantity('BB', 10)
        self.q.set_quantity('CC', 10)
        self.q.set_quantity('D', 10)
        self.q.set_quantity('E', 10)
        self.q.set_t(1)
        self.q.set_quantity('AA', 9)
        self.q.set_quantity('BB', 9)
        self.q.set_quantity('F', 1)
        self.q.set_t(2)
        self.q.set_quantity('CC', 9)
        self.q.set_quantity('D', 9)
        self.q.set_quantity('G', 1)
        self.q.set_t(3)
        self.q.set_quantity('AA', 8)
        self.q.set_quantity('F', 0)
        self.q.set_quantity('H', 1)

    def tearDown(self):
        del self.p
        del self.q

    def test_init(self):
        self.assertEqual(5, len(self.p.get_items()))
        self.assertEqual(1, len(self.p.get_times()))

        self.assertEqual(8, len(self.q.get_items()))
        self.assertEqual(4, len(self.q.get_times()))

        q = Population(population=self.q)
        self.assertEqual(8, len(q.get_items()))
        self.assertEqual(4, len(q.get_times()))

    def test_str(self):
        self.assertIsInstance(self.p.__str__(), str)

    def test_expand_population(self):
        old_shape = self.q._population.shape
        self.q._expand_population((0, 0))
        self.assertEqual(old_shape, self.q._population.shape)
        new_shape = self.q._population.shape[0] + 100, self.q._population.shape[1] + 200
        self.q._expand_population((100, 200))
        self.assertEqual(new_shape, self.q._population.shape)

    def test_get_times(self):
        self.assertEqual([0, 1, 2, 3], self.q.get_times())

    def test_get_slice_by_items(self):
        smiles = "O=C(O)CCOC(=O)O"
        actual = self.p.get_slice_by_items([smiles])
        self.assertEqual(25, actual.get_quantity(smiles))  # Test that we have a Population, and that value is correct

        actual = self.p.get_slice_by_items([smiles, "CN"])
        self.assertEqual(25, actual.get_quantity(smiles))
        self.assertEqual(50, actual.get_quantity("CN"))
        self.assertEqual(50, self.p.get_quantity("C=CC(=C)N"))
        self.assertEqual(0, actual.get_quantity("C=CC(=C)N"))  # not in slice, so 0
        self.assertEqual(2, actual._population.shape[1])

    def test_get_slice_by_time(self):
        actual = self.p.get_slice_by_time([0])
        self.assertEqual(1, actual._population.shape[0])
        self.assertEqual(325, actual.get_population_size())  # Test that we have a Population, and that value is correct

        self.assertEqual(10, self.q.get_slice_by_time([0]).get_quantity("AA"))
        self.assertEqual(50, self.q.get_slice_by_time([0]).get_population_size())  # Test that we have a Population, and that value is correct
        self.assertEqual(48, self.q.get_slice_by_time([2]).get_population_size())  # Test that we have a Population, and that value is correct
        self.assertEqual(47, self.q.get_slice_by_time([3]).get_population_size())  # Test that we have a Population, and that value is correct
        actual = self.q.get_slice_by_time([1, 3])
        self.assertEqual(2, actual._population.shape[0])
        self.assertEqual([1, 3], actual._t)
        self.assertEqual(8, actual.get_quantity("AA"))
        self.assertEqual(1, actual.get_quantity("H"))
        self.assertEqual(1, actual.get_quantity("G"))

        # At slice t=3:
        #  <setvalue item='AA' quantity='8' />
        #  <setvalue item='BB' quantity='9' />
        #  <setvalue item='CC' quantity='9' />
        #  <setvalue item='D' quantity='9' />
        #  <setvalue item='E' quantity='10' />
        #  <setvalue item='F' quantity='0' />
        #  <setvalue item='G' quantity='1' />
        #  <setvalue item='H' quantity='1' />
        actual = self.q.get_slice_by_time([3])
        self.assertEqual(1, actual._population.shape[0])
        self.assertEqual(8, actual.get_quantity("AA"))
        self.assertEqual(0, actual.get_quantity("F"))
        self.assertEqual(1, actual.get_quantity("G"))
        self.assertEqual(actual.get_population_size(), self.q.get_population_size())
        self.assertEqual(actual.get_population_size(), len(actual.get_population()))
        self.assertEqual(actual.get_population(), self.q.get_population())

    def test_get_last_slice(self):
        actual = self.q.get_last_slice()
        self.assertEqual(1, actual._population.shape[0])
        self.assertEqual(8, actual.get_quantity("AA"))
        self.assertEqual(0, actual.get_quantity("F"))
        self.assertEqual(1, actual.get_quantity("G"))

    def test_get_changing_items(self):
        actual = self.q.get_changing_items()
        self.assertEqual(7, len(actual))

    def test_get_unchanging_items(self):
        actual = self.q.get_unchanging_items()
        self.assertEqual(1, len(actual))

    def test_get_population_size(self):
        self.assertEqual(325, self.p.get_population_size())

    def test_population_set_quantity(self):
        smarts = 'O=C(O)CCOC(=O)O'
        self.assertEqual(25, self.p.get_quantity(smarts))

        self.p.set_quantity(smarts, 15)
        self.assertEqual(15, self.p.get_quantity(smarts))

        self.p.set_quantity(smarts, 0)
        self.assertEqual(0, self.p.get_quantity(smarts))

        self.p.set_quantity('CCC', 25)  # new molecule
        self.assertEqual(25, self.p.get_quantity('CCC'))

    def test_set_t(self):
        smarts = 'O=C(O)CCOC(=O)O'
        self.assertEqual(25, self.p.get_quantity(smarts))
        self.assertEqual(1, len(self.p.get_times()))
        self.p.set_t(1.25)
        self.assertEqual(25, self.p.get_quantity(smarts))  # quantities carry through to new t
        self.assertEqual(2, len(self.p.get_times()))
        self.p.set_quantity(smarts, 20)
        self.assertEqual(20, self.p.get_quantity(smarts))
        self.assertEqual(25, self.p.get_slice_by_time([0]).get_quantity(smarts))

    def test_get_population(self):
        p = Population()
        p.set_quantity('a', 5)
        p.set_quantity('b', 2)
        self.assertEqual(['a', 'a', 'a', 'a', 'a', 'b', 'b'], p.get_population())
        self.assertEqual(p.get_population_size(), len(p.get_population()))
        self.assertEqual(self.q.get_population_size(), len(self.q.get_population()))

    def test_regression_initial_population(self):
        actual = self.q.get_slice_by_time([0])
        self.assertEqual(8, len(actual.get_items()))

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
