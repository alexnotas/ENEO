import math
import unittest

from src import utils


class TestUtils(unittest.TestCase):
    def test_km_m_roundtrip(self):
        km = 12.345
        meters = utils.km_to_m(km)
        self.assertAlmostEqual(utils.m_to_km(meters), km, places=9)

    def test_convert_energy_j_to_mt(self):
        one_mt = 4.184e15
        self.assertAlmostEqual(utils.convert_energy_j_to_mt(one_mt), 1.0, places=9)

    def test_curvature_adjustment_factor_zero_distance(self):
        factor = utils.curvature_adjustment_factor(0.0, 1000.0)
        self.assertAlmostEqual(factor, 1.0, places=9)

    def test_compute_scaling_factor_zero(self):
        self.assertEqual(utils.compute_scaling_factor(0.0), 1.0)


if __name__ == '__main__':
    unittest.main()
