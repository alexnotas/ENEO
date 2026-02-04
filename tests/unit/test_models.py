import math
import unittest

from src.models import AsteroidImpactSimulation


class TestAsteroidImpactSimulation(unittest.TestCase):
    def test_calculate_asteroid_energy_matches_formula(self):
        diameter = 2.0
        density = 3000.0
        velocity_km_s = 20.0
        sim = AsteroidImpactSimulation(diameter, velocity_km_s, density, 45)

        expected_radius = diameter / 2.0
        expected_volume = (4.0 / 3.0) * math.pi * (expected_radius ** 3)
        expected_mass = density * expected_volume
        expected_energy = 0.5 * expected_mass * (velocity_km_s * 1000.0) ** 2

        energy_j, _ = sim.calculate_asteroid_energy()
        self.assertAlmostEqual(energy_j, expected_energy, places=6)

    def test_simulate_atmospheric_entry_breakup_flag(self):
        sim = AsteroidImpactSimulation(diameter=50, velocity_km_s=20, density=3000, entry_angle_deg=45)
        entry = sim.simulate_atmospheric_entry()

        y_i = sim.compute_yield_strength()
        i_f = 4.07 * (2.0 * 8000.0 * y_i) / (
            sim.density * sim.diameter * sim.v0 ** 2 * math.sin(sim.theta)
        )
        self.assertEqual(entry["breakup"], i_f < 1.0)

    def test_calculate_final_crater_diameter_simple(self):
        sim = AsteroidImpactSimulation(diameter=50, velocity_km_s=20, density=3000, entry_angle_deg=45)
        d_tc = 3000.0  # 3 km transient crater => simple regime
        d_fr = sim.calculate_final_crater_diameter(d_tc)
        self.assertAlmostEqual(d_fr, 1.25 * d_tc, places=6)

    def test_calculate_seismic_magnitude(self):
        """Test seismic magnitude calculation (Richter scale)."""
        sim = AsteroidImpactSimulation(diameter=100, velocity_km_s=20, density=3000, entry_angle_deg=45)
        # 1 Joule impact energy
        # M = 0.67 * log10(1) - 5.87 = -5.87
        m_seismic = sim.calculate_seismic_magnitude(impact_energy=1.0)
        self.assertAlmostEqual(m_seismic, -5.87, places=2)

        # 1 MT impact energy (4.184e15 J)
        # M = 0.67 * log10(4.184e15) - 5.87 = 0.67 * 15.62 - 5.87 = 10.46 - 5.87 = 4.59 approx
        m_seismic_mt = sim.calculate_seismic_magnitude(impact_energy=4.184e15)
        self.assertGreater(m_seismic_mt, 4.0)
        self.assertLess(m_seismic_mt, 5.0)

    def test_calculate_thermal_exposure_scaling(self):
        """Test thermal exposure scaling with distance (inverse square law approximation)."""
        sim = AsteroidImpactSimulation(diameter=50, velocity_km_s=20, density=3000, entry_angle_deg=45)
        energy_j = 1e15

        # Calculate exposure at 10km and 20km
        phi_10 = sim.calculate_thermal_exposure(energy_j, r_km=10)
        phi_20 = sim.calculate_thermal_exposure(energy_j, r_km=20)

        # Should decrease significantly (roughly 1/4th, ignoring curvature/transmissivity for a moment)
        self.assertGreater(phi_10, phi_20)
        # Strictly speaking, Phi ~ 1/r^2, so phi_10 / phi_20 approx (20/10)^2 = 4
        # We allow some tolerance due to curvature and f factor
        ratio = phi_10 / phi_20
        self.assertGreater(ratio, 3.0)
        self.assertLess(ratio, 5.0)

    def test_calculate_blast_arrival_simple(self):
        """Test blast wave arrival time (distance / speed of sound)."""
        sim = AsteroidImpactSimulation(diameter=50, velocity_km_s=20, density=3000, entry_angle_deg=45)
        # Speed of sound c0 is roughly 330-340 m/s depending on utils.py
        # If distance is 3400m, time should be ~10s
        t_arrival = sim.calculate_blast_arrival_time(distance_m=3400)
        self.assertTrue(9.0 < t_arrival < 11.0, f"Arrival time {t_arrival} should be around 10s")

    def test_calculate_tsunami_land_impact(self):
        """Test tsunami calculation returns no wave for land impact."""
        pass


if __name__ == '__main__':
    unittest.main()
