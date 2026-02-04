import unittest

from src.results import run_simulation_full


class TestResultsPipeline(unittest.TestCase):
    def test_run_simulation_full_returns_expected_structure(self):
        results_text, results_data = run_simulation_full(
            diameter=50,
            density=3000,
            velocity_km_s=20,
            entry_angle=45,
            r_distance=10,
            lat=None,
            lon=None,
        )

        self.assertIsInstance(results_text, str)
        self.assertIsInstance(results_data, dict)
        self.assertIn("input_parameters", results_data)
        self.assertIn("atmospheric_entry", results_data)
        self.assertIn("energy", results_data)
        self.assertIn("vulnerability_analysis", results_data)

        event_type = results_data["atmospheric_entry"]["event_type"]
        self.assertIn(event_type, {"airburst", "ground impact", "intact"})

        zones = results_data["vulnerability_analysis"]["zones"]
        self.assertTrue(all(zone["start_distance"] <= zone["end_distance"] for zone in zones))


if __name__ == "__main__":
    unittest.main()
