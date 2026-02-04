import unittest
from unittest import mock

from src import population_calculator
from src import gdp_calculator


class TestPopulationCalculator(unittest.TestCase):
    def test_population_in_zones_without_raster(self):
        vulnerability_zones = [
            {"threshold": 1.0, "start_distance": 0.0, "end_distance": 1.0}
        ]

        with mock.patch.object(population_calculator, "RASTER", None):
            result = population_calculator.calculate_population_in_zones(
                lat=0.0,
                lon=0.0,
                vulnerability_zones=vulnerability_zones,
            )

        self.assertEqual(result["total_casualties"], 0)
        self.assertEqual(result["note"], "Error: Population data not available")
        self.assertIn("Population raster data is not loaded", result["warnings"])


class TestGDPCalculator(unittest.TestCase):
    def test_normalize_country_name_mappings(self):
        self.assertEqual(
            gdp_calculator.normalize_country_name("United States"),
            "united states of america",
        )
        self.assertEqual(
            gdp_calculator.normalize_country_name("Czechia"),
            "czech republic",
        )

    def test_calculate_economic_damage_without_data(self):
        countries = [{"name": "Testland", "fid": 1, "total_casualties": 10}]
        with mock.patch.object(gdp_calculator, "load_gdp_data", return_value=None), \
             mock.patch.object(gdp_calculator, "merged_df", None):
            result = gdp_calculator.calculate_economic_damage(countries)

        self.assertEqual(result["total_economic_damage"], 0)
        self.assertEqual(result["countries_without_data"], 1)
        self.assertIn("GDP data could not be loaded", result["note"])


if __name__ == "__main__":
    unittest.main()
