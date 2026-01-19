import unittest
import sys
import os
import math

# Add project root to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.models import AsteroidImpactSimulation
from src.utils import convert_energy_j_to_mt

class TestPhysicsScenarios(unittest.TestCase):
    
    def test_scenario_a_airburst(self):
        print("\n=== Scenario A (Airburst) Results ===")
        print("Parameters: D=50m, v=22km/s, rho=3100kg/m3, angle=55, dist=2km")
        
        diameter = 50
        velocity = 22
        density = 3100
        angle = 55
        distance_km = 2
        
        sim = AsteroidImpactSimulation(diameter, velocity, density, angle)
        
        # 1. Initial KE
        ke_j, ke_mt = sim.calculate_asteroid_energy()
        print(f"Initial kinetic energy: {ke_j:.2e} J")
        
        # Run Sim
        res = sim.simulate_atmospheric_entry()
        
        # 2. Break-up Altitude
        if res['z_star']:
            print(f"Break-up altitude: {res['z_star']:.0f} m")
        else:
            print("Break-up altitude: None")
        
        # 3. Airburst Altitude
        if res['airburst_altitude']:
            print(f"Airburst altitude: {res['airburst_altitude']:.0f} m")
        else:
            print("Airburst altitude: None")
        
        # 4. Speed after burst
        print(f"Speed after burst: {res['post_breakup_velocity']/1000:.2f} km/s")
        
        # 5. Overpressure & Wind
        # Calculate Airburst Energy components
        mass = density * (4.0/3.0) * math.pi * ((diameter/2)**3)
        ke_initial = 0.5 * mass * (sim.v0**2)
        ke_post = 0.5 * mass * (res['post_breakup_velocity']**2)
        ke_internal = ke_initial - ke_post
        
        # Main energy used for airburst intensity
        burst_energy_main = max(ke_post, ke_internal)
        
        # Standard ENEO calculation
        op_main = sim.calculate_overpressure_airburst(
            distance_km * 1000, 
            res['airburst_altitude'], 
            burst_energy_main, 
            res['z_star']
        )
        print(f"Overpressure: {op_main:.2f} Pa")
        
        # Check alternative energy source for range explanation if needed
        op_alt = sim.calculate_overpressure_airburst(
            distance_km * 1000, 
            res['airburst_altitude'], 
            ke_post, 
            res['z_star']
        )
        if abs(op_main - op_alt) > 100:
             print(f"(Alternative Overpressure calc: {op_alt:.2f} Pa)")
        
        wind = sim.calculate_peak_wind_velocity(op_main)
        print(f"Wind speed: {wind:.2f} m/s")

    def test_scenario_b_ground(self):
        print("\n=== Scenario B (Ground Impact) Results ===")
        print("Parameters: D=250m, v=27km/s, rho=3100kg/m3, angle=35, dist=10km")
        
        diameter = 250
        velocity = 27
        density = 3100
        angle = 35
        distance_km = 10
        
        sim = AsteroidImpactSimulation(diameter, velocity, density, angle)
        
        # 1. Initial KE
        ke_j, ke_mt = sim.calculate_asteroid_energy()
        print(f"Initial kinetic energy: {ke_j:.2e} J")
        
        # Run Sim
        res = sim.simulate_atmospheric_entry()
        print(f"Break-up altitude: {res['z_star']:.0f} m")
        print(f"Impact speed: {res['post_breakup_velocity']/1000:.2f} km/s")
        
        # Impact Energy
        imp_e, _ = sim.calculate_impact_energy(res['post_breakup_velocity'])
        print(f"Energy released on impact: {imp_e:.2e} J")
        
        # Crater
        dtc = sim.calculate_transient_crater_diameter(res['post_breakup_velocity'], 2600, angle)
        print(f"Transient crater diameter: {dtc/1000:.2f} km")
        
        dfr = sim.calculate_final_crater_diameter(dtc)
        print(f"Final crater diameter: {dfr/1000:.2f} km")
        
        # Thermal
        therm = sim.calculate_thermal_exposure(imp_e, distance_km)
        print(f"Thermal fluence: {therm:.2e} J/m2")
        
        # Seismic
        mag = sim.calculate_seismic_magnitude(imp_e)
        print(f"Seismic magnitude (Richter): {mag:.2f}")
        
        # Ejecta
        ejecta = sim.calculate_ejecta_thickness(dtc, distance_km)
        print(f"Ejecta thickness: {ejecta:.2f} m")
        
        # Overpressure
        op = sim.calculate_overpressure_ground_new(distance_km * 1000, imp_e)
        print(f"Overpressure: {op:.2f} Pa")
        
        # Wind
        wind = sim.calculate_peak_wind_velocity(op)
        print(f"Wind speed: {wind:.2f} m/s")

if __name__ == '__main__':
    unittest.main()
