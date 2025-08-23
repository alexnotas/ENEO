import math
import numpy as np # Ensure numpy is imported if not already at the top
from models import AsteroidImpactSimulation
from utils import (
    km_to_m, m_to_km, convert_energy_j_to_mt,
    g_E_atmos, C_D, rho0, rho_target, get_ocean_depth_from_geotiff, WATER_DENSITY_CONSTANT, ETOPO_FILE_PATH
)
from vulnerability_models import (
    fun_CraterVulnerability, fun_SeismicVulnerability, fun_OverpressureVulnerability,
    fun_ThermRadVulnerability, fun_HighWindVulnerability, fun_EjectaBlanketVulnerability
)
from thresholds import (
    SELECTED_THERMAL_THRESHOLDS, SEISMIC_THRESHOLDS, EJECTA_THRESHOLDS, BLAST_THRESHOLDS,
    EF_WIND_THRESHOLDS, TSUNAMI_AMPLITUDE_THRESHOLDS,
    get_wind_damage_category, get_ejecta_damage_category, get_thermal_damage_category # Added new imports
)
import math # Ensure math is imported for airburst energy calculation and tsunami calcs

def collect_simulation_results(sim: AsteroidImpactSimulation, entry_results,
                               initial_energy_joules, initial_energy_megatons,
                               specific_energy_joules, specific_energy_type,
                               r_distance, tsunami_data=None): # Added tsunami_data
    """
    Aggregates all simulation results into a structured dictionary.

    This function takes the simulation object, atmospheric entry results, energy calculations,
    and other data points to create a comprehensive output dictionary that summarizes
    the entire impact event.

    Args:
        sim (AsteroidImpactSimulation): The simulation object containing asteroid parameters.
        entry_results (dict): A dictionary with the results of the atmospheric entry simulation.
        initial_energy_joules (float): The initial kinetic energy of the asteroid in Joules.
        initial_energy_megatons (float): The initial kinetic energy of the asteroid in Megatons of TNT.
        specific_energy_joules (float): The event-specific energy (impact or airburst) in Joules.
        specific_energy_type (str): A string indicating the type of energy ('Impact', 'Airburst', 'Initial').
        r_distance (float): The distance from the impact/airburst point for which some effects are calculated.
        tsunami_data (dict, optional): Tsunami simulation results, if applicable. Defaults to None.

    Returns:
        dict: A dictionary containing the structured results of the simulation.
    """
    # Initialize the results dictionary with the core information.
    results = {
        # Store the initial input parameters for reference and debugging.
        'input_parameters': {
            'diameter': sim.diameter,
            'density': sim.density,
            'velocity': sim.velocity_km_s,
            'entry_angle': sim.entry_angle_deg,
            'distance': r_distance
        },
        # Store the detailed results from the atmospheric entry phase.
        'atmospheric_entry': {
            'breakup': entry_results['breakup'],
            'breakup_altitude': entry_results.get('z_star'),
            'airburst_altitude': entry_results.get('airburst_altitude'),
            'velocity_at_breakup': entry_results.get('v_breakup'),
            'final_velocity': entry_results['post_breakup_velocity'],
            'dispersion_length': entry_results.get('dispersion_length'),
            'event_type': entry_results['event_type']
        },
        # Store the energy calculations, including initial and event-specific values.
        'energy': {
            'initial_energy_joules': initial_energy_joules,
            'initial_energy_megatons': initial_energy_megatons,
            'specific_energy_joules': specific_energy_joules,
            'specific_energy_megatons': convert_energy_j_to_mt(specific_energy_joules),
            'specific_energy_type': specific_energy_type # e.g., "Impact", "Airburst", "Initial"
        }
    }
    # If the event is a "ground impact", calculate and add crater formation details.
    if entry_results["event_type"] == "ground impact":
        # Retrieve the surface impact velocity and dispersion length from entry results.
        v_surface = entry_results['post_breakup_velocity']
        l_dispersion = entry_results.get('dispersion_length', 0)
        # Calculate the transient crater diameter, which is the initial crater size before collapse.
        D_tc = sim.calculate_transient_crater_diameter(v_surface, rho_target, sim.entry_angle_deg)
        # Calculate the final crater diameter after gravitational collapse and modification.
        D_fr = sim.calculate_final_crater_diameter(D_tc)
        # Calculate the final crater depth.
        depth = sim.calculate_crater_depth(D_tc)
        # Structure the crater results.
        crater_results = {
            "crater_formation": {
                "transient_diameter": D_tc,
                "final_diameter": D_fr,
                "depth": depth
            }
        }
        # If the dispersion length of fragments is greater than the transient crater diameter,
        # it's likely that a field of smaller craters would form rather than a single large one.
        if l_dispersion >= D_tc:
            crater_results["crater_formation"]["note"] = "Crater field likely due to large dispersion." # Example note
        # Add the crater results to the main results dictionary.
        results['crater'] = crater_results
    
    # Add Tsunami results to the main dictionary if they have been calculated and passed.
    if tsunami_data:
        results['tsunami'] = tsunami_data

    # The vulnerability analysis calculates danger zones based on a combined vulnerability model.
    # It defines zones where the population is expected to experience a certain level of harm.
    # Define the vulnerability thresholds to analyze, from 100% down to 5%.
    thresholds_to_apply_to_population = [round(x, 2) for x in np.arange(1.0, 0.04, -0.05)] # e.g., [1.0, 0.95, 0.90, ..., 0.05]
    vulnerability_zones = []
    current_start_distance = 0.0 # This tracks the start distance for the current vulnerability band.
    
    # The lowest threshold we are applying to the population
    lowest_applied_threshold = min(thresholds_to_apply_to_population) if thresholds_to_apply_to_population else 0.0

    # Iterate through the defined thresholds to create concentric danger zones.
    for applied_threshold_for_zone in thresholds_to_apply_to_population:
        # Determine the actual vulnerability level that defines the outer boundary of this zone.
        actual_vuln_level_defining_this_band_outer_edge = applied_threshold_for_zone
        
        # Apply the 'F - 0.025' logic for all zones EXCEPT the lowest one.
        # This logic helps define the boundary between vulnerability bands.
        # The lowest zone (e.g., 5%) will use its own value as the boundary.
        if applied_threshold_for_zone > lowest_applied_threshold:
            actual_vuln_level_defining_this_band_outer_edge = applied_threshold_for_zone - 0.025
        # Ensure the calculated edge is not below a very small positive number,
        # especially if applied_threshold_for_zone - 0.025 results in zero or negative.
        # However, find_vulnerability_distance should handle r_min correctly.
        # We also need to ensure that actual_vuln_level_defining_this_band_outer_edge
        # does not become less than the next lower applied_threshold_for_zone.
        # For instance, for 0.10 zone, edge is 0.075. For 0.05 zone, edge is 0.05.
        # This logic seems fine as is.

        # Ensure the target vulnerability level is not below a practical minimum (e.g., 0.001)
        # This is important if applied_threshold_for_zone - 0.025 could result in a very low or negative number.
        # Given the lowest applied_threshold_for_zone is 0.05, and it uses 0.05 directly,
        # the next lowest (0.10) would use 0.10 - 0.025 = 0.075, which is fine.
        actual_vuln_level_defining_this_band_outer_edge = max(0.001, actual_vuln_level_defining_this_band_outer_edge)

        # Find the maximum distance from the impact point where the vulnerability
        # is at least at the level of the calculated outer edge.
        end_distance_for_this_band = find_vulnerability_distance(sim, actual_vuln_level_defining_this_band_outer_edge, entry_results)
        
        # Skip creating a zone if it's negligibly small and starts at the impact point.
        # This avoids creating tiny, meaningless zones like "0.0 km to 0.01 km".
        if abs(current_start_distance) < 1e-9 and end_distance_for_this_band <= 0.01:
            # Do not add this zone. current_start_distance remains 0.0, so the next
            # significant zone will correctly start from 0.0.
            continue

        # If a valid, non-overlapping zone is found, add it to the list.
        if end_distance_for_this_band > current_start_distance:
            vulnerability_zones.append({
                "threshold": applied_threshold_for_zone, # This is the factor applied to population
                "start_distance": current_start_distance,
                "end_distance": end_distance_for_this_band
            })
            # The end of the current zone becomes the start of the next one.
            current_start_distance = end_distance_for_this_band # Update for the start of the next band
    
    # Add the calculated vulnerability zones to the main results dictionary.
    results["vulnerability_analysis"] = {
        "zones": vulnerability_zones
    }
    
    # Return the complete, structured results.
    return results


def run_simulation_full(diameter, density, velocity_km_s, entry_angle, r_distance, lat=None, lon=None): # Added lat, lon
    """
    Runs a full asteroid impact simulation from start to finish.

    This function orchestrates the entire simulation process, including:
    1. Initializing the asteroid simulation object.
    2. Calculating the initial kinetic energy.
    3. Simulating the atmospheric entry and determining the event type (ground impact vs. airburst).
    4. Calculating event-specific energy (impact or airburst energy).
    5. Checking for transitional events (where small changes in input cause large changes in output).
    6. Calculating detailed effects for cratering, thermal radiation, seismic activity, ejecta, and airblast.
    7. Collecting all results into a structured format.

    Args:
        diameter (float): Diameter of the asteroid in meters.
        density (float): Density of the asteroid in kg/m^3.
        velocity_km_s (float): Velocity of the asteroid in km/s.
        entry_angle (float): Entry angle of the asteroid in degrees from horizontal.
        r_distance (float): The specific distance from the impact point (in km) for which to calculate effects.
        lat (float, optional): Latitude of the impact point. Used for tsunami calculations. Defaults to None.
        lon (float, optional): Longitude of the impact point. Used for tsunami calculations. Defaults to None.

    Returns:
        tuple: A tuple containing:
            - dict: A structured dictionary of all simulation results.
            - str: A formatted string summarizing the impact energy.
            - str: A formatted string summarizing the atmospheric entry results.
            - str: A formatted string summarizing crater and melt details.
            - str: A formatted string summarizing thermal radiation effects.
            - str: A formatted string summarizing seismic effects.
            - str: A formatted string summarizing ejecta effects.
            - str: A formatted string summarizing airblast effects.
            - str: A formatted string summarizing tsunami effects (if applicable).
    """
    # Step 1: Initialize the simulation object with the asteroid's physical properties.
    sim = AsteroidImpactSimulation(diameter, velocity_km_s, density, entry_angle)
    # Step 2: Calculate the initial kinetic energy before atmospheric entry.
    initial_energy_joules, initial_energy_megatons = sim.calculate_asteroid_energy()
    impact_text = f"Kinetic energy (before entry): {initial_energy_joules:.2e} J\nEquivalent: {initial_energy_megatons:.2f} MT\n\n"
    # Step 3: Simulate the asteroid's passage through the atmosphere.
    entry_results = sim.simulate_atmospheric_entry()
    atm_text = "Atmospheric Entry Results:\n"
    # Based on whether the asteroid breaks up, format the summary text.
    if not entry_results["breakup"]:
        atm_text += "Asteroid remains intact during entry.\n"
        atm_text += f"I_f = {entry_results['I_f']:.3f}\n"
        atm_text += f"Final impact velocity ≈ {m_to_km(entry_results['post_breakup_velocity']):.2f} km/s\n"
    else:
        atm_text += "Asteroid breaks up in the atmosphere.\n"
        atm_text += f"I_f = {entry_results['I_f']:.3f}\n"
        atm_text += f"Breakup altitude (z*): {m_to_km(entry_results['z_star']):.2f} km\n"
        atm_text += f"Velocity at breakup: {m_to_km(entry_results['v_breakup']):.2f} km/s\n"
        atm_text += f"Dispersion length: {m_to_km(entry_results['dispersion_length']):.2f} km\n"
        # Further detail depends on whether it's an airburst or ground impact.
        if entry_results["event_type"] == "airburst":
            atm_text += f"Airburst altitude (z_b): {m_to_km(entry_results['airburst_altitude']):.2f} km\n"
            atm_text += f"Residual velocity at airburst: {m_to_km(entry_results['post_breakup_velocity']):.2f} km/s\n"
        else:
            atm_text += "Fragments reach ground (crater-forming event).\n"
            atm_text += f"Surface impact velocity: {m_to_km(entry_results['post_breakup_velocity']):.2f} km/s\n"
            atm_text += f"Pancake factor at ground: {entry_results['pancake_factor_ground']:.2f}\n"
    # Calculate the percentage of velocity lost during atmospheric entry.
    perc_reduction = 100.0 * (sim.v0 - entry_results["post_breakup_velocity"]) / sim.v0
    atm_text += f"\nVelocity reduction: {perc_reduction:.2f}%\n"
    # v_swarm = entry_results["post_breakup_velocity"] # Not directly used for energy type decision here
    # v_terminal = math.sqrt((density * (sim.diameter/2) * g_E_atmos) / (3 * C_D * rho0)) # Not directly used
    # This is the key velocity for determining subsequent effects (impact energy, etc.).
    v_impact_final = entry_results["post_breakup_velocity"] # This is the key velocity for impact/airburst

    # Step 4: Determine the specific energy for the event (impact or airburst).
    # This energy value will be used for calculating the magnitude of various effects.
    specific_energy_joules_for_event = initial_energy_joules # Default to initial energy
    specific_energy_type_for_event = "Initial" # Default type

    # If it's a ground impact (either by fragments or an intact object), calculate impact energy.
    if entry_results["event_type"] == "ground impact" or entry_results["event_type"] == "intact": # MODIFIED to include intact
        imp_energy_joules, imp_energy_MT = sim.calculate_impact_energy(v_impact_final) # Use v_impact_final
        if entry_results["event_type"] == "ground impact":
            impact_text += f"Impact Energy (Fragments): {imp_energy_joules:.2e} J ({imp_energy_MT:.2f} MT)\n\n"
        else: # intact
            impact_text += f"Impact Energy (Intact Object): {imp_energy_joules:.2e} J ({imp_energy_MT:.2f} MT)\n\n"
        # This becomes the primary energy for subsequent calculations.
        specific_energy_joules_for_event = imp_energy_joules
        specific_energy_type_for_event = "Impact"
    # If it's an airburst, calculate the energy released in the air.
    elif entry_results["event_type"] == "airburst":
        # Calculate airburst energy (consistent with logic in airblast_text section)
        mass = sim.density * (4.0/3.0) * math.pi * ((sim.diameter/2)**3)
        # Use sim.v0 (initial velocity in m/s) for KE_initial_for_airburst
        KE_initial_for_airburst = 0.5 * mass * (sim.v0)**2
        # entry_results["post_breakup_velocity"] is residual velocity at airburst altitude
        KE_post_for_airburst = 0.5 * mass * (entry_results["post_breakup_velocity"])**2
        KE_internal_for_airburst = KE_initial_for_airburst - KE_post_for_airburst
        # The airburst energy is the greater of the remaining kinetic energy or the energy lost.
        airburst_energy_joules = max(KE_post_for_airburst, KE_internal_for_airburst)
        airburst_energy_MT = convert_energy_j_to_mt(airburst_energy_joules)
        impact_text += f"Airburst Energy: {airburst_energy_joules:.2e} J ({airburst_energy_MT:.2f} MT)\n\n"
        # This becomes the primary energy for subsequent calculations.
        specific_energy_joules_for_event = airburst_energy_joules
        specific_energy_type_for_event = "Airburst"
    else: # "intact" or other cases
        # Use v_impact_final for intact cases as well
        impact_text += f"Estimated surface impact velocity: {m_to_km(v_impact_final):.2f} km/s\n"
        impact_text += "No specific impact or airburst energy computed for this event type.\n\n"
        # specific_energy_joules_for_event remains initial_energy_joules
        # specific_energy_type_for_event remains "Initial"

    # Step 5: Check for "transitional" events. These are cases where a small change
    # in the asteroid's diameter could lead to a different event type (e.g., ground impact to airburst)
    # and potentially more severe overpressure effects.
    transitional = False
    if entry_results["event_type"] == "ground impact":
        # Use v_impact_final here
        imp_energy, _ = sim.calculate_impact_energy(v_impact_final)
        candidate_distances = [1, 5, 10]
        candidate_deltas = [5, 10, 20, 30, 40]
        # Test what happens if the diameter were slightly smaller.
        for delta in candidate_deltas:
            if sim.diameter - delta <= 0:
                continue
            temp_sim = AsteroidImpactSimulation(sim.diameter - delta, sim.velocity_km_s, sim.density, sim.entry_angle_deg)
            temp_entry = temp_sim.simulate_atmospheric_entry()
            for d in candidate_distances:
                D_test = km_to_m(d)
                # Calculate overpressure for the original asteroid.
                p_original_test = sim.calculate_overpressure_ground_new(D_test, imp_energy)
                # Calculate overpressure for the slightly smaller asteroid.
                if temp_entry["event_type"] == "ground impact":
                    v_terminal_temp = math.sqrt((density * (temp_sim.diameter/2) * g_E_atmos) / (3 * C_D * rho0))
                    v_swarm_temp = temp_entry["post_breakup_velocity"]
                    v_surface_temp = max(v_swarm_temp, v_terminal_temp)
                    candidate_imp_energy, _ = temp_sim.calculate_impact_energy(v_surface_temp)
                    p_candidate = temp_sim.calculate_overpressure_ground_new(D_test, candidate_imp_energy)
                elif temp_entry["event_type"] == "airburst":
                    mass = temp_sim.density * (4.0/3.0) * math.pi * ((temp_sim.diameter/2)**3)
                    KE_initial = 0.5 * mass * (temp_sim.v0**2)
                    KE_post = 0.5 * mass * (temp_entry["post_breakup_velocity"]**2)
                    KE_internal = KE_initial - KE_post
                    airburst_energy = max(KE_post, KE_internal)
                    p_candidate = temp_sim.calculate_overpressure_airburst(D_test, temp_entry.get("airburst_altitude", 0), airburst_energy, temp_entry["z_star"])
                else:
                    continue
                # If the smaller asteroid produces a higher overpressure, it's a transitional event.
                if p_candidate > p_original_test:
                    transitional = True
                    break
            if transitional:
                break
    if transitional:
        impact_text += "Note: Transition region detected due to sensitivity to small diameter changes.\n\n"
    # Step 6: Calculate detailed effects based on the event type.
    # Crater & Melt (only for ground impact events)
    crater_text = ""
    D_fr = 0 # Initialize D_fr for potential use in ejecta if crater doesn't form but event type implies it
    D_tc = 0 # Initialize D_tc
    if entry_results["event_type"] == "ground impact":
        # Use v_impact_final and sim.entry_angle_deg
        D_tc = sim.calculate_transient_crater_diameter(v_impact_final, rho_target, sim.entry_angle_deg)
        D_tc_km = m_to_km(D_tc)
        crater_type = "simple" if D_tc_km <= 3.2 else "complex"
        D_fr = sim.calculate_final_crater_diameter(D_tc)
        final_depth = sim.calculate_crater_depth(D_tc)
        transient_depth = D_tc / (2 * math.sqrt(2))
        V_br = sim.calculate_breccia_volume(D_fr)
        # Use v_impact_final and sim.entry_angle_deg
        imp_energy, _ = sim.calculate_impact_energy(v_impact_final)
        V_m = sim.calculate_melt_volume(imp_energy, sim.entry_angle_deg)
        t_m = sim.calculate_melt_sheet_thickness(V_m, D_tc)
        crater_text += "Crater & Melt Details:\n"
        crater_text += f"Crater type: {crater_type.capitalize()}\n"
        crater_text += f"Transient crater diameter: {D_tc_km:.2f} km\n"
        crater_text += f"Final crater diameter: {m_to_km(D_fr):.2f} km\n"
        crater_text += f"Transient crater depth: {transient_depth:.2f} m\n"
        crater_text += f"Final crater depth: {final_depth:.2f} m\n"
        crater_text += f"Breccia volume: {V_br:.2e} m³\n"
        crater_text += f"Impact melt volume: {V_m:.2e} m³\n"
        crater_text += f"Melt sheet thickness: {t_m:.2e} m\n"
    else:
        crater_text += "Airburst event: no crater formed.\n"
    # Thermal Radiation
    thermal_text = ""
    if entry_results["event_type"] == "ground impact":
        # Use v_impact_final
        # A fireball and thermal radiation are only produced if the impact velocity is high enough.
        if m_to_km(v_impact_final) < 15.0:
            thermal_text += "No fireball is created, therefore, there is no thermal radiation damage.\n"
            phi_ground = 0.0
        else:
            imp_energy, imp_energy_MT = sim.calculate_impact_energy(v_impact_final)
            # Calculate thermal exposure at the specified distance.
            phi_ground = sim.calculate_thermal_exposure(imp_energy, r_distance)
            thermal_text += "Thermal Radiation (Ground Impact):\n"
            thermal_text += f"Calculated thermal exposure at {r_distance:.2f} km: {phi_ground:.2e} J/m²\n"
            thermal_cat_ground = get_thermal_damage_category(phi_ground)
            thermal_text += f"Thermal Effect Category at {r_distance:.2f} km: {thermal_cat_ground}\n"
            
            # Add fireball characteristics for ground impact
            R_f_ground = sim.calculate_fireball_radius(imp_energy)
            T_t_ground = sim.calculate_time_of_max_radiation(R_f_ground, v_impact_final)
            tau_t_ground = sim.calculate_irradiation_duration(imp_energy, R_f_ground)
            
            thermal_text += f"Fireball radius: {R_f_ground:.2f} m\n"
            thermal_text += f"Time of maximum radiation: {T_t_ground:.2f} s\n"
            thermal_text += f"Irradiation duration: {tau_t_ground:.2f} s\n"
            
            # Calculate danger zones for different levels of thermal damage.
            previous_bound = 0.0
            thermal_zones = []
            
            for desc, phi_1Mt in SELECTED_THERMAL_THRESHOLDS:
                # Scale the standard threshold to the energy of this specific impact.
                scaled_threshold = sim.calculate_ignition_exposure(imp_energy_MT, phi_1Mt * 1e6)
                
                # Find the maximum distance at which this threshold is met or exceeded.
                # This is done using a binary search for efficiency.
                max_dist = 0.0
                left, right = 0.01, 1000000.0
                while right - left > 0.01:
                    mid = (left + right) / 2
                    current_phi = sim.calculate_thermal_exposure(imp_energy, mid)
                    if current_phi >= scaled_threshold:
                        max_dist = mid
                        left = mid
                    else:
                        right = mid
                
                if max_dist > previous_bound:
                    thermal_zones.append((desc, previous_bound, max_dist))
                    previous_bound = max_dist
        
            thermal_text += "\nThermal Danger Zones:\n"
            for desc, lower, upper in thermal_zones:
                thermal_text += f"{desc}: {lower:.2f} km - {upper:.2f} km\n"
            
    else: # Airburst
        v_airburst = entry_results["v_breakup"] # Velocity at breakup for airburst
        z_b = entry_results["airburst_altitude"]
        D_m = km_to_m(r_distance)

        # Similar to ground impact, a fireball requires a minimum velocity.
        if m_to_km(v_airburst) < 15.0:
            thermal_text += "No fireball is created, therefore, there is no thermal radiation damage.\n"
            phi = 0.0
        else:
            # --- Airburst thermal radiation (unchanged) ---
            mass = sim.density * (4.0/3.0) * math.pi * ((sim.diameter/2)**3)
            KE_initial_for_airburst_thermal = 0.5 * mass * (sim.v0)**2
            KE_post_for_airburst_thermal = 0.5 * mass * (entry_results["post_breakup_velocity"])**2
            KE_internal_for_airburst_thermal = KE_initial_for_airburst_thermal - KE_post_for_airburst_thermal
            airburst_energy_joules_for_thermal = max(KE_post_for_airburst_thermal, KE_internal_for_airburst_thermal)

            # Calculate thermal flux, considering both distance and airburst altitude.
            phi = sim.calculate_airburst_thermal_flux(airburst_energy_joules_for_thermal, z_b, D_m)
            thermal_text += "Thermal Radiation (Airburst):\n"
            thermal_text += f"Calculated thermal flux density at {r_distance:.2f} km: {phi:.2e} J/m²\n"
            thermal_cat_airburst = get_thermal_damage_category(phi)
            thermal_text += f"Thermal Effect Category at {r_distance:.2f} km: {thermal_cat_airburst}\n"
            
            # Airburst calculations use airburst_energy_joules_for_thermal and breakup velocity (v_airburst)
            R_f_airburst = sim.calculate_fireball_radius(airburst_energy_joules_for_thermal)
            T_t_airburst = sim.calculate_time_of_max_radiation(R_f_airburst, v_airburst)
            tau_t_airburst = sim.calculate_irradiation_duration(airburst_energy_joules_for_thermal, R_f_airburst)
            
            thermal_text += f"Fireball radius: {R_f_airburst:.2f} m\n"
            thermal_text += f"Time of maximum radiation: {T_t_airburst:.2f} s\n"
            thermal_text += f"Irradiation duration: {tau_t_airburst:.2f} s\n"
            
            E_MT = convert_energy_j_to_mt(airburst_energy_joules_for_thermal)
            # Helper function to find the max distance for a given thermal threshold in an airburst.
            def find_airburst_thermal_max_distance(threshold):
                low = 0.01
                high = 1e8
                tol = 0.01
                best = low
                # Binary search for the distance.
                for _ in range(100):
                    mid = (low + high) / 2.0
                    # Critical fix: Convert km to m before passing to flux calculation
                    current_phi = sim.calculate_airburst_thermal_flux(airburst_energy_joules_for_thermal, z_b, km_to_m(mid))
                    if current_phi >= threshold:
                        best = mid
                        low = mid
                    else:
                        high = mid
                # Also need km_to_m here in the final check
                return best if sim.calculate_airburst_thermal_flux(airburst_energy_joules_for_thermal, z_b, km_to_m(best)) >= threshold else 0.0

            previous_bound = 0.0
            thermal_zones = []
            
            for desc, phi_1Mt in SELECTED_THERMAL_THRESHOLDS:
                scaled_threshold = sim.calculate_ignition_exposure(E_MT, phi_1Mt * 1e6)
                max_dist = find_airburst_thermal_max_distance(scaled_threshold)
                if max_dist > previous_bound:
                    # Remove the m_to_km conversion - distances are already in km
                    thermal_zones.append((desc, previous_bound, max_dist))
                    previous_bound = max_dist
            thermal_text += "\nThermal Danger Zones:\n"
            for desc, lower, upper in thermal_zones:
                thermal_text += f"{desc}: {lower:.2f} km - {upper:.2f} km\n"
    # Seismic Effects (only for ground impact events)
    seismic_text = ""
    if entry_results["event_type"] == "ground impact":
        # Use v_impact_final
        imp_energy, _ = sim.calculate_impact_energy(v_impact_final)
        # Calculate the Richter magnitude of the seismic event.
        M = sim.calculate_seismic_magnitude(imp_energy)
        # Calculate the effective magnitude at the specified distance, accounting for attenuation.
        M_eff = sim.calculate_effective_seismic_magnitude(M, r_distance)
        T_s = sim.calculate_seismic_arrival_time(r_distance)
        seismic_text += "Seismic Effects:\n"
        seismic_text += f"Seismic magnitude: {M:.2f}\n"
        seismic_text += f"Effective magnitude at {r_distance:.2f} km: {M_eff:.2f}\n"
        seismic_text += f"Modified Mercalli Intensity at {r_distance:.2f} km: {sim.map_magnitude_to_mmi(M_eff)}\n"
        seismic_text += f"Seismic arrival time at {r_distance:.2f} km: {T_s:.2f} s\n"
        seismic_text += f"\nSeismic Danger Zones:\n"
        # Helper function to find the max distance for a given seismic magnitude threshold.
        def find_seismic_max_distance(threshold):
            low = 0.01
            high = 20000
            tol = 0.01
            best = low
            # Binary search for the distance.
            while high - low > tol:
                mid = (low + high) / 2
                current_meff = sim.calculate_effective_seismic_magnitude(M, mid)
                if current_meff >= threshold:
                    best = mid
                    low = mid
                else:
                    high = mid
            return best
        previous_max = 0.0
        for desc, thresh in SEISMIC_THRESHOLDS:
            current_max = find_seismic_max_distance(thresh)
            if current_max == 0.01 and sim.calculate_effective_seismic_magnitude(M, 0.01) < thresh:
                zone_str = "None"
            elif current_max <= previous_max:
                zone_str = "None"
            else:
                zone_str = f"{previous_max:.2f}-{current_max:.2f} km"
                previous_max = current_max

            seismic_text += f"{desc}: {zone_str}\n"
    else:
        seismic_text += "Airburst event: no seismic effects.\n"
    # Ejecta (only for ground impact events)
    ejecta_text = ""
    if entry_results["event_type"] == "ground impact":
        # Calculate the thickness of the ejecta blanket at the specified distance.
        t_e = sim.calculate_ejecta_thickness(D_tc, r_distance)
        # Calculate the mean diameter of ejecta fragments at that distance.
        d_mean = sim.calculate_mean_fragment_diameter(D_fr, r_distance)
        T_ejecta = sim.calculate_ejecta_arrival_time(r_distance)
        ejecta_text += "Ejecta Deposit:\n"
        ejecta_text += f"Ejecta thickness at {r_distance:.2f} km: {t_e:.4f} m\n"
        ejecta_cat = get_ejecta_damage_category(t_e)
        ejecta_text += f"Ejecta Category at {r_distance:.2f} km: {ejecta_cat}\n"
        ejecta_text += f"Mean fragment diameter at {r_distance:.2f} km: {d_mean:.2f} m\n"
        ejecta_text += f"Ejecta arrival time at {r_distance:.2f} km: {T_ejecta:.2f} s\n" if T_ejecta is not None else f"Ejecta arrival time at {r_distance:.2f} km: N/A\n"
        ejecta_text += f"\nEjecta Danger Zones:\n"
        # Helper function to find the max distance for a given ejecta thickness threshold.
        def find_ejecta_max_distance(threshold):
            low = 0.01
            high = 20000
            tol = 0.01
            best = low
            # Binary search for the distance.
            while high - low > tol:
                mid = (low + high) / 2
                current_thickness = sim.calculate_ejecta_thickness(D_tc, mid)
                if current_thickness >= threshold:
                    best = mid
                    low = mid
                else:
                    high = mid
            return best
        previous_max = 0.0
        for desc, thresh in EJECTA_THRESHOLDS:
            current_max = find_ejecta_max_distance(thresh)
            if current_max == 0.01 and sim.calculate_ejecta_thickness(D_tc, 0.01) < thresh:
                zone_str = "None"
            elif current_max <= previous_max:
                zone_str = "None"
            else:
                zone_str = f"{previous_max:.2f}-{current_max:.2f} km"
                previous_max = current_max
            ejecta_text += f"{desc}: {zone_str}\n"

    else:
        ejecta_text += "Airburst event: no ejecta computed.\n"
    # Airblast Effects (calculated for both ground impact and airburst)
    airblast_text = ""
    D_m = km_to_m(r_distance)
    if entry_results["event_type"] == "ground impact":
        # Use v_impact_final
        imp_energy, _ = sim.calculate_impact_energy(v_impact_final)
        # Calculate overpressure. If it's a transitional event, special calculations may apply.
        if transitional:
            # Note: The logic for transitional overpressure seems to be missing here,
            # it falls through to the standard calculation.
            p_overpressure = sim.calculate_overpressure_ground_new(D_m, imp_energy)
        else:
            p_overpressure = sim.calculate_overpressure_ground_new(D_m, imp_energy)
        wind_velocity = sim.calculate_peak_wind_velocity(p_overpressure)
        damage_zone = sim.calculate_damage_category(p_overpressure)
        airblast_text += "Ground Impact Air Blast:\n"
        airblast_text += f"Overpressure: {p_overpressure:.2f} Pa\n"
        airblast_text += f"Damage Category: {damage_zone}\n"
        I, intensity_db = sim.calculate_sound_intensity(p_overpressure, wind_velocity)
        airblast_text += f"Sound Intensity: {I:.2e} J/m²\n"
        airblast_text += f"SPL: {intensity_db:.2f} dB\n"
    elif entry_results["event_type"] == "airburst":
        # Calculate airburst energy.
        mass = density * (4.0/3.0) * math.pi * ((sim.diameter/2)**3)
        KE_initial = 0.5 * mass * (sim.v0**2)
        KE_post = 0.5 * mass * (entry_results["post_breakup_velocity"]**2)
        KE_internal = KE_initial - KE_post
        airburst_energy = max(KE_post, KE_internal)
        # Calculate overpressure from the airburst.
        p_overpressure = sim.calculate_overpressure_airburst(D_m, entry_results.get("airburst_altitude", 0), airburst_energy, entry_results["z_star"])
        wind_velocity = sim.calculate_peak_wind_velocity(p_overpressure)
        damage_zone = sim.calculate_damage_category(p_overpressure)
        airblast_text += "Airburst Blast Wave:\n"
        airblast_text += f"Overpressure: {p_overpressure:.2f} Pa\n"
        airblast_text += f"Damage Category: {damage_zone}\n"
        I, intensity_db = sim.calculate_sound_intensity(p_overpressure, wind_velocity)
        airblast_text += f"Sound Intensity: {I:.2e} J/m²\n"
        airblast_text += f"SPL: {intensity_db:.2f} dB\n"

    # Calculate the arrival time of the blast wave.
    T_b = sim.calculate_blast_arrival_time(D_m, burst_altitude_m=entry_results.get("airburst_altitude"))
    airblast_text += f"\nBlast arrival time: {T_b:.2f} s\n"
    
    # Add blast zone calculations to airblast text
    airblast_text += "\nBlast Danger Zones:\n"
    # Helper function to find the max distance for a given overpressure threshold.
    def find_blast_max_distance(threshold):
        low = 0.01
        high = 1e8 # Search up to a very large distance
        # Binary search for the distance.
        for _ in range(100): # 100 iterations for precision
            mid = (low + high) / 2.0
            mid_m = km_to_m(mid)
            # Calculate overpressure based on event type.
            if entry_results["event_type"] == "ground impact":
                imp_energy, _ = sim.calculate_impact_energy(v_impact_final)
                current_p = sim.calculate_overpressure_ground_new(mid_m, imp_energy)
            elif entry_results["event_type"] == "airburst":
                mass = density * (4.0/3.0) * math.pi * ((sim.diameter/2)**3)
                KE_initial = 0.5 * mass * (sim.v0**2)
                KE_post = 0.5 * mass * (entry_results["post_breakup_velocity"]**2)
                KE_internal = KE_initial - KE_post
                airburst_energy = max(KE_post, KE_internal)
                current_p = sim.calculate_overpressure_airburst(mid_m, entry_results.get("airburst_altitude", 0), airburst_energy, entry_results["z_star"])
            else:
                current_p = 0 # No overpressure for other event types.

            if current_p >= threshold:
                low = mid
            else:
                high = mid
        return low # Return the max distance found

    previous_max = 0.0
    for desc, thresh in BLAST_THRESHOLDS:
        current_max = find_blast_max_distance(thresh)

        zone_text_to_add = ""

        if current_max == 0.0: 
            zone_text_to_add = f"{desc}: None\n"
            # previous_max is not updated
        elif current_max <= previous_max: 
            zone_text_to_add = f"{desc}: None\n"
            # previous_max is not updated
        else: # current_max > previous_max AND current_max > 0.0
            zone_text_to_add = f"{desc}: {previous_max:.2f}-{current_max:.2f} km\n"
            
            # If this zone is "0.00-0.01 km", don't update previous_max to 0.01.
            # This ensures the next significant zone starts from 0.00.
            if not (previous_max == 0.0 and current_max == 0.01):
                previous_max = current_max
            # else: previous_max remains 0.0 (or its prior value if not the very first zone)
        
        airblast_text += zone_text_to_add

    # --- NEW: Wind Effects (EF Scale) Section ---
    wind_effects_content = "" # Changed variable name to avoid conflict with section header text
    
    # Add Peak Wind Velocity calculated at r_distance (from Airblast section)
    # This assumes 'wind_velocity' variable from the Airblast section is correctly scoped and holds the value at r_distance.
    if 'wind_velocity' in locals() or 'wind_velocity' in globals():
        wind_effects_content += f"Peak wind velocity at {r_distance:.2f} km: {wind_velocity:.2f} m/s\n"
        wind_cat = get_wind_damage_category(wind_velocity)
        wind_effects_content += f"EF Scale Category at {r_distance:.2f} km: {wind_cat}\n"
    else:
        wind_effects_content += f"Peak wind velocity at {r_distance:.2f} km: Not calculated (likely no significant overpressure).\n"
        wind_effects_content += f"EF Scale Category at {r_distance:.2f} km: Below EF0\n"

    # Helper function to find max distance for a given wind speed threshold
    def find_wind_ef_scale_max_distance(wind_speed_threshold_mps):
        low = 0.01  # km
        high = 20000 # km
        tol = 0.01   # km, Tolerance for binary search
        best_distance_km = 0.0

        energy_for_calc = 0
        if entry_results["event_type"] == "airburst":
            # Recalculate airburst_energy consistently (or use 'airburst_energy' variable if defined in outer scope)
            mass = sim.density * (4.0/3.0) * math.pi * ((sim.diameter/2)**3)
            KE_initial = 0.5 * mass * (sim.v0**2)
            KE_post = 0.5 * mass * (entry_results["post_breakup_velocity"]**2)
            KE_internal = KE_initial - KE_post
            energy_for_calc = max(KE_post, KE_internal)
        else: # ground impact
            # This relies on 'imp_energy' being available from the scope of run_simulation_full,
            # which should be the case if the Airblast section (where it's defined for ground impacts) precedes this.
            energy_for_calc = imp_energy 

        # Binary search for distance
        for _ in range(100): # Max iterations
            mid_km = (low + high) / 2.0
            if mid_km == 0: # Avoid issues with zero distance
                current_wind_mps = 0
            else:
                current_pressure_pa = 0
                if entry_results["event_type"] == "airburst":
                    current_pressure_pa = sim.calculate_overpressure_airburst(
                        km_to_m(mid_km),
                        entry_results["airburst_altitude"],
                        energy_for_calc, # Use the determined energy
                        entry_results["z_star"]
                    )
                else: # ground impact
                    current_pressure_pa = sim.calculate_overpressure_ground_new(
                        km_to_m(mid_km),
                        energy_for_calc # Use the determined energy
                    )
                current_wind_mps = sim.calculate_peak_wind_velocity(current_pressure_pa)

            if current_wind_mps >= wind_speed_threshold_mps:
                best_distance_km = mid_km
                low = mid_km
            else:
                high = mid_km
            
            if high - low < tol:
                break
        
        # Final check for the 'best_distance_km'
        if best_distance_km > 0:
            final_pressure_pa = 0
            if entry_results["event_type"] == "airburst":
                final_pressure_pa = sim.calculate_overpressure_airburst(km_to_m(best_distance_km), entry_results["airburst_altitude"], energy_for_calc, entry_results["z_star"])
            else: # ground impact
                final_pressure_pa = sim.calculate_overpressure_ground_new(km_to_m(best_distance_km), energy_for_calc)
            final_wind_mps = sim.calculate_peak_wind_velocity(final_pressure_pa)
            
            if final_wind_mps >= wind_speed_threshold_mps:
                return best_distance_km
            else: 
                # Check if wind speed at 0.01km (min distance) meets threshold
                pressure_at_min_dist = 0
                if entry_results["event_type"] == "airburst":
                    pressure_at_min_dist = sim.calculate_overpressure_airburst(km_to_m(0.01), entry_results["airburst_altitude"], energy_for_calc, entry_results["z_star"])
                else: # ground impact
                    pressure_at_min_dist = sim.calculate_overpressure_ground_new(km_to_m(0.01), energy_for_calc)
                wind_at_min_dist = sim.calculate_peak_wind_velocity(pressure_at_min_dist)
                if wind_at_min_dist >= wind_speed_threshold_mps:
                    return 0.01 
                return 0.0 
        return 0.0 # If best_distance_km remained 0 or initial checks failed

    # Sort EF_WIND_THRESHOLDS by severity (descending min_mps) to process strongest winds first
    sorted_ef_thresholds_desc = sorted(EF_WIND_THRESHOLDS, key=lambda x: x[1][0], reverse=True)
    
    previous_ef_max_km = 0.0
    
    for ef_desc, (min_mps, max_mps) in sorted_ef_thresholds_desc:
        current_ef_max_km = find_wind_ef_scale_max_distance(min_mps)
        
        zone_text_to_add = ""
        
        # Format the wind speed range string
        wind_range_str = f"{min_mps}-{max_mps} m/s"
        if max_mps == float('inf'):
            wind_range_str = f">{min_mps} m/s"

        if current_ef_max_km == 0.0: # Threshold not met at all
            zone_text_to_add = f"{ef_desc} ({wind_range_str}): None\n"
        elif current_ef_max_km <= previous_ef_max_km: # This EF scale's extent is within or equal to the previous (more severe) one
            zone_text_to_add = f"{ef_desc} ({wind_range_str}): None\n"
        else: # A new, wider band for this EF scale
            zone_text_to_add = f"{ef_desc} ({wind_range_str}): {previous_ef_max_km:.2f}-{current_ef_max_km:.2f} km\n"
            # Update previous_ef_max_km, ensuring not to set it to a tiny value if it's the first small zone
            # This logic is generally fine: if the first zone is 0.00-0.01, previous_ef_max_km becomes 0.01.
            previous_ef_max_km = current_ef_max_km
        
        wind_effects_content += zone_text_to_add
        
    # --- Tsunami Effects ---
    tsunami_text = ""
    tsunami_results_data = None # To store data for collect_simulation_results
    if (entry_results["event_type"] == "ground impact" or entry_results["event_type"] == "intact") and lat is not None and lon is not None:
        v_surface_for_tsunami = entry_results['post_breakup_velocity'] # Velocity at surface
        
        # Call the tsunami calculation method from the simulation object
        tsunami_calc = sim.calculate_tsunami_effects(v_surface_for_tsunami, lat, lon)
        tsunami_results_data = tsunami_calc # Store for results_data

        if tsunami_calc["error_message"]:
            tsunami_text += f"Tsunami calculation error: {tsunami_calc['error_message']}\n"
        elif tsunami_calc["is_on_land"]:
            tsunami_text += "Impact on land: No tsunami generated.\n"
            if tsunami_calc.get('ocean_depth') is not None: # ocean_depth might be 0 if on land
                 tsunami_text += f"Note: Ocean depth at location reported as {tsunami_calc['ocean_depth']:.2f} m.\n"
        else:
            tsunami_text += f"Ocean depth at impact: {tsunami_calc['ocean_depth']:.2f} m\n"
            tsunami_text += f"Transient cavity diameter in water: {tsunami_calc['transient_cavity_diameter_water']:.2f} m\n"
            tsunami_text += f"Max amplitude at source: {tsunami_calc['max_amplitude_at_source']:.2f} m\n"

            if tsunami_calc['max_amplitude_at_source'] > 0 and tsunami_calc['transient_cavity_diameter_water'] > 0:
                amplitude_at_r_distance = sim.calculate_tsunami_amplitude_at_distance(
                    tsunami_calc['max_amplitude_at_source'],
                    tsunami_calc['transient_cavity_diameter_water'],
                    km_to_m(r_distance) # r_distance is the user input distance in km
                )
                tsunami_text += f"Amplitude at {r_distance:.2f} km: {amplitude_at_r_distance:.2f} m\n"
                if tsunami_results_data: # Should always be true here as it's tsunami_calc
                    tsunami_results_data['amplitude_at_r_distance'] = amplitude_at_r_distance

            tsunami_text += "\nTsunami Danger Zones:\n"
            
            def find_tsunami_max_distance(amplitude_threshold_m):
                if tsunami_calc['max_amplitude_at_source'] <= amplitude_threshold_m:
                    return 0.0
                
                low_km = 0.01
                high_km = 20000.0 
                tol_km = 0.01
                best_distance_km = 0.0

                for _ in range(100):
                    if high_km - low_km < tol_km:
                        break
                    mid_km = (low_km + high_km) / 2
                    if mid_km == 0: # Avoid division by zero if d_tc_water is 0 or mid_km is 0
                        current_amplitude_at_mid = tsunami_calc['max_amplitude_at_source']
                    else:
                        current_amplitude_at_mid = sim.calculate_tsunami_amplitude_at_distance(
                            tsunami_calc['max_amplitude_at_source'],
                            tsunami_calc['transient_cavity_diameter_water'],
                            km_to_m(mid_km)
                        )
                    if current_amplitude_at_mid >= amplitude_threshold_m:
                        best_distance_km = mid_km
                        low_km = mid_km
                    else:
                        high_km = mid_km
                
                # Final check
                if best_distance_km > 0:
                    final_amp_at_best = sim.calculate_tsunami_amplitude_at_distance(
                        tsunami_calc['max_amplitude_at_source'],
                        tsunami_calc['transient_cavity_diameter_water'],
                        km_to_m(best_distance_km)
                    )
                    if final_amp_at_best < amplitude_threshold_m:
                        # If it undershoots significantly, it might mean the threshold is never truly met beyond source effects
                        # For simplicity, we trust the binary search's best_distance_km if it's > 0
                        pass
                return best_distance_km

            sorted_tsunami_thresholds = sorted(TSUNAMI_AMPLITUDE_THRESHOLDS, key=lambda x: x[1], reverse=True)
            
            if 'zones' not in tsunami_results_data:
                tsunami_results_data['zones'] = []

            previous_max_km_for_tsunami = 0.0 # Initialize previous max distance

            for desc, amp_m_threshold in sorted_tsunami_thresholds:
                current_zone_outer_boundary_km = find_tsunami_max_distance(amp_m_threshold)
                zone_text_to_add = ""

                # Check if the current zone extends beyond the previous one and is significant
                if current_zone_outer_boundary_km > previous_max_km_for_tsunami and current_zone_outer_boundary_km > 0.01:
                    zone_start_km = previous_max_km_for_tsunami
                    zone_end_km = current_zone_outer_boundary_km
                    
                    zone_text_to_add = f"{desc}: {zone_start_km:.2f}-{zone_end_km:.2f} km\n"
                    
                    tsunami_results_data['zones'].append({
                        'description': desc,
                        'amplitude_threshold_m': amp_m_threshold,
                        'start_distance_km': zone_start_km,
                        'end_distance_km': zone_end_km, # This is the outer radius of the band
                    })
                    previous_max_km_for_tsunami = zone_end_km # Update for the next band
                else:
                    # This threshold does not create a new, wider band or is insignificant
                    zone_text_to_add = f"{desc}: None\n"
                    # Do not add to tsunami_results_data['zones'] as it's not a distinct new band for visualization data
                
                tsunami_text += zone_text_to_add
            # DEBUG: Print the tsunami_results_data for verification
            # print(f"DEBUG Tsunami Results Data: {tsunami_results_data}")
    elif lat is None or lon is None:
        tsunami_text += "Tsunami calculation skipped: Impact location (latitude/longitude) not provided.\n"
    else: # Airburst or other non-surface impact event
        tsunami_text += "Tsunami calculation not applicable for airburst events.\n"


    # --- Vulnerability Models Section for results_text ---
    vuln_model_text_lines = ["=== Vulnerability Models ==="]

    if entry_results["event_type"] == "ground impact":
        # Ensure p_overpressure, wind_velocity, phi_ground, M_eff, t_e are the values at r_distance
        # These should be available from their respective calculation blocks above.

        v_pressure = sim.calculate_overpressure_vulnerability(p_overpressure) # p_overpressure at r_distance
        vuln_model_text_lines.append(f"Overpressure Vulnerability: {v_pressure:.4f}")

        v_wind = sim.calculate_wind_vulnerability(wind_velocity) # wind_velocity at r_distance
        vuln_model_text_lines.append(f"Wind Vulnerability: {v_wind:.4f}")

        v_thermal = sim.calculate_thermal_vulnerability(phi_ground) # phi_ground at r_distance
        vuln_model_text_lines.append(f"Thermal Vulnerability: {v_thermal:.4f}")

        v_seismic = sim.calculate_seismic_vulnerability(M_eff) # M_eff at r_distance
        vuln_model_text_lines.append(f"Seismic Vulnerability: {v_seismic:.4f}")

        # Ensure D_tc is calculated if not already available for t_e calculation context
        # D_tc was calculated in the "Crater & Melt" section
        v_ejecta = sim.calculate_ejecta_vulnerability(t_e) # t_e at r_distance
        vuln_model_text_lines.append(f"Ejecta Vulnerability: {v_ejecta:.4f}")
        
        combined_vuln = 1.0 - (
            (1.0 - v_pressure) * (1.0 - v_wind) * (1.0 - v_thermal) *
            (1.0 - v_seismic) * (1.0 - v_ejecta)
        )
        vuln_model_text_lines.append(f"Combined Vulnerability: {combined_vuln:.4f}")

    elif entry_results["event_type"] == "airburst":
        # Ensure p_overpressure, wind_velocity, phi are the values at r_distance

        v_pressure = sim.calculate_overpressure_vulnerability(p_overpressure) # p_overpressure at r_distance
        vuln_model_text_lines.append(f"Overpressure Vulnerability: {v_pressure:.4f}")

        v_wind = sim.calculate_wind_vulnerability(wind_velocity) # wind_velocity at r_distance
        vuln_model_text_lines.append(f"Wind Vulnerability: {v_wind:.4f}")

        v_thermal = sim.calculate_thermal_vulnerability(phi) # phi at r_distance (for airburst)
        vuln_model_text_lines.append(f"Thermal Vulnerability: {v_thermal:.4f}")

        vuln_model_text_lines.append("Seismic Vulnerability: N/A")
        vuln_model_text_lines.append("Ejecta Vulnerability: N/A")
        
        combined_vuln = 1.0 - (
            (1.0 - v_pressure) * (1.0 - v_wind) * (1.0 - v_thermal)
        )
        vuln_model_text_lines.append(f"Combined Vulnerability: {combined_vuln:.4f}")
    else: # E.g., "intact" event type if it doesn't fit ground impact or airburst logic for these effects
        vuln_model_text_lines.append("Overpressure Vulnerability: N/A")
        vuln_model_text_lines.append("Wind Vulnerability: N/A")
        vuln_model_text_lines.append("Thermal Vulnerability: N/A")
        vuln_model_text_lines.append("Seismic Vulnerability: N/A")
        vuln_model_text_lines.append("Ejecta Vulnerability: N/A")
        vuln_model_text_lines.append("Combined Vulnerability: N/A")

    vuln_text_for_output = "\n" + "\n".join(vuln_model_text_lines) + "\n"
    # Replace the old placeholder vuln_text (which was just "\n")
    full_results = ("\n=== Atmospheric Entry ===\n" + atm_text +
                    "\n=== Impact & Energy ===\n" + impact_text +
                    "\n=== Crater & Melt ===\n" + crater_text +
                    "\n=== Thermal Radiation ===\n" + thermal_text +
                    "\n=== Seismic Effects ===\n" + seismic_text +
                    "\n=== Ejecta ===\n" + ejecta_text +
                    "\n=== Airblast ===\n" + airblast_text +
                    "\n=== Wind Effects ===\n" + wind_effects_content +
                    "\n=== Tsunami Effects ===\n" + tsunami_text + # Added Tsunami
                    vuln_text_for_output) # USE THE GENERATED vuln_text_for_output HERE
    
    # Add a debug print statement here to check the generated wind_effects_content
    # print(f"DEBUG Wind Effects Content:\n{wind_effects_content}")
    # print(f"DEBUG Vulnerability Models Content:\n{vuln_text_for_output}") # DEBUG LINE


    return full_results, collect_simulation_results(sim, entry_results, initial_energy_joules, initial_energy_megatons, specific_energy_joules_for_event, specific_energy_type_for_event, r_distance, tsunami_results_data)

def find_vulnerability_distance(sim, threshold, entry_results, r_min=0.01, r_max=20000.0, tol=0.01):
    """
    Determines the maximum distance at which the combined vulnerability
    equals or exceeds a specified threshold using a point calculation.

    This function employs a binary search to efficiently find the distance.
    It calculates the combined vulnerability from various effects (overpressure,
    wind, thermal, seismic, ejecta) at different distances, considering whether
    the event is a ground impact or an airburst. This version calculates
    vulnerability at a single point for a given distance, not averaged over a ring.

    Args:
        sim (AsteroidImpactSimulation): The simulation object.
        threshold (float): The target vulnerability threshold (0.0 to 1.0).
        entry_results (dict): Results from the atmospheric entry simulation.
        r_min (float, optional): Minimum search distance in kilometers. Defaults to 0.01.
        r_max (float, optional): Maximum search distance in kilometers. Defaults to 20000.0.
        tol (float, optional): Tolerance for the binary search convergence in kilometers. Defaults to 0.01.

    Returns:
        float: The maximum distance (in km) where combined vulnerability >= threshold.
               Returns 0.0 if the threshold is not met even at the closest distances
               or if the calculated distance is effectively zero.
    """
    def calculate_vulnerability_at_distance(r_km):
        """Helper to calculate combined vulnerability at a specific distance r_km."""
        if r_km <= 0: return 1.0 # Assume max vulnerability at or before impact point
        
        D_m_local = km_to_m(r_km)
        point_vulnerability = 0.0 # Combined vulnerability at this distance
        
        if entry_results["event_type"] == "ground impact":
            v_surface_local = entry_results['post_breakup_velocity']
            imp_energy_local, _ = sim.calculate_impact_energy(v_surface_local) 
            
            p_overpressure_local = sim.calculate_overpressure_ground_new(D_m_local, imp_energy_local)
            v_pressure = sim.calculate_overpressure_vulnerability(p_overpressure_local)
            
            wind_velocity_local = sim.calculate_peak_wind_velocity(p_overpressure_local)
            v_wind = sim.calculate_wind_vulnerability(wind_velocity_local)
            
            phi_ground_local = sim.calculate_thermal_exposure(imp_energy_local, r_km)
            v_thermal = sim.calculate_thermal_vulnerability(phi_ground_local)
            
            M_local = sim.calculate_seismic_magnitude(imp_energy_local)
            M_eff_local = sim.calculate_effective_seismic_magnitude(M_local, r_km)
            v_seismic = sim.calculate_seismic_vulnerability(M_eff_local)
            
            D_tc_local = sim.calculate_transient_crater_diameter(v_surface_local, rho_target, sim.entry_angle_deg)
            t_e_local = sim.calculate_ejecta_thickness(D_tc_local, r_km) if D_tc_local > 0 else 0
            v_ejecta = sim.calculate_ejecta_vulnerability(t_e_local)
            
            point_vulnerability = 1.0 - (
                (1.0 - v_pressure) * (1.0 - v_wind) * (1.0 - v_thermal) *
                (1.0 - v_seismic) * (1.0 - v_ejecta)
            )
            
        elif entry_results["event_type"] == "airburst":
            z_b_local = entry_results.get("airburst_altitude", 0)
            mass_local = sim.density * (4.0/3.0) * math.pi * ((sim.diameter/2)**3)
            KE_initial_local = 0.5 * mass_local * (sim.v0**2)
            KE_post_local = 0.5 * mass_local * (entry_results["post_breakup_velocity"]**2) 
            airburst_energy_local = max(KE_post_local, KE_initial_local - KE_post_local) 
            
            p_overpressure_local = sim.calculate_overpressure_airburst(D_m_local, z_b_local, airburst_energy_local, entry_results.get("z_star",0))
            v_pressure = sim.calculate_overpressure_vulnerability(p_overpressure_local)
            
            wind_velocity_local = sim.calculate_peak_wind_velocity(p_overpressure_local)
            v_wind = sim.calculate_wind_vulnerability(wind_velocity_local)
            
            phi_airburst_local = sim.calculate_airburst_thermal_flux(airburst_energy_local, z_b_local, D_m_local)
            v_thermal = sim.calculate_thermal_vulnerability(phi_airburst_local)
            
            point_vulnerability = 1.0 - (
                (1.0 - v_pressure) * (1.0 - v_wind) * (1.0 - v_thermal)
            )
        return point_vulnerability

    left_km, right_km = r_min, r_max
    best_distance_km = 0.0 

    if calculate_vulnerability_at_distance(r_min) < threshold:
        return 0.0 

    best_distance_km = r_min 

    for _ in range(100): 
        mid_km = (left_km + right_km) / 2
        if mid_km <= 0: 
            break 
        if calculate_vulnerability_at_distance(mid_km) >= threshold:
            best_distance_km = mid_km 
            left_km = mid_km 
        else:
            right_km = mid_km 
        
        if (right_km - left_km) < tol:
            break
            
    if calculate_vulnerability_at_distance(best_distance_km) >= threshold:
        return best_distance_km
    else:
        return 0.0