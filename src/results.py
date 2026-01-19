"""
ENEO Asteroid Impact Simulation - Results Processing Module

This module orchestrates the complete asteroid impact simulation pipeline and formats
results for presentation. It coordinates atmospheric entry calculations, impact effects
modeling, and vulnerability assessments to produce comprehensive impact reports.

Key Functions:
- collect_simulation_results(): Structures simulation data into organized result sets
- run_simulation_full(): Executes complete impact simulation with all effects
- find_vulnerability_distance(): Determines damage zone boundaries using binary search
- find_specific_vulnerability_distance(): Calculates specific hazard distances (e.g., thermal, seismic)

The module integrates multiple physics models to calculate crater formation, thermal
radiation, seismic effects, airblast damage, ejecta distribution, wind effects,
tsunami generation, and population vulnerability across distance zones.

Author: Alexandros Notas
Institution: National Technical University of Athens
Date: July 2025
"""

import math
import numpy as np # Ensure numpy is imported if not already at the top
from src.models import AsteroidImpactSimulation
from src.utils import (
    km_to_m, m_to_km, convert_energy_j_to_mt,
    g_E_atmos, C_D, rho0, rho_target, get_ocean_depth_from_geotiff, WATER_DENSITY_CONSTANT, ETOPO_FILE_PATH
)
from src.vulnerability_models import (
    fun_CraterVulnerability, fun_SeismicVulnerability, fun_OverpressureVulnerability,
    fun_ThermRadVulnerability, fun_HighWindVulnerability, fun_EjectaBlanketVulnerability
)
from src.thresholds import (
    get_selected_thermal_thresholds, get_seismic_thresholds, get_ejecta_thresholds, get_blast_thresholds,
    get_ef_wind_thresholds, get_tsunami_amplitude_thresholds,
    get_wind_damage_category, get_ejecta_damage_category, get_thermal_damage_category # Added new imports
)
from src.translation_utils import get_translation
import math # Ensure math is imported for airburst energy calculation and tsunami calcs

def collect_simulation_results(sim: AsteroidImpactSimulation, entry_results,
                               initial_energy_joules, initial_energy_megatons,
                               specific_energy_joules, specific_energy_type,
                               r_distance, tsunami_data=None): # Added tsunami_data
    """
    Organize simulation results into structured data format for API responses and analysis.
    
    This function takes raw simulation outputs and packages them into a comprehensive
    result dictionary containing input parameters, atmospheric entry details, energy
    calculations, crater formation data, tsunami effects, and vulnerability zones.
    
    Args:
        sim: AsteroidImpactSimulation object with impact parameters
        entry_results: Atmospheric entry simulation results dictionary
        initial_energy_joules: Original kinetic energy before atmospheric entry
        initial_energy_megatons: Initial energy in MT TNT equivalent
        specific_energy_joules: Energy used for damage calculations (impact/airburst/initial)
        specific_energy_type: Type of energy used ("Impact", "Airburst", or "Initial")
        r_distance: Distance from impact point for effect calculations (km)
        tsunami_data: Optional tsunami calculation results
        
    Returns:
        dict: Structured results containing all simulation outputs organized by category
    """
    # Package basic simulation parameters and atmospheric entry results
    results = {
        'input_parameters': {
            'diameter': sim.diameter,
            'density': sim.density,
            'velocity': sim.velocity_km_s,
            'entry_angle': sim.entry_angle_deg,
            'distance': r_distance
        },
        'atmospheric_entry': {
            'breakup': entry_results['breakup'],
            'breakup_altitude': entry_results.get('z_star'),
            'airburst_altitude': entry_results.get('airburst_altitude'),
            'velocity_at_breakup': entry_results.get('v_breakup'),
            'final_velocity': entry_results['post_breakup_velocity'],
            'dispersion_length': entry_results.get('dispersion_length'),
            'event_type': entry_results['event_type']
        },
        'energy': {
            'initial_energy_joules': initial_energy_joules,
            'initial_energy_megatons': initial_energy_megatons,
            'specific_energy_joules': specific_energy_joules,
            'specific_energy_megatons': convert_energy_j_to_mt(specific_energy_joules),
            'specific_energy_type': specific_energy_type # e.g., "Impact", "Airburst", "Initial"
        }
    }
    
    # Add crater formation results for ground impacts
    if entry_results["event_type"] == "ground impact":
        v_surface = entry_results['post_breakup_velocity']
        l_dispersion = entry_results.get('dispersion_length', 0)
        
        # Calculate crater dimensions using scaling laws
        D_tc = sim.calculate_transient_crater_diameter(v_surface, rho_target, sim.entry_angle_deg)
        D_fr = sim.calculate_final_crater_diameter(D_tc)
        depth = sim.calculate_crater_depth(D_tc)
        
        crater_results = {
            "crater_formation": {
                "transient_diameter": D_tc,
                "final_diameter": D_fr,
                "depth": depth
            }
        }
        
        # Add note for crater field formation when dispersion is large
        if l_dispersion >= D_tc:
            crater_results["crater_formation"]["note"] = "Crater field likely due to large dispersion." # Example note
        results['crater'] = crater_results
    
    # Include tsunami results if ocean impact occurred
    if tsunami_data:
        results['tsunami'] = tsunami_data

    # Generate vulnerability zones for population impact assessment
    # Create graduated vulnerability thresholds from 100% down to 5% in 5% increments
    thresholds_to_apply_to_population = [round(x, 2) for x in np.arange(1.0, 0.04, -0.05)] # e.g., [1.0, 0.95, 0.90, ..., 0.05]
    vulnerability_zones = []
    current_start_distance = 0.0 # Starting point for current vulnerability zone
    
    # Find the minimum threshold for boundary condition handling
    lowest_applied_threshold = min(thresholds_to_apply_to_population) if thresholds_to_apply_to_population else 0.0

    # Process each vulnerability threshold to create concentric damage zones
    for applied_threshold_for_zone in thresholds_to_apply_to_population:
        # Determine boundary vulnerability level for this zone's outer edge
        # Apply offset logic: most zones use threshold - 2.5%, lowest zone uses its own value
        actual_vuln_level_defining_this_band_outer_edge = applied_threshold_for_zone
        
        # Apply the offset for all zones except the minimum threshold zone
        if applied_threshold_for_zone > lowest_applied_threshold:
            actual_vuln_level_defining_this_band_outer_edge = applied_threshold_for_zone - 0.025
        
        # Ensure boundary level doesn't drop below practical minimum
        actual_vuln_level_defining_this_band_outer_edge = max(0.001, actual_vuln_level_defining_this_band_outer_edge)

        # Find maximum distance where vulnerability meets or exceeds the boundary level
        end_distance_for_this_band = find_vulnerability_distance(sim, actual_vuln_level_defining_this_band_outer_edge, entry_results)
        
        # Skip negligible zones (very small distances that don't represent meaningful danger areas)
        if abs(current_start_distance) < 1e-9 and end_distance_for_this_band <= 0.01:
            # Skip this insignificant zone, keep current_start_distance for next meaningful zone
            continue

        # Add zone if it represents a meaningful distance range
        if end_distance_for_this_band > current_start_distance:
            vulnerability_zones.append({
                "threshold": applied_threshold_for_zone, # Population vulnerability factor applied
                "start_distance": current_start_distance,
                "end_distance": end_distance_for_this_band
            })
            current_start_distance = end_distance_for_this_band # Set start for next zone
    
    # Package vulnerability analysis results
    results["vulnerability_analysis"] = {
        "zones": vulnerability_zones
    }
    
    return results


def run_simulation_full(diameter, density, velocity_km_s, entry_angle, r_distance, lat=None, lon=None): # Added lat, lon
    """
    Execute complete asteroid impact simulation with all physics models and effects.
    
    This is the main simulation orchestrator that runs the full impact analysis pipeline:
    1. Atmospheric entry and breakup modeling
    2. Impact energy calculations (initial, impact, or airburst)
    3. Crater formation and melt production
    4. Thermal radiation and fireball effects
    5. Seismic wave generation
    6. Ejecta blanket distribution
    7. Airblast and overpressure calculations
    8. Wind damage assessment (EF scale)
    9. Tsunami generation (for ocean impacts)
    10. Population vulnerability modeling
    
    Args:
        diameter: Asteroid diameter in meters
        density: Asteroid density in kg/m³
        velocity_km_s: Initial velocity in km/s
        entry_angle: Entry angle in degrees from horizontal
        r_distance: Distance from impact for effect calculations in km
        lat: Impact latitude (required for tsunami calculations)
        lon: Impact longitude (required for tsunami calculations)
        
    Returns:
        tuple: (formatted_text_results, structured_data_results)
            - formatted_text_results: Human-readable simulation report
            - structured_data_results: Machine-readable data for APIs/visualization
    """
    # Initialize simulation object and calculate basic energetics
    sim = AsteroidImpactSimulation(diameter, velocity_km_s, density, entry_angle)
    initial_energy_joules, initial_energy_megatons = sim.calculate_asteroid_energy()
    impact_text = f"Kinetic energy (before entry): {initial_energy_joules:.2e} J\nEquivalent: {initial_energy_megatons:.2f} MT\n\n"
    
    # Simulate atmospheric entry to determine event type (intact, breakup, airburst, ground impact)
    entry_results = sim.simulate_atmospheric_entry()
    atm_text = "Atmospheric Entry Results:\n"
    # Format atmospheric entry results based on breakup behavior
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
        if entry_results["event_type"] == "airburst":
            atm_text += f"Airburst altitude (z_b): {m_to_km(entry_results['airburst_altitude']):.2f} km\n"
            atm_text += f"Residual velocity at airburst: {m_to_km(entry_results['post_breakup_velocity']):.2f} km/s\n"
        else:
            atm_text += "Fragments reach ground (crater-forming event).\n"
            atm_text += f"Surface impact velocity: {m_to_km(entry_results['post_breakup_velocity']):.2f} km/s\n"
            atm_text += f"Pancake factor at ground: {entry_results['pancake_factor_ground']:.2f}\n"
    
    # Calculate velocity reduction due to atmospheric effects
    perc_reduction = 100.0 * (sim.v0 - entry_results["post_breakup_velocity"]) / sim.v0
    atm_text += f"\nVelocity reduction: {perc_reduction:.2f}%\n"
    
    # Determine which energy value to use for damage calculations
    v_impact_final = entry_results["post_breakup_velocity"] # Key velocity for all calculations
    specific_energy_joules_for_event = initial_energy_joules # Default to initial energy
    specific_energy_type_for_event = "Initial" # Default type

    # Calculate event-specific energy based on impact type
    if entry_results["event_type"] == "ground impact" or entry_results["event_type"] == "intact": # MODIFIED to include intact
        # Use final impact velocity to calculate surface impact energy
        imp_energy_joules, imp_energy_MT = sim.calculate_impact_energy(v_impact_final) # Use v_impact_final
        if entry_results["event_type"] == "ground impact":
            impact_text += f"Impact Energy (Fragments): {imp_energy_joules:.2e} J ({imp_energy_MT:.2f} MT)\n\n"
        else: # intact
            impact_text += f"Impact Energy (Intact Object): {imp_energy_joules:.2e} J ({imp_energy_MT:.2f} MT)\n\n"
        specific_energy_joules_for_event = imp_energy_joules
        specific_energy_type_for_event = "Impact"
    elif entry_results["event_type"] == "airburst":
        # Calculate airburst energy from kinetic energy difference
        mass = sim.density * (4.0/3.0) * math.pi * ((sim.diameter/2)**3)
        # Energy before fragmentation
        KE_initial_for_airburst = 0.5 * mass * (sim.v0)**2
        # Residual kinetic energy at airburst altitude  
        KE_post_for_airburst = 0.5 * mass * (entry_results["post_breakup_velocity"])**2
        # Internal energy released during fragmentation
        KE_internal_for_airburst = KE_initial_for_airburst - KE_post_for_airburst
        # Use the larger of residual or internal energy for airburst effects
        airburst_energy_joules = max(KE_post_for_airburst, KE_internal_for_airburst)
        airburst_energy_MT = convert_energy_j_to_mt(airburst_energy_joules)
        impact_text += f"Airburst Energy: {airburst_energy_joules:.2e} J ({airburst_energy_MT:.2f} MT)\n\n"
        specific_energy_joules_for_event = airburst_energy_joules
        specific_energy_type_for_event = "Airburst"
    else: # "intact" or other cases
        # Fall back to velocity reporting for undefined event types
        impact_text += f"Estimated surface impact velocity: {m_to_km(v_impact_final):.2f} km/s\n"
        impact_text += "No specific impact or airburst energy computed for this event type.\n\n"
        # Energy values remain at initial defaults

    # Transitional Event Check
    transitional = False
    if entry_results["event_type"] == "ground impact":
        # Use v_impact_final here
        imp_energy, _ = sim.calculate_impact_energy(v_impact_final)
        candidate_distances = [1, 5, 10]
        candidate_deltas = [5, 10, 20, 30, 40]
        for delta in candidate_deltas:
            if sim.diameter - delta <= 0:
                continue
            temp_sim = AsteroidImpactSimulation(sim.diameter - delta, sim.velocity_km_s, sim.density, sim.entry_angle_deg)
            temp_entry = temp_sim.simulate_atmospheric_entry()
            for d in candidate_distances:
                D_test = km_to_m(d)
                p_original_test = sim.calculate_overpressure_ground_new(D_test, imp_energy)
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
                if p_candidate > p_original_test:
                    transitional = True
                    break
            if transitional:
                break
    if transitional:
        impact_text += "Note: Transition region detected due to sensitivity to small diameter changes.\n\n"
    # Crater & Melt (if ground impact)
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
        if m_to_km(v_impact_final) < 15.0:
            thermal_text += f"{get_translation('effectMessages.noFireball', 'No fireball is created, therefore, there is no thermal radiation damage.')}\n"
            phi_ground = 0.0
        else:
            imp_energy, imp_energy_MT = sim.calculate_impact_energy(v_impact_final)
            phi_ground = sim.calculate_thermal_exposure(imp_energy, r_distance)
            thermal_text += f"{get_translation('effectSections.thermalRadiationGroundImpact', 'Thermal Radiation (Ground Impact)')}:\n"
            thermal_text += f"{get_translation('effectDescriptions.calculatedThermalExposure', 'Calculated thermal exposure at')} {r_distance:.2f} km: {phi_ground:.2e} J/m²\n"
            thermal_cat_ground = get_thermal_damage_category(phi_ground)
            thermal_text += f"{get_translation('effectDescriptions.thermalEffectCategory', 'Thermal Effect Category at')} {r_distance:.2f} km: {thermal_cat_ground}\n"
            
            # Add fireball characteristics for ground impact
            R_f_ground = sim.calculate_fireball_radius(imp_energy)
            T_t_ground = sim.calculate_time_of_max_radiation(R_f_ground, v_impact_final)
            tau_t_ground = sim.calculate_irradiation_duration(imp_energy, R_f_ground)
            
            thermal_text += f"{get_translation('effectDescriptions.fireballRadius', 'Fireball radius')}: {R_f_ground:.2f} m\n"
            thermal_text += f"{get_translation('effectDescriptions.timeOfMaxRadiation', 'Time of maximum radiation')}: {T_t_ground:.2f} s\n"
            thermal_text += f"{get_translation('effectDescriptions.irradiationDuration', 'Irradiation duration')}: {tau_t_ground:.2f} s\n"
            
            # Calculate scaled thresholds and find maximum distances
            previous_bound = 0.0
            thermal_zones = []
            
            for desc, phi_1Mt in get_selected_thermal_thresholds():
                # Calculate scaled threshold for this impact energy
                scaled_threshold = sim.calculate_ignition_exposure(imp_energy_MT, phi_1Mt * 1e6)
                
                # Find maximum distance for this threshold
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
        
            for desc, lower, upper in thermal_zones:
                thermal_text += f"{desc}: {lower:.2f} km - {upper:.2f} km\n"
            
    else: # Airburst
        v_airburst = entry_results.get("v_breakup") # Safely get velocity
        z_b = entry_results.get("airburst_altitude") # Safely get altitude
        D_m = km_to_m(r_distance)

        if v_airburst is None or z_b is None:
            thermal_text += f"{get_translation('effectMessages.burnedUpCompletely', 'Asteroid likely burned up completely')}; {get_translation('effectMessages.noThermalRadiation', 'no significant thermal radiation at the surface.')}\n"
            phi = 0.0
        elif m_to_km(v_airburst) < 15.0:
            thermal_text += f"{get_translation('effectMessages.noFireball', 'No fireball is created, therefore, there is no thermal radiation damage.')}\n"
            phi = 0.0
        else:
            # --- Airburst thermal radiation (unchanged) ---
            mass = sim.density * (4.0/3.0) * math.pi * ((sim.diameter/2)**3)
            KE_initial_for_airburst_thermal = 0.5 * mass * (sim.v0)**2
            KE_post_for_airburst_thermal = 0.5 * mass * (entry_results["post_breakup_velocity"])**2
            KE_internal_for_airburst_thermal = KE_initial_for_airburst_thermal - KE_post_for_airburst_thermal
            airburst_energy_joules_for_thermal = max(KE_post_for_airburst_thermal, KE_internal_for_airburst_thermal)

            phi = sim.calculate_airburst_thermal_flux(airburst_energy_joules_for_thermal, z_b, D_m)
            thermal_text += f"{get_translation('effectSections.thermalRadiationAirburst', 'Thermal Radiation (Airburst)')}:\n"
            thermal_text += f"{get_translation('effectDescriptions.calculatedThermalFlux', 'Calculated thermal flux density at')} {r_distance:.2f} km: {phi:.2e} J/m²\n"
            thermal_cat_airburst = get_thermal_damage_category(phi)
            thermal_text += f"{get_translation('effectDescriptions.thermalEffectCategory', 'Thermal Effect Category at')} {r_distance:.2f} km: {thermal_cat_airburst}\n"
            
            # Airburst calculations use airburst_energy_joules_for_thermal and breakup velocity (v_airburst)
            R_f_airburst = sim.calculate_fireball_radius(airburst_energy_joules_for_thermal)
            T_t_airburst = sim.calculate_time_of_max_radiation(R_f_airburst, v_airburst)
            tau_t_airburst = sim.calculate_irradiation_duration(airburst_energy_joules_for_thermal, R_f_airburst)
            
            thermal_text += f"{get_translation('effectDescriptions.fireballRadius', 'Fireball radius')}: {R_f_airburst:.2f} m\n"
            thermal_text += f"{get_translation('effectDescriptions.timeOfMaxRadiation', 'Time of maximum radiation')}: {T_t_airburst:.2f} s\n"
            thermal_text += f"{get_translation('effectDescriptions.irradiationDuration', 'Irradiation duration')}: {tau_t_airburst:.2f} s\n"
            
            E_MT = convert_energy_j_to_mt(airburst_energy_joules_for_thermal)
            def find_airburst_thermal_max_distance(threshold):
                low = 0.01
                high = 1e8
                tol = 0.01
                best = low
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
            
            for desc, phi_1Mt in get_selected_thermal_thresholds():
                scaled_threshold = sim.calculate_ignition_exposure(E_MT, phi_1Mt * 1e6)
                max_dist = find_airburst_thermal_max_distance(scaled_threshold)
                if max_dist > previous_bound:
                    # Remove the m_to_km conversion - distances are already in km
                    thermal_zones.append((desc, previous_bound, max_dist))
                    previous_bound = max_dist
            for desc, lower, upper in thermal_zones:
                thermal_text += f"{desc}: {lower:.2f} km - {upper:.2f} km\n"
    # Seismic Effects
    seismic_text = ""
    if entry_results["event_type"] == "ground impact":
        # Use v_impact_final
        imp_energy, _ = sim.calculate_impact_energy(v_impact_final)
        M = sim.calculate_seismic_magnitude(imp_energy)
        M_eff = sim.calculate_effective_seismic_magnitude(M, r_distance)
        T_s = sim.calculate_seismic_arrival_time(r_distance)
        seismic_text += f"{get_translation('effectSections.seismicEffects', 'Seismic Effects')}:\n"
        seismic_text += f"{get_translation('effectDescriptions.seismicMagnitude', 'Seismic magnitude')}: {M:.2f}\n"
        seismic_text += f"{get_translation('effectDescriptions.effectiveMagnitude', 'Effective magnitude at')} {r_distance:.2f} km: {M_eff:.2f}\n"
        seismic_text += f"{get_translation('effectDescriptions.modifiedMercalliIntensity', 'Modified Mercalli Intensity at')} {r_distance:.2f} km: {sim.map_magnitude_to_mmi(M_eff)}\n"
        seismic_text += f"{get_translation('effectDescriptions.seismicArrivalTime', 'Seismic arrival time at')} {r_distance:.2f} km: {T_s:.2f} s\n"
        def find_seismic_max_distance(threshold):
            low = 0.01
            high = 20000
            tol = 0.01
            best = low
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
        for desc, thresh in get_seismic_thresholds():
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
        seismic_text += f"{get_translation('effectMessages.noSeismicAirburst', 'Airburst event: no seismic effects.')}\n"
    # Ejecta
    ejecta_text = ""
    if entry_results["event_type"] == "ground impact":
        t_e = sim.calculate_ejecta_thickness(D_tc, r_distance)
        d_mean = sim.calculate_mean_fragment_diameter(D_fr, r_distance)
        T_ejecta = sim.calculate_ejecta_arrival_time(r_distance)
        ejecta_text += f"{get_translation('effectSections.ejectaDeposit', 'Ejecta Deposit')}:\n"
        ejecta_text += f"{get_translation('effectDescriptions.ejectaThickness', 'Ejecta thickness')} {get_translation('translations.at', 'at')} {r_distance:.2f} km: {t_e:.4f} m\n"
        ejecta_cat = get_ejecta_damage_category(t_e)
        ejecta_text += f"{get_translation('effectDescriptions.ejectaDamageCategory', 'Ejecta Damage Category')} {get_translation('translations.at', 'at')} {r_distance:.2f} km: {ejecta_cat}\n"
        ejecta_text += f"{get_translation('effectDescriptions.meanFragmentDiameter', 'Mean fragment diameter at')} {r_distance:.2f} km: {d_mean:.2f} m\n"
        ejecta_text += f"{get_translation('effectDescriptions.ejectaArrivalTime', 'Ejecta arrival time at')} {r_distance:.2f} km: {T_ejecta:.2f} s\n" if T_ejecta is not None else f"{get_translation('effectDescriptions.ejectaArrivalTime', 'Ejecta arrival time at')} {r_distance:.2f} km: N/A\n"
        def find_ejecta_max_distance(threshold):
            low = 0.01
            high = 20000
            tol = 0.01
            best = low
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
        for desc, thresh in get_ejecta_thresholds():
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
        ejecta_text += f"{get_translation('effectMessages.noEjectaAirburst', 'Airburst event: no ejecta.')}\n"
    # Airblast Effects
    airblast_text = ""
    D_m = km_to_m(r_distance)
    if entry_results["event_type"] == "ground impact":
        # Use v_impact_final
        imp_energy, _ = sim.calculate_impact_energy(v_impact_final)
        if transitional:
            mass = sim.density * (4.0/3.0) * math.pi * ((sim.diameter/2)**3)
            KE_initial = 0.5 * mass * (sim.v0)**2
            KE_post = 0.5 * mass * (entry_results["post_breakup_velocity"])**2
            KE_internal = KE_initial - KE_post
            airburst_energy = max(KE_post, KE_internal)
            p_overpressure = sim.calculate_overpressure_airburst(D_m, entry_results.get("airburst_altitude", 0), airburst_energy, entry_results["z_star"])
        else:
            p_overpressure = sim.calculate_overpressure_ground_new(D_m, imp_energy)
        wind_velocity = sim.calculate_peak_wind_velocity(p_overpressure)
        damage_zone = sim.calculate_damage_category(p_overpressure)
        airblast_text += f"{get_translation('effectDescriptions.groundImpactAirBlast', 'Ground Impact Air Blast')}:\n"
        airblast_text += f"{get_translation('effectDescriptions.overpressure', 'Overpressure')}: {p_overpressure:.2f} Pa\n"
        airblast_text += f"{get_translation('effectDescriptions.damageCategory', 'Damage Category')}: {damage_zone}\n"
        I, intensity_db = sim.calculate_sound_intensity(p_overpressure, wind_velocity)
        airblast_text += f"{get_translation('effectDescriptions.soundIntensity', 'Sound Intensity')}: {I:.2e} J/m²\n"
        airblast_text += f"{get_translation('effectDescriptions.spl', 'SPL')}: {intensity_db:.2f} dB\n"
    elif entry_results["event_type"] == "airburst":
        mass = density * (4.0/3.0) * math.pi * ((sim.diameter/2)**3)
        KE_initial = 0.5 * mass * (sim.v0)**2
        KE_post = 0.5 * mass * (entry_results["post_breakup_velocity"])**2
        KE_internal = KE_initial - KE_post
        airburst_energy = max(KE_post, KE_internal)
        p_overpressure = sim.calculate_overpressure_airburst(D_m,
                                                             entry_results["airburst_altitude"],
                                                             airburst_energy,
                                                             entry_results["z_star"])
        airblast_text += f"{get_translation('effectDescriptions.airburstAirBlast', 'Airburst Air Blast')}:\n"
        if r_distance > 3 * m_to_km(entry_results["airburst_altitude"]):
            airblast_text += f"{get_translation('effectDescriptions.overpressure', 'Overpressure')}: {p_overpressure:.2f} Pa\n"
        else:
            p_range_min = p_overpressure
            p_range_max = p_overpressure * 2
            airblast_text += f"{get_translation('effectDescriptions.overpressure', 'Overpressure')}: {p_range_min:.2f} - {p_range_max:.2f} Pa\n"
        wind_velocity = sim.calculate_peak_wind_velocity(p_overpressure)
        damage_zone = sim.calculate_damage_category(p_overpressure)
        airblast_text += f"{get_translation('effectDescriptions.damageCategory', 'Damage Category')}: {damage_zone}\n"
        I, intensity_db = sim.calculate_sound_intensity(p_overpressure, wind_velocity)
        airblast_text += f"{get_translation('effectDescriptions.soundIntensity', 'Sound Intensity')}: {I:.2e} J/m²\n"
        airblast_text += f"{get_translation('effectDescriptions.spl', 'SPL')}: {intensity_db:.2f} dB\n"
    T_b = sim.calculate_blast_arrival_time(D_m, burst_altitude_m=entry_results.get("airburst_altitude"))
    airblast_text += f"\n{get_translation('effectDescriptions.blastArrivalTime', 'Blast arrival time')}: {T_b:.2f} s\n"
    
    # Add blast zone calculations to airblast text
    def find_blast_max_distance(threshold):
        low = 0.01
        high = 20000
        tol = 0.01
        best = low
        # Ensure airburst_energy or imp_energy are correctly scoped or passed if not already
        # Assuming they are available in the outer scope of run_simulation_full
        
        # Determine which energy to use based on event type for clarity
        energy_for_calc = 0
        if entry_results["event_type"] == "airburst":
            # Ensure airburst_energy is defined; it should be from earlier in run_simulation_full
            mass = sim.density * (4.0/3.0) * math.pi * ((sim.diameter/2)**3)
            KE_initial = 0.5 * mass * (sim.v0**2)
            KE_post = 0.5 * mass * (entry_results["post_breakup_velocity"]**2)
            KE_internal = KE_initial - KE_post
            energy_for_calc = max(KE_post, KE_internal)
        else: # ground impact
            # Ensure imp_energy is defined; it should be from earlier in run_simulation_full
            v_surface_for_blast = entry_results['post_breakup_velocity'] # Or however v_surface is determined for blast energy
            # Recalculate or retrieve imp_energy if it's not in the immediate scope
            # For simplicity, assuming imp_energy (from earlier ground impact calculations) is accessible
            # If not, it needs to be recalculated or passed.
            energy_for_calc = imp_energy # Assuming imp_energy is in scope

        while high - low > tol:
            mid = (low + high) / 2
            current_pressure = 0
            if entry_results["event_type"] == "airburst":
                current_pressure = sim.calculate_overpressure_airburst(
                    km_to_m(mid),
                    entry_results["airburst_altitude"],
                    energy_for_calc, # Use defined energy_for_calc
                    entry_results["z_star"]
                )
            else: # ground impact
                current_pressure = sim.calculate_overpressure_ground_new(
                    km_to_m(mid),
                    energy_for_calc # Use defined energy_for_calc
                )
            if current_pressure >= threshold:
                best = mid
                low = mid
            else:
                high = mid
        
        # Final check to ensure the pressure at 'best' actually meets the threshold
        final_pressure_at_best = 0
        if entry_results["event_type"] == "airburst":
            final_pressure_at_best = sim.calculate_overpressure_airburst(km_to_m(best), entry_results["airburst_altitude"], energy_for_calc, entry_results["z_star"])
        else: # ground impact
            final_pressure_at_best = sim.calculate_overpressure_ground_new(km_to_m(best), energy_for_calc)

        if final_pressure_at_best >= threshold:
            return best
        else: # If even at 0.01km, the threshold isn't met
            return 0.0

    previous_max = 0.0
    for desc, thresh in get_blast_thresholds():
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
        wind_effects_content += f"{get_translation('effectDescriptions.windVelocity', 'Wind velocity')} {get_translation('translations.at', 'at')} {r_distance:.2f} km: {wind_velocity:.2f} m/s\n"
        wind_cat = get_wind_damage_category(wind_velocity)
        wind_effects_content += f"{get_translation('effectDescriptions.efScaleCategory', 'EF Scale Category')} {get_translation('translations.at', 'at')} {r_distance:.2f} km: {wind_cat}\n"
    else:
        wind_effects_content += f"{get_translation('effectDescriptions.windVelocity', 'Wind velocity')} {get_translation('translations.at', 'at')} {r_distance:.2f} km: Not calculated (likely no significant overpressure).\n"
        wind_effects_content += f"{get_translation('effectDescriptions.efScaleCategory', 'EF Scale Category')} {get_translation('translations.at', 'at')} {r_distance:.2f} km: Below EF0\n"

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
    sorted_ef_thresholds_desc = sorted(get_ef_wind_thresholds(), key=lambda x: x[1][0], reverse=True)
    
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
            error_prefix = get_translation('tsunamiResults.calculationError', 'Tsunami calculation error')
            tsunami_text += f"{error_prefix}: {tsunami_calc['error_message']}\n"
        elif tsunami_calc["is_on_land"]:
            tsunami_text += get_translation('tsunamiResults.noTsunamiOnLand', 'Impact on land: No tsunami generated.') + "\n"
            if tsunami_calc.get('ocean_depth') is not None: # ocean_depth might be 0 if on land
                 depth_value = f"{tsunami_calc['ocean_depth']:.2f}"
                 note_template = get_translation('tsunamiResults.oceanDepthNote', 'Note: Ocean depth at location reported as {depth} m.')
                 tsunami_text += note_template.format(depth=depth_value) + "\n"
        else:
            ocean_depth_template = get_translation('tsunamiResults.oceanDepthAtImpact', 'Ocean depth at impact: {depth} m')
            tsunami_text += ocean_depth_template.format(depth=f"{tsunami_calc['ocean_depth']:.2f}") + "\n"
            cavity_template = get_translation('tsunamiResults.transientCavityDiameterWater', 'Transient cavity diameter in water: {diameter} m')
            tsunami_text += cavity_template.format(diameter=f"{tsunami_calc['transient_cavity_diameter_water']:.2f}") + "\n"
            max_amp_template = get_translation('tsunamiResults.maxAmplitudeAtSource', 'Max amplitude at source: {amplitude} m')
            tsunami_text += max_amp_template.format(amplitude=f"{tsunami_calc['max_amplitude_at_source']:.2f}") + "\n"

            if tsunami_calc['max_amplitude_at_source'] > 0 and tsunami_calc['transient_cavity_diameter_water'] > 0:
                amplitude_at_r_distance = sim.calculate_tsunami_amplitude_at_distance(
                    tsunami_calc['max_amplitude_at_source'],
                    tsunami_calc['transient_cavity_diameter_water'],
                    km_to_m(r_distance) # r_distance is the user input distance in km
                )
                amplitude_label = get_translation('tsunamiResults.amplitudeLabel', 'Amplitude')
                tsunami_text += f"{amplitude_label} {get_translation('translations.at', 'at')} {r_distance:.2f} km: {amplitude_at_r_distance:.2f} m\n"
                if tsunami_results_data: # Should always be true here as it's tsunami_calc
                    tsunami_results_data['amplitude_at_r_distance'] = amplitude_at_r_distance

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

            sorted_tsunami_thresholds = sorted(get_tsunami_amplitude_thresholds(), key=lambda x: x[1], reverse=True)
            
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
        vuln_model_text_lines.append(f"{get_translation('vulnerabilityLabels.overpressureVulnerability', 'Overpressure Vulnerability')}: {v_pressure:.4f}")

        v_wind = sim.calculate_wind_vulnerability(wind_velocity) # wind_velocity at r_distance
        vuln_model_text_lines.append(f"{get_translation('vulnerabilityLabels.windVulnerability', 'Wind Vulnerability')}: {v_wind:.4f}")

        v_thermal = sim.calculate_thermal_vulnerability(phi_ground) # phi_ground at r_distance
        vuln_model_text_lines.append(f"{get_translation('vulnerabilityLabels.thermalVulnerability', 'Thermal Vulnerability')}: {v_thermal:.4f}")

        v_seismic = sim.calculate_seismic_vulnerability(M_eff) # M_eff at r_distance
        vuln_model_text_lines.append(f"{get_translation('vulnerabilityLabels.seismicVulnerability', 'Seismic Vulnerability')}: {v_seismic:.4f}")

        # Ensure D_tc is calculated if not already available for t_e calculation context
        # D_tc was calculated in the "Crater & Melt" section
        v_ejecta = sim.calculate_ejecta_vulnerability(t_e) # t_e at r_distance
        vuln_model_text_lines.append(f"{get_translation('vulnerabilityLabels.ejectaVulnerability', 'Ejecta Vulnerability')}: {v_ejecta:.4f}")
        
        combined_vuln = 1.0 - (
            (1.0 - v_pressure) * (1.0 - v_wind) * (1.0 - v_thermal) *
            (1.0 - v_seismic) * (1.0 - v_ejecta)
        )
        vuln_model_text_lines.append(f"{get_translation('vulnerabilityLabels.combinedVulnerability', 'Combined Vulnerability')}: {combined_vuln:.4f}")

    elif entry_results["event_type"] == "airburst":
        # Ensure p_overpressure, wind_velocity, phi are the values at r_distance

        v_pressure = sim.calculate_overpressure_vulnerability(p_overpressure) # p_overpressure at r_distance
        vuln_model_text_lines.append(f"{get_translation('vulnerabilityLabels.overpressureVulnerability', 'Overpressure Vulnerability')}: {v_pressure:.4f}")

        v_wind = sim.calculate_wind_vulnerability(wind_velocity) # wind_velocity at r_distance
        vuln_model_text_lines.append(f"{get_translation('vulnerabilityLabels.windVulnerability', 'Wind Vulnerability')}: {v_wind:.4f}")

        v_thermal = sim.calculate_thermal_vulnerability(phi) # phi at r_distance (for airburst)
        vuln_model_text_lines.append(f"{get_translation('vulnerabilityLabels.thermalVulnerability', 'Thermal Vulnerability')}: {v_thermal:.4f}")

        vuln_model_text_lines.append(f"{get_translation('vulnerabilityLabels.seismicVulnerability', 'Seismic Vulnerability')}: N/A")
        vuln_model_text_lines.append(f"{get_translation('vulnerabilityLabels.ejectaVulnerability', 'Ejecta Vulnerability')}: N/A")
        
        combined_vuln = 1.0 - (
            (1.0 - v_pressure) * (1.0 - v_wind) * (1.0 - v_thermal)
        )
        vuln_model_text_lines.append(f"{get_translation('vulnerabilityLabels.combinedVulnerability', 'Combined Vulnerability')}: {combined_vuln:.4f}")
    else: # E.g., "intact" event type if it doesn't fit ground impact or airburst logic for these effects
        vuln_model_text_lines.append(f"{get_translation('vulnerabilityLabels.overpressureVulnerability', 'Overpressure Vulnerability')}: N/A")
        vuln_model_text_lines.append(f"{get_translation('vulnerabilityLabels.windVulnerability', 'Wind Vulnerability')}: N/A")
        vuln_model_text_lines.append(f"{get_translation('vulnerabilityLabels.thermalVulnerability', 'Thermal Vulnerability')}: N/A")
        vuln_model_text_lines.append(f"{get_translation('vulnerabilityLabels.seismicVulnerability', 'Seismic Vulnerability')}: N/A")
        vuln_model_text_lines.append(f"{get_translation('vulnerabilityLabels.ejectaVulnerability', 'Ejecta Vulnerability')}: N/A")
        vuln_model_text_lines.append(f"{get_translation('vulnerabilityLabels.combinedVulnerability', 'Combined Vulnerability')}: N/A")

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
    Find maximum distance where combined vulnerability meets or exceeds threshold.
    
    Uses binary search to efficiently determine the outer boundary of damage zones.
    Calculates combined vulnerability from multiple effects (overpressure, wind,
    thermal, seismic, ejecta) at each test distance until threshold is reached.
    
    Args:
        sim: AsteroidImpactSimulation object
        threshold: Target vulnerability level (0.0 to 1.0)
        entry_results: Atmospheric entry simulation results
        r_min: Minimum search distance in km (default 0.01)
        r_max: Maximum search distance in km (default 20000)
        tol: Search tolerance in km (default 0.01)
        
    Returns:
        float: Maximum distance in km where vulnerability >= threshold,
               or 0.0 if threshold not met
    """
    def calculate_vulnerability_at_distance(r_km):
        """Calculate combined vulnerability from all effects at specified distance."""
        if r_km <= 0: return 1.0 # Maximum vulnerability at impact point
        
        D_m_local = km_to_m(r_km)
        point_vulnerability = 0.0 # Initialize combined vulnerability
        
        # Calculate vulnerability for ground impact events
        if entry_results["event_type"] == "ground impact":
            v_surface_local = entry_results['post_breakup_velocity']
            imp_energy_local, _ = sim.calculate_impact_energy(v_surface_local) 
            
            # Overpressure vulnerability from blast wave
            p_overpressure_local = sim.calculate_overpressure_ground_new(D_m_local, imp_energy_local)
            v_pressure = sim.calculate_overpressure_vulnerability(p_overpressure_local)
            
            # Wind vulnerability from blast-induced winds
            wind_velocity_local = sim.calculate_peak_wind_velocity(p_overpressure_local)
            v_wind = sim.calculate_wind_vulnerability(wind_velocity_local)
            
            # Thermal vulnerability from fireball radiation
            phi_ground_local = sim.calculate_thermal_exposure(imp_energy_local, r_km)
            v_thermal = sim.calculate_thermal_vulnerability(phi_ground_local)
            
            # Seismic vulnerability from ground shaking
            M_local = sim.calculate_seismic_magnitude(imp_energy_local)
            M_eff_local = sim.calculate_effective_seismic_magnitude(M_local, r_km)
            v_seismic = sim.calculate_seismic_vulnerability(M_eff_local)
            
            # Ejecta vulnerability from debris loading
            D_tc_local = sim.calculate_transient_crater_diameter(v_surface_local, rho_target, sim.entry_angle_deg)
            t_e_local = sim.calculate_ejecta_thickness(D_tc_local, r_km) if D_tc_local > 0 else 0
            v_ejecta = sim.calculate_ejecta_vulnerability(t_e_local)
            
            # Combine all vulnerabilities using survival probability method
            point_vulnerability = 1.0 - (
                (1.0 - v_pressure) * (1.0 - v_wind) * (1.0 - v_thermal) *
                (1.0 - v_seismic) * (1.0 - v_ejecta)
            )
            
        # Calculate vulnerability for airburst events  
        elif entry_results["event_type"] == "airburst":
            z_b_local = entry_results.get("airburst_altitude", 0)
            mass_local = sim.density * (4.0/3.0) * math.pi * ((sim.diameter/2)**3)
            KE_initial_local = 0.5 * mass_local * (sim.v0**2)
            KE_post_local = 0.5 * mass_local * (entry_results["post_breakup_velocity"]**2) 
            airburst_energy_local = max(KE_post_local, KE_initial_local - KE_post_local) 
            
            # Overpressure vulnerability from airburst blast
            p_overpressure_local = sim.calculate_overpressure_airburst(D_m_local, z_b_local, airburst_energy_local, entry_results.get("z_star",0))
            v_pressure = sim.calculate_overpressure_vulnerability(p_overpressure_local)
            
            # Wind vulnerability from blast winds
            wind_velocity_local = sim.calculate_peak_wind_velocity(p_overpressure_local)
            v_wind = sim.calculate_wind_vulnerability(wind_velocity_local)
            
            # Thermal vulnerability from airburst fireball
            phi_airburst_local = sim.calculate_airburst_thermal_flux(airburst_energy_local, z_b_local, D_m_local)
            v_thermal = sim.calculate_thermal_vulnerability(phi_airburst_local)
            
            # Combine airburst vulnerabilities (no seismic or ejecta effects)
            point_vulnerability = 1.0 - (
                (1.0 - v_pressure) * (1.0 - v_wind) * (1.0 - v_thermal)
            )
        return point_vulnerability

    # Perform binary search to find threshold distance
    left_km, right_km = r_min, r_max
    best_distance_km = 0.0 

    # Check if threshold is met at minimum distance
    if calculate_vulnerability_at_distance(r_min) < threshold:
        return 0.0 

    best_distance_km = r_min 

    # Binary search for maximum distance meeting threshold
    for _ in range(100): # Maximum iterations to prevent infinite loops
        mid_km = (left_km + right_km) / 2
        if mid_km <= 0: 
            break 
        if calculate_vulnerability_at_distance(mid_km) >= threshold:
            best_distance_km = mid_km 
            left_km = mid_km 
        else:
            right_km = mid_km 
        
        # Check convergence tolerance
        if (right_km - left_km) < tol:
            break
            
    # Verify final result meets threshold requirement
    if calculate_vulnerability_at_distance(best_distance_km) >= threshold:
        return best_distance_km
    else:
        return 0.0

def find_specific_vulnerability_distance(sim, entry_results, vuln_type, threshold, r_min=0.01, r_max=20000.0, tol=0.01):
    """
    Calculates the maximum distance at which a specific vulnerability type 
    (e.g., thermal, seismic) meets or exceeds a defined threshold.
    
    This function uses a binary search algorithm to efficiently find this boundary 
    distance. Beyond this distance, the vulnerability's intensity drops below the 
    given threshold. It accounts for different impact scenarios (ground impacts vs. airbursts) 
    and calculates relevant vulnerability metrics for each distance evaluated.
    
    The binary search converges quickly, making it suitable for dynamic calculations.
    
    Args:
        sim (AsteroidImpactSimulation): An initialized simulation object.
        entry_results (dict): Results from the atmospheric entry simulation, including
            event_type, post_breakup_velocity, airburst_altitude, v_breakup, and z_star.
        vuln_type (str): The type of vulnerability to assess (e.g., 'thermal', 'overpressure').
        threshold (float): The vulnerability threshold value (0.0 to 1.0).
        r_min (float, optional): Minimum search distance in km (default: 0.01 km).
        r_max (float, optional): Maximum search distance in km (default: 20000.0 km).
        tol (float, optional): Tolerance for search convergence in km (default: 0.01 km).
    
    Returns:
        float: The maximum distance (km) where vulnerability is ≥ threshold. 
               Returns r_min if the threshold isn't met even at close distances.
    
    Mathematical Approach:
        The binary search maintains these conditions:
        - The left boundary (left) of the search interval represents a distance where 
          vulnerability(distance) ≥ threshold.
        - The right boundary (right) represents a distance where vulnerability(distance) < threshold.
        
        The search stops when the interval (right - left) is smaller than the tolerance.
    
    Details on vulnerability calculation:
        - For Ground Impacts: Uses models based on impact energy and crater formation.
        - For Airbursts: Employs models based on burst energy and atmospheric blast wave propagation.
        - Returns 0.0 if a vulnerability type isn't applicable (e.g., ejecta from an airburst).
    """
    from src.utils import km_to_m, rho_target
    
    def calculate_specific_vulnerability_at_distance(r_km, vuln_type):
        """Calculate vulnerability value at specific distance for given hazard type."""
        D_m = km_to_m(r_km)  # Convert distance to meters
        
        # Ground impact vulnerability calculations
        if entry_results["event_type"] == "ground impact":
            v_surface = entry_results['post_breakup_velocity']
            imp_energy, _ = sim.calculate_impact_energy(v_surface)
            
            if vuln_type == 'thermal':
                phi_ground = sim.calculate_thermal_exposure(imp_energy, r_km)
                return sim.calculate_thermal_vulnerability(phi_ground)
            
            elif vuln_type == 'overpressure':
                p_overpressure = sim.calculate_overpressure_ground_new(D_m, imp_energy)
                return sim.calculate_overpressure_vulnerability(p_overpressure)
            
            elif vuln_type == 'wind':
                p_overpressure = sim.calculate_overpressure_ground_new(D_m, imp_energy)
                wind_velocity = sim.calculate_peak_wind_velocity(p_overpressure)
                return sim.calculate_wind_vulnerability(wind_velocity)
            
            elif vuln_type == 'seismic':
                M = sim.calculate_seismic_magnitude(imp_energy)
                M_eff = sim.calculate_effective_seismic_magnitude(M, r_km)
                return sim.calculate_seismic_vulnerability(M_eff)
            
            elif vuln_type == 'ejecta':
                D_tc = sim.calculate_transient_crater_diameter(v_surface, rho_target, sim.entry_angle_deg)
                t_e = sim.calculate_ejecta_thickness(D_tc, r_km)
                return sim.calculate_ejecta_vulnerability(t_e)
                
        # Airburst vulnerability calculations
        elif entry_results["event_type"] == "airburst":
            z_b = entry_results["airburst_altitude"]
            mass = sim.density * (4.0/3.0) * 3.14159 * ((sim.diameter/2)**3)
            KE_initial = 0.5 * mass * (sim.v0**2)
            KE_post = 0.5 * mass * (entry_results["post_breakup_velocity"]**2)
            KE_internal = KE_initial - KE_post
            airburst_energy = max(KE_post, KE_internal)
            
            if vuln_type == 'thermal':
                phi = sim.calculate_airburst_thermal_flux(airburst_energy, z_b, D_m)
                return sim.calculate_thermal_vulnerability(phi)
            
            elif vuln_type == 'overpressure':
                p_overpressure = sim.calculate_overpressure_airburst(D_m, z_b, airburst_energy, entry_results["z_star"])
                return sim.calculate_overpressure_vulnerability(p_overpressure)
            
            elif vuln_type == 'wind':
                p_overpressure = sim.calculate_overpressure_airburst(D_m, z_b, airburst_energy, entry_results["z_star"])
                wind_velocity = sim.calculate_peak_wind_velocity(p_overpressure)
                return sim.calculate_wind_vulnerability(wind_velocity)
            
            elif vuln_type == 'seismic':
                return 0.0  # No seismic effects for airburst events
            
            elif vuln_type == 'ejecta':
                return 0.0  # No ejecta for airburst events
        
        return 0.0  # Default case

    # Binary search implementation for efficient distance finding
    left, right = r_min, r_max
    best_distance = r_min
    
    # Converge on the boundary where vulnerability meets threshold
    while right - left > tol:
        mid = (left + right) / 2
        if calculate_specific_vulnerability_at_distance(mid, vuln_type) >= threshold:
            best_distance = mid
            left = mid  # Vulnerability still meets threshold, search further out
        else:
            right = mid  # Vulnerability below threshold, search closer in
    
    return best_distance