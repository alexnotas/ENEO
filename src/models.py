"""
ENEO Asteroid Impact Simulation - Core Physics Models

This module contains the main AsteroidImpactSimulation class that models the complete
physics of asteroid impact events from atmospheric entry through ground effects.
"""

import math
import numpy as N
from src.utils import (
    H, rho0, C_D, g_E_atmos, g_E_crater, fp_limit, rho_target, P0, c0, 
    acoustic_efficiency, R_EARTH,
    km_to_m, m_to_km, convert_energy_j_to_mt, compute_scaling_factor,
    curvature_adjustment_factor, get_ocean_depth_from_geotiff, WATER_DENSITY_CONSTANT
)
from src.thresholds import get_blast_thresholds, SEISMIC_VULNERABILITY_THRESHOLD, EJECTA_VULNERABILITY_THRESHOLD

from src.vulnerability_models import (
    fun_CraterVulnerability, fun_SeismicVulnerability,
    fun_OverpressureVulnerability, fun_ThermRadVulnerability,
    fun_HighWindVulnerability, fun_EjectaBlanketVulnerability
)
from src.translation_utils import get_translation

class AsteroidImpactSimulation:
    """
    Complete NEO impact simulation.
    """
    
    def __init__(self, diameter, velocity_km_s, density=3000, entry_angle_deg=55):
        self.diameter = diameter
        self.velocity_km_s = velocity_km_s
        self.density = density
        self.entry_angle_deg = entry_angle_deg
        self.v0 = km_to_m(velocity_km_s)
        self.theta = math.radians(entry_angle_deg)
    
    # Energy Calculations
    def calculate_asteroid_energy(self):
        """Calculate initial kinetic energy of asteroid before atmospheric entry."""
        radius = self.diameter / 2.0
        volume = (4.0 / 3.0) * math.pi * (radius ** 3)
        mass = self.density * volume
        kinetic_energy = 0.5 * mass * (self.v0 ** 2)
        energy_mt = convert_energy_j_to_mt(kinetic_energy)
        return kinetic_energy, energy_mt

    def calculate_impact_energy(self, v_impact):
        """Calculate impact energy using actual impact velocity."""
        radius = self.diameter / 2.0
        volume = (4.0/3.0) * math.pi * (radius ** 3)
        mass = self.density * volume
        impact_energy = 0.5 * mass * (v_impact ** 2)
        energy_mt = convert_energy_j_to_mt(impact_energy)
        return impact_energy, energy_mt

    # Atmospheric Entry Calculations
    def compute_yield_strength(self):
        """Calculate asteroid material yield strength based on density."""
        return 10 ** (2.107 + 0.0624 * math.sqrt(self.density))
    
    def v_before_breakup(self, z):
        """Calculate asteroid velocity at altitude z before potential breakup."""
        rho_z = rho0 * math.exp(-z/H)
        return self.v0 * math.exp(- (3 * C_D * H * rho_z) / (4 * self.density * self.diameter * math.sin(self.theta)))

    def pancake_factor(self, z, z_star, l, L0):
        """Calculate pancaking factor for flattened asteroid fragments."""
        return math.sqrt(1 + (2 * H / l)**2 * (math.exp((z_star - z)/(2 * H)) - 1)**2)
    
    def compute_residual_velocity(self, z, v_z_star, z_star, L0, l):
        """Calculate residual velocity after atmospheric passage using numerical integration."""
        N = 100
        dz = (z_star - z) / N
        integral = 0.0
        
        for i in range(N + 1):
            z_prime = z + i * dz
            weight = 1.0 if i not in (0, N) else 0.5
            integral += weight * math.exp((z_star - z_prime)/H) * (self.pancake_factor(z_prime, z_star, l, L0)**2) * dz
        
        exponent = - (3/4) * (C_D * (rho0 * math.exp(-z_star/H))) / (self.density * L0 * math.sin(self.theta)) * integral
        return v_z_star * math.exp(exponent)
    
    def simulate_atmospheric_entry(self):
        """Simulate complete atmospheric entry process including potential breakup."""
        Y_i = self.compute_yield_strength()
        v_i = self.v0
        I_f = 4.07 * (C_D * H * Y_i) / (self.density * self.diameter * v_i**2 * math.sin(self.theta))
        
        # Check if asteroid survives atmospheric passage intact
        if I_f >= 1.0: 
            v_at_surface_intact = self.v_before_breakup(0)
            return {
                "breakup": False,
                "I_f": I_f,
                "z_star": None,
                "airburst_altitude": None,
                "v_breakup": None,
                "dispersion_length": None,
                "post_breakup_velocity": v_at_surface_intact,
                "pancake_factor_ground": 1.0,
                "event_type": "intact"
            }
        
        # Breakup scenario
        z_star = -H * (math.log(Y_i / (rho0 * v_i**2)) + 1.308 - 0.314 * I_f - 1.303 * math.sqrt(1 - I_f))
        v_z_star = self.v_before_breakup(z_star)
        rho_z_star = rho0 * math.exp(-z_star/H)
        l = self.diameter * math.sin(self.theta) * math.sqrt(self.density / (C_D * rho_z_star))
        
        alpha = math.sqrt(fp_limit**2 - 1)
        z_b = z_star - 2 * H * math.log(1 + (l / (2 * H)) * alpha)
        event_type = "airburst" if z_b >= 0 else "ground impact"
        z_b = z_b if z_b >= 0 else 0
        
        z_eval = z_b if event_type == "airburst" else 0
        v_post = self.compute_residual_velocity(z_eval, v_z_star, z_star, self.diameter, l)
        fp_ground = self.pancake_factor(0, z_star, l, self.diameter)
        
        return {
            "breakup": True,
            "I_f": I_f,
            "z_star": z_star,
            "v_breakup": v_z_star,
            "dispersion_length": l,
            "airburst_altitude": z_b,
            "post_breakup_velocity": v_post,
            "event_type": event_type,
            "pancake_factor_ground": fp_ground
        }
    
    # Crater and Melt Calculations
    def calculate_transient_crater_diameter(self, v_i, rho_t, theta_deg, g=g_E_crater):
        """Calculate initial crater diameter immediately after impact."""
        theta = math.radians(theta_deg)
        if rho_t == WATER_DENSITY_CONSTANT:
            coefficient = 1.365
        else:
            coefficient = 1.161
        return coefficient * ((self.density / rho_t)**(1/3)) * (self.diameter**0.78) * (v_i**0.44) * (g**(-0.22)) * (math.sin(theta)**(1/3))

    def calculate_final_crater_diameter(self, D_tc):
        """Calculate final crater diameter after gravitational collapse."""
        D_tc_km = D_tc / 1000.0
        if D_tc_km <= 3.2:
            return 1.25 * D_tc
        else:
            D_fr_km = 1.17 * (D_tc_km ** 1.13) / (3.2 ** 0.13)
            return D_fr_km * 1000.0
    
    def calculate_crater_depth(self, D_tc):
        """Calculate final crater depth based on transient crater size."""
        D_tc_km = D_tc / 1000.0
        if D_tc_km <= 3.2:
            return D_tc / (2 * math.sqrt(2))
        else:
            D_fr = self.calculate_final_crater_diameter(D_tc)
            D_fr_km = D_fr / 1000.0
            d_fr_km = 0.294 * (D_fr_km ** 0.301)
            return d_fr_km * 1000.0

    def calculate_breccia_volume(self, D_fr):
        """Calculate volume of fractured rock (breccia) produced by impact."""
        return 0.032 * (D_fr ** 3)

    def calculate_melt_volume(self, impact_energy, theta_deg):
        """Calculate volume of impact melt produced during crater formation."""
        return 8.9e-12 * impact_energy * math.sin(math.radians(theta_deg))

    def calculate_melt_sheet_thickness(self, V_m, D_tc):
        """Calculate thickness of impact melt sheet in crater."""
        return 4 * V_m / (math.pi * (D_tc ** 2))
    
    # Seismic and Blast Calculations
    def calculate_seismic_magnitude(self, impact_energy):
        """Calculate seismic magnitude generated by impact energy."""
        return 0.67 * math.log10(impact_energy) - 5.87

    def calculate_effective_seismic_magnitude(self, M, r_km):
        """Calculate effective seismic magnitude at specific distance."""
        if r_km < 60:
            return M - 0.0238 * r_km
        elif r_km < 700:
            return M - 0.0048 * r_km - 1.1644
        else:
            Delta = r_km / 6371.0
            return M - 1.66 * math.log10(Delta) - 6.399

    def calculate_seismic_arrival_time(self, r_km):
        """Calculate time for seismic waves to reach specific distance."""
        return r_km / 5.0

    def map_magnitude_to_mmi(self, M_eff):
        """Convert effective seismic magnitude to Modified Mercalli Intensity."""
        if M_eff < 1: return "Not felt"
        elif M_eff < 2: return "I"
        elif M_eff < 3: return "I-II"
        elif M_eff < 4: return "III-IV"
        elif M_eff < 5: return "IV-V"
        elif M_eff < 6: return "VI-VII"
        elif M_eff < 7: return "VII-VIII"
        elif M_eff < 8: return "IX-X"
        elif M_eff < 9: return "X-XI"
        else: return "XII"

    def calculate_fireball_radius(self, impact_energy):
        """Calculate maximum radius of impact fireball."""
        return 0.002 * (impact_energy ** (1/3))

    def calculate_time_of_max_radiation(self, R_f, v_impact):
        """Calculate time when thermal radiation peaks."""
        return R_f / v_impact

    def calculate_irradiation_duration(self, impact_energy, R_f, T_star=3000):
        """Calculate duration of significant thermal radiation."""
        sigma = 5.67e-8
        eta = 3e-3
        return (eta * impact_energy) / (2 * math.pi * (R_f ** 2) * sigma * (T_star ** 4))

    def calculate_thermal_exposure(self, impact_energy, r_km, eta=3e-3):
        """Calculate thermal exposure including proper scaling and curvature effects."""
        r = km_to_m(r_km)
        Phi = (eta * impact_energy) / (2 * math.pi * (r ** 2))
        R_f = self.calculate_fireball_radius(impact_energy)
        f = curvature_adjustment_factor(r, R_f)
        Phi *= f
        return Phi

    def calculate_airburst_thermal_flux(self, airburst_energy, z_b, D):
        """Calculate thermal flux density for airburst events."""
        D_los = math.sqrt(z_b**2 + D**2)
        if D_los == 0:
            return float('inf')

        eta_airburst = 0.007
        phi = (eta_airburst * airburst_energy) / (2 * math.pi * D_los**2)
        return phi

    def calculate_ignition_exposure(self, E_MT, phi_1Mt):
        """Calculate scaled ignition exposure threshold."""
        return phi_1Mt * (E_MT ** (1/6))

    def calculate_ejecta_thickness(self, D_tc, r_km):
        """Calculate thickness of ejecta blanket at specific distance."""
        r = km_to_m(r_km)
        return (D_tc ** 4) / (112 * (r ** 3))

    def calculate_mean_fragment_diameter(self, D_fr, r_km, alpha=2.65):
        """Calculate average size of ejecta fragments at distance."""
        D_fr_km = D_fr / 1000.0
        d_c = 2400 * ((D_fr_km / 2) ** (-1.62))
        return d_c * ((D_fr_km / (2 * r_km)) ** alpha)

    def calculate_ejecta_arrival_time(self, r_distance_km):
        """Calculate time for ejecta to reach specific distance."""
        if r_distance_km >= 10000:
            return None
        
        r_m = km_to_m(r_distance_km)
        Delta = r_m / 6371000.0
        tan_half_Delta = math.tan(Delta / 2.0)
        
        v_e_sq = (2 * g_E_atmos * 6371000.0 * tan_half_Delta) / (1 + tan_half_Delta)
        v_e = math.sqrt(v_e_sq)
        
        ratio = v_e_sq / (g_E_atmos * 6371000.0)
        e_sq = 0.5 * (((ratio - 1)**2) + 1)
        e = -math.sqrt(e_sq) if ratio <= 1 else math.sqrt(e_sq)
        a = v_e_sq / (2 * g_E_atmos * (1 - e**2))
        
        term1 = 2 * math.atan(math.sqrt((1 - e) / (1 + e)) * math.tan(Delta / 4.0))
        term2 = (e * math.sqrt(1 - e**2) * math.sin(Delta / 2.0)) / (1 + e * math.cos(Delta / 2.0))
        T_e = (2 * a**1.5) / (6371000.0 * math.sqrt(g_E_atmos)) * (term1 - term2)
        return T_e

    def calculate_peak_wind_velocity(self, p, P0=P0, c0=c0):
        """Calculate peak wind velocity behind blast wave."""
        return (5 * p / (7 * P0)) * c0 / math.sqrt(1 + (6 * p)/(7 * P0))

    def calculate_blast_arrival_time(self, distance_m, burst_altitude_m=None):
        """Calculate time for blast wave to reach target distance."""
        if burst_altitude_m is not None:
            effective_distance = math.sqrt(distance_m**2 + burst_altitude_m**2)
            return effective_distance / c0
        else:
            return distance_m / c0

    def calculate_sound_intensity(self, p_peak, u_peak, I_ref=1e-12):
        """Calculate sound intensity and decibel level from blast parameters."""
        I = acoustic_efficiency * (p_peak * u_peak) / 2.0
        intensity_db = 10 * math.log10(I / I_ref)
        return I, intensity_db

    # Ground Impact Overpressure
    def calculate_overpressure_ground_new(self, D, impact_energy):
        """Calculate overpressure for ground impact events using scaled distance."""
        E_kt = impact_energy / 4.184e12
        scaling = E_kt ** (1/3) if E_kt > 0 else 1.0
        D1 = D / scaling if D > 0 else 1e-6
        
        p_x = 75000.0
        D_x = 290.0
        
        p_D = (p_x * D_x) / (4 * D1) * (1 + 3 * ((D_x / D1) ** 1.3))
        return p_D

    # Airburst Overpressure with Interpolation
    def calculate_overpressure_airburst(self, D, z_b, impact_energy, z_star):
        """
        Calculate overpressure for airburst events with complex interpolation depending on 
        asteroid size and burst altitude.
        """
        E_kt = impact_energy / 4.184e12
        scaling = E_kt ** (1/3) if E_kt > 0 else 1.0
        D1 = D / scaling if D > 0 else 1e-6
        z_b1 = z_b / scaling if scaling > 0 else z_b

        if z_b1 == 0:
            return self.calculate_overpressure_ground_new(D, impact_energy)

        p_x = 75000.0

        # Small asteroid regime (< 25m diameter)
        if self.diameter < 25:
            p_special = 3.14e11 * ((D1**2 + z_b1**2) ** (-1.3)) + 1.8e7 * ((D1**2 + z_b1**2) ** (-0.565))
            D_x = 289.0 + 0.65 * z_b1
            p_mach = (p_x * D_x) / (4 * D1) * (1 + 3 * ((D_x / D1) ** 1.3)) if D1 > 0 else float('inf')
            
            if z_b1 < 550:
                return p_mach
            else:
                return p_special
        
        # Large asteroid with high breakup altitude regime
        if self.diameter > 25 and z_star >= 45000:
            p_0 = 3.14e11 * (z_b1**(-2.6)) + 1.8e7 * (z_b1**(-1.13)) if z_b1 > 0 else float('inf')
            beta = 34.87 * (z_b1**(-1.73)) if z_b1 > 0 else float('inf')
            p_exp = p_0 * math.exp(-beta * D1)
            p_special = 3.14e11 * ((D1**2 + z_b1**2) ** (-1.3)) + 1.8e7 * ((D1**2 + z_b1**2) ** (-0.565))
            
            if z_b1 >= 550:
                return max(p_exp, p_special)
            else: # z_b1 < 550
                D_x = 289.0 + 0.65 * z_b1
                p_mach = (p_x * D_x) / (4 * D1) * (1 + 3 * ((D_x / D1) ** 1.3)) if D1 > 0 else float('inf')
                z_b_thresh = 550.0
                band_z = 0.1 * z_b_thresh
                base_value = max(p_exp, p_special)

                # Distance interpolation
                denominator_D_m1 = 1.2 * (550 - z_b1)
                D_m1 = (550 * z_b1) / denominator_D_m1 if denominator_D_m1 != 0 else float('inf')
                
                if D_m1 == float('inf'):
                    final_value = base_value
                else:
                    band_D = 0.15 * D_m1
                    if D1 < D_m1 - band_D:
                        final_value = base_value
                    elif D1 > D_m1 + band_D:
                        final_value = p_mach
                    else:
                        denominator_t_D = 2 * band_D
                        t_D = (D1 - (D_m1 - band_D)) / denominator_t_D if denominator_t_D != 0 else 0.5
                        
                        if self.diameter > 45:
                            final_value = (1 - t_D**2) * base_value + t_D**2 * p_mach
                        else:
                            final_value = (1 - t_D) * base_value + t_D * p_mach
                return final_value
        else: # Medium asteroids or low breakup altitude
            if self.diameter < 45 and z_star >= 46000:
                return 3.14e11 * ((D1**2 + z_b1**2) ** (-1.3)) + 1.8e7 * ((D1**2 + z_b1**2) ** (-0.565))

            p_0 = 3.14e11 * (z_b1**(-2.6)) + 1.8e7 * (z_b1**(-1.13)) if z_b1 > 0 else float('inf')
            beta = 34.87 * (z_b1**(-1.73)) if z_b1 > 0 else float('inf')
            p_exp = p_0 * math.exp(-beta * D1)
            D_x = 289.0 + 0.65 * z_b1
            p_mach = (p_x * D_x) / (4 * D1) * (1 + 3 * ((D_x / D1) ** 1.3)) if D1 > 0 else float('inf')
            z_b_thresh = 550.0
            band_z = 0.1 * z_b_thresh
            
            base_value = p_exp

            if z_b1 >= z_b_thresh and z_b1 < z_b_thresh + band_z:
                t_z = (z_b1 - z_b_thresh) / band_z
                base_value = (1 - t_z) * p_mach + t_z * p_exp
            
            if z_b1 < z_b_thresh:
                denominator_D_m1 = 1.2 * (550 - z_b1)
                D_m1 = (550 * z_b1) / denominator_D_m1 if denominator_D_m1 != 0 else float('inf')
                if D_m1 == float('inf'):
                    final_value = base_value 
                else:
                    band_D = 0.15 * D_m1
                    if D1 < D_m1 - band_D:
                        final_value = base_value 
                    elif D1 > D_m1 + band_D:
                        final_value = p_mach     
                    else:
                        denominator_t_D = 2 * band_D
                        t_D = (D1 - (D_m1 - band_D)) / denominator_t_D if denominator_t_D != 0 else 0.5
                        
                        if self.diameter > 45:
                            final_value = (1 - t_D**2) * base_value + t_D**2 * p_mach
                        else:
                            final_value = (1 - t_D) * base_value + t_D * p_mach
            else:
                final_value = base_value
            
            return final_value

    def calculate_damage_category(self, overpressure_pa: float) -> str:
        """Map overpressure value to descriptive damage category."""
        for desc, threshold_value in get_blast_thresholds():
            if overpressure_pa >= threshold_value:
                return desc
        return get_translation("thresholds.fallbackMessages.lightStructuralDamage", "Light or no structural damage")

    def calculate_overpressure_vulnerability(self, press):
        return fun_OverpressureVulnerability(press)

    def calculate_wind_vulnerability(self, v_wind):
        return fun_HighWindVulnerability(v_wind)

    def calculate_thermal_vulnerability(self, phi):
        return fun_ThermRadVulnerability(phi)

    def calculate_seismic_vulnerability(self, eff_mag):
        return fun_SeismicVulnerability(eff_mag)

    def calculate_ejecta_vulnerability(self, t_ejecta):
        return fun_EjectaBlanketVulnerability(t_ejecta)

    def calculate_crater_vulnerability(self, dist_m, final_crater_diam):
        return fun_CraterVulnerability(dist_m, final_crater_diam)

    # Tsunami Calculations
    def calculate_tsunami_effects(self, impact_velocity_at_surface, latitude, longitude):
        """Calculate tsunami effects for ocean impacts."""
        tsunami_results = {
            "ocean_depth": 0.0,
            "transient_cavity_diameter_water": 0.0,
            "max_amplitude_at_source": 0.0,
            "cavity_radius_water": 0.0,
            "is_on_land": True,
            "error_message": None
        }

        # Validate input coordinates
        if latitude is None or longitude is None:
            tsunami_results["error_message"] = "Latitude or longitude not provided for tsunami calculation."
            return tsunami_results # Cannot calculate without location

        # Determine water depth at impact location
        ocean_depth = get_ocean_depth_from_geotiff(latitude, longitude)

        if ocean_depth is None: # Critical error reading GeoTIFF
            tsunami_results["error_message"] = "Could not determine ocean depth due to map data error."
            return tsunami_results
        
        tsunami_results["ocean_depth"] = ocean_depth

        # Check if impact occurs on land
        if ocean_depth <= 0:
            tsunami_results["is_on_land"] = True
            return tsunami_results  # No tsunami for land impacts
        
        tsunami_results["is_on_land"] = False

        # Calculate transient cavity diameter in water using water density
        d_tc_water = self.calculate_transient_crater_diameter(
            v_i=impact_velocity_at_surface,
            rho_t=WATER_DENSITY_CONSTANT, # Target is water
            theta_deg=self.entry_angle_deg, # Use original impact angle
            g=g_E_crater
        )
        tsunami_results["transient_cavity_diameter_water"] = d_tc_water

        # Calculate tsunami wave parameters
        if d_tc_water <= 0:
            # No significant cavity means no substantial wave generation
            return tsunami_results

        # Calculate maximum wave amplitude (limited by water depth)
        max_amplitude_at_source = min(0.14 * d_tc_water, ocean_depth)
        cavity_radius_water = d_tc_water / 2.0

        tsunami_results["max_amplitude_at_source"] = max_amplitude_at_source
        tsunami_results["cavity_radius_water"] = cavity_radius_water
        
        return tsunami_results

    def calculate_tsunami_amplitude_at_distance(self, max_amplitude_at_source, d_tc_water, distance_m):
        """Calculate tsunami wave amplitude at specific distance from impact."""
        # Validate input parameters
        if max_amplitude_at_source <= 0 or d_tc_water <= 0:
            return 0.0  # No wave if no source amplitude or cavity
        if distance_m <= 0:
            return max_amplitude_at_source # At impact site

        cavity_radius_water = d_tc_water / 2.0

        # Inside cavity - maintain maximum amplitude
        if distance_m <= cavity_radius_water:
            return max_amplitude_at_source
        else: 
            # Outside cavity - amplitude decreases with distance
            return max_amplitude_at_source * (d_tc_water / (2 * distance_m))
