"""
ENEO Asteroid Impact Simulation - Core Physics Models

This module contains the main AsteroidImpactSimulation class that models the complete
physics of asteroid impact events from atmospheric entry through ground effects.
Handles atmospheric fragmentation, crater formation, blast waves, thermal radiation,
seismic effects, ejecta distribution, and tsunami generation.

Key Features:
- Atmospheric entry and breakup modeling with pancaking effects
- Ground and airburst impact scenarios with different overpressure calculations
- Crater formation and thermal effects for land and water impacts
- Seismic wave propagation and magnitude calculations
- Ejecta blanket distribution and fragment size modeling
- Tsunami generation and wave amplitude calculations for ocean impacts
- Comprehensive vulnerability assessment integration

Dependencies:
- Physical constants and utility functions from utils module
- Damage thresholds and categories from thresholds module
- Vulnerability calculation functions from vulnerability_models module
"""

import math
import numpy as N
from utils import (
    H, rho0, C_D, g_E_atmos, g_E_crater, fp_limit, rho_target, P0, c0, 
    acoustic_efficiency, R_EARTH,
    km_to_m, m_to_km, convert_energy_j_to_mt, compute_scaling_factor,
    curvature_adjustment_factor, get_ocean_depth_from_geotiff, WATER_DENSITY_CONSTANT
)
from thresholds import (
    # OVERPRESSURE_VULNERABILITY_THRESHOLDS, # Now handled in vulnerability_models.py
    # WIND_VULNERABILITY_THRESHOLDS,         # Now handled in vulnerability_models.py
    # THERMAL_VULNERABILITY_THRESHOLD,     # Now handled in vulnerability_models.py
    SEISMIC_VULNERABILITY_THRESHOLD, # Retained if still used by fun_SeismicVulnerability directly or for other logic
    EJECTA_VULNERABILITY_THRESHOLD,  # Retained if still used by fun_EjectaBlanketVulnerability directly or for other logic
    get_blast_thresholds
)
# Ensure all necessary vulnerability functions are imported
from vulnerability_models import (
    fun_CraterVulnerability, fun_SeismicVulnerability,
    fun_OverpressureVulnerability, fun_ThermRadVulnerability,
    fun_HighWindVulnerability, fun_EjectaBlanketVulnerability
)
from translation_utils import get_translation

class AsteroidImpactSimulation:
    """
    Complete asteroid impact simulation modeling atmospheric entry through ground effects.
    
    This class simulates the entire impact sequence: atmospheric entry with potential
    fragmentation, surface impact or airburst, crater formation, blast wave propagation,
    thermal radiation, seismic effects, ejecta distribution, and tsunami generation
    for ocean impacts.
    
    Attributes:
        diameter: Asteroid diameter in meters
        velocity_km_s: Initial velocity in km/s
        density: Asteroid density in kg/m³
        entry_angle_deg: Entry angle from horizontal in degrees
        v0: Initial velocity converted to m/s
        theta: Entry angle in radians
    """
    
    def __init__(self, diameter, velocity_km_s, density=3000, entry_angle_deg=55):
        """
        Initialize asteroid impact simulation with physical parameters.
        
        Args:
            diameter: Asteroid diameter in meters
            velocity_km_s: Initial velocity in km/s
            density: Asteroid density in kg/m³ (default: 3000 for typical rocky asteroid)
            entry_angle_deg: Entry angle from horizontal in degrees (default: 55)
        """
        self.diameter = diameter
        self.velocity_km_s = velocity_km_s
        self.density = density
        self.entry_angle_deg = entry_angle_deg
        self.v0 = km_to_m(velocity_km_s)  # Convert to m/s for calculations
        self.theta = math.radians(entry_angle_deg)  # Convert to radians
        # self.latitude = None # Store if needed globally, or pass to methods
        # self.longitude = None
    
    # Energy Calculations
    def calculate_asteroid_energy(self):
        """
        Calculate initial kinetic energy of asteroid before atmospheric entry.
        
        Uses spherical volume assumption to determine mass from diameter and density,
        then calculates kinetic energy from mass and initial velocity.
        
        Returns:
            tuple: (kinetic_energy_joules, energy_megatons)
        """
        radius = self.diameter / 2.0
        volume = (4.0 / 3.0) * math.pi * (radius ** 3)  # Spherical volume
        mass = self.density * volume  # Mass from density and volume
        kinetic_energy = 0.5 * mass * (self.v0 ** 2)  # Kinetic energy formula
        energy_mt = convert_energy_j_to_mt(kinetic_energy)  # Convert to megatons
        return kinetic_energy, energy_mt

    def calculate_impact_energy(self, v_impact):
        """
        Calculate impact energy using actual impact velocity.
        
        Similar to initial energy calculation but uses the velocity at impact
        (which may be reduced due to atmospheric deceleration).
        
        Args:
            v_impact: Impact velocity in m/s
            
        Returns:
            tuple: (impact_energy_joules, energy_megatons)
        """
        radius = self.diameter / 2.0
        volume = (4.0/3.0) * math.pi * (radius ** 3)
        mass = self.density * volume
        impact_energy = 0.5 * mass * (v_impact ** 2)
        energy_mt = convert_energy_j_to_mt(impact_energy)
        return impact_energy, energy_mt

    # Atmospheric Entry Calculations
    def compute_yield_strength(self):
        """
        Calculate asteroid material yield strength based on density.
        
        Uses empirical relationship between material density and mechanical strength.
        Higher density materials generally have higher yield strengths.
        
        Returns:
            float: Yield strength in Pa
        """
        return 10 ** (2.107 + 0.0624 * math.sqrt(self.density))
    
    def v_before_breakup(self, z):
        """
        Calculate asteroid velocity at altitude z before potential breakup.
        
        Models velocity reduction due to atmospheric drag using exponential
        atmospheric density profile and drag equation.
        
        Args:
            z: Altitude in meters
            
        Returns:
            float: Velocity in m/s at altitude z
        """
        rho_z = rho0 * math.exp(-z/H)  # Atmospheric density at altitude z
        return self.v0 * math.exp(- (3 * C_D * H * rho_z) / (4 * self.density * self.diameter * math.sin(self.theta)))

    def pancake_factor(self, z, z_star, l, L0):
        """
        Calculate pancaking factor for flattened asteroid fragments.
        
        Models how asteroid fragments flatten during atmospheric passage,
        affecting their aerodynamic properties and drag.
        
        Args:
            z: Current altitude in meters
            z_star: Breakup altitude in meters
            l: Dispersion length in meters
            L0: Initial diameter in meters
            
        Returns:
            float: Pancaking factor (dimensionless)
        """
        return math.sqrt(1 + (2 * H / l)**2 * (math.exp((z_star - z)/(2 * H)) - 1)**2)
    
    def compute_residual_velocity(self, z, v_z_star, z_star, L0, l):
        """
        Calculate residual velocity after atmospheric passage using numerical integration.
        
        Integrates the effects of atmospheric drag and pancaking on velocity reduction
        from breakup altitude to final altitude (surface or airburst height).
        
        Args:
            z: Final altitude in meters
            v_z_star: Velocity at breakup altitude in m/s
            z_star: Breakup altitude in meters
            L0: Initial diameter in meters
            l: Dispersion length in meters
            
        Returns:
            float: Residual velocity in m/s
        """
        N = 100  # Number of integration steps
        dz = (z_star - z) / N  # Step size
        integral = 0.0
        
        # Numerical integration using trapezoidal rule
        for i in range(N + 1):
            z_prime = z + i * dz
            weight = 1.0 if i not in (0, N) else 0.5  # Trapezoidal weights
            integral += weight * math.exp((z_star - z_prime)/H) * (self.pancake_factor(z_prime, z_star, l, L0)**2) * dz
        
        # Calculate velocity reduction exponent
        exponent = - (3/4) * (C_D * (rho0 * math.exp(-z_star/H))) / (self.density * L0 * math.sin(self.theta)) * integral
        return v_z_star * math.exp(exponent)
    
    def simulate_atmospheric_entry(self):
        """
        Simulate complete atmospheric entry process including potential breakup.
        
        Determines if asteroid breaks up during atmospheric passage and calculates
        resulting velocities, altitudes, and event characteristics. Handles both
        intact passage and fragmentation scenarios.
        
        Returns:
            dict: Atmospheric entry results containing:
                - breakup: Whether asteroid fragments
                - event_type: "intact", "airburst", or "ground impact"
                - velocities and altitudes for various stages
                - dispersion and pancaking parameters
        """
        # Calculate material properties and breakup criteria
        Y_i = self.compute_yield_strength()
        v_i = self.v0
        I_f = 4.07 * (C_D * H * Y_i) / (self.density * self.diameter * v_i**2 * math.sin(self.theta))
        
        # Check if asteroid survives atmospheric passage intact
        if I_f >= 1.0: # Asteroid does not break up
            # Calculate surface velocity for intact object
            v_at_surface_intact = self.v_before_breakup(0) # Velocity at z=0 (surface)
            return {
                "breakup": False,
                "I_f": I_f,
                "z_star": None, # No breakup altitude
                "airburst_altitude": None, # No airburst
                "v_breakup": None, # No breakup velocity
                "dispersion_length": None, # No dispersion
                "post_breakup_velocity": v_at_surface_intact, # Velocity at surface for intact object
                "pancake_factor_ground": 1.0, # No pancaking if no breakup
                "event_type": "intact"
            }
        
        # Breakup scenario - calculate fragmentation dynamics
        z_star = -H * (math.log(Y_i / (rho0 * v_i**2)) + 1.308 - 0.314 * I_f - 1.303 * math.sqrt(1 - I_f))
        v_z_star = self.v_before_breakup(z_star)  # Velocity at breakup
        rho_z_star = rho0 * math.exp(-z_star/H)   # Atmospheric density at breakup
        l = self.diameter * math.sin(self.theta) * math.sqrt(self.density / (C_D * rho_z_star))  # Dispersion length
        
        # Calculate airburst altitude and determine event type
        alpha = math.sqrt(fp_limit**2 - 1)
        z_b = z_star - 2 * H * math.log(1 + (l / (2 * H)) * alpha)
        event_type = "airburst" if z_b >= 0 else "ground impact"
        z_b = z_b if z_b >= 0 else 0  # Ensure non-negative altitude
        
        # Calculate final parameters
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
        """
        Calculate initial crater diameter immediately after impact.
        
        Uses scaling laws that account for projectile and target properties,
        impact velocity, angle, and gravitational field. Different coefficients
        are used for water vs. land impacts.
        
        Args:
            v_i: Impact velocity in m/s
            rho_t: Target density in kg/m³
            theta_deg: Impact angle in degrees
            g: Gravitational acceleration (default: Earth surface gravity)
            
        Returns:
            float: Transient crater diameter in meters
        """
        theta = math.radians(theta_deg)
        
        # Select coefficient based on target material
        if rho_t == WATER_DENSITY_CONSTANT:
            coefficient = 1.365  # Formula for water impacts
        else:
            coefficient = 1.161   # Standard formula for land impacts
        
        return coefficient * ((self.density / rho_t)**(1/3)) * (self.diameter**0.78) * (v_i**0.44) * (g**(-0.22)) * (math.sin(theta)**(1/3))

    def calculate_final_crater_diameter(self, D_tc):
        """
        Calculate final crater diameter after gravitational collapse.
        
        Large craters undergo gravitational modification that increases their
        final diameter beyond the initial transient crater size.
        
        Args:
            D_tc: Transient crater diameter in meters
            
        Returns:
            float: Final crater diameter in meters
        """
        D_tc_km = D_tc / 1000.0
        
        # Small craters maintain simple scaling
        if D_tc_km <= 3.2:
            return 1.25 * D_tc
        else:
            # Large craters undergo complex collapse
            D_fr_km = 1.17 * (D_tc_km ** 1.13) / (3.2 ** 0.13)
            return D_fr_km * 1000.0
    
    def calculate_crater_depth(self, D_tc):
        """
        Calculate final crater depth based on transient crater size.
        
        Depth-to-diameter ratios depend on crater size, with different
        relationships for simple vs. complex craters.
        
        Args:
            D_tc: Transient crater diameter in meters
            
        Returns:
            float: Final crater depth in meters
        """
        D_tc_km = D_tc / 1000.0
        
        # Simple craters have constant depth-to-diameter ratio
        if D_tc_km <= 3.2:
            return D_tc / (2 * math.sqrt(2))
        else:
            # Complex craters are shallower relative to diameter
            D_fr = self.calculate_final_crater_diameter(D_tc)
            D_fr_km = D_fr / 1000.0
            d_fr_km = 0.294 * (D_fr_km ** 0.301)
            return d_fr_km * 1000.0

    def calculate_breccia_volume(self, D_fr):
        """
        Calculate volume of fractured rock (breccia) produced by impact.
        
        Breccia forms from shock-fractured target material around the crater.
        Volume scales with final crater diameter.
        
        Args:
            D_fr: Final crater diameter in meters
            
        Returns:
            float: Breccia volume in cubic meters
        """
        return 0.032 * (D_fr ** 3)

    def calculate_melt_volume(self, impact_energy, theta_deg):
        """
        Calculate volume of impact melt produced during crater formation.
        
        High-energy impacts generate molten rock that pools in the crater.
        Melt volume depends on impact energy and angle.
        
        Args:
            impact_energy: Impact energy in Joules
            theta_deg: Impact angle in degrees
            
        Returns:
            float: Melt volume in cubic meters
        """
        return 8.9e-12 * impact_energy * math.sin(math.radians(theta_deg))

    def calculate_melt_sheet_thickness(self, V_m, D_tc):
        """
        Calculate thickness of impact melt sheet in crater.
        
        Assumes melt distributes evenly across crater floor area.
        
        Args:
            V_m: Melt volume in cubic meters
            D_tc: Transient crater diameter in meters
            
        Returns:
            float: Melt sheet thickness in meters
        """
        return 4 * V_m / (math.pi * (D_tc ** 2))
    
    # Seismic and Blast Calculations
    def calculate_seismic_magnitude(self, impact_energy):
        """
        Calculate seismic magnitude generated by impact energy.
        
        Converts impact energy to equivalent earthquake magnitude using
        empirical relationship between energy release and seismic magnitude.
        
        Args:
            impact_energy: Impact energy in Joules
            
        Returns:
            float: Seismic magnitude (Richter scale)
        """
        return 0.67 * math.log10(impact_energy) - 5.87

    def calculate_effective_seismic_magnitude(self, M, r_km):
        """
        Calculate effective seismic magnitude at specific distance.
        
        Accounts for seismic wave attenuation with distance using different
        attenuation relationships for near, intermediate, and far distances.
        
        Args:
            M: Source seismic magnitude
            r_km: Distance from impact in kilometers
            
        Returns:
            float: Effective magnitude felt at distance
        """
        if r_km < 60:
            # Near field - linear attenuation
            return M - 0.0238 * r_km
        elif r_km < 700:
            # Intermediate field - modified linear attenuation
            return M - 0.0048 * r_km - 1.1644
        else:
            # Far field - geometric spreading dominates
            Delta = r_km / 6371.0  # Angular distance in radians
            return M - 1.66 * math.log10(Delta) - 6.399

    def calculate_seismic_arrival_time(self, r_km):
        """
        Calculate time for seismic waves to reach specific distance.
        
        Uses average seismic wave velocity for travel time estimation.
        
        Args:
            r_km: Distance from impact in kilometers
            
        Returns:
            float: Arrival time in seconds
        """
        return r_km / 5.0  # Assumes 5 km/s average seismic velocity

    def map_magnitude_to_mmi(self, M_eff):
        """
        Convert effective seismic magnitude to Modified Mercalli Intensity.
        
        Maps numerical magnitude values to descriptive intensity scales
        that describe expected damage and human perception.
        
        Args:
            M_eff: Effective seismic magnitude
            
        Returns:
            str: Modified Mercalli Intensity description
        """
        if M_eff < 1:
            return "Not felt"
        elif M_eff < 2:
            return "I"      # Barely perceptible
        elif M_eff < 3:
            return "I-II"   # Felt by few
        elif M_eff < 4:
            return "III-IV" # Felt by many, hanging objects sway
        elif M_eff < 5:
            return "IV-V"   # Felt by all, dishes rattle
        elif M_eff < 6:
            return "VI-VII" # Damage to poorly built structures
        elif M_eff < 7:
            return "VII-VIII" # Considerable damage to ordinary buildings
        elif M_eff < 8:
            return "IX-X"   # Buildings shifted off foundations
        elif M_eff < 9:
            return "X-XI"   # Most structures destroyed
        else:
            return "XII"    # Total destruction

    def calculate_fireball_radius(self, impact_energy):
        """
        Calculate maximum radius of impact fireball.
        
        The fireball is the hot, luminous gas cloud that forms during impact.
        Radius scales with impact energy.
        
        Args:
            impact_energy: Impact energy in Joules
            
        Returns:
            float: Fireball radius in meters
        """
        return 0.002 * (impact_energy ** (1/3))

    def calculate_time_of_max_radiation(self, R_f, v_impact):
        """
        Calculate time when thermal radiation peaks.
        
        Maximum radiation occurs when fireball reaches peak expansion.
        
        Args:
            R_f: Fireball radius in meters
            v_impact: Impact velocity in m/s
            
        Returns:
            float: Time of maximum radiation in seconds
        """
        return R_f / v_impact

    def calculate_irradiation_duration(self, impact_energy, R_f, T_star=3000):
        """
        Calculate duration of significant thermal radiation.
        
        Time period during which thermal radiation can cause ignition
        or burns, based on fireball temperature and energy content.
        
        Args:
            impact_energy: Impact energy in Joules
            R_f: Fireball radius in meters
            T_star: Fireball temperature in Kelvin (default: 3000K)
            
        Returns:
            float: Irradiation duration in seconds
        """
        sigma = 5.67e-8  # Stefan-Boltzmann constant
        eta = 3e-3       # Thermal efficiency
        return (eta * impact_energy) / (2 * math.pi * (R_f ** 2) * sigma * (T_star ** 4))

    def calculate_thermal_exposure(self, impact_energy, r_km, eta=3e-3):
        """
        Calculate thermal exposure including proper scaling and curvature effects.
        
        Thermal exposure determines potential for ignition and burns at distance.
        Includes geometric spreading and Earth curvature corrections.
        
        Args:
            impact_energy: Impact energy in Joules
            r_km: Distance from impact in kilometers
            eta: Thermal efficiency (default: 0.003)
            
        Returns:
            float: Thermal exposure in J/m²
        """
        r = km_to_m(r_km)
        
        # Initial thermal exposure calculation with geometric spreading
        Phi = (eta * impact_energy) / (2 * math.pi * (r ** 2))
        
        # Apply Earth curvature adjustment for long distances
        R_f = self.calculate_fireball_radius(impact_energy)
        f = curvature_adjustment_factor(r, R_f)
        Phi *= f
        return Phi

    def calculate_airburst_thermal_flux(self, airburst_energy, z_b, D):
        """
        Calculate thermal flux density for airburst events.
        
        Airbursts have different thermal characteristics than ground impacts,
        with energy distributed over a larger area but potentially less
        atmospheric absorption.
        
        Args:
            airburst_energy: Energy of airburst in Joules
            z_b: Burst altitude in meters
            D: Ground distance from burst point in meters
            
        Returns:
            float: Thermal flux density in J/m²
        """
        # Calculate line-of-sight distance to burst point
        D_los = math.sqrt(z_b**2 + D**2)
        if D_los == 0:
            return float('inf')  # At the point of burst

        eta_airburst = 0.007  # Specific thermal efficiency for airbursts
        phi = (eta_airburst * airburst_energy) / (2 * math.pi * D_los**2)
        
        # Note: Curvature adjustment not applied for airbursts
        return phi

    def calculate_ignition_exposure(self, E_MT, phi_1Mt):
        """
        Calculate scaled ignition exposure threshold.
        
        Ignition thresholds scale with impact energy according to
        established scaling laws for thermal radiation.
        
        Args:
            E_MT: Impact energy in megatons
            phi_1Mt: Reference ignition threshold for 1 MT impact
            
        Returns:
            float: Scaled ignition threshold in J/m²
        """
        return phi_1Mt * (E_MT ** (1/6))

    def calculate_ejecta_thickness(self, D_tc, r_km):
        """
        Calculate thickness of ejecta blanket at specific distance.
        
        Impact ejecta forms a blanket of debris around the crater,
        with thickness decreasing rapidly with distance.
        
        Args:
            D_tc: Transient crater diameter in meters
            r_km: Distance from crater center in kilometers
            
        Returns:
            float: Ejecta thickness in meters
        """
        r = km_to_m(r_km)
        return (D_tc ** 4) / (112 * (r ** 3))

    def calculate_mean_fragment_diameter(self, D_fr, r_km, alpha=2.65):
        """
        Calculate average size of ejecta fragments at distance.
        
        Fragment size depends on crater size and distance, with larger
        fragments landing closer to the crater.
        
        Args:
            D_fr: Final crater diameter in meters
            r_km: Distance from crater in kilometers
            alpha: Fragmentation scaling exponent (default: 2.65)
            
        Returns:
            float: Mean fragment diameter in meters
        """
        D_fr_km = D_fr / 1000.0
        d_c = 2400 * ((D_fr_km / 2) ** (-1.62))  # Critical fragment size
        return d_c * ((D_fr_km / (2 * r_km)) ** alpha)

    def calculate_ejecta_arrival_time(self, r_distance_km):
        """
        Calculate time for ejecta to reach specific distance.
        
        Uses ballistic trajectory calculation for ejecta flight time,
        accounting for Earth's gravity and curvature.
        
        Args:
            r_distance_km: Distance from crater in kilometers
            
        Returns:
            float: Arrival time in seconds, or None if beyond range
        """
        if r_distance_km >= 10000:
            return None  # Beyond reasonable ejecta range
        
        r_m = km_to_m(r_distance_km)
        Delta = r_m / 6371000.0  # Angular distance
        tan_half_Delta = math.tan(Delta / 2.0)
        
        # Calculate escape velocity needed
        v_e_sq = (2 * g_E_atmos * 6371000.0 * tan_half_Delta) / (1 + tan_half_Delta)
        v_e = math.sqrt(v_e_sq)
        
        # Calculate orbital parameters
        ratio = v_e_sq / (g_E_atmos * 6371000.0)
        e_sq = 0.5 * (((ratio - 1)**2) + 1)
        e = -math.sqrt(e_sq) if ratio <= 1 else math.sqrt(e_sq)
        a = v_e_sq / (2 * g_E_atmos * (1 - e**2))
        
        # Calculate flight time using orbital mechanics
        term1 = 2 * math.atan(math.sqrt((1 - e) / (1 + e)) * math.tan(Delta / 4.0))
        term2 = (e * math.sqrt(1 - e**2) * math.sin(Delta / 2.0)) / (1 + e * math.cos(Delta / 2.0))
        T_e = (2 * a**1.5) / (6371000.0 * math.sqrt(g_E_atmos)) * (term1 - term2)
        return T_e

    def calculate_peak_wind_velocity(self, p, P0=P0, c0=c0):
        """
        Calculate peak wind velocity behind blast wave.
        
        Converts overpressure to wind velocity using gas dynamics
        relationships for shock wave propagation.
        
        Args:
            p: Overpressure in Pa
            P0: Ambient pressure (default from utils)
            c0: Sound speed (default from utils)
            
        Returns:
            float: Peak wind velocity in m/s
        """
        return (5 * p / (7 * P0)) * c0 / math.sqrt(1 + (6 * p)/(7 * P0))

    def calculate_blast_arrival_time(self, distance_m, burst_altitude_m=None):
        """
        Calculate time for blast wave to reach target distance.
        
        For airbursts, uses slant distance to burst point.
        For ground impacts, uses surface distance.
        
        Args:
            distance_m: Surface distance in meters
            burst_altitude_m: Burst altitude for airbursts (None for ground impact)
            
        Returns:
            float: Blast arrival time in seconds
        """
        if burst_altitude_m is not None:
            # Calculate slant distance for airburst
            effective_distance = math.sqrt(distance_m**2 + burst_altitude_m**2)
            return effective_distance / c0
        else:
            # Direct surface distance for ground impact
            return distance_m / c0

    def calculate_sound_intensity(self, p_peak, u_peak, I_ref=1e-12):
        """
        Calculate sound intensity and decibel level from blast parameters.
        
        Converts blast pressure and wind velocity to acoustic intensity
        and sound pressure level for hearing damage assessment.
        
        Args:
            p_peak: Peak overpressure in Pa
            u_peak: Peak wind velocity in m/s
            I_ref: Reference intensity in W/m² (default: 1e-12)
            
        Returns:
            tuple: (intensity_W/m², intensity_dB)
        """
        I = acoustic_efficiency * (p_peak * u_peak) / 2.0
        intensity_db = 10 * math.log10(I / I_ref)
        return I, intensity_db

    # Ground Impact Overpressure
    def calculate_overpressure_ground_new(self, D, impact_energy):
        """
        Calculate overpressure for ground impact events.
        
        Uses scaled distance approach to determine blast overpressure
        at specified distance from ground impact crater.
        
        Args:
            D: Distance from impact in meters
            impact_energy: Impact energy in Joules
            
        Returns:
            float: Overpressure in Pa
        """
        E_kt = impact_energy / 4.184e12  # Convert to kilotons TNT equivalent
        scaling = E_kt ** (1/3) if E_kt > 0 else 1.0  # Cube root scaling
        D1 = D / scaling if D > 0 else 1e-6  # Scaled distance, avoid division by zero
        
        # Ground impact overpressure parameters
        p_x = 75000.0  # Reference pressure
        D_x = 290.0    # Reference distance
        
        # Calculate overpressure using empirical formula
        p_D = (p_x * D_x) / (4 * D1) * (1 + 3 * ((D_x / D1) ** 1.3))
        return p_D

    # Airburst Overpressure with Interpolation
    def calculate_overpressure_airburst(self, D, z_b, impact_energy, z_star):
        """
        Calculate overpressure for airburst events with complex interpolation.
        
        This method handles the complex physics of airburst overpressure generation,
        including different regimes based on asteroid size, burst altitude, and
        distance. Uses multiple interpolation schemes to transition smoothly
        between different pressure calculation methods.
        
        Key features:
        - Different formulations for small vs. large asteroids
        - Mach vs. exponential pressure decay regimes
        - Smooth interpolation between calculation methods
        - Special handling for high-altitude vs. low-altitude bursts
        
        Args:
            D: Ground distance from burst point in meters
            z_b: Burst altitude in meters
            impact_energy: Impact energy in Joules
            z_star: Atmospheric breakup altitude in meters
            
        Returns:
            float: Overpressure in Pa
        """
        # Convert to scaled parameters
        E_kt = impact_energy / 4.184e12  # Convert to kilotons TNT equivalent
        scaling = E_kt ** (1/3) if E_kt > 0 else 1.0
        D1 = D / scaling if D > 0 else 1e-6 # Avoid division by zero if D is 0
        z_b1 = z_b / scaling if scaling > 0 else z_b # Avoid division by zero if scaling is 0

        # Handle ground impact case after scaling
        if z_b1 == 0: # Ground impact case after scaling
            return self.calculate_overpressure_ground_new(D, impact_energy)

        p_x = 75000.0  # Reference pressure parameter

        # Small asteroid regime (< 25m diameter)
        if self.diameter < 25:
            # Calculate both pressure formulations
            p_special = 3.14e11 * ((D1**2 + z_b1**2) ** (-1.3)) + 1.8e7 * ((D1**2 + z_b1**2) ** (-0.565))
            D_x = 289.0 + 0.65 * z_b1
            p_mach = (p_x * D_x) / (4 * D1) * (1 + 3 * ((D_x / D1) ** 1.3)) if D1 > 0 else float('inf')
            
            # Choose appropriate formulation based on altitude
            if z_b1 < 550:
                return p_mach  # Low altitude - use Mach regime
            else:
                return p_special  # High altitude - use special formula
            p_special = 3.14e11 * ((D1**2 + z_b1**2) ** (-1.3)) + 1.8e7 * ((D1**2 + z_b1**2) ** (-0.565))
            D_x = 289.0 + 0.65 * z_b1
            p_mach = (p_x * D_x) / (4 * D1) * (1 + 3 * ((D_x / D1) ** 1.3)) if D1 > 0 else float('inf')
            if z_b1 < 550:
                return p_mach
            else:
                return p_special
        
        # Large asteroid with high breakup altitude regime
        if self.diameter > 25 and z_star >= 45000:
            # Calculate pressure components
            p_0 = 3.14e11 * (z_b1**(-2.6)) + 1.8e7 * (z_b1**(-1.13)) if z_b1 > 0 else float('inf')
            beta = 34.87 * (z_b1**(-1.73)) if z_b1 > 0 else float('inf')
            p_exp = p_0 * math.exp(-beta * D1)  # Exponential decay regime
            p_special = 3.14e11 * ((D1**2 + z_b1**2) ** (-1.3)) + 1.8e7 * ((D1**2 + z_b1**2) ** (-0.565))
            
            # High altitude burst - use maximum of exponential and special formulas
            if z_b1 >= 550:
                return max(p_exp, p_special)
            else: # z_b1 < 550 - Apply complex interpolation
                # Calculate Mach regime pressure
                D_x = 289.0 + 0.65 * z_b1
                p_mach = (p_x * D_x) / (4 * D1) * (1 + 3 * ((D_x / D1) ** 1.3)) if D1 > 0 else float('inf')
                z_b_thresh = 550.0
                band_z = 0.1 * z_b_thresh # 10% interpolation band for altitude
                base_value = max(p_exp, p_special)

                # Altitude-based interpolation (this block may be unreachable due to condition)
                if z_b1 > z_b_thresh: 
                    if z_b1 < z_b_thresh + band_z:
                        t_z = (z_b1 - z_b_thresh) / band_z
                        interp_z = (1 - t_z) * p_mach + t_z * base_value
                        base_value = interp_z
                
                # Distance-based interpolation for low altitude
                if z_b1 < z_b_thresh:
                    # Calculate critical distance for interpolation
                    denominator_D_m1 = 1.2 * (550 - z_b1)
                    D_m1 = (550 * z_b1) / denominator_D_m1 if denominator_D_m1 != 0 else float('inf')
                    
                    if D_m1 == float('inf'):
                        final_value = base_value
                    else:
                        band_D = 0.15 * D_m1 # 15% interpolation band for distance
                        if D1 < D_m1 - band_D:
                            final_value = base_value  # Far field - use base formula
                        elif D1 > D_m1 + band_D:
                            final_value = p_mach      # Near field - use Mach formula
                        else:
                            # Interpolation region - smooth transition
                            denominator_t_D = 2 * band_D
                            if denominator_t_D == 0:
                                t_D = 0.5 # Handle edge case
                            else:
                                t_D = (D1 - (D_m1 - band_D)) / denominator_t_D
                            
                            # Different interpolation for very large asteroids
                            if self.diameter > 45:
                                final_value = (1 - t_D**2) * base_value + t_D**2 * p_mach
                            else:
                                final_value = (1 - t_D) * base_value + t_D * p_mach
                else: 
                    # High altitude case (unreachable due to outer condition)
                    final_value = base_value
                
                return final_value
        else: # Medium asteroids or low breakup altitude regime
            # Special case for medium asteroids with very high breakup altitude
            if self.diameter < 45 and z_star >= 46000:
                return 3.14e11 * ((D1**2 + z_b1**2) ** (-1.3)) + 1.8e7 * ((D1**2 + z_b1**2) ** (-0.565))

            # Standard calculation for medium asteroids
            p_0 = 3.14e11 * (z_b1**(-2.6)) + 1.8e7 * (z_b1**(-1.13)) if z_b1 > 0 else float('inf')
            beta = 34.87 * (z_b1**(-1.73)) if z_b1 > 0 else float('inf')
            p_exp = p_0 * math.exp(-beta * D1)  # Exponential decay regime
            D_x = 289.0 + 0.65 * z_b1
            p_mach = (p_x * D_x) / (4 * D1) * (1 + 3 * ((D_x / D1) ** 1.3)) if D1 > 0 else float('inf')
            z_b_thresh = 550.0
            band_z = 0.1 * z_b_thresh # 10% interpolation band for altitude
            
            base_value = p_exp # Use exponential decay as base

            # Apply altitude-based interpolation near threshold
            if z_b1 >= z_b_thresh and z_b1 < z_b_thresh + band_z: # e.g. 550 <= z_b1 < 605
                t_z = (z_b1 - z_b_thresh) / band_z
                base_value = (1 - t_z) * p_mach + t_z * p_exp
            
            # Distance-based interpolation for low altitude bursts
            if z_b1 < z_b_thresh: # z_b1 < 550
                denominator_D_m1 = 1.2 * (550 - z_b1)
                D_m1 = (550 * z_b1) / denominator_D_m1 if denominator_D_m1 != 0 else float('inf')
                if D_m1 == float('inf'):
                    final_value = base_value 
                else:
                    band_D = 0.15 * D_m1 # 15% interpolation band for distance
                    if D1 < D_m1 - band_D:
                        final_value = base_value 
                    elif D1 > D_m1 + band_D:
                        final_value = p_mach     
                    else:
                        # Smooth interpolation between regimes
                        denominator_t_D = 2 * band_D
                        if denominator_t_D == 0:
                             t_D = 0.5
                        else:
                            t_D = (D1 - (D_m1 - band_D)) / denominator_t_D
                        
                        # Different interpolation for very large asteroids
                        if self.diameter > 45:
                            final_value = (1 - t_D**2) * base_value + t_D**2 * p_mach
                        else:
                            final_value = (1 - t_D) * base_value + t_D * p_mach
            else: # z_b1 >= z_b_thresh (z_b1 >= 550)
                final_value = base_value
            
            return final_value

    def calculate_damage_category(self, overpressure_pa: float) -> str:
        """
        Map overpressure value to descriptive damage category.
        
        Converts numerical overpressure to human-readable damage description
        based on established damage thresholds.
        
        Args:
            overpressure_pa: Overpressure in Pascal
            
        Returns:
            str: Damage category description
        """
        for desc, threshold_value in get_blast_thresholds():
            if overpressure_pa >= threshold_value:
                return desc
        return get_translation("thresholds.fallbackMessages.lightStructuralDamage", "Light or no structural damage")

    def calculate_overpressure_vulnerability(self, press):
        """
        Calculate population vulnerability from overpressure exposure.
        
        Wrapper method that calls the detailed vulnerability calculation
        from the vulnerability_models module.
        
        Args:
            press: Overpressure in Pa
            
        Returns:
            float: Vulnerability fraction (0-1)
        """
        return fun_OverpressureVulnerability(press)

    def calculate_wind_vulnerability(self, v_wind):
        """
        Calculate population vulnerability from high wind exposure.
        
        Args:
            v_wind: Wind velocity in m/s
            
        Returns:
            float: Vulnerability fraction (0-1)
        """
        return fun_HighWindVulnerability(v_wind)

    def calculate_thermal_vulnerability(self, phi):
        """
        Calculate population vulnerability from thermal radiation exposure.
        
        Args:
            phi: Thermal exposure in J/m²
            
        Returns:
            float: Vulnerability fraction (0-1)
        """
        return fun_ThermRadVulnerability(phi)

    def calculate_seismic_vulnerability(self, eff_mag):
        """
        Calculate population vulnerability from seismic effects.
        
        Args:
            eff_mag: Effective seismic magnitude
            
        Returns:
            float: Vulnerability fraction (0-1)
        """
        return fun_SeismicVulnerability(eff_mag)

    def calculate_ejecta_vulnerability(self, t_ejecta): # t_ejecta in meters
        """
        Calculate population vulnerability from ejecta blanket burial.
        
        Args:
            t_ejecta: Ejecta thickness in meters
            
        Returns:
            float: Vulnerability fraction (0-1)
        """
        return fun_EjectaBlanketVulnerability(t_ejecta)

    def calculate_crater_vulnerability(self, dist_m, final_crater_diam):
        """
        Calculate population vulnerability from crater formation.
        
        Args:
            dist_m: Distance from crater center in meters
            final_crater_diam: Final crater diameter in meters
            
        Returns:
            float: Vulnerability fraction (0-1)
        """
        return fun_CraterVulnerability(dist_m, final_crater_diam)

    # Tsunami Calculations
    def calculate_tsunami_effects(self, impact_velocity_at_surface, latitude, longitude):
        """
        Calculate tsunami effects for ocean impacts.
        
        Determines if impact occurs in water and calculates resulting tsunami
        parameters including wave amplitude and cavity dimensions. Returns
        detailed tsunami characteristics or indicates land impact.
        
        Args:
            impact_velocity_at_surface: Impact velocity in m/s
            latitude: Impact latitude in decimal degrees
            longitude: Impact longitude in decimal degrees
            
        Returns:
            dict: Tsunami parameters including:
                - ocean_depth: Water depth at impact site
                - transient_cavity_diameter_water: Initial cavity size
                - max_amplitude_at_source: Peak wave height at impact
                - cavity_radius_water: Cavity radius
                - is_on_land: Whether impact occurs on land
                - error_message: Any calculation errors
        """
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
        """
        Calculate tsunami wave amplitude at specific distance from impact.
        
        Models tsunami wave propagation with amplitude decay as waves spread
        outward from the impact site. Uses simple geometric spreading model.
        
        Args:
            max_amplitude_at_source: Peak wave height at impact site in meters
            d_tc_water: Transient cavity diameter in water in meters
            distance_m: Distance from impact site in meters
            
        Returns:
            float: Wave amplitude at specified distance in meters
        """
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