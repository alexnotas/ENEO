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
    DAMAGE_CATEGORIES
)
# Ensure all necessary vulnerability functions are imported
from vulnerability_models import (
    fun_CraterVulnerability, fun_SeismicVulnerability,
    fun_OverpressureVulnerability, fun_ThermRadVulnerability,
    fun_HighWindVulnerability, fun_EjectaBlanketVulnerability
)

class AsteroidImpactSimulation:
    """
    A class to simulate the effects of an asteroid impact, from atmospheric entry to ground effects.

    This class models the physical processes involved in an asteroid impact event, including:
    - Energy calculation based on asteroid properties.
    - Atmospheric entry dynamics, including breakup and airburst events.
    - Crater formation and melt production.
    - Seismic, thermal, and blast wave effects.
    - Tsunami generation for impacts in water.

    Attributes:
        diameter (float): The diameter of the asteroid in meters.
        velocity_km_s (float): The entry velocity of the asteroid in kilometers per second.
        density (float): The density of the asteroid in kg/m³.
        entry_angle_deg (float): The entry angle of the asteroid in degrees from the horizontal.
        v0 (float): The initial velocity in m/s.
        theta (float): The entry angle in radians.
    """
    def __init__(self, diameter, velocity_km_s, density=3000, entry_angle_deg=55):
        """
        Initializes the AsteroidImpactSimulation with asteroid parameters.

        Args:
            diameter (float): Asteroid diameter in meters.
            velocity_km_s (float): Asteroid velocity in km/s.
            density (float, optional): Asteroid density in kg/m³. Defaults to 3000.
            entry_angle_deg (float, optional): Entry angle in degrees. Defaults to 55.
        """
        self.diameter = diameter
        self.velocity_km_s = velocity_km_s
        self.density = density
        self.entry_angle_deg = entry_angle_deg
        self.v0 = km_to_m(velocity_km_s)
        self.theta = math.radians(entry_angle_deg)
        # self.latitude = None # Store if needed globally, or pass to methods
        # self.longitude = None
    
    # Energy Calculations
    def calculate_asteroid_energy(self):
        """
        Calculates the initial kinetic energy of the asteroid.

        Returns:
            tuple: A tuple containing:
                - kinetic_energy (float): The kinetic energy in Joules.
                - energy_mt (float): The kinetic energy in megatons of TNT.
        """
        radius = self.diameter / 2.0
        volume = (4.0 / 3.0) * math.pi * (radius ** 3)
        mass = self.density * volume
        kinetic_energy = 0.5 * mass * (self.v0 ** 2)
        energy_mt = convert_energy_j_to_mt(kinetic_energy)
        return kinetic_energy, energy_mt

    def calculate_impact_energy(self, v_impact):
        """
        Calculates the kinetic energy at impact.

        Args:
            v_impact (float): The velocity of the asteroid at impact in m/s.

        Returns:
            tuple: A tuple containing:
                - impact_energy (float): The impact energy in Joules.
                - energy_mt (float): The impact energy in megatons of TNT.
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
        Computes the yield strength of the asteroid based on its density.

        Returns:
            float: The yield strength in Pascals.
        """
        return 10 ** (2.107 + 0.0624 * math.sqrt(self.density))
    
    def v_before_breakup(self, z):
        """
        Calculates the velocity of the asteroid at a given altitude before breakup.

        Args:
            z (float): The altitude in meters.

        Returns:
            float: The velocity at altitude z in m/s.
        """
        rho_z = rho0 * math.exp(-z/H)
        return self.v0 * math.exp(- (3 * C_D * H * rho_z) / (4 * self.density * self.diameter * math.sin(self.theta)))

    def pancake_factor(self, z, z_star, l, L0):
        """
        Calculates the pancake factor, which represents the flattening of the asteroid after breakup.

        Args:
            z (float): The current altitude in meters.
            z_star (float): The breakup altitude in meters.
            l (float): The dispersion length scale.
            L0 (float): The initial diameter of the asteroid.

        Returns:
            float: The pancake factor.
        """
        return math.sqrt(1 + (2 * H / l)**2 * (math.exp((z_star - z)/(2 * H)) - 1)**2)
    
    def compute_residual_velocity(self, z, v_z_star, z_star, L0, l):
        """
        Computes the residual velocity of the fragmented asteroid cloud.

        Args:
            z (float): The altitude at which to calculate the velocity.
            v_z_star (float): The velocity at the breakup altitude.
            z_star (float): The breakup altitude.
            L0 (float): The initial diameter of the asteroid.
            l (float): The dispersion length scale.

        Returns:
            float: The residual velocity in m/s.
        """
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
        """
        Simulates the atmospheric entry of the asteroid, determining if it breaks up.

        Returns:
            dict: A dictionary containing the results of the simulation, including
                  breakup status, altitude, velocity, and event type (intact, airburst, or ground impact).
        """
        Y_i = self.compute_yield_strength()
        v_i = self.v0
        # Fragmentation index I_f determines if the asteroid survives entry intact.
        I_f = 4.07 * (C_D * H * Y_i) / (self.density * self.diameter * v_i**2 * math.sin(self.theta))
        
        if I_f >= 1.0: # Asteroid does not break up
            # For an intact object, calculate its velocity at the surface (z=0)
            # due to atmospheric drag, even if it doesn't fragment.
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
        # Breakup scenario
        # Calculate the altitude where breakup begins (z_star).
        z_star = -H * (math.log(Y_i / (rho0 * v_i**2)) + 1.308 - 0.314 * I_f - 1.303 * math.sqrt(1 - I_f))
        v_z_star = self.v_before_breakup(z_star)
        rho_z_star = rho0 * math.exp(-z_star/H)
        
        # Calculate dispersion length scale 'l'.
        l = self.diameter * math.sin(self.theta) * math.sqrt(self.density / (C_D * rho_z_star))
        alpha = math.sqrt(fp_limit**2 - 1)
        
        # Calculate the altitude of maximum "pancaking" (z_b), which is the airburst altitude.
        z_b = z_star - 2 * H * math.log(1 + (l / (2 * H)) * alpha)
        
        # Determine if it's an airburst or ground impact based on z_b.
        event_type = "airburst" if z_b >= 0 else "ground impact"
        z_b = z_b if z_b >= 0 else 0 # Airburst altitude cannot be negative.
        
        # The evaluation altitude is the airburst altitude or ground level.
        z_eval = z_b if event_type == "airburst" else 0
        
        # Compute the final velocity of the fragment cloud.
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
        Calculates the transient crater diameter.

        Args:
            v_i (float): Impact velocity in m/s.
            rho_t (float): Target density in kg/m³.
            theta_deg (float): Impact angle in degrees.
            g (float, optional): Gravitational acceleration. Defaults to g_E_crater.

        Returns:
            float: The transient crater diameter in meters.
        """
        theta = math.radians(theta_deg)
        # Determine the coefficient based on the target material
        # Use 1.365 for water, 1.161 for land (original formula)
        if rho_t == WATER_DENSITY_CONSTANT:
            coefficient = 1.365
        else:
            coefficient = 1.161
        
        return coefficient * ((self.density / rho_t)**(1/3)) * (self.diameter**0.78) * (v_i**0.44) * (g**(-0.22)) * (math.sin(theta)**(1/3))

    def calculate_final_crater_diameter(self, D_tc):
        """
        Calculates the final crater diameter from the transient crater diameter.

        Args:
            D_tc (float): Transient crater diameter in meters.

        Returns:
            float: The final crater diameter in meters.
        """
        D_tc_km = D_tc / 1000.0
        if D_tc_km <= 3.2: # Simple crater regime
            return 1.25 * D_tc
        else: # Complex crater regime
            D_fr_km = 1.17 * (D_tc_km ** 1.13) / (3.2 ** 0.13)
            return D_fr_km * 1000.0
    
    def calculate_crater_depth(self, D_tc):
        """
        Calculates the final crater depth from the transient crater diameter.

        Args:
            D_tc (float): Transient crater diameter in meters.

        Returns:
            float: The final crater depth in meters.
        """
        D_tc_km = D_tc / 1000.0
        if D_tc_km <= 3.2: # Simple crater
            return D_tc / (2 * math.sqrt(2))
        else: # Complex crater
            D_fr = self.calculate_final_crater_diameter(D_tc)
            D_fr_km = D_fr / 1000.0
            d_fr_km = 0.294 * (D_fr_km ** 0.301)
            return d_fr_km * 1000.0

    def calculate_breccia_volume(self, D_fr):
        """
        Calculates the volume of the breccia lens.

        Args:
            D_fr (float): Final crater diameter in meters.

        Returns:
            float: The volume of breccia in cubic meters.
        """
        return 0.032 * (D_fr ** 3)

    def calculate_melt_volume(self, impact_energy, theta_deg):
        """
        Calculates the volume of impact melt.

        Args:
            impact_energy (float): Impact energy in Joules.
            theta_deg (float): Impact angle in degrees.

        Returns:
            float: The volume of melt in cubic meters.
        """
        return 8.9e-12 * impact_energy * math.sin(math.radians(theta_deg))

    def calculate_melt_sheet_thickness(self, V_m, D_tc):
        """
        Calculates the thickness of the melt sheet in the crater.

        Args:
            V_m (float): Volume of melt in cubic meters.
            D_tc (float): Transient crater diameter in meters.

        Returns:
            float: The thickness of the melt sheet in meters.
        """
        return 4 * V_m / (math.pi * (D_tc ** 2))
    
    # Seismic and Blast Calculations
    def calculate_seismic_magnitude(self, impact_energy):
        """
        Calculates the Richter magnitude of the impact-induced earthquake.

        Args:
            impact_energy (float): Impact energy in Joules.

        Returns:
            float: The seismic magnitude (Richter scale).
        """
        return 0.67 * math.log10(impact_energy) - 5.87

    def calculate_effective_seismic_magnitude(self, M, r_km):
        """
        Calculates the effective seismic magnitude at a given distance.

        Args:
            M (float): The initial seismic magnitude.
            r_km (float): The distance from the impact in kilometers.

        Returns:
            float: The effective seismic magnitude at the given distance.
        """
        if r_km < 60:
            return M - 0.0238 * r_km
        elif r_km < 700:
            return M - 0.0048 * r_km - 1.1644
        else:
            Delta = r_km / 6371.0
            return M - 1.66 * math.log10(Delta) - 6.399

    def calculate_seismic_arrival_time(self, r_km):
        """
        Calculates the arrival time of seismic waves.

        Args:
            r_km (float): The distance from the impact in kilometers.

        Returns:
            float: The arrival time in seconds.
        """
        return r_km / 5.0

    def map_magnitude_to_mmi(self, M_eff):
        """
        Maps effective seismic magnitude to the Modified Mercalli Intensity (MMI) scale.

        Args:
            M_eff (float): The effective seismic magnitude.

        Returns:
            str: The MMI scale rating.
        """
        if M_eff < 1:
            return "Not felt"
        elif M_eff < 2:
            return "I"
        elif M_eff < 3:
            return "I-II"
        elif M_eff < 4:
            return "III-IV"
        elif M_eff < 5:
            return "IV-V"
        elif M_eff < 6:
            return "VI-VII"
        elif M_eff < 7:
            return "VII-VIII"
        elif M_eff < 8:
            return "IX-X"
        elif M_eff < 9:
            return "X-XI"
        else:
            return "XII"

    def calculate_fireball_radius(self, impact_energy):
        """
        Calculates the radius of the fireball.

        Args:
            impact_energy (float): Impact energy in Joules.

        Returns:
            float: The fireball radius in meters.
        """
        return 0.002 * (impact_energy ** (1/3))

    def calculate_time_of_max_radiation(self, R_f, v_impact):
        """
        Calculates the time to reach maximum radiation from the fireball.

        Args:
            R_f (float): Fireball radius in meters.
            v_impact (float): Impact velocity in m/s.

        Returns:
            float: The time of maximum radiation in seconds.
        """
        return R_f / v_impact

    def calculate_irradiation_duration(self, impact_energy, R_f, T_star=3000):
        """
        Calculates the duration of thermal irradiation.

        Args:
            impact_energy (float): Impact energy in Joules.
            R_f (float): Fireball radius in meters.
            T_star (float, optional): Effective temperature of the fireball. Defaults to 3000 K.

        Returns:
            float: The duration of irradiation in seconds.
        """
        sigma = 5.67e-8
        eta = 3e-3
        return (eta * impact_energy) / (2 * math.pi * (R_f ** 2) * sigma * (T_star ** 4))

    def calculate_thermal_exposure(self, impact_energy, r_km, eta=3e-3):
        """Calculate thermal exposure including proper scaling and curvature effects"""

        r = km_to_m(r_km)
        # Initial thermal exposure calculation
        Phi = (eta * impact_energy) / (2 * math.pi * (r ** 2))
        # Calculate fireball radius and apply curvature adjustment
        R_f = self.calculate_fireball_radius(impact_energy)
        f = curvature_adjustment_factor(r, R_f)
        Phi *= f
        return Phi

    def calculate_airburst_thermal_flux(self, airburst_energy, z_b, D):
        """
        Calculate the thermal energy flux density (J/m²) at a target distance for an airburst event.
        This model uses a thermal efficiency (eta) of 0.007 for airbursts.
        The flux is calculated based on the energy released and the slant distance to the target.

        Args:
            airburst_energy (float): The energy of the airburst in Joules.
            z_b (float): The burst altitude in meters.
            D (float): The ground distance from the burst in meters.

        Returns:
            float: The thermal energy flux density (phi) in J/m².
        """
        # Slant distance from the airburst point to the target on the ground.
        D_los = math.sqrt(z_b**2 + D**2)
        if D_los == 0:
            return float('inf')  # At the point of burst

        eta_airburst = 0.007  # Specific thermal efficiency for airbursts
        phi = (eta_airburst * airburst_energy) / (2 * math.pi * D_los**2)
        
        # Note: Curvature adjustment is not applied here as in the ground impact model.
        return phi

    def calculate_ignition_exposure(self, E_MT, phi_1Mt):
        """
        Calculate scaled ignition exposure threshold.

        Args:
            E_MT (float): Impact energy in megatons of TNT.
            phi_1Mt (float): The baseline thermal exposure for a 1 Mt event.

        Returns:
            float: The scaled ignition exposure threshold.
        """
        # Use proper scaling law for thermal radiation
        return phi_1Mt * (E_MT ** (1/6))

    def calculate_ejecta_thickness(self, D_tc, r_km):
        """
        Calculates the thickness of the ejecta blanket at a given distance.

        Args:
            D_tc (float): Transient crater diameter in meters.
            r_km (float): Distance from the crater rim in kilometers.

        Returns:
            float: The ejecta thickness in meters.
        """
        r = km_to_m(r_km)
        return (D_tc ** 4) / (112 * (r ** 3))

    def calculate_mean_fragment_diameter(self, D_fr, r_km, alpha=2.65):
        """
        Calculates the mean diameter of ejecta fragments.

        Args:
            D_fr (float): Final crater diameter in meters.
            r_km (float): Distance from the impact in kilometers.
            alpha (float, optional): A coefficient for fragment size distribution. Defaults to 2.65.

        Returns:
            float: The mean fragment diameter in meters.
        """
        D_fr_km = D_fr / 1000.0
        d_c = 2400 * ((D_fr_km / 2) ** (-1.62))
        return d_c * ((D_fr_km / (2 * r_km)) ** alpha)

    def calculate_ejecta_arrival_time(self, r_distance_km):
        """
        Calculates the arrival time of the ejecta.

        Args:
            r_distance_km (float): The distance from the impact in kilometers.

        Returns:
            float or None: The arrival time in seconds, or None if the distance is too great.
        """
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
        """
        Calculates the peak wind velocity from the overpressure.

        Args:
            p (float): Peak overpressure in Pascals.
            P0 (float, optional): Ambient atmospheric pressure. Defaults to P0.
            c0 (float, optional): Speed of sound. Defaults to c0.

        Returns:
            float: The peak wind velocity in m/s.
        """
        return (5 * p / (7 * P0)) * c0 / math.sqrt(1 + (6 * p)/(7 * P0))

    def calculate_blast_arrival_time(self, distance_m, burst_altitude_m=None):
        """
        Calculates the arrival time of the blast wave.

        Args:
            distance_m (float): The distance from the impact in meters.
            burst_altitude_m (float, optional): The altitude of the airburst in meters. Defaults to None.

        Returns:
            float: The arrival time in seconds.
        """
        if burst_altitude_m is not None:
            effective_distance = math.sqrt(distance_m**2 + burst_altitude_m**2)
            return effective_distance / c0
        else:
            return distance_m / c0

    def calculate_sound_intensity(self, p_peak, u_peak, I_ref=1e-12):
        """
        Calculates the sound intensity in decibels.

        Args:
            p_peak (float): Peak overpressure in Pascals.
            u_peak (float): Peak wind velocity in m/s.
            I_ref (float, optional): Reference sound intensity. Defaults to 1e-12 W/m².

        Returns:
            tuple: A tuple containing:
                - I (float): The sound intensity in W/m².
                - intensity_db (float): The sound intensity in decibels.
        """
        I = acoustic_efficiency * (p_peak * u_peak) / 2.0
        intensity_db = 10 * math.log10(I / I_ref)
        return I, intensity_db

    # Ground Impact Overpressure
    def calculate_overpressure_ground_new(self, D, impact_energy):
        """
        Calculates the peak overpressure for a ground impact.

        Args:
            D (float): The distance from the impact in meters.
            impact_energy (float): The impact energy in Joules.

        Returns:
            float: The peak overpressure in Pascals.
        """
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
        Calculates the peak overpressure for an airburst, with interpolation between different regimes.

        Args:
            D (float): The ground distance from the burst in meters.
            z_b (float): The burst altitude in meters.
            impact_energy (float): The impact energy in Joules.
            z_star (float): The breakup altitude in meters.

        Returns:
            float: The peak overpressure in Pascals.
        """
        # Convert energy to kilotons for scaling laws
        E_kt = impact_energy / 4.184e12
        scaling = E_kt ** (1/3) if E_kt > 0 else 1.0
        
        # Scaled distance and altitude
        D1 = D / scaling if D > 0 else 1e-6 # Avoid division by zero if D is 0
        z_b1 = z_b / scaling if scaling > 0 else z_b # Avoid division by zero if scaling is 0

        if z_b1 == 0: # Ground impact case after scaling
            return self.calculate_overpressure_ground_new(D, impact_energy)

        p_x = 75000.0 # Reference pressure

        # Regime for smaller asteroids (< 25m diameter)
        if self.diameter < 25:
            p_special = 3.14e11 * ((D1**2 + z_b1**2) ** (-1.3)) + 1.8e7 * ((D1**2 + z_b1**2) ** (-0.565))
            D_x = 289.0 + 0.65 * z_b1
            p_mach = (p_x * D_x) / (4 * D1) * (1 + 3 * ((D_x / D1) ** 1.3)) if D1 > 0 else float('inf')
            if z_b1 < 550:
                return p_mach
            else:
                return p_special
        
        # Regime for larger asteroids (> 25m) that break up at high altitude (>= 45km)
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
                band_z = 0.1 * z_b_thresh # This is a 10% band for z_b1
                base_value = max(p_exp, p_special)

                # This condition (z_b1 > z_b_thresh) is z_b1 > 550.
                # Since we are in the 'else' branch where z_b1 < 550, this block is not entered.
                # Keeping it as per original logic, though it seems unreachable here.
                if z_b1 > z_b_thresh: 
                    if z_b1 < z_b_thresh + band_z:
                        t_z = (z_b1 - z_b_thresh) / band_z
                        interp_z = (1 - t_z) * p_mach + t_z * base_value
                        base_value = interp_z
                
                # This condition (z_b1 < z_b_thresh) is z_b1 < 550.
                # This is always true in this 'else' block.
                if z_b1 < z_b_thresh:
                    # D_m1 calculation, ensuring (550 - z_b1) is not zero
                    denominator_D_m1 = 1.2 * (550 - z_b1)
                    D_m1 = (550 * z_b1) / denominator_D_m1 if denominator_D_m1 != 0 else float('inf')
                    
                    if D_m1 == float('inf'):
                        final_value = base_value
                    else:
                        band_D = 0.15 * D_m1 # This is the 15% interpolation band for D1
                        if D1 < D_m1 - band_D:
                            final_value = base_value
                        elif D1 > D_m1 + band_D:
                            final_value = p_mach
                        else:
                            # Ensure 2 * band_D is not zero
                            denominator_t_D = 2 * band_D
                            if denominator_t_D == 0: # Should not happen if D_m1 is finite and > 0
                                t_D = 0.5 # Or handle as an edge case
                            else:
                                t_D = (D1 - (D_m1 - band_D)) / denominator_t_D
                            

                            if self.diameter > 45:
                                final_value = (1 - t_D**2) * base_value + t_D**2 * p_mach
                            else:
                                final_value = (1 - t_D) * base_value + t_D * p_mach
                else: 
                    # This 'else' (z_b1 >= z_b_thresh) is unreachable due to outer 'else' (z_b1 < 550).
                    # Keeping it as per original logic.
                    final_value = base_value
                
                return final_value
        else: # Other conditions: e.g., larger asteroids breaking up at lower altitudes.
            if self.diameter < 45 and z_star >= 46000:
                return 3.14e11 * ((D1**2 + z_b1**2) ** (-1.3)) + 1.8e7 * ((D1**2 + z_b1**2) ** (-0.565))

            p_0 = 3.14e11 * (z_b1**(-2.6)) + 1.8e7 * (z_b1**(-1.13)) if z_b1 > 0 else float('inf')
            beta = 34.87 * (z_b1**(-1.73)) if z_b1 > 0 else float('inf')
            p_exp = p_0 * math.exp(-beta * D1)
            D_x = 289.0 + 0.65 * z_b1
            p_mach = (p_x * D_x) / (4 * D1) * (1 + 3 * ((D_x / D1) ** 1.3)) if D1 > 0 else float('inf')
            z_b_thresh = 550.0
            band_z = 0.1 * z_b_thresh # This is a 10% band for z_b1
            
            base_value = p_exp # Simplified from: p_exp if z_b1 >= z_b_thresh else p_exp

            # Interpolation near the z_b_thresh boundary
            if z_b1 >= z_b_thresh and z_b1 < z_b_thresh + band_z: # e.g. 550 <= z_b1 < 605
                t_z = (z_b1 - z_b_thresh) / band_z
                base_value = (1 - t_z) * p_mach + t_z * p_exp
            
            if z_b1 < z_b_thresh: # z_b1 < 550
                denominator_D_m1 = 1.2 * (550 - z_b1)
                D_m1 = (550 * z_b1) / denominator_D_m1 if denominator_D_m1 != 0 else float('inf')
                if D_m1 == float('inf'):
                    final_value = base_value 
                else:
                    band_D = 0.15 * D_m1 # This is the 15% interpolation band for D1
                    if D1 < D_m1 - band_D:
                        final_value = base_value
                    elif D1 > D_m1 + band_D:
                        final_value = p_mach
                    else:
                        denominator_t_D = 2 * band_D
                        if denominator_t_D == 0:
                             t_D = 0.5
                        else:
                            t_D = (D1 - (D_m1 - band_D)) / denominator_t_D
                        
                        # Interpolation based on diameter
                        if self.diameter > 45:
                            final_value = (1 - t_D**2) * base_value + t_D**2 * p_mach
                        else:
                            final_value = (1 - t_D) * base_value + t_D * p_mach
            else: # z_b1 >= z_b_thresh (z_b1 >= 550)
                final_value = base_value
            
            return final_value

    def calculate_damage_category(self, overpressure_pa: float) -> str:
        """
        Maps a peak overpressure value to a structural damage category.

        Args:
            overpressure_pa (float): The peak overpressure in Pascals.

        Returns:
            str: A description of the expected structural damage.
        """
        for desc, threshold_value in DAMAGE_CATEGORIES:
            if overpressure_pa >= threshold_value:
                return desc
        return "Light or no structural damage"

    def calculate_overpressure_vulnerability(self, press):
        """
        Calculates the vulnerability to overpressure.

        Args:
            press (float): The peak overpressure in Pascals.

        Returns:
            float: The vulnerability score.
        """
        return fun_OverpressureVulnerability(press)

    def calculate_wind_vulnerability(self, v_wind):
        """
        Calculates the vulnerability to high winds.

        Args:
            v_wind (float): The peak wind velocity in m/s.

        Returns:
            float: The vulnerability score.
        """
        return fun_HighWindVulnerability(v_wind)

    def calculate_thermal_vulnerability(self, phi):
        """
        Calculates the vulnerability to thermal radiation.

        Args:
            phi (float): The thermal exposure in J/m².

        Returns:
            float: The vulnerability score.
        """
        return fun_ThermRadVulnerability(phi)

    def calculate_seismic_vulnerability(self, eff_mag):
        """
        Calculates the vulnerability to seismic effects.

        Args:
            eff_mag (float): The effective seismic magnitude.

        Returns:
            float: The vulnerability score.
        """
        return fun_SeismicVulnerability(eff_mag)

    def calculate_ejecta_vulnerability(self, t_ejecta): # t_ejecta in meters
        """
        Calculates the vulnerability to the ejecta blanket.

        Args:
            t_ejecta (float): The thickness of the ejecta blanket in meters.

        Returns:
            float: The vulnerability score.
        """
        return fun_EjectaBlanketVulnerability(t_ejecta)

    def calculate_crater_vulnerability(self, dist_m, final_crater_diam):
        """
        Calculates the vulnerability to cratering effects.

        Args:
            dist_m (float): The distance from the impact in meters.
            final_crater_diam (float): The final crater diameter in meters.

        Returns:
            float: The vulnerability score.
        """
        # pixel_len is no longer used by the simplified fun_CraterVulnerability
        return fun_CraterVulnerability(dist_m, final_crater_diam)

    # Tsunami Calculations
    def calculate_tsunami_effects(self, impact_velocity_at_surface, latitude, longitude):
        """
        Calculates tsunami effects if the impact is in water.
        impact_velocity_at_surface: velocity in m/s at the surface.
        latitude, longitude: coordinates of the impact site.

        Returns:
            dict: A dictionary with tsunami parameters. Returns is_on_land=True
                  if the impact is on land, or an error message if data is missing.
        """
        tsunami_results = {
            "ocean_depth": 0.0,
            "transient_cavity_diameter_water": 0.0,
            "max_amplitude_at_source": 0.0,
            "cavity_radius_water": 0.0,
            "is_on_land": True,
            "error_message": None
        }

        if latitude is None or longitude is None:
            tsunami_results["error_message"] = "Latitude or longitude not provided for tsunami calculation."
            return tsunami_results # Cannot calculate without location

        # Get ocean depth from the GeoTIFF data file.
        ocean_depth = get_ocean_depth_from_geotiff(latitude, longitude)

        if ocean_depth is None: # Critical error reading GeoTIFF
            tsunami_results["error_message"] = "Could not determine ocean depth due to map data error."
            return tsunami_results
        
        tsunami_results["ocean_depth"] = ocean_depth

        if ocean_depth <= 0:
            tsunami_results["is_on_land"] = True
            # No error message needed if on land, it's a valid scenario
            return tsunami_results
        
        tsunami_results["is_on_land"] = False

        # Calculate transient cavity diameter in water
        # Using the existing method, but with water density as target density
        # g_E_crater is used from utils (9.8)
        d_tc_water = self.calculate_transient_crater_diameter(
            v_i=impact_velocity_at_surface,
            rho_t=WATER_DENSITY_CONSTANT, # Target is water
            theta_deg=self.entry_angle_deg, # Use the original impact angle
            g=g_E_crater
        )
        tsunami_results["transient_cavity_diameter_water"] = d_tc_water

        if d_tc_water <= 0:
            # This case implies no significant cavity, so no wave.
            # max_amplitude_at_source will remain 0.0
            return tsunami_results

        # The initial wave amplitude is limited by the ocean depth.
        max_amplitude_at_source = min(0.14 * d_tc_water, ocean_depth)
        cavity_radius_water = d_tc_water / 2.0

        tsunami_results["max_amplitude_at_source"] = max_amplitude_at_source
        tsunami_results["cavity_radius_water"] = cavity_radius_water
        
        return tsunami_results

    def calculate_tsunami_amplitude_at_distance(self, max_amplitude_at_source, d_tc_water, distance_m):
        """
        Calculates the tsunami wave amplitude at a specific distance from the source.
        Assumes d_tc_water > 0 and ocean_depth > 0 and max_amplitude_at_source has been calculated.
        """
        if max_amplitude_at_source <= 0 or d_tc_water <= 0:
            return 0.0
        if distance_m <= 0:
             # Should not happen if called with r_distance > 0
            return max_amplitude_at_source # Or raise error

        cavity_radius_water = d_tc_water / 2.0

        # Within the initial cavity, amplitude is at its maximum.
        if distance_m <= cavity_radius_water:
            return max_amplitude_at_source
        else: # Outside the cavity, amplitude decays.
            # Ensure distance_m is not zero to prevent division by zero, though d_tc_water could also be zero
            if distance_m == 0: return max_amplitude_at_source # technically at source
            return max_amplitude_at_source * (d_tc_water / (2 * distance_m))