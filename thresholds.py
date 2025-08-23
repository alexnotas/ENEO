"""
Defines physical constants and thresholds for various impact-related effects.

This module contains predefined lists and dictionaries that store the thresholds
for damage caused by airblast (overpressure), seismic events, ejecta deposition,
and thermal radiation. These values are used throughout the simulation to
categorize the severity of different physical phenomena and to calculate their
impact on infrastructure and population.

The module also provides helper functions to map a calculated physical value
(e.g., wind speed, ejecta thickness) to a human-readable damage category.
"""
import math

# Airblast/Overpressure Thresholds (in Pascals)
# These values represent the peak overpressure required to cause specific types
# of structural damage. The list is ordered from most severe to least severe damage.
BLAST_THRESHOLDS = [
    ("Highway girder bridges will collapse; Vehicles will be largely displaced and grossly distorted and will require rebuilding before use", 426000),
    ("Multistory steel-framed office-type buildings will suffer extreme frame distortion, incipient collapse", 297000),
    ("Highway truss bridges will suffer substantial distortion of bracing;may collapse", 100000),
    ("Multistory wall-bearing buildings will experience severe damage;may collapse", 42600),
    ("Wood frame buildings will experience severe damage;may collapse", 26800),
    ("Glass windows shatter", 6900)
]

# A duplicate of BLAST_THRESHOLDS, likely for specific damage assessment logic.
# It's recommended to consolidate these if their purpose is identical.
DAMAGE_CATEGORIES = [
    ("Highway girder bridges will collapse; Vehicles will be largely displaced and grossly distorted and will require rebuilding before use", 426000),
    ("Multistory steel-framed office-type buildings will suffer extreme frame distortion, incipient collapse", 297000),
    ("Highway truss bridges will suffer substantial distortion of bracing;may collapse", 100000),
    ("Multistory wall-bearing buildings will experience severe damage;may collapse", 42600),
    ("Wood frame buildings will experience severe damage;may collapse", 26800),
    ("Glass windows shatter", 6900)
]

# Seismic Thresholds (Richter scale magnitude)
# Defines categories for earthquake intensity based on the Richter scale.
# The list is ordered from highest magnitude to lowest.
SEISMIC_THRESHOLDS = [
    ("10+ Richter", 10.0),
    ("9-10 Richter", 9.0),
    ("8-9 Richter", 8.0),
    ("7-8 Richter", 7.0),
    ("6-7 Richter", 6.0),
    ("5-6 Richter", 5.0),
    ("4-5 Richter", 4.0),
    ("3-4 Richter", 3.0),
    ("2-3 Richter", 2.0),
    ("1-2 Richter", 1.0),
    ("0-1 Richter", 0.0)
]

# Ejecta Thresholds (thickness in meters)
# Defines categories for the severity of ejecta deposition based on its thickness.
# The list is ordered from thickest to thinnest deposition.
EJECTA_THRESHOLDS = [
    (">100 m", 100),
    (">50 m", 50),
    (">25 m", 25),
    (">10 m", 10),
    (">1 m", 1),
    (">0.1 m", 0.1)
]

# Thermal Radiation Thresholds (normalized to energy flux in MJ/m² for a 1 Mt event)
# These values represent the thermal energy required to cause specific effects,
# such as ignition of materials or burns.
THERMAL_THRESHOLDS = [
    ("Clothing ignition", 1.0),
    ("Plywood ignition", 0.67),
    ("Grass ignition", 0.38),
    ("Newspaper ignition", 0.33),
    ("Deciduous trees ignition", 0.25),
    ("Third degree burns", 0.42),
    ("Second degree burns", 0.25),
    ("First degree burns", 0.13)
]

# A curated list of thermal thresholds used for primary results reporting.
# This list focuses on key indicators of damage and human harm.
SELECTED_THERMAL_THRESHOLDS = [
    ("Clothing ignition", 1.0),
    ("Third degree burns", 0.42),
    ("Deciduous trees ignition & Second degree burns", 0.25),
    ("First degree burns", 0.13)
]

# Vulnerability Thresholds
# These constants define the lower and upper bounds for calculating a vulnerability
# score (from 0 to 1) for different physical effects. A score of 0 means no
# effect, while a score of 1 means complete destruction or maximum effect.

# Overpressure Vulnerability (in Pascals)
# Defines the pressure range over which vulnerability increases from 0 to 1.
OVERPRESSURE_VULNERABILITY_THRESHOLDS = {
    "no_effect": 150000.0,  # Pa, below this threshold, vulnerability is 0
    "complete_destruction": 900000.0  # Pa, above this threshold, vulnerability is 1
}

# Thermal Radiation Vulnerability (in Watts/m²)
# Defines the minimum thermal flux below which vulnerability is considered 0.
THERMAL_VULNERABILITY_THRESHOLD = 85000.0  # W/m², below this threshold, vulnerability is 0

# High Wind Vulnerability (in m/s)
# Defines the wind speed range over which vulnerability increases from 0 to 1.
WIND_VULNERABILITY_THRESHOLDS = {
    "no_effect": 10.0,  # m/s, below this threshold, vulnerability is 0
    "complete_destruction": 250.0  # m/s, above this threshold, vulnerability is 1
}

# Seismic Vulnerability (Richter scale magnitude)
# Defines the minimum earthquake magnitude below which vulnerability is 0.
SEISMIC_VULNERABILITY_THRESHOLD = 4.5  # Richter, below this threshold, vulnerability is 0

# Ejecta Blanket Vulnerability (thickness in meters)
# Defines the minimum ejecta thickness below which vulnerability is 0.
EJECTA_VULNERABILITY_THRESHOLD = 0.001  # m, below this threshold, vulnerability is 0

# Enhanced Fujita Wind Scale Thresholds (in m/s)
# Maps wind speeds to the Enhanced Fujita (EF) scale for tornado/wind damage.
# Format: ("Description", (min_wind_speed_mps, max_wind_speed_mps))
# The list is ordered from most severe (EF5) to least severe (EF0).
EF_WIND_THRESHOLDS = [
    ("EF5 - Incredible Damage (Severe general destruction)", (89, float('inf'))),  # Using float('inf') for >89 m/s
    ("EF4 - Devastating Damage (Institutional buildings severely damaged, all family residence walls collapse)", (74, 89)),
    ("EF3 - Severe Damage (Metal towers collapse, most walls in homes collapse)", (60, 74)),
    ("EF2 - Considerable Damage (Trees debark, wooden towers break, homes shift off foundation)", (49, 60)),
    ("EF1 - Moderate Damage (Tree trunks snap, windows break, facade tears off)", (38, 49)),
    ("EF0 - Light Damage (Large tree branches broken, strip mall roofs begin to uplift)", (29, 38))
]

# Tsunami Amplitude Thresholds (in meters)
# Defines minimum wave amplitudes for categorizing tsunami danger levels.
# The simulation reports the distance at which the wave amplitude exceeds these values.
TSUNAMI_AMPLITUDE_THRESHOLDS = [
    (">100m", 100.0),
    (">10m", 10.0),
    (">1m", 1.0)
]

# Helper functions to map effect values to categories

def get_wind_damage_category(wind_speed_mps):
    """
    Maps a given wind speed to its corresponding Enhanced Fujita (EF) scale category.

    Args:
        wind_speed_mps (float): The wind speed in meters per second.

    Returns:
        str: A string describing the EF scale damage category (e.g., "EF5 - Incredible Damage").
             Returns a message for light or no damage if the speed is below EF0.
    """
    # EF_WIND_THRESHOLDS is assumed to be sorted from EF5 (most severe) down to EF0
    for category_name, (lower_bound_mps, _) in EF_WIND_THRESHOLDS:
        if wind_speed_mps >= lower_bound_mps:
            # Returns the descriptive part like "EF5 - Incredible Damage"
            return category_name.split('(')[0].strip()
    return "Light or no wind damage (Below EF0)"

def get_ejecta_damage_category(thickness_m):
    """
    Maps a given ejecta thickness to a descriptive damage category.

    Args:
        thickness_m (float): The ejecta thickness in meters.

    Returns:
        str: A string describing the ejecta category (e.g., ">10 m").
             Returns a message for negligible ejecta if the thickness is below the lowest threshold.
    """
    # EJECTA_THRESHOLDS is assumed to be sorted from thickest to thinnest
    for category_name, threshold_val_m in EJECTA_THRESHOLDS:
        if thickness_m >= threshold_val_m:
            return category_name
    return "Negligible ejecta (less than 0.1 m)"

def get_thermal_damage_category(phi_J_m2):
    """
    Maps a thermal energy flux to its most severe corresponding damage category.

    Args:
        phi_J_m2 (float): The thermal energy flux in Joules per square meter (J/m²).

    Returns:
        str: A string describing the most severe thermal effect (e.g., "Clothing ignition").
             Returns a message for no significant effects if the flux is below all thresholds.
    """
    # SELECTED_THERMAL_THRESHOLDS are (description, threshold_MJ_m2_for_1MT)
    # We assume the critical flux for an effect is constant.
    sorted_thermal_thresholds = sorted(SELECTED_THERMAL_THRESHOLDS, key=lambda x: x[1], reverse=True)
    phi_MJ_m2 = phi_J_m2 / 1e6  # Convert Joules/m² to MegaJoules/m² for comparison.

    for category_name, threshold_MJ_m2 in sorted_thermal_thresholds:
        if phi_MJ_m2 >= threshold_MJ_m2:
            return category_name
    return "No significant thermal effects"