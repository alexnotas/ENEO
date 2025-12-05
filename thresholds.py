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
from translation_utils import get_translation

# Airblast/Overpressure Thresholds (in Pascals)
# These values represent the peak overpressure required to cause specific types
# of structural damage. The list is ordered from most severe to least severe damage.
def get_blast_thresholds():
    """Get localized blast threshold descriptions."""
    return [
        (get_translation("thresholds.blast.highwayGirderBridges", "Highway girder bridges will collapse; Vehicles will be largely displaced and grossly distorted and will require rebuilding before use"), 426000),
        (get_translation("thresholds.blast.multistorySteelFramed", "Multistory steel-framed office-type buildings will suffer extreme frame distortion, incipient collapse"), 297000),
        (get_translation("thresholds.blast.highwayTrussBridges", "Highway truss bridges will suffer substantial distortion of bracing;may collapse"), 100000),
        (get_translation("thresholds.blast.multistoryWallBearing", "Multistory wall-bearing buildings will experience severe damage;may collapse"), 42600),
        (get_translation("thresholds.blast.woodFrameBuildings", "Wood frame buildings will experience severe damage;may collapse"), 26800),
        (get_translation("thresholds.blast.glassWindows", "Glass windows shatter"), 6900)
    ]

# Dynamic threshold lists that update when language changes
def get_blast_thresholds_static():
    """Get blast thresholds - use this for static access."""
    return get_blast_thresholds()

BLAST_THRESHOLDS = get_blast_thresholds()

# A duplicate of BLAST_THRESHOLDS, likely for specific damage assessment logic.
# It's recommended to consolidate these if their purpose is identical.
def get_damage_categories():
    """Get damage categories - same as blast thresholds."""
    return get_blast_thresholds()

DAMAGE_CATEGORIES = get_blast_thresholds()

# Seismic Thresholds (Richter scale magnitude)
# Defines categories for earthquake intensity based on the Richter scale.
# The list is ordered from highest magnitude to lowest.
def get_seismic_thresholds():
    """Get localized seismic threshold descriptions."""
    return [
        (get_translation("thresholds.seismic.richter10Plus", "10+ Richter"), 10.0),
        (get_translation("thresholds.seismic.richter9to10", "9-10 Richter"), 9.0),
        (get_translation("thresholds.seismic.richter8to9", "8-9 Richter"), 8.0),
        (get_translation("thresholds.seismic.richter7to8", "7-8 Richter"), 7.0),
        (get_translation("thresholds.seismic.richter6to7", "6-7 Richter"), 6.0),
        (get_translation("thresholds.seismic.richter5to6", "5-6 Richter"), 5.0),
        (get_translation("thresholds.seismic.richter4to5", "4-5 Richter"), 4.0),
        (get_translation("thresholds.seismic.richter3to4", "3-4 Richter"), 3.0),
        (get_translation("thresholds.seismic.richter2to3", "2-3 Richter"), 2.0),
        (get_translation("thresholds.seismic.richter1to2", "1-2 Richter"), 1.0),
        (get_translation("thresholds.seismic.richter0to1", "0-1 Richter"), 0.0)
    ]

def get_seismic_thresholds_static():
    """Get seismic thresholds - use this for static access."""
    return get_seismic_thresholds()

SEISMIC_THRESHOLDS = get_seismic_thresholds()

# Ejecta Thresholds (thickness in meters)
# Defines categories for the severity of ejecta deposition based on its thickness.
# The list is ordered from thickest to thinnest deposition.
def get_ejecta_thresholds():
    """Get localized ejecta threshold descriptions."""
    return [
        (get_translation("thresholds.ejecta.over100m", ">100 m"), 100),
        (get_translation("thresholds.ejecta.over50m", ">50 m"), 50),
        (get_translation("thresholds.ejecta.over25m", ">25 m"), 25),
        (get_translation("thresholds.ejecta.over10m", ">10 m"), 10),
        (get_translation("thresholds.ejecta.over1m", ">1 m"), 1),
        (get_translation("thresholds.ejecta.over01m", ">0.1 m"), 0.1)
    ]

def get_ejecta_thresholds_static():
    """Get ejecta thresholds - use this for static access."""
    return get_ejecta_thresholds()

EJECTA_THRESHOLDS = get_ejecta_thresholds()

# Thermal Radiation Thresholds (normalized to energy flux in MJ/m² for a 1 Mt event)
# These values represent the thermal energy required to cause specific effects,
# such as ignition of materials or burns.
def get_thermal_thresholds():
    """Get localized thermal threshold descriptions."""
    return [
        (get_translation("thresholds.thermal.clothingIgnition", "Clothing ignition"), 1.0),
        (get_translation("thresholds.thermal.plywoodIgnition", "Plywood ignition"), 0.67),
        (get_translation("thresholds.thermal.grassIgnition", "Grass ignition"), 0.38),
        (get_translation("thresholds.thermal.newspaperIgnition", "Newspaper ignition"), 0.33),
        (get_translation("thresholds.thermal.deciduousTreesIgnition", "Deciduous trees ignition"), 0.25),
        (get_translation("thresholds.thermal.thirdDegreeBurns", "Third degree burns"), 0.42),
        (get_translation("thresholds.thermal.secondDegreeBurns", "Second degree burns"), 0.25),
        (get_translation("thresholds.thermal.firstDegreeBurns", "First degree burns"), 0.13)
    ]

def get_thermal_thresholds_static():
    """Get thermal thresholds - use this for static access."""
    return get_thermal_thresholds()

THERMAL_THRESHOLDS = get_thermal_thresholds()

# A curated list of thermal thresholds used for primary results reporting.
# This list focuses on key indicators of damage and human harm.
def get_selected_thermal_thresholds():
    """Get localized selected thermal threshold descriptions."""
    return [
        (get_translation("thresholds.thermal.clothingIgnition", "Clothing ignition"), 1.0),
        (get_translation("thresholds.thermal.thirdDegreeBurns", "Third degree burns"), 0.42),
        (get_translation("thresholds.thermal.deciduousTreesAndSecondDegree", "Deciduous trees ignition & Second degree burns"), 0.25),
        (get_translation("thresholds.thermal.firstDegreeBurns", "First degree burns"), 0.13)
    ]

def get_selected_thermal_thresholds_static():
    """Get selected thermal thresholds - use this for static access."""
    return get_selected_thermal_thresholds()

SELECTED_THERMAL_THRESHOLDS = get_selected_thermal_thresholds()

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
def get_ef_wind_thresholds():
    """Get localized Enhanced Fujita wind threshold descriptions."""
    return [
        (get_translation("thresholds.wind.ef5IncredibleDamage", "EF5 - Incredible Damage (Severe general destruction)"), (89, float('inf'))),  # Using float('inf') for >89 m/s
        (get_translation("thresholds.wind.ef4DevastatingDamage", "EF4 - Devastating Damage (Institutional buildings severely damaged, all family residence walls collapse)"), (74, 89)),
        (get_translation("thresholds.wind.ef3SevereDamage", "EF3 - Severe Damage (Metal towers collapse, most walls in homes collapse)"), (60, 74)),
        (get_translation("thresholds.wind.ef2ConsiderableDamage", "EF2 - Considerable Damage (Trees debark, wooden towers break, homes shift off foundation)"), (49, 60)),
        (get_translation("thresholds.wind.ef1ModerateDamage", "EF1 - Moderate Damage (Tree trunks snap, windows break, facade tears off)"), (38, 49)),
        (get_translation("thresholds.wind.ef0LightDamage", "EF0 - Light Damage (Large tree branches broken, strip mall roofs begin to uplift)"), (29, 38))
    ]

def get_ef_wind_thresholds_static():
    """Get Enhanced Fujita wind thresholds - use this for static access."""
    return get_ef_wind_thresholds()

EF_WIND_THRESHOLDS = get_ef_wind_thresholds()

# Tsunami Amplitude Thresholds (in meters)
# Defines minimum wave amplitudes for categorizing tsunami danger levels.
# The simulation reports the distance at which the wave amplitude exceeds these values.
def get_tsunami_amplitude_thresholds():
    """Get localized tsunami amplitude threshold descriptions."""
    return [
        (get_translation("thresholds.tsunami.over100m", ">100m"), 100.0),
        (get_translation("thresholds.tsunami.over10m", ">10m"), 10.0),
        (get_translation("thresholds.tsunami.over1m", ">1m"), 1.0)
    ]

def get_tsunami_amplitude_thresholds_static():
    """Get tsunami amplitude thresholds - use this for static access."""
    return get_tsunami_amplitude_thresholds()

TSUNAMI_AMPLITUDE_THRESHOLDS = get_tsunami_amplitude_thresholds()

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
    for category_name, (lower_bound_mps, _) in get_ef_wind_thresholds():
        if wind_speed_mps >= lower_bound_mps:
            # Returns the descriptive part like "EF5 - Incredible Damage"
            return category_name.split('(')[0].strip()
    return get_translation("thresholds.fallbackMessages.lightOrNoWind", "Light or no wind damage (Below EF0)")

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
    for category_name, threshold_val_m in get_ejecta_thresholds():
        if thickness_m >= threshold_val_m:
            return category_name
    return get_translation("thresholds.fallbackMessages.negligibleEjecta", "Negligible ejecta (less than 0.1 m)")

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
    sorted_thermal_thresholds = sorted(get_selected_thermal_thresholds(), key=lambda x: x[1], reverse=True)
    phi_MJ_m2 = phi_J_m2 / 1e6  # Convert Joules/m² to MegaJoules/m² for comparison.

    for category_name, threshold_MJ_m2 in sorted_thermal_thresholds:
        if phi_MJ_m2 >= threshold_MJ_m2:
            return category_name
    return get_translation("thresholds.fallbackMessages.noSignificantThermal", "No significant thermal effects")