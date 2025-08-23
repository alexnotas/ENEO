import math

# Airblast/Overpressure Thresholds (Pascal)
BLAST_THRESHOLDS = [
    ("Highway girder bridges will collapse; Vehicles will be largely displaced and grossly distorted and will require rebuilding before use", 426000),
    ("Multistory steel-framed office-type buildings will suffer extreme frame distortion, incipient collapse", 297000),
    ("Highway truss bridges will suffer substantial distortion of bracing;may collapse", 100000),
    ("Multistory wall-bearing buildings will experience severe damage;may collapse", 42600),
    ("Wood frame buildings will experience severe damage;may collapse", 26800),
    ("Glass windows shatter", 6900)
]

DAMAGE_CATEGORIES = [
    ("Highway girder bridges will collapse; Vehicles will be largely displaced and grossly distorted and will require rebuilding before use", 426000),
    ("Multistory steel-framed office-type buildings will suffer extreme frame distortion, incipient collapse", 297000),
    ("Highway truss bridges will suffer substantial distortion of bracing;may collapse", 100000),
    ("Multistory wall-bearing buildings will experience severe damage;may collapse", 42600),
    ("Wood frame buildings will experience severe damage;may collapse", 26800),
    ("Glass windows shatter", 6900)
]

# Seismic Thresholds (Richter scale)
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

# Ejecta Thresholds (meters)
EJECTA_THRESHOLDS = [
    (">100 m", 100),
    (">50 m", 50),
    (">25 m", 25),
    (">10 m", 10),
    (">1 m", 1),
    (">0.1 m", 0.1)
]

# Thermal Radiation Thresholds (MJ/m² for 1 Mt event)
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

SELECTED_THERMAL_THRESHOLDS = [
    ("Clothing ignition", 1.0),
    ("Third degree burns", 0.42),
    ("Deciduous trees ignition & Second degree burns", 0.25),
    ("First degree burns", 0.13)
]

# Vulnerability Thresholds

# Overpressure Vulnerability
OVERPRESSURE_VULNERABILITY_THRESHOLDS = {
    "no_effect": 150000.0,  # Pa, below this threshold, vulnerability is 0
    "complete_destruction": 900000.0  # Pa, above this threshold, vulnerability is 1
}

# Thermal Radiation Vulnerability
THERMAL_VULNERABILITY_THRESHOLD = 85000.0  # W/m², below this threshold, vulnerability is 0

# High Wind Vulnerability
WIND_VULNERABILITY_THRESHOLDS = {
    "no_effect": 10.0,  # m/s, below this threshold, vulnerability is 0
    "complete_destruction": 250.0  # m/s, above this threshold, vulnerability is 1
}

# Seismic Vulnerability
SEISMIC_VULNERABILITY_THRESHOLD = 4.5  # Richter, below this threshold, vulnerability is 0

# Ejecta Blanket Vulnerability
EJECTA_VULNERABILITY_THRESHOLD = 0.001  # m, below this threshold, vulnerability is 0

# Enhanced Fujita Wind Scale Thresholds (m/s)
# Format: ("Description", (min_wind_speed_mps, max_wind_speed_mps))
# Max wind speed for EF5 is effectively unbounded (>89 m/s).
# Ordered from highest wind speed (most severe) to lowest.
EF_WIND_THRESHOLDS = [
    ("EF5 - Incredible Damage (Severe general destruction)", (89, float('inf'))),  # Using float('inf') for >89 m/s
    ("EF4 - Devastating Damage (Institutional buildings severely damaged, all family residence walls collapse)", (74, 89)),
    ("EF3 - Severe Damage (Metal towers collapse, most walls in homes collapse)", (60, 74)),
    ("EF2 - Considerable Damage (Trees debark, wooden towers break, homes shift off foundation)", (49, 60)),
    ("EF1 - Moderate Damage (Tree trunks snap, windows break, facade tears off)", (38, 49)),
    ("EF0 - Light Damage (Large tree branches broken, strip mall roofs begin to uplift)", (29, 38))
]

# Tsunami Amplitude Thresholds (meters)
# Description: (Display Name, Amplitude in meters)
# These define the minimum amplitude for a zone to be considered "dangerous" at that level.
# The simulation will report the distance up to which the wave EXCEEDS this amplitude.
TSUNAMI_AMPLITUDE_THRESHOLDS = [
    (">100m", 100.0),
    (">10m", 10.0),
    (">1m", 1.0)
]

# Helper functions to map effect values to categories

def get_wind_damage_category(wind_speed_mps):
    """Maps wind speed in m/s to an EF scale category description."""
    # EF_WIND_THRESHOLDS is assumed to be sorted from EF5 (most severe) down to EF0
    for category_name, (lower_bound_mps, _) in EF_WIND_THRESHOLDS:
        if wind_speed_mps >= lower_bound_mps:
            # Returns the descriptive part like "EF5 - Incredible Damage"
            return category_name.split('(')[0].strip()
    return "Light or no wind damage (Below EF0)"

def get_ejecta_damage_category(thickness_m):
    """Maps ejecta thickness in meters to a category description."""
    # EJECTA_THRESHOLDS is assumed to be sorted from thickest to thinnest
    for category_name, threshold_val_m in EJECTA_THRESHOLDS:
        if thickness_m >= threshold_val_m:
            return category_name
    return "Negligible ejecta (less than 0.1 m)"

def get_thermal_damage_category(phi_J_m2):
    """Maps thermal flux in J/m^2 to a thermal effect category description."""
    # SELECTED_THERMAL_THRESHOLDS are (description, threshold_MJ_m2_for_1MT)
    # We assume the critical flux for an effect is constant.
    # Sort by threshold value (MJ/m^2) descending to get the most severe match first.
    sorted_thermal_thresholds = sorted(SELECTED_THERMAL_THRESHOLDS, key=lambda x: x[1], reverse=True)
    phi_MJ_m2 = phi_J_m2 / 1e6  # Convert actual flux to MJ/m^2 for comparison

    for category_name, threshold_MJ_m2 in sorted_thermal_thresholds:
        if phi_MJ_m2 >= threshold_MJ_m2:
            return category_name
    return "No significant thermal effects"