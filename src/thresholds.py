"""
Physical constants and damage thresholds for impact effects.
"""
from src.translation_utils import get_translation

# ==========================================
# Airblast / Overpressure (Pascals)
# ==========================================
def get_blast_thresholds():
    return [
        (get_translation("thresholds.blast.highwayGirderBridges", "Highway girder bridges collapse"), 426000),
        (get_translation("thresholds.blast.multistorySteelFramed", "Multistory steel-framed buildings collapse"), 297000),
        (get_translation("thresholds.blast.highwayTrussBridges", "Highway truss bridges distortion"), 100000),
        (get_translation("thresholds.blast.multistoryWallBearing", "Multistory wall-bearing buildings severe damage"), 42600),
        (get_translation("thresholds.blast.woodFrameBuildings", "Wood frame buildings severe damage"), 26800),
        (get_translation("thresholds.blast.glassWindows", "Glass windows shatter"), 6900)
    ]

DAMAGE_CATEGORIES = get_blast_thresholds()

# ==========================================
# Seismic (Richter)
# ==========================================
def get_seismic_thresholds():
    return [
        (get_translation(f"thresholds.seismic.richter{i}to{i+1}", f"{i}-{i+1} Richter"), float(i)) 
        for i in range(10, -1, -1)
    ]

SEISMIC_THRESHOLDS = get_seismic_thresholds()

# ==========================================
# Ejecta (Thickness in meters)
# ==========================================
def get_ejecta_thresholds():
    return [
        (get_translation("thresholds.ejecta.over100m", ">100 m"), 100),
        (get_translation("thresholds.ejecta.over50m", ">50 m"), 50),
        (get_translation("thresholds.ejecta.over25m", ">25 m"), 25),
        (get_translation("thresholds.ejecta.over10m", ">10 m"), 10),
        (get_translation("thresholds.ejecta.over1m", ">1 m"), 1),
        (get_translation("thresholds.ejecta.over01m", ">0.1 m"), 0.1)
    ]

EJECTA_THRESHOLDS = get_ejecta_thresholds()

# ==========================================
# Thermal Radiation (MJ/m²)
# ==========================================
def get_thermal_thresholds():
    return [
        (get_translation("thresholds.thermal.clothingIgnition", "Clothing ignition"), 1.0),
        (get_translation("thresholds.thermal.plywoodIgnition", "Plywood ignition"), 0.67),
        (get_translation("thresholds.thermal.thirdDegreeBurns", "Third degree burns"), 0.42),
        (get_translation("thresholds.thermal.grassIgnition", "Grass ignition"), 0.38),
        (get_translation("thresholds.thermal.newspaperIgnition", "Newspaper ignition"), 0.33),
        (get_translation("thresholds.thermal.deciduousTreesIgnition", "Deciduous trees ignition"), 0.25),
        (get_translation("thresholds.thermal.secondDegreeBurns", "Second degree burns"), 0.25),
        (get_translation("thresholds.thermal.firstDegreeBurns", "First degree burns"), 0.13)
    ]

THERMAL_THRESHOLDS = get_thermal_thresholds()

def get_selected_thermal_thresholds():
    """Subset of thermal thresholds for primary reporting."""
    return [
        (get_translation("thresholds.thermal.clothingIgnition", "Clothing ignition"), 1.0),
        (get_translation("thresholds.thermal.thirdDegreeBurns", "Third degree burns"), 0.42),
        (get_translation("thresholds.thermal.deciduousTreesAndSecondDegree", "Deciduous trees & 2nd deg burns"), 0.25),
        (get_translation("thresholds.thermal.firstDegreeBurns", "First degree burns"), 0.13)
    ]

SELECTED_THERMAL_THRESHOLDS = get_selected_thermal_thresholds()

# ==========================================
# Vulnerability Cutoffs (0.0 to 1.0 scaling)
# ==========================================
OVERPRESSURE_VULNERABILITY_THRESHOLDS = { "no_effect": 150000.0, "complete_destruction": 900000.0 }
THERMAL_VULNERABILITY_THRESHOLD = 85000.0  # W/m²
WIND_VULNERABILITY_THRESHOLDS = { "no_effect": 10.0, "complete_destruction": 250.0 } # m/s
SEISMIC_VULNERABILITY_THRESHOLD = 4.5  # Richter
EJECTA_VULNERABILITY_THRESHOLD = 0.001 # m

# ==========================================
# Wind Scale (Enhanced Fujita)
# ==========================================
def get_ef_wind_thresholds():
    return [
        (get_translation("thresholds.wind.ef5", "EF5 - Incredible Damage"), (89, float('inf'))),
        (get_translation("thresholds.wind.ef4", "EF4 - Devastating Damage"), (74, 89)),
        (get_translation("thresholds.wind.ef3", "EF3 - Severe Damage"), (60, 74)),
        (get_translation("thresholds.wind.ef2", "EF2 - Considerable Damage"), (49, 60)),
        (get_translation("thresholds.wind.ef1", "EF1 - Moderate Damage"), (38, 49)),
        (get_translation("thresholds.wind.ef0", "EF0 - Light Damage"), (29, 38))
    ]

EF_WIND_THRESHOLDS = get_ef_wind_thresholds()

# ==========================================
# Tsunami (Amplitude in m)
# ==========================================
def get_tsunami_amplitude_thresholds():
    return [
        (get_translation("thresholds.tsunami.over100m", ">100m"), 100.0),
        (get_translation("thresholds.tsunami.over10m", ">10m"), 10.0),
        (get_translation("thresholds.tsunami.over1m", ">1m"), 1.0)
    ]

TSUNAMI_AMPLITUDE_THRESHOLDS = get_tsunami_amplitude_thresholds()

# ==========================================
# Helpers
# ==========================================
def get_wind_damage_category(wind_speed_mps):
    for category_name, (lower_bound, _) in get_ef_wind_thresholds():
        if wind_speed_mps >= lower_bound:
            return category_name.split('(')[0].strip()
    return get_translation("thresholds.fallbackMessages.lightOrNoWind", "Light/No Wind Damage (<EF0)")

def get_ejecta_damage_category(thickness_m):
    for category_name, threshold in get_ejecta_thresholds():
        if thickness_m >= threshold:
            return category_name
    return get_translation("thresholds.fallbackMessages.negligibleEjecta", "Negligible (<0.1m)")

def get_thermal_damage_category(phi_J_m2):
    phi_MJ_m2 = phi_J_m2 / 1e6
    # Sort specifically for this logic check
    sorted_thresholds = sorted(get_selected_thermal_thresholds(), key=lambda x: x[1], reverse=True)
    
    for category_name, threshold in sorted_thresholds:
        if phi_MJ_m2 >= threshold:
            return category_name
    return get_translation("thresholds.fallbackMessages.noSignificantThermal", "No significant thermal effects")
