"""
A collection of utility functions and physical constants for asteroid impact simulations.

This module provides a centralized place for various helper functions and constants
that are used across the different calculation modules of the application. It includes:
- Physical and atmospheric constants (e.g., Earth's gravity, atmospheric density).
- Unit conversion functions (e.g., kilometers to meters).
- Mathematical formulas for specific physical calculations (e.g., overpressure,
  curvature adjustment).
- Data lookup functions, such as retrieving ocean depth from a GeoTIFF file.

The constants and functions are designed to be imported and used by other modules
like `models.py` to ensure consistency and avoid code duplication.
"""
import math
import rasterio

# =============================================================================
# GLOBAL PHYSICAL AND ENVIRONMENTAL CONSTANTS
# =============================================================================

# Atmospheric constants
H = 8000.0  # Atmospheric scale height in meters
rho0 = 1.0  # Atmospheric density at sea level (kg/m^3) - simplified value
C_D = 2.0  # Drag coefficient for a spherical object
g_E_atmos = 9.81  # Standard gravity for atmospheric calculations (m/s^2)
P0 = 101325.0  # Standard atmospheric pressure at sea level (Pascals)
c0 = 343.0  # Speed of sound in air at sea level (m/s)

# Asteroid and target properties
g_E_crater = 9.8  # Gravity used for crater calculations (m/s^2)
fp_limit = 7.0  # A dimensionless parameter used in crater scaling laws
rho_target = 2600.0  # Assumed density of the target rock (kg/m^3)

# Energy and effects constants
acoustic_efficiency = 1e-8  # Fraction of impact energy converted to acoustic energy

# Earth-related constants
R_EARTH = 6371000.0  # Mean radius of the Earth in meters

# Simulation thresholds
BURST_ALTITUDE_THRESHOLD = 1000.0  # Altitude (m) below which an airburst is treated differently for reflection

# File paths and data-specific constants
ETOPO_FILE_PATH = "/home/debian/flaskapp/maps/ETOPO_2022_v1_60s_N90W180_surface.tif"
WATER_DENSITY_CONSTANT = 1000  # Density of water in kg/m^3

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def km_to_m(km):
    """Converts kilometers to meters."""
    return km * 1000.0

def m_to_km(m):
    """Converts meters to kilometers."""
    return m / 1000.0

def convert_energy_j_to_mt(energy_j):
    """Converts energy from Joules to Megatons of TNT."""
    return energy_j / 4.184e15

def compute_scaling_factor(energy_kt):
    """
    Computes a scaling factor based on the cube root of energy in kilotons.
    Used in various empirical formulas for impact effects.
    """
    return energy_kt ** (1/3) if energy_kt > 0 else 1.0

def common_p0(r1, z_b1):
    """
    Calculates a baseline overpressure value based on scaled distance and burst height.
    This is an empirical formula used in airblast calculations.

    Args:
        r1 (float): Scaled ground distance.
        z_b1 (float): Scaled burst height.

    Returns:
        float: The calculated baseline overpressure.
    """
    return 3.14e11 * ((r1**2 + z_b1**2) ** (-1.3)) + 1.8e7 * ((r1**2 + z_b1**2) ** (-0.565))

def curvature_adjustment_factor(r_distance_m, R_f):
    """
    Calculates an adjustment factor to account for the Earth's curvature, which
    can shield locations from the direct line of sight of an airburst.

    Args:
        r_distance_m (float): The ground distance from the impact in meters.
        R_f (float): The radius of the fireball in meters.

    Returns:
        float: A factor between 0 and 1, where 0 means fully shielded and 1
               means fully exposed. Returns 0 if the fireball radius is non-positive.
    """
    if R_f <= 0:
        return 0.0
    Delta = r_distance_m / R_EARTH
    h = (1 - math.cos(Delta)) * R_EARTH
    if h >= R_f:
        return 0.0
    delta = math.acos(h / R_f)
    f = (2 / math.pi) * (delta - (h / R_f) * math.sin(delta))
    return f

# =============================================================================
# OCEAN DEPTH AND TSUNAMI HELPERS
# =============================================================================

def get_ocean_depth_from_geotiff(lat, lon, file_path=ETOPO_FILE_PATH):
    """
    Reads the ocean depth from a GeoTIFF file (like ETOPO1) for a given latitude and longitude.

    The ETOPO dataset uses negative values for bathymetry (ocean depth) and
    positive values for topography (land elevation). This function returns the
    depth as a positive number if the location is in the ocean, and 0 if it's on
    land. It includes error handling for out-of-bounds coordinates and file I/O issues.

    Args:
        lat (float): Latitude of the point.
        lon (float): Longitude of the point.
        file_path (str, optional): The path to the GeoTIFF file. Defaults to ETOPO_FILE_PATH.

    Returns:
        float or None: The depth of the ocean in meters (as a positive value),
                       0 if the location is on land, or None if a critical error
                       occurs (e.g., file not found, coordinates out of bounds).
    """
    try:
        with rasterio.open(file_path) as src:
            if not (src.bounds.left <= lon <= src.bounds.right and \
                    src.bounds.bottom <= lat <= src.bounds.top):
                print(f"Warning: Coordinates ({lat}, {lon}) are outside the map bounds.")
                return None  # Indicate out of bounds

            # Get the pixel coordinates for the given lat/lon
            row, col = src.index(lon, lat)
            # Read the single pixel value from the first band
            elevation = src.read(1)[row, col].item()
            
            # Check if the pixel value is the designated 'no data' value
            nodata_val = src.nodatavals[0]
            if nodata_val is not None and elevation == nodata_val:
                print(f"Warning: Impact site at ({lat:.2f}, {lon:.2f}) is in an area with no data. Assuming land.")
                return 0  # No data, assume land for tsunami purposes

            # Negative elevation indicates ocean depth
            if elevation < 0:
                ocean_depth = -elevation
                return ocean_depth
            else:
                # Positive elevation indicates land
                return 0  # On land
    except IndexError:
        # This can happen if coordinates are valid but fall on the exact edge
        print(f"Warning: Could not retrieve data for coordinates ({lat:.2f}, {lon:.2f}). Assuming land.")
        return 0  # Error, assume land
    except rasterio.errors.RasterioIOError as e:
        print(f"Error opening or reading GeoTIFF for ocean depth: {e}")
        return None  # Critical error, indicates a problem with the file itself
    except Exception as e:
        print(f"An unexpected error occurred while getting ocean depth: {e}")
        return None  # Critical error for any other unexpected issue
