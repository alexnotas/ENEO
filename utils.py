"""
ENEO Asteroid Impact Simulation - Utility Functions and Constants Module

This module provides essential utility functions and physical constants used throughout
the asteroid impact simulation system. It includes:

1. Physical and atmospheric constants for impact calculations
2. Unit conversion utilities (distance, energy)
3. Geometric and mathematical helper functions
4. Ocean depth retrieval from geospatial data
5. Scaling factors for impact effects

The constants and functions are based on established atmospheric models,
impact physics research, and geophysical standards.

Author: Alexandros Notas
Institution: National Technical University of Athens
Date: July 2025
"""

import math
import rasterio # Add this import for geospatial data processing

# =============================================================================
# PHYSICAL AND ATMOSPHERIC CONSTANTS
# =============================================================================

# Atmospheric scale height (meters) - characteristic height over which atmospheric density decreases by factor e
H = 8000.0

# Sea level atmospheric density (kg/m³) - reference density for atmospheric calculations
rho0 = 1.0

# Drag coefficient (dimensionless) - aerodynamic drag coefficient for spherical objects
C_D = 2.0

# Earth's gravitational acceleration in atmosphere (m/s²) - standard gravity for atmospheric entry
g_E_atmos = 9.81

# Earth's gravitational acceleration for crater formation (m/s²) - gravity used in crater scaling laws
g_E_crater = 9.8

# Pancake factor limit (dimensionless) - maximum flattening ratio for fragmenting meteoroids
fp_limit = 7.0

# Target rock density (kg/m³) - typical density of crustal rock material
rho_target = 2600.0

# Standard atmospheric pressure at sea level (Pascal) - reference pressure for blast calculations
P0 = 101325.0

# Speed of sound in air at sea level (m/s) - used for acoustic and blast wave calculations
c0 = 343.0

# Acoustic efficiency factor (dimensionless) - fraction of impact energy converted to acoustic waves
acoustic_efficiency = 1e-8

# Earth's mean radius (meters) - used for great circle distance calculations and geometric corrections
R_EARTH = 6371000.0

# Burst altitude threshold (meters) - critical altitude distinguishing ground vs. airburst effects
BURST_ALTITUDE_THRESHOLD = 1000.0

# =============================================================================
# GEOSPATIAL DATA PATHS AND CONSTANTS
# =============================================================================

# Path to ETOPO elevation/bathymetry data file - global topographic dataset
ETOPO_FILE_PATH = "/home/debian/flaskapp/maps/ETOPO_2022_v1_60s_N90W180_surface.tif"

# Water density constant (kg/m³) - standard density of seawater for tsunami calculations
WATER_DENSITY_CONSTANT = 1000  # kg/m^3
H = 8000.0
rho0 = 1.0
C_D = 2.0
g_E_atmos = 9.81
g_E_crater = 9.8
fp_limit = 7.0
rho_target = 2600.0
P0 = 101325.0
c0 = 343.0
acoustic_efficiency = 1e-8
R_EARTH = 6371000.0
BURST_ALTITUDE_THRESHOLD = 1000.0

ETOPO_FILE_PATH = "/home/debian/flaskapp/maps/ETOPO_2022_v1_60s_N90W180_surface.tif"
WATER_DENSITY_CONSTANT = 1000  # kg/m^3

# Transition settings for altitude burst vs. regular reflection
BURST_ALTITUDE_THRESHOLD = 1000.0  # meters

# =============================================================================
# UNIT CONVERSION UTILITIES
# =============================================================================

def km_to_m(km):
    """
    Convert kilometers to meters.
    
    Parameters
    ----------
    km : float
        Distance in kilometers
        
    Returns
    -------
    float
        Distance in meters
    """
    return km * 1000.0

def m_to_km(m):
    """
    Convert meters to kilometers.
    
    Parameters
    ----------
    m : float
        Distance in meters
        
    Returns
    -------
    float
        Distance in kilometers
    """
    return m / 1000.0

def convert_energy_j_to_mt(energy_j):
    """
    Convert energy from Joules to Megatons of TNT equivalent.
    
    Uses the standard conversion factor where 1 MT TNT = 4.184 × 10^15 Joules.
    This is the internationally accepted TNT equivalent for nuclear yield measurements.
    
    Parameters
    ----------
    energy_j : float
        Energy in Joules
        
    Returns
    -------
    float
        Energy in Megatons TNT equivalent
    """
    return energy_j / 4.184e15

# =============================================================================
# IMPACT SCALING AND GEOMETRIC FUNCTIONS
# =============================================================================

def compute_scaling_factor(energy_kt):
    """
    Compute cube-root scaling factor for impact effects.
    
    Many impact effects (crater size, blast radius, etc.) scale with the cube root
    of the impact energy. This is based on dimensional analysis and similarity
    scaling laws in explosion physics.
    
    Parameters
    ----------
    energy_kt : float
        Impact energy in kilotons TNT equivalent
        
    Returns
    -------
    float
        Scaling factor (energy_kt^(1/3)) or 1.0 if energy <= 0
    """
    return energy_kt ** (1/3) if energy_kt > 0 else 1.0

def common_p0(r1, z_b1):
    """
    Calculate reference overpressure using empirical blast scaling laws.
    
    This function implements a two-term empirical fit for blast overpressure
    as a function of scaled distance. The formula combines near-field and
    far-field blast decay behavior.
    
    Parameters
    ----------
    r1 : float
        Horizontal distance from blast center (meters)
    z_b1 : float
        Burst altitude above surface (meters)
        
    Returns
    -------
    float
        Reference overpressure in Pascal
        
    Mathematical Form
    -----------------
    P = 3.14×10^11 × (r²+z²)^(-1.3) + 1.8×10^7 × (r²+z²)^(-0.565)
    
    The two terms represent different physical regimes:
    - First term: near-field shock behavior (steep decay)
    - Second term: far-field acoustic behavior (gradual decay)
    """
    return 3.14e11 * ((r1**2 + z_b1**2) ** (-1.3)) + 1.8e7 * ((r1**2 + z_b1**2) ** (-0.565))

def curvature_adjustment_factor(r_distance_m, R_f):
    """
    Calculate geometric adjustment factor for Earth's curvature effects on thermal radiation.
    
    When thermal radiation travels long distances, Earth's curvature can block part of the
    fireball from view. This function calculates the fraction of the fireball visible
    from a given distance, accounting for the spherical geometry of Earth.
    
    Parameters
    ----------
    r_distance_m : float
        Horizontal distance from impact point (meters)
    R_f : float
        Fireball radius (meters)
        
    Returns
    -------
    float
        Visibility fraction (0.0 = completely blocked, 1.0 = fully visible)
        
    Algorithm
    ---------
    1. Calculate angular distance on Earth's surface
    2. Determine horizon height at the observation point
    3. Check if fireball extends above the horizon
    4. Calculate visible fraction using geometric integration
    
    Mathematical Basis
    ------------------
    Uses spherical geometry to determine the intersection of the fireball
    sphere with the line of sight from the observer, accounting for Earth's
    curvature blocking the lower portion of the fireball.
    """
    # Return zero if fireball radius is invalid
    if R_f <= 0:
        return 0.0
    
    # Calculate angular distance on Earth's surface (radians)
    Delta = r_distance_m / R_EARTH
    
    # Calculate height of horizon due to Earth's curvature
    h = (1 - math.cos(Delta)) * R_EARTH
    
    # If horizon is higher than fireball, no thermal radiation reaches observer
    if h >= R_f:
        return 0.0
    
    # Calculate the angle subtended by the visible portion of the fireball
    delta = math.acos(h / R_f)
    
    # Calculate the fraction of fireball visible (geometric integration result)
    # Formula derived from integration over the visible circular segment
    f = (2 / math.pi) * (delta - (h / R_f) * math.sin(delta))
    return f

# =============================================================================
# GEOSPATIAL DATA PROCESSING FUNCTIONS
# =============================================================================

def get_ocean_depth_from_geotiff(lat, lon, file_path=ETOPO_FILE_PATH):
    """
    Retrieve ocean depth from ETOPO bathymetry/elevation data at specified coordinates.
    
    This function reads elevation data from a GeoTIFF file and converts it to ocean depth
    for tsunami calculations. The ETOPO dataset uses the convention where negative values
    represent depths below sea level and positive values represent elevations above sea level.
    
    Parameters
    ----------
    lat : float
        Latitude in decimal degrees (-90 to 90)
    lon : float
        Longitude in decimal degrees (-180 to 180)
    file_path : str, optional
        Path to the ETOPO GeoTIFF file (default: ETOPO_FILE_PATH)
        
    Returns
    -------
    float or None
        Ocean depth in meters (positive value) if location is underwater
        0.0 if location is on land or no data available
        None if coordinates are out of bounds or critical error occurs
        
    Data Convention
    ---------------
    ETOPO data format:
    - Negative values: Below sea level (ocean depth)
    - Positive values: Above sea level (land elevation)
    - NoData values: Areas with no elevation data
    
    Error Handling
    --------------
    The function implements comprehensive error handling for:
    - Out-of-bounds coordinates
    - Missing or corrupted data files
    - NoData values in the dataset
    - Coordinate transformation errors
    """
    try:
        # Open the GeoTIFF file using rasterio for geospatial data access
        with rasterio.open(file_path) as src:
            # Check if coordinates fall within the dataset bounds
            if not (src.bounds.left <= lon <= src.bounds.right and \
                    src.bounds.bottom <= lat <= src.bounds.top):
                print(f"Warning: Coordinates ({lat}, {lon}) are outside the map bounds.")
                return None # Indicate out of bounds

            # Convert geographic coordinates to raster indices
            row, col = src.index(lon, lat)
            
            # Read elevation value from the raster at the specified location
            elevation = src.read(1)[row, col].item()
            
            # Check for NoData values in the dataset
            nodata_val = src.nodatavals[0]
            if nodata_val is not None and elevation == nodata_val:
                print(f"Warning: Impact site at ({lat:.2f}, {lon:.2f}) is in an area with no data. Assuming land.")
                return 0 # No data, assume land for tsunami purposes

            # Convert elevation to ocean depth based on ETOPO convention
            if elevation < 0:
                # Negative elevation = below sea level = ocean depth
                ocean_depth = -elevation
                return ocean_depth
            else:
                # Positive elevation = above sea level = land
                # print(f"Info: Impact site at ({lat:.2f}, {lon:.2f}) is on land (elevation: {elevation:.2f} m).")
                return 0 # On land
                
    except IndexError:
        # Handle coordinate transformation errors or invalid indices
        print(f"Warning: Could not retrieve data for coordinates ({lat:.2f}, {lon:.2f}). Assuming land.")
        return 0 # Error, assume land
    except rasterio.errors.RasterioIOError as e:
        # Handle file I/O errors (missing file, corrupted data, etc.)
        print(f"Error opening or reading GeoTIFF for ocean depth: {e}")
        return None # Critical error
    except Exception as e:
        # Handle any other unexpected errors
        print(f"An unexpected error occurred while getting ocean depth: {e}")
        return None # Critical error