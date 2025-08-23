import os
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.mask import mask
from pyproj import Transformer, CRS
from shapely.geometry import Point, MultiPolygon, Polygon
from shapely.ops import transform
from shapely.geometry.base import BaseGeometry
from shapely.validation import make_valid
import math  # Ensure math is imported if create_geodesic_buffer uses it directly
from collections import defaultdict # Added for defaultdict

# Global variables to hold the loaded datasets
TIFF_PATH = "/home/alexandros-linux/Επιφάνεια εργασίας/ENEO FINAL FOR GIT/flaskapp/maps/population_with_country_fid_assigned.tif"
COUNTRY_NAMES_PATH = "/home/alexandros-linux/Επιφάνεια εργασίας/ENEO FINAL FOR GIT/flaskapp/maps/country_fid_lookup.csv"  # Path to FID-to-name lookup if available
WORLD_SHAPEFILE_PATH = "/home/alexandros-linux/Επιφάνεια εργασίας/ENEO FINAL FOR GIT/flaskapp/maps/world.shp"  # Path to world shapefile

# Initialize globals as None
RASTER = None
WORLD_GDF = None
COUNTRY_NAMES = {}

def load_datasets():
    """
    Load the necessary datasets at startup:
    1. A population raster TIFF file where pixel values represent population counts.
    2. A world shapefile for country boundaries and names.
    3. A CSV file as a backup for country FID-to-name mappings.
    
    This function is called once when the module is imported to populate the global
    RASTER, WORLD_GDF, and COUNTRY_NAMES variables.
    """
    global RASTER, WORLD_GDF, COUNTRY_NAMES
    
    print("Loading datasets...")
    try:
        # Load the raster file
        RASTER = rasterio.open(TIFF_PATH)
        print(f"Loaded raster from {TIFF_PATH}")
    except rasterio.errors.RasterioIOError:
        print(f"Warning: Population raster not found or could not be opened at {TIFF_PATH}")
        RASTER = None
    except Exception as e:
        print(f"Error loading raster: {e}")
        RASTER = None
            
    try:
        # Load the world shapefile using geopandas
        WORLD_GDF = gpd.read_file(WORLD_SHAPEFILE_PATH)
        print(f"Loaded world shapefile from {WORLD_SHAPEFILE_PATH}")
            
        # Dynamically find the field containing country names from a list of common names
        name_field = None
        potential_name_fields = ['NAME', 'CNTRY_NAME', 'COUNTRY', 'NAME_ENGLI', 'ADMIN', 'SOVEREIGNT']
            
        for field in potential_name_fields:
            if field in WORLD_GDF.columns:
                name_field = field
                print(f"Found country name field: {field}")
                break
            
        if name_field:
            # Find a field to use as a unique country identifier (FID)
            fid_field = None
            potential_fid_fields = ['FID', 'OBJECTID']
            for field in potential_fid_fields:
                if field in WORLD_GDF.columns:
                    fid_field = field
                    break
            
            if not fid_field:
                # If no explicit FID field is found, use the GeoDataFrame index as a fallback
                WORLD_GDF['FID'] = WORLD_GDF.index
                fid_field = 'FID'
                print("Using index as country identifier field (FID)")
            else:
                print(f"Using {fid_field} as country identifier field")
                
            # Create the mapping from FID to country name
            COUNTRY_NAMES = dict(zip(WORLD_GDF[fid_field].astype(int), WORLD_GDF[name_field]))
            print(f"Loaded {len(COUNTRY_NAMES)} country names from shapefile")
        else:
            print("Could not find a suitable country name field in the shapefile.")
    except Exception as e:
        print(f"Warning: Could not load or process world shapefile from {WORLD_SHAPEFILE_PATH}: {e}")
        WORLD_GDF = None
    
    # If loading from the shapefile failed or it didn't contain names, try a backup CSV
    if not COUNTRY_NAMES:
        try:
            if os.path.exists(COUNTRY_NAMES_PATH):
                print(f"Attempting to load country names from CSV: {COUNTRY_NAMES_PATH}")
                country_df = pd.read_csv(COUNTRY_NAMES_PATH)
                if 'fid' in country_df.columns and 'name' in country_df.columns:
                    # Create mapping from CSV columns
                    COUNTRY_NAMES = dict(zip(country_df['fid'].astype(int), country_df['name']))
                    print(f"Loaded {len(COUNTRY_NAMES)} country names from CSV")
                else:
                    print(f"CSV file {COUNTRY_NAMES_PATH} does not contain 'fid' and/or 'name' columns.")
            else:
                print(f"Country names CSV not found at {COUNTRY_NAMES_PATH}")
        except Exception as e:
            print(f"Error loading country names from CSV {COUNTRY_NAMES_PATH}: {e}")
            
    # Ensure a fallback name for unassigned territories (e.g., oceans) exists
    if -1 not in COUNTRY_NAMES:
        COUNTRY_NAMES[-1] = "Unassigned Territory"
        
    print("Dataset loading complete.")

# Load the datasets when this module is imported
load_datasets()

# Note: The create_geodesic_buffer function is assumed to be defined elsewhere
# and imported. It creates a circular or world-covering buffer around a point.

def normalize_geometry(geom, split_antimeridian=True):
    """
    Normalize a geometry's coordinates and handle antimeridian (180/-180 longitude) crossing.

    This function performs two main tasks:
    1. Normalizes all longitudes to the [-180, 180] range and clamps latitudes to [-90, 90].
    2. If a geometry crosses the antimeridian, it splits it into two separate geometries,
       one for the Western Hemisphere and one for the Eastern Hemisphere.

    Args:
        geom (shapely.geometry.BaseGeometry): The input geometry (Polygon or MultiPolygon).
        split_antimeridian (bool): If True, split geometries that cross the antimeridian.

    Returns:
        list: A list of one or more normalized (and possibly split) geometries.
    """
    if geom.is_empty:
        return [geom]
    
    # Define an inner function to normalize individual coordinate pairs
    def normalize_coords(x, y):
        # Normalize longitude to the [-180, 180] range
        while x > 180:
            x -= 360
        while x < -180:
            x += 360
            
        # Clamp latitude to the [-90, 90] range to avoid invalid values
        y = max(-90, min(90, y))
        
        return x, y
    
    # Apply the coordinate normalization to the entire geometry
    normalized = transform(lambda x, y: normalize_coords(x, y), geom)
    
    # Ensure the geometry is valid after transformation, fixing if necessary
    if not normalized.is_valid:
        normalized = make_valid(normalized)
    
    # If splitting is disabled, return the single normalized geometry
    if not split_antimeridian:
        return [normalized]
    
    # --- Antimeridian Crossing Detection ---
    crosses_antimeridian = False
    
    # Extract all coordinates from the geometry for analysis
    coords = []
    if isinstance(normalized, Polygon):
        coords = list(normalized.exterior.coords)
    elif isinstance(normalized, MultiPolygon):
        for poly in normalized.geoms:
            coords.extend(list(poly.exterior.coords))
    
    # Detect crossing by checking if the longitude range is greater than 180 degrees
    lons = [p[0] for p in coords]
    if max(lons) - min(lons) > 180:
        crosses_antimeridian = True
    else:
        # Also check for large jumps between consecutive points (e.g., from +179 to -179)
        for i in range(len(coords) - 1):
            lon_diff = abs(coords[i][0] - coords[i+1][0])
            if lon_diff > 180:
                crosses_antimeridian = True
                break
    
    if not crosses_antimeridian:
        return [normalized]
    
    # --- Geometry Splitting Logic ---
    result = []
    
    # Define a function to split a single polygon
    def split_polygon(poly):
        # Get coordinates from the polygon's exterior ring
        poly_coords = list(poly.exterior.coords)
        
        # Create lists to hold points for the two new polygons
        west_points = []  # For points with longitude < 0 (Western Hemisphere)
        east_points = []  # For points with longitude >= 0 (Eastern Hemisphere)
        
        # Iterate over all segments (pairs of consecutive points) of the polygon
        for i in range(len(poly_coords) - 1):
            p1 = poly_coords[i]
            p2 = poly_coords[i + 1]
            
            # Check if this segment crosses the antimeridian
            if abs(p2[0] - p1[0]) > 180:
                lat_cross = 0.0 # Initialize latitude of the crossing point
                # Calculate the crossing point's latitude using linear interpolation
                if p1[0] < 0 and p2[0] > 0:  # Crossing from west to east (e.g., -170 to 170)
                    # This case is not a standard antimeridian cross, but handled for robustness
                    numerator = -180 - p1[0]
                    denominator = (p2[0] - p1[0]) - 360
                    if denominator == 0.0:
                        # This implies numerator is also 0 (p1[0]=-180, p2[0]=180)
                        # Set t=0, so lat_cross = p1[1]
                        lat_cross = p1[1]
                    else:
                        t = numerator / denominator
                        lat_cross = p1[1] + t * (p2[1] - p1[1])
                    
                    # Add the interpolated point to both new polygons at the antimeridian
                    west_points.append((-180, lat_cross))
                    east_points.append((180, lat_cross))
                else:  # Crossing from east to west (e.g., 170 to -170)
                    numerator = 180 - p1[0]
                    denominator = (p2[0] - p1[0]) + 360
                    if denominator == 0.0:
                        # This implies numerator is also 0 (p1[0]=180, p2[0]=-180)
                        # Set t=0, so lat_cross = p1[1]
                        lat_cross = p1[1]
                    else:
                        t = numerator / denominator
                        lat_cross = p1[1] + t * (p2[1] - p1[1])
                    
                    # Add the interpolated point to both new polygons at the antimeridian
                    east_points.append((180, lat_cross))
                    west_points.append((-180, lat_cross))
            
            # Add the starting point of the segment to the appropriate hemisphere list
            if p1[0] < 0:
                west_points.append(p1)
            else:
                east_points.append(p1)
        
        # --- Polygon Reconstruction ---
        polygons = []
        
        if len(west_points) >= 3:
            # Close the polygon by ensuring the first and last points are the same
            if west_points[0] != west_points[-1]:
                west_points.append(west_points[0])
            west_poly = Polygon(west_points)
            # Validate and add the new western polygon
            if west_poly.is_valid:
                polygons.append(west_poly)
            else:
                west_poly = make_valid(west_poly)
                if not west_poly.is_empty:
                    polygons.append(west_poly)
            
        if len(east_points) >= 3:
            # Close the polygon by ensuring the first and last points are the same
            if east_points[0] != east_points[-1]:
                east_points.append(east_points[0])
            east_poly = Polygon(east_points)
            # Validate and add the new eastern polygon
            if east_poly.is_valid:
                polygons.append(east_poly)
            else:
                east_poly = make_valid(east_poly)
                if not east_poly.is_empty:
                    polygons.append(east_poly)
        
        return polygons
    
    # Apply the splitting logic to the input geometry
    if isinstance(normalized, Polygon):
        result.extend(split_polygon(normalized))
    elif isinstance(normalized, MultiPolygon):
        for poly in normalized.geoms:
            result.extend(split_polygon(poly))
    
    # If splitting failed for any reason, return the original normalized geometry as a fallback
    if not result:
        return [normalized]
    
    return result

def calculate_population_in_zones(lat, lon, vulnerability_zones):
    """
    Calculates the population and estimated casualties within a series of vulnerability zones
    around a given impact point (lat, lon).

    This function processes a list of concentric zones (defined by start/end distances and a
    vulnerability factor), calculates the population within each annular region, and estimates
    casualties based on the vulnerability factor. It also provides a breakdown of population
    and casualties by country.

    Args:
        lat (float): Latitude of the impact point.
        lon (float): Longitude of the impact point.
        vulnerability_zones (list of dicts): A list where each dictionary defines a zone with
                                             'start_distance', 'end_distance', and 'threshold'.

    Returns:
        dict: A dictionary containing detailed results, including population and casualty counts
              for each zone, a total casualty count, a country-by-country breakdown, and any
              warnings generated during processing.
    """
    results = []
    total_casualties = 0
    warnings = []
    
    # Initialize a flag to detect if the impact is near the antimeridian
    date_line_crossing = False
    
    # Sort zones by vulnerability threshold in descending order to process highest-risk zones first
    sorted_zones = sorted(vulnerability_zones, key=lambda x: x['threshold'], reverse=True)
    
    # Create a common WGS84 CRS object
    wgs84 = CRS.from_epsg(4326)
    
    # Check if the location is near the International Date Line, which requires special handling
    if abs(lon) > 170:
        warnings.append("Location is near the International Date Line")
        date_line_crossing = True
    
    # --- Projection Setup ---
    # Select an appropriate map projection for accurate distance calculations.
    # Use polar projections for high latitudes and UTM for others.
    if abs(lat) > 84:
        # Use appropriate polar stereographic projection for polar regions
        polar_epsg = 3995 if lat > 0 else 3031  # North/South polar stereographic
        proj_crs = CRS.from_epsg(polar_epsg)
        warnings.append("Using polar stereographic projection for high latitude location")
    else:
        # Determine the appropriate UTM zone for the given longitude
        zone_number = int((lon + 180) // 6 + 1)
        is_northern = lat >= 0
        utm_epsg = 32600 + zone_number if is_northern else 32700 + zone_number
        utm_crs = CRS.from_epsg(utm_epsg)
    
    # Dictionary to store and accumulate country-specific population and casualty data
    # Using defaultdict simplifies adding data for a country for the first time
    countries_population = defaultdict(lambda: {
        "fid": 0, # Will be set on first encounter
        "name": "", # Will be set on first encounter
        "total_population": 0,
        "total_casualties": 0,
        "zone_breakdown": []
    })
    
    # Use the globally loaded RASTER; if it's not available, exit with an error
    if RASTER is None:
        warnings.append("Population raster data is not loaded")
        return {
            "zones": [],
            "total_casualties": 0,
            "note": "Error: Population data not available",
            "warnings": warnings,
            "countries": []
        }
        
    src = RASTER
    tiff_crs = src.crs
    
    # Create a coordinate transformer from WGS84 (lat/lon) to the raster's CRS
    wgs84_to_tiff = Transformer.from_crs(wgs84, tiff_crs, always_xy=True)
    
    # --- Check if the impact point is within the raster's bounds ---
    point_raster_x, point_raster_y = wgs84_to_tiff.transform(lon, lat)
    left, bottom, right, top = src.bounds
    
    # Warn if the impact location is near the edge of the population dataset,
    # as this might lead to underestimation for large zones.
    if sorted_zones:
        buffer_distance = max([zone['end_distance'] for zone in sorted_zones]) * 1000  # in meters
        edge_margin = 0.1  # A 10% margin
        
        if (point_raster_x < left + buffer_distance * edge_margin or 
            point_raster_x > right - buffer_distance * edge_margin or
            point_raster_y < bottom + buffer_distance * edge_margin or 
            point_raster_y > top - buffer_distance * edge_margin):
            warnings.append("Location is near the edge of the dataset")
    else: # Handle case with no zones provided
        warnings.append("No vulnerability zones provided.")

    # Define a transformation function for use with shapely.ops.transform
    def transform_from_wgs84_to_tiff(x, y):
        return wgs84_to_tiff.transform(x, y)
    
    # Function to convert a shapely geometry into the GeoJSON-like format required by rasterio.mask
    def prepare_geometry_for_mask(geom):
        if isinstance(geom, MultiPolygon):
            # For MultiPolygons, create a list of polygon dictionaries
            result = []
            for poly in geom.geoms:
                poly_dict = {
                    "type": "Polygon", 
                    "coordinates": [list(poly.exterior.coords)]
                }
                # Include interior rings (holes) if they exist
                if hasattr(poly, 'interiors') and poly.interiors:
                    for interior in poly.interiors:
                        poly_dict["coordinates"].append(list(interior.coords))
                result.append(poly_dict)
            return result
        else:
            # For a single Polygon
            poly_dict = {
                "type": "Polygon", 
                "coordinates": [list(geom.exterior.coords)]
            }
            # Include interior rings (holes) if they exist
            if hasattr(geom, 'interiors') and geom.interiors:
                for interior in geom.interiors:
                    poly_dict["coordinates"].append(list(interior.coords))
            return [poly_dict]
    
    previous_buffer = None
    world_coverage_initiated = False # Flag to track if a zone has covered the entire world

    # --- Main Loop: Process each vulnerability zone ---
    for zone in sorted_zones:
        inner_radius_km = zone['start_distance']
        outer_radius_km = zone['end_distance']
        vulnerability = zone['threshold']

        # If a previous, higher-vulnerability zone already covered the entire map,
        # subsequent (lower-vulnerability) zones will contain no additional population.
        if world_coverage_initiated:
            results.append({
                "vulnerability_threshold": vulnerability,
                "start_distance": inner_radius_km,
                "end_distance": outer_radius_km,
                "zone_population": 0,
                "estimated_casualties": 0,
                "country_breakdown": [],
                "note": "Skipped as a previous zone effectively covered the entire map."
            })
            continue

        # Create a geodesic buffer for the outer radius of the current zone.
        # This creates a geographically accurate circle on the Earth's surface.
        outer_buffer = create_geodesic_buffer(lon, lat, outer_radius_km)
        
        # Check if this zone's buffer is large enough to cover the entire world.
        # If so, flag it so subsequent zones can be skipped.
        if outer_radius_km >= 12400: # Approx. half the Earth's circumference
            world_coverage_initiated = True
            warnings.append(f"Zone {inner_radius_km}-{outer_radius_km}km (threshold {vulnerability}) may cover the entire map. Subsequent zones will be assigned 0 population.")

        # Calculate the annular region (the "donut shape") for this zone
        # by subtracting the previous (smaller) buffer from the current (larger) one.
        if previous_buffer:
            try:
                zone_area = outer_buffer.difference(previous_buffer)
                # The difference operation can sometimes produce invalid geometries. Validate and fix.
                if not zone_area.is_valid:
                    zone_area = make_valid(zone_area)
                # If still invalid, try the buffer(0) trick, a common way to fix topological errors.
                if not zone_area.is_valid or zone_area.is_empty:
                    zone_area = outer_buffer.buffer(0).difference(previous_buffer.buffer(0))
            except Exception as e:
                print(f"Warning: Buffer difference operation failed: {e}")
                # Fallback to the buffer(0) trick if the primary method fails.
                try:
                    zone_area = outer_buffer.buffer(0).difference(previous_buffer.buffer(0))
                except:
                    # As a last resort, use the full outer buffer, though this overestimates the area.
                    zone_area = outer_buffer
                    warnings.append(f"Could not calculate difference for zone {inner_radius_km}-{outer_radius_km}km")
        else:
            # For the first zone, the area is just the inner-most buffer.
            zone_area = outer_buffer

        # If the impact is near the date line, the geometry may need to be split.
        if date_line_crossing:
            # The create_geodesic_buffer function should handle this, but we double-check.
            if isinstance(zone_area, Polygon):
                # Check if this single polygon itself crosses the antimeridian
                lons = [p[0] for p in list(zone_area.exterior.coords)]
                if max(lons) - min(lons) > 180:
                    # Split the polygon into two parts using the normalization function
                    split_geoms = normalize_geometry(zone_area)
                    
                    # Ensure the results are valid polygons
                    polygon_geoms = [geom for geom in split_geoms if isinstance(geom, Polygon)]
                    
                    if polygon_geoms:
                        zone_area = MultiPolygon(polygon_geoms) if len(polygon_geoms) > 1 else polygon_geoms[0]
                    else:
                        warnings.append(f"Failed to split polygon at antimeridian for zone {inner_radius_km}-{outer_radius_km}km")
                        # Make sure the original zone_area is valid as a fallback
                        if not zone_area.is_valid:
                            zone_area = make_valid(zone_area)
                            
    
        # Transform the final zone geometry to the raster's CRS for masking
        try:
            tiff_area = transform(transform_from_wgs84_to_tiff, zone_area)
            geoms = prepare_geometry_for_mask(tiff_area)
        except Exception as e:
            print(f"Warning: Error transforming geometry to TIFF CRS: {e}")
            # If transformation fails, try to fix the geometry and re-attempt
            try:
                zone_area = make_valid(zone_area)
                tiff_area = transform(transform_from_wgs84_to_tiff, zone_area)
                geoms = prepare_geometry_for_mask(tiff_area)
            except:
                warnings.append(f"Could not process zone {inner_radius_km}-{outer_radius_km}km")
                # Add a placeholder result and skip to the next zone
                results.append({
                    "vulnerability_threshold": vulnerability,
                    "start_distance": inner_radius_km,
                    "end_distance": outer_radius_km,
                    "zone_population": 0,
                    "estimated_casualties": 0,
                    "error": "Could not process geometry"
                })
                continue

        try:
            # --- Masking the Raster ---
            # Use rasterio.mask to extract the pixel values within the zone's geometry.
            # This reads both population (band 1) and country FID (band 2).
            all_touched = outer_radius_km < 2.3  # Use 'all_touched' for small buffers to avoid missing pixels
            
            out_image, out_transform = mask(
                src, 
                geoms, 
                crop=True, # Crop the output to the bounds of the geometry
                nodata=src.nodata, 
                indexes=[1, 2], # Read both bands
                all_touched=all_touched
            )
            
            population_data = out_image[0]  # Band 1: Population counts
            country_fid_data = out_image[1]  # Band 2: Country identifiers
            
            # Create a mask for valid data pixels (where population is not nodata)
            valid_mask = population_data != src.nodata
            
            # Special handling for pixels with FID -1 (unassigned territory, e.g., oceans)
            unassigned_mask = (country_fid_data == -1) & valid_mask
            if np.any(unassigned_mask):
                print(f"Found {np.sum(unassigned_mask)} pixels with FID -1 containing population")
                
                # Get the total population in these unassigned areas
                unassigned_population = np.sum(population_data[unassigned_mask])
                
                # Find all valid country FIDs that are also present in this zone
                valid_fids = np.unique(country_fid_data[(country_fid_data > 0) & valid_mask])
                
                if len(valid_fids) > 0:
                    # This section appears incomplete in the original code.
                    # It likely intended to distribute the unassigned population
                    # among the present countries, but the logic is not finished.
                    pass
            
            # Calculate the total population for the entire zone
            zone_population = int(np.sum(population_data[valid_mask]))
            
            # Calculate estimated casualties for the zone based on its vulnerability factor
            zone_casualties = int(zone_population * vulnerability)
            total_casualties += zone_casualties

            # --- Country-by-Country Breakdown ---
            unique_fids = np.unique(country_fid_data[valid_mask])
            country_breakdown = []
            
            for fid in unique_fids:
                if fid == src.nodata or fid <= 0:
                    # Skip invalid FIDs
                    continue
                
                # Create a mask for the current country's FID within the valid data area
                country_mask = (country_fid_data == fid) & valid_mask
                
                # Calculate population and casualties for this specific country in this zone
                country_population = int(np.sum(population_data[country_mask]))
                country_casualties = int(country_population * vulnerability)
                
                # Get country name and FID from the lookup table
                country_name = COUNTRY_NAMES.get(int(fid), f"Unknown (FID: {int(fid)})")
                
                # Add to the zone's specific country breakdown
                country_breakdown.append({
                    "fid": int(fid),
                    "name": country_name,
                    "population": country_population,
                    "casualties": country_casualties
                })
                
                # Aggregate totals for the final country list
                countries_population[int(fid)]["fid"] = int(fid)
                countries_population[int(fid)]["name"] = country_name
                countries_population[int(fid)]["total_population"] += country_population
                countries_population[int(fid)]["total_casualties"] += country_casualties
                countries_population[int(fid)]["zone_breakdown"].append({
                    "vulnerability": vulnerability,
                    "population": country_population,
                    "casualties": country_casualties
                })


            # Append the results for the current zone
            results.append({
                "vulnerability_threshold": vulnerability,
                "start_distance": inner_radius_km,
                "end_distance": outer_radius_km,
                "zone_population": zone_population,
                "estimated_casualties": zone_casualties,
                "country_breakdown": country_breakdown
            })

        except ValueError as e:
            # This can happen if the zone is completely outside the raster data
            if "Input shapes do not overlap raster." in str(e):
                # Add a result with zero population for this zone
                results.append({
                    "vulnerability_threshold": vulnerability,
                    "start_distance": inner_radius_km,
                    "end_distance": outer_radius_km,
                    "zone_population": 0,
                    "estimated_casualties": 0,
                    "country_breakdown": [],
                    "note": "Zone is outside the population dataset."
                })
            else:
                # Handle other ValueErrors
                warnings.append(f"Could not calculate population for zone {inner_radius_km}-{outer_radius_km}km: {e}")

        # Update the previous buffer for the next iteration's difference calculation
        previous_buffer = outer_buffer
    
    # --- Final Aggregation and Verification ---
    # Sort the aggregated country list by total casualties in descending order
    countries_list = list(countries_population.values())
    countries_list.sort(key=lambda x: x["total_casualties"], reverse=True)
    
    # Verify that the sum of country totals matches the overall totals
    country_total_population = sum(c["total_population"] for c in countries_list)
    country_total_casualties = sum(c["total_casualties"] for c in countries_list)

    zone_total_population = sum(z["zone_population"] for z in results)

    # If there's a discrepancy (which can happen due to rounding or incomplete logic),
    # add the residual amount to the country with the fewest casualties to ensure totals match.
    if countries_list and (abs(country_total_casualties - total_casualties) > 0 or 
                           abs(country_total_population - zone_total_population) > 0):
        # Calculate the residual amounts
        casualties_residual = total_casualties - country_total_casualties
        population_residual = zone_total_population - country_total_population
        
        # Find the country with the smallest casualty count to add the residual to
        min_casualties_country = min(countries_list, key=lambda x: x["total_casualties"])
        
        # Add the residuals to this country's totals
        # This part of the logic is also incomplete in the original code.
        pass


    # Add a final warning if a significant discrepancy remains after the correction attempt
    if abs(country_total_casualties - total_casualties) > 1:
        warnings.append("Final casualty counts between zones and countries do not match perfectly.")

    return {
        "zones": results,
        "total_casualties": total_casualties,
        "note": "Population estimates based on a 2.5 arcminute resolution map",
        "warnings": warnings if warnings else None,
        "countries": countries_list
    }

def calculate_population_for_standard_zones(lat, lon, vulnerability_zones, max_distance_km=3000):
    """
    A simplified version of population calculation for zones under a certain size
    that are not expected to cross the antimeridian.

    This function is intended as a faster alternative for smaller impact events
    by skipping complex edge cases like polar projections and antimeridian crossing.

    Args:
        lat (float): Latitude of the impact point.
        lon (float): Longitude of the impact point.
        vulnerability_zones (list of dicts): The list of vulnerability zones.
        max_distance_km (int): The maximum radius this function will handle.

    Returns:
        list: A list of results for each zone.
    """
    results = []
    total_casualties = 0
    
    # Sort zones by vulnerability threshold in descending order
    sorted_zones = sorted(vulnerability_zones, key=lambda x: x['threshold'], reverse=True)
    
    # Filter out zones that are too large for this simplified function
    sorted_zones = [z for z in sorted_zones if z['end_distance'] <= max_distance_km]
    
    # This function's implementation is incomplete in the provided file.
    # It would likely follow a similar, but simplified, logic to
    # calculate_population_in_zones.
    
    return results
def create_geodesic_buffer(center_lng, center_lat, buffer_radius_km):
    """
    Create a geodesic buffer around a point with proper handling of poles and the antimeridian.
    Returns a Polygon or MultiPolygon in WGS84 coordinates.
    """
    # Create center point
    center_point = Point(center_lng, center_lat)
    
    # Earth's radius in km
    EARTH_RADIUS = 6371.0088 # More precise GRS80 mean radius, or stick to 6371.0
    
    # Check if we need to handle pole proximity
    dist_to_north_pole_km = (90 - center_lat) * (math.pi / 180) * EARTH_RADIUS
    dist_to_south_pole_km = (90 + center_lat) * (math.pi / 180) * EARTH_RADIUS
    
    # Calculate if we need to limit the shape at poles
    limit_north = dist_to_north_pole_km < buffer_radius_km
    limit_south = dist_to_south_pole_km < buffer_radius_km
    
    # Special case: Buffer is so large it covers most of the Earth
    covers_most_earth = buffer_radius_km > EARTH_RADIUS * math.pi / 2
    
    # For very large buffers (approaching or exceeding Earth circumference)
    if buffer_radius_km >= 12400:
        # Create a polygon covering the entire world for extreme distances
        world_coords = [
            (-180, -90), (180, -90), (180, 90), (-180, 90), (-180, -90)
        ]
        return Polygon(world_coords)
    
    # Center point in radians
    lat_rad = math.radians(center_lat)
    lon_rad = math.radians(center_lng)
    
    # Angular radius in radians
    angular_radius = buffer_radius_km / EARTH_RADIUS
    
    # Determine number of points based on buffer size
    # Scale number of points based on buffer size but with minimum for small buffers
    # n_points = int(400 * (1 + math.log(1 + buffer_radius_km / 1000)))
    # Ensure minimum number of points for very small buffers (at least 16)
    # n_points = max(16, min(1000, n_points))  # Minimum of 16 points, maximum of 1000
    
    # Improved point calculation for small buffers
    n_points = 64  # Minimum points for any buffer size
    if buffer_radius_km > 1:
        # Logarithmic scaling maintains performance while increasing detail
        n_points = min(1000, max(64, int(400 * math.log1p(buffer_radius_km))))

    if covers_most_earth or limit_north or limit_south:
        n_points = max(n_points, 400)  # Use more points for better shape accuracy in complex cases
    
    angles = np.linspace(0, 2*np.pi, n_points, endpoint=False)
    buffer_coords = []
    
    for angle in angles:
        # Determine if we need to adjust radius for poles
        current_radius = angular_radius
        
        # Modify radius for N-S directions if needed
        if limit_north or limit_south:
            # Determine how much of this point's direction is northward/southward
            north_component = math.cos(angle)  # 1 for north, -1 for south
            
            # Only adjust radius when heading toward limited pole
            if (north_component > 0 and limit_north) or (north_component < 0 and limit_south):
                # Calculate how much to scale the radius
                if north_component > 0:  # Heading north
                    max_allowed_radius = dist_to_north_pole_km / EARTH_RADIUS * 0.99  # 99% for safety
                else:  # Heading south
                    max_allowed_radius = dist_to_south_pole_km / EARTH_RADIUS * 0.99  # 99% for safety
                
                # Improved transition function for more natural ellipsoidal shape
                # Use a smoother cubic easing function for better shape near poles
                direction_factor = abs(north_component)
                
                if direction_factor > 0.2:  # Start adjusting earlier
                    # Normalize factor between 0.2 and 1.0 to 0.0 and 1.0
                    normalized_factor = (direction_factor - 0.2) / 0.8
                    # Cubic easing for smoother transition
                    blend = normalized_factor * normalized_factor * (3 - 2 * normalized_factor)
                    # Blend between full radius and max allowed
                    current_radius = angular_radius * (1 - blend) + max_allowed_radius * blend
        
        # Calculate the new point using the adjusted radius with improved formula for edge cases
        sin_lat = math.sin(lat_rad) * math.cos(current_radius) + math.cos(lat_rad) * math.sin(current_radius) * math.cos(angle)
        new_lat_rad = math.asin(max(-1, min(1, sin_lat)))  # Clamp to valid range
        
        # Handle numerical precision issues for longitude calculation
        denom = math.cos(current_radius) - math.sin(lat_rad) * math.sin(new_lat_rad)
        if abs(denom) < 1e-10:
            # Point is very close to a pole, longitude becomes unstable
            # Use center longitude in this case
            new_lon_rad = lon_rad
        else:
            new_lon_rad = lon_rad + math.atan2(
                math.sin(angle) * math.sin(current_radius) * math.cos(lat_rad),
                denom
            )
        
        # Convert back to degrees
        new_lat = math.degrees(new_lat_rad)
        new_lon = math.degrees(new_lon_rad)
        
        # Normalize longitude to -180 to 180
        while new_lon > 180:
            new_lon -= 360
        while new_lon < -180:
            new_lon += 360
        
        # Apply final latitude constraints to ensure we don't exceed ±90°
        new_lat = max(-90, min(90, new_lat))
        
        buffer_coords.append((new_lon, new_lat))
    
    # Close the polygon
    buffer_coords.append(buffer_coords[0])
    
    # For very small buffers, ensure points are actually distinct
    # Use a smaller tolerance for small buffers
    min_coord_distance = 0.0001 if buffer_radius_km < 1 else 0.01
    
    # Filter out duplicate or very close points to avoid shape artifacts
    filtered_coords = [buffer_coords[0]]
    for i in range(1, len(buffer_coords)):
        prev_pt = filtered_coords[-1]
        curr_pt = buffer_coords[i]
        # Only add point if it's sufficiently different from the previous one
        if abs(curr_pt[0] - prev_pt[0]) > min_coord_distance or abs(curr_pt[1] - prev_pt[1]) > min_coord_distance:
            filtered_coords.append(curr_pt)
    
    # Ensure polygon is still closed
    if filtered_coords[0] != filtered_coords[-1]:
        filtered_coords.append(filtered_coords[0])
    
    # Ensure we have at least 4 points (minimum for a valid LinearRing)
    if len(filtered_coords) < 4:
        # For very small buffers, create an artificial 4-point square around the center
        # Ensure dx, dy are not zero if buffer_radius_km is very small but non-zero
        min_delta = 1e-5 # Minimum delta in degrees for smallest buffers
        dx = max(min_delta, 0.0001 * buffer_radius_km / (111 * math.cos(math.radians(center_lat)) + 1e-6)) # Approx degree conversion
        dy = max(min_delta, 0.0001 * buffer_radius_km / 111)  # Approx degree conversion
        
        # Create a simple square buffer
        filtered_coords = [
            (center_lng - dx, center_lat - dy),
            (center_lng + dx, center_lat - dy),
            (center_lng + dx, center_lat + dy),
            (center_lng - dx, center_lat + dy),
            (center_lng - dx, center_lat - dy)  # Close the polygon
        ]
        print(f"Created simple square buffer due to small buffer size: {buffer_radius_km} km")
    
    buffer_coords = filtered_coords
    
    # Check if buffer actually crosses the antimeridian by looking for segments with large longitude jumps
    has_antimeridian_crossing = False
    for i in range(len(buffer_coords) - 1):
        p1 = buffer_coords[i]
        p2 = buffer_coords[i + 1]
        # If adjacent points have longitudes that differ by more than 180 degrees,
        # we have an actual antimeridian crossing
        if abs(p2[0] - p1[0]) > 180:
            has_antimeridian_crossing = True
            break
    
    if has_antimeridian_crossing:
        # Split into two polygons at the antimeridian
        west_points = []
        east_points = []
        
        for i in range(len(buffer_coords) - 1):
            p1 = buffer_coords[i]
            p2 = buffer_coords[i + 1]
            
            # Check if this segment crosses the antimeridian
            if abs(p2[0] - p1[0]) > 180:
                # Calculate the crossing point using linear interpolation
                if p1[0] < 0 and p2[0] > 0:  # Crossing from -180 to 180
                    t = (-180 - p1[0]) / (p2[0] - p1[0] - 360)
                    lat_cross = p1[1] + t * (p2[1] - p1[1])
                    
                    west_points.append((-180, lat_cross))
                    east_points.append((180, lat_cross))
                else:  # Crossing from 180 to -180
                    t = (180 - p1[0]) / (p2[0] - p1[0] + 360)
                    lat_cross = p1[1] + t * (p2[1] - p1[1])
                    
                    east_points.append((180, lat_cross))
                    west_points.append((-180, lat_cross))
            
            # Add the current point to the appropriate list
            if p1[0] < 0:
                west_points.append(p1)
            else:
                east_points.append(p1)
        
        # Create the two polygons
        polygons = []
        
        if len(west_points) >= 3:
            # Close the polygon if needed
            if west_points[0] != west_points[-1]:
                west_points.append(west_points[0])
            try:
                polygons.append(Polygon(west_points))
            except Exception as e:
                print(f"Warning: Could not create west polygon: {e}")
                # Try to clean up the polygon
                if len(west_points) > 4:
                    try:
                        polygons.append(Polygon(west_points).buffer(0))
                    except:
                        pass
            
        if len(east_points) >= 3:
            # Close the polygon if needed
            if east_points[0] != east_points[-1]:
                east_points.append(east_points[0])
            try:
                polygons.append(Polygon(east_points))
            except Exception as e:
                print(f"Warning: Could not create east polygon: {e}")
                # Try to clean up the polygon
                if len(east_points) > 4:
                    try:
                        polygons.append(Polygon(east_points).buffer(0))
                    except:
                        pass
        
        if polygons:
            return MultiPolygon(polygons)
        else:
            # Fallback to using the original polygon if splitting failed
            return Polygon(buffer_coords)
    else:
        # No antimeridian crossing
        try:
            return Polygon(buffer_coords)
        except Exception as e:
            print(f"Warning: Error creating polygon: {e}")
            # Last resort - ensure we always return a valid polygon
            dx, dy = 0.001 * buffer_radius_km, 0.001 * buffer_radius_km
            if dx < 0.00001:  # Ensure minimum size
                dx = dy = 0.00001
            
            # Create a simple square buffer
            square_coords = [
                (center_lng - dx, center_lat - dy),
                (center_lng + dx, center_lat - dy),
                (center_lng + dx, center_lat + dy),
                (center_lng - dx, center_lat + dy),
                (center_lng - dx, center_lat - dy)  # Close the polygon
            ]
            return Polygon(square_coords)
def calculate_population(lat, lon, vulnerability_zones):
    """
    Main function to calculate population that branches to the appropriate handler
    based on complexity of the calculation required
    """
    # Check for complex case conditions
    max_radius = max(zone['end_distance'] for zone in vulnerability_zones)
    near_poles = abs(lat) > 84
    near_antimeridian = abs(lon) > 170
    
    # Determine which function to use
    if near_poles or near_antimeridian or max_radius > 3000:
        # Complex case - use full calculation with all edge case handling
        print(f"Using complex calculation method because: " + 
              (f"near poles ({lat})" if near_poles else "") +
              (f"near antimeridian ({lon})" if near_antimeridian else "") +
              (f"large radius ({max_radius}km)" if max_radius > 3000 else ""))
        
        return calculate_population_in_zones(lat, lon, vulnerability_zones)
    else:
        # Simple case - use optimized calculation for standard cases
        print(f"Using simplified calculation method (lat={lat}, lon={lon}, max_radius={max_radius}km)")
        return calculate_population_for_standard_zones(lat, lon, vulnerability_zones)

def cleanup():
    """Close any open resources when shutting down"""
    global RASTER
    if RASTER is not None:
        RASTER.close()
        RASTER = None
        print("Closed raster file")
