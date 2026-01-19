"""
ENEO Asteroid Impact Simulation - Population Impact Calculator Module

This module calculates population casualties within asteroid impact damage zones using
high-resolution global population and country boundary datasets. It handles complex
geographical scenarios including polar regions, antimeridian crossings, and very
large impact zones that span continents.

Key Features:
- Global population raster processing with 2.5 arcminute resolution
- Country-specific casualty attribution using FID mappings
- Geodesic buffer creation for accurate distance calculations
- Antimeridian and polar region handling for worldwide coverage
- Optimized calculations for standard vs. complex geographical scenarios
- Vulnerability zone processing with graduated casualty thresholds

Data Sources:
- Population: Global gridded population with country assignments
- Boundaries: World shapefile with country names and identifiers

Author: Alexandros Notas
Institution: National Technical University of Athens
Date: July 2025
"""

import rasterio
from rasterio.mask import mask
from pyproj import Transformer, CRS
from shapely.geometry import Point, MultiPolygon, Polygon
from shapely.ops import transform
import numpy as np
import os
import geopandas as gpd
import pandas as pd
import datetime
from shapely.geometry.base import BaseGeometry
from shapely.validation import make_valid
import math # Ensure math is imported if create_geodesic_buffer uses it directly
from collections import defaultdict # Added for defaultdict

# =============================================================================
# GLOBAL DATA PATHS AND CONFIGURATION
# =============================================================================

# File paths for population and country data (Use your local path here)
# dynamic path resolution
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
MAPS_DIR = os.path.join(BASE_DIR, 'maps')

TIFF_PATH = os.path.join(MAPS_DIR, "population_with_country_fid_assigned.tif")
COUNTRY_NAMES_PATH = os.path.join(MAPS_DIR, "country_fid_lookup.csv")
WORLD_SHAPEFILE_PATH = os.path.join(MAPS_DIR, "world.shp")


# Global variables to hold loaded datasets (populated at module import)
RASTER = None  # Global population raster with country FID assignments
WORLD_GDF = None  # World countries geodataframe
COUNTRY_NAMES = {}  # Mapping from FID to country names

def load_datasets():
    """
    Load global population and country datasets at module startup.
    
    Loads and prepares:
    1. Population raster with embedded country FID assignments
    2. World shapefile for country boundary and name information
    3. Country name lookup mapping from FID to readable names
    
    Implements fallback loading from CSV if shapefile processing fails.
    Populates global variables RASTER, WORLD_GDF, and COUNTRY_NAMES.
    """
    global RASTER, WORLD_GDF, COUNTRY_NAMES
    
    print("Loading datasets...")
    
    # Load population raster with country assignments
    try:
        # Open raster file containing population data and country FIDs
        RASTER = rasterio.open(TIFF_PATH)
        print(f"Loaded raster from {TIFF_PATH}")
    except rasterio.errors.RasterioIOError:
        print(f"Warning: Population raster not found or could not be opened at {TIFF_PATH}")
        RASTER = None
    except Exception as e:
        print(f"Error loading raster: {e}")
        RASTER = None
            
    # Load world shapefile for country information
    try:
        # Load country boundaries and attributes
        WORLD_GDF = gpd.read_file(WORLD_SHAPEFILE_PATH)
        print(f"Loaded world shapefile from {WORLD_SHAPEFILE_PATH}")
            
        # Attempt to identify country name field from common possibilities
        name_field = None
        potential_name_fields = ['NAME', 'CNTRY_NAME', 'COUNTRY', 'NAME_ENGLI', 'ADMIN', 'SOVEREIGNT']
            
        # Find available country name field in shapefile attributes
        for field in potential_name_fields:
            if field in WORLD_GDF.columns:
                name_field = field
                print(f"Found country name field: {field}")
                break
            
        if name_field:
            # Locate country identifier field (FID or equivalent)
            fid_field = None
            potential_fid_fields = ['FID', 'OBJECTID']
            for field in potential_fid_fields:
                if field in WORLD_GDF.columns:
                    fid_field = field
                    break
            
            # Use row index as FID if no explicit identifier found
            if not fid_field:
                WORLD_GDF['FID'] = WORLD_GDF.index
                fid_field = 'FID'
                print("Using index as country identifier field (FID)")
            else:
                print(f"Using {fid_field} as country identifier field")
                
            # Create FID to country name mapping for fast lookup
            COUNTRY_NAMES = dict(zip(WORLD_GDF[fid_field].astype(int), WORLD_GDF[name_field]))
            print(f"Loaded {len(COUNTRY_NAMES)} country names from shapefile")
        else:
            print("Could not find a suitable country name field in the shapefile.")
    except Exception as e:
        print(f"Warning: Could not load or process world shapefile from {WORLD_SHAPEFILE_PATH}: {e}")
        WORLD_GDF = None
    
    # Fallback to CSV lookup if shapefile method failed
    if not COUNTRY_NAMES:
        try:
            if os.path.exists(COUNTRY_NAMES_PATH):
                print(f"Attempting to load country names from CSV: {COUNTRY_NAMES_PATH}")
                country_df = pd.read_csv(COUNTRY_NAMES_PATH)
                if 'fid' in country_df.columns and 'name' in country_df.columns:
                    COUNTRY_NAMES = dict(zip(country_df['fid'].astype(int), country_df['name']))
                    print(f"Loaded {len(COUNTRY_NAMES)} country names from CSV")
                else:
                    print(f"CSV file {COUNTRY_NAMES_PATH} does not contain 'fid' and/or 'name' columns.")
            else:
                print(f"Country names CSV not found at {COUNTRY_NAMES_PATH}")
        except Exception as e:
            print(f"Error loading country names from CSV {COUNTRY_NAMES_PATH}: {e}")
            
    # Ensure fallback name exists for unassigned territories
    if -1 not in COUNTRY_NAMES:
        COUNTRY_NAMES[-1] = "Unassigned Territory"
        
    print("Dataset loading complete.")

# Initialize datasets when module is imported
load_datasets()

def normalize_geometry(geom, split_antimeridian=True):
    """
    Normalize geometry coordinates and handle antimeridian crossings.
    
    Ensures geometry coordinates are within standard ranges and optionally
    splits geometries that cross the International Date Line into separate
    polygons for proper processing.
    
    Args:
        geom: Shapely geometry to normalize
        split_antimeridian: Whether to split at 180/-180 longitude line
        
    Returns:
        list: List of normalized geometries (single item if no splitting)
    """
    if geom.is_empty:
        return [geom]
    
    # Normalize coordinates to standard geographic ranges
    def normalize_coords(x, y):
        # Wrap longitude to [-180, 180] range
        while x > 180:
            x -= 360
        while x < -180:
            x += 360
            
        # Clamp latitude to valid [-90, 90] range
        y = max(-90, min(90, y))
        
        return x, y
    
    # Apply coordinate normalization to all geometry points
    normalized = transform(lambda x, y: normalize_coords(x, y), geom)
    
    # Ensure geometry validity after coordinate changes
    if not normalized.is_valid:
        normalized = make_valid(normalized)
    
    # Return early if antimeridian splitting not requested
    if not split_antimeridian:
        return [normalized]
    
    # Detect antimeridian crossing by analyzing coordinate distribution
    crosses_antimeridian = False
    
    # Extract coordinates from geometry for analysis
    coords = []
    if isinstance(normalized, Polygon):
        coords = list(normalized.exterior.coords)
    elif isinstance(normalized, MultiPolygon):
        for poly in normalized.geoms:
            coords.extend(list(poly.exterior.coords))
    
    # Check for antimeridian crossing using longitude range analysis
    lons = [p[0] for p in coords]
    if max(lons) - min(lons) > 180:
        crosses_antimeridian = True
    else:
        # Also check for large jumps between consecutive coordinate pairs
        for i in range(len(coords) - 1):
            lon_diff = abs(coords[i][0] - coords[i+1][0])
            if lon_diff > 180:
                crosses_antimeridian = True
                break
    
    if not crosses_antimeridian:
        return [normalized]
    
    # Split geometry into eastern and western hemispheres
    result = []
    
    # Split function for a polygon
    def split_polygon(poly):
        # Get coordinates from the polygon
        poly_coords = list(poly.exterior.coords)
        
        # Create two hemispheres of the polygon
        west_points = []  # For points with longitude < 0
        east_points = []  # For points with longitude > 0
        
        # Process all segments (pairs of consecutive points)
        for i in range(len(poly_coords) - 1):
            p1 = poly_coords[i]
            p2 = poly_coords[i + 1]
            
            # Check if this segment crosses the antimeridian
            if abs(p2[0] - p1[0]) > 180:
                lat_cross = 0.0 # Initialize lat_cross
                # Calculate the crossing point using linear interpolation
                if p1[0] < 0 and p2[0] > 0:  # Crossing from -180 to 180
                    numerator = -180 - p1[0]
                    denominator = (p2[0] - p1[0]) - 360
                    if denominator == 0.0:
                        # This implies numerator is also 0 (p1[0]=-180, p2[0]=180)
                        # Set t=0, so lat_cross = p1[1]
                        lat_cross = p1[1]
                    else:
                        t = numerator / denominator
                        lat_cross = p1[1] + t * (p2[1] - p1[1])
                    
                    west_points.append((-180, lat_cross))
                    east_points.append((180, lat_cross))
                else:  # Crossing from 180 to -180
                    numerator = 180 - p1[0]
                    denominator = (p2[0] - p1[0]) + 360
                    if denominator == 0.0:
                        # This implies numerator is also 0 (p1[0]=180, p2[0]=-180)
                        # Set t=0, so lat_cross = p1[1]
                        lat_cross = p1[1]
                    else:
                        t = numerator / denominator
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
            west_poly = Polygon(west_points)
            if west_poly.is_valid:
                polygons.append(west_poly)
            else:
                west_poly = make_valid(west_poly)
                if not west_poly.is_empty:
                    polygons.append(west_poly)
            
        if len(east_points) >= 3:
            # Close the polygon if needed
            if east_points[0] != east_points[-1]:
                east_points.append(east_points[0])
            east_poly = Polygon(east_points)
            if east_poly.is_valid:
                polygons.append(east_poly)
            else:
                east_poly = make_valid(east_poly)
                if not east_poly.is_empty:
                    polygons.append(east_poly)
        
        return polygons
    
    # Split the geometry
    if isinstance(normalized, Polygon):
        result.extend(split_polygon(normalized))
    elif isinstance(normalized, MultiPolygon):
        for poly in normalized.geoms:
            result.extend(split_polygon(poly))
    
    # If splitting failed, return the original normalized geometry
    if not result:
        return [normalized]
    
    return result

def calculate_population_in_zones(lat, lon, vulnerability_zones):
    results = []
    total_casualties = 0
    warnings = []
    
    # Initialize crossing variables at the beginning
    date_line_crossing = False
    
    # Sort zones by vulnerability threshold in descending order
    sorted_zones = sorted(vulnerability_zones, key=lambda x: x['threshold'], reverse=True)
    
    # Create common transformers early
    wgs84 = CRS.from_epsg(4326)
    
    # Check for international date line crossing
    if abs(lon) > 170:
        warnings.append("Location is near the International Date Line")
        date_line_crossing = True
    
    # Get UTM or polar projection as appropriate
    if abs(lat) > 84:
        # Use appropriate polar stereographic projection instead of UTM
        polar_epsg = 3995 if lat > 0 else 3031  # North/South polar stereographic
        proj_crs = CRS.from_epsg(polar_epsg)
        # transformer_to_proj = Transformer.from_crs(wgs84, proj_crs, always_xy=True) # Not directly used later
        # transformer_from_proj = Transformer.from_crs(proj_crs, wgs84, always_xy=True) # Not directly used later
        warnings.append("Using polar stereographic projection for high latitude location")
    else:
        # Get UTM zone
        zone_number = int((lon + 180) // 6 + 1)
        is_northern = lat >= 0
        utm_epsg = 32600 + zone_number if is_northern else 32700 + zone_number
        
        # Create coordinate transformers
        utm_crs = CRS.from_epsg(utm_epsg)
        # transformer_to_proj = Transformer.from_crs(wgs84, utm_crs, always_xy=True) # Not directly used later
        # transformer_from_proj = Transformer.from_crs(utm_crs, wgs84, always_xy=True) # Not directly used later
    
    # Dictionary to store country-specific population data
    # Using defaultdict to simplify accumulation
    countries_population = defaultdict(lambda: {
        "fid": 0, # Will be set on first encounter
        "name": "", # Will be set on first encounter
        "total_population": 0,
        "total_casualties": 0,
        "zone_breakdown": []
    })
    
    # Use the globally loaded RASTER instead of opening it every time
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
    
    # Create WGS84 to TIFF transformer once
    wgs84_to_tiff = Transformer.from_crs(wgs84, tiff_crs, always_xy=True)
    
    # Check raster bounds
    point_raster_x, point_raster_y = wgs84_to_tiff.transform(lon, lat)
    left, bottom, right, top = src.bounds
    
    # Check if point is outside or near edge of dataset
    # Ensure sorted_zones is not empty before accessing 'end_distance'
    if sorted_zones:
        buffer_distance = max([zone['end_distance'] for zone in sorted_zones]) * 1000  # in meters
        edge_margin = 0.1  # fraction of buffer to consider "near edge"
        
        if (point_raster_x < left + buffer_distance * edge_margin or 
            point_raster_x > right - buffer_distance * edge_margin or
            point_raster_y < bottom + buffer_distance * edge_margin or 
            point_raster_y > top - buffer_distance * edge_margin):
            warnings.append("Location is near the edge of the dataset")
    else: # Handle case with no zones
        warnings.append("No vulnerability zones provided.")


    # Set up transformer from WGS84 to tiff CRS (used for transforming zone_area)
    # This transformer is directly from WGS84 to TIFF CRS for the zone_area geometry
    # The 'proj_crs' (UTM or Polar) is not directly used for transforming the final zone_area to TIFF CRS
    # as zone_area is always in WGS84 after create_geodesic_buffer.
    # So, wgs84_to_tiff is the correct transformer here.
    
    # Pre-define transformation functions
    # def transform_to_tiff(x, y): # This was defined but not used with proj_crs -> tiff_crs
    #     return transformer_to_tiff.transform(x, y)
        
    def transform_from_wgs84_to_tiff(x, y): # This is the key one for masking
        return wgs84_to_tiff.transform(x, y)
    
    # Function to extract geometry details for masking
    def prepare_geometry_for_mask(geom):
        if isinstance(geom, MultiPolygon):
            result = []
            for poly in geom.geoms:
                poly_dict = {
                    "type": "Polygon", 
                    "coordinates": [list(poly.exterior.coords)]
                }
                if hasattr(poly, 'interiors') and poly.interiors:
                    for interior in poly.interiors:
                        poly_dict["coordinates"].append(list(interior.coords))
                result.append(poly_dict)
            return result
        else:
            poly_dict = {
                "type": "Polygon", 
                "coordinates": [list(geom.exterior.coords)]
            }
            if hasattr(geom, 'interiors') and geom.interiors:
                for interior in geom.interiors:
                    poly_dict["coordinates"].append(list(interior.coords))
            return [poly_dict]
    
    previous_buffer = None
    world_coverage_initiated = False # Flag to track if a zone has covered the world

    for zone in sorted_zones:
        inner_radius_km = zone['start_distance']
        outer_radius_km = zone['end_distance']
        vulnerability = zone['threshold']

        if world_coverage_initiated:
            # A previous zone already covered the entire map.
            # Subsequent zones get 0 population/casualties.
            results.append({
                "vulnerability_threshold": vulnerability,
                "start_distance": inner_radius_km,
                "end_distance": outer_radius_km,
                "zone_population": 0,
                "estimated_casualties": 0,
                "country_breakdown": [],
                "note": "Skipped as a previous zone effectively covered the entire map."
            })
            # previous_buffer remains from the last processed zone.
            # If the last processed zone was the one covering the world, 
            # previous_buffer will be that world-covering buffer.
            continue

        # Use the improved geodesic buffer creation method
        outer_buffer = create_geodesic_buffer(lon, lat, outer_radius_km)
        
        # Check if this zone's outer buffer covers the world
        # Using 20000km as the threshold, consistent with create_geodesic_buffer
        if outer_radius_km >= 12400:
            world_coverage_initiated = True
            # This current zone will be processed. Subsequent zones will be skipped.
            warnings.append(f"Zone {inner_radius_km}-{outer_radius_km}km (threshold {vulnerability}) may cover the entire map. Subsequent zones will be assigned 0 population.")

        # Calculate the zone area (difference between outer and previous buffer)
        if previous_buffer:
            # Handle any issues with buffer difference operations
            try:
                zone_area = outer_buffer.difference(previous_buffer)
                # Sometimes the difference operation can fail or produce invalid geometry
                if not zone_area.is_valid:
                    zone_area = make_valid(zone_area)
                # If still invalid or empty, try buffer(0) trick
                if not zone_area.is_valid or zone_area.is_empty:
                    zone_area = outer_buffer.buffer(0).difference(previous_buffer.buffer(0))
            except Exception as e:
                print(f"Warning: Buffer difference operation failed: {e}")
                # Fallback: Use buffer(0) to fix topology issues
                try:
                    zone_area = outer_buffer.buffer(0).difference(previous_buffer.buffer(0))
                except:
                    # Last resort: just use the outer buffer
                    zone_area = outer_buffer
                    warnings.append(f"Could not calculate difference for zone {inner_radius_km}-{outer_radius_km}km")
        else:
            zone_area = outer_buffer

        # Handle date line crossing if needed
        if date_line_crossing:
            # The create_geodesic_buffer function already handles the antimeridian
            # but we'll ensure the geometry is split properly for masking
            if isinstance(zone_area, Polygon):
                # Check if this single polygon crosses the antimeridian
                lons = [p[0] for p in list(zone_area.exterior.coords)]
                if max(lons) - min(lons) > 180:
                    # Split at antimeridian
                    split_geoms = normalize_geometry(zone_area)
                    
                    # Filter out any non-Polygon geometries
                    polygon_geoms = [geom for geom in split_geoms if isinstance(geom, Polygon)]
                    
                    if polygon_geoms:
                        zone_area = MultiPolygon(polygon_geoms) if len(polygon_geoms) > 1 else polygon_geoms[0]
                    else:
                        # If no valid polygons remain, keep the original zone_area
                        warnings.append(f"Failed to split polygon at antimeridian for zone {inner_radius_km}-{outer_radius_km}km")
                        # Make sure zone_area is valid
                        if not zone_area.is_valid:
                            zone_area = make_valid(zone_area)
    
        # Transform to TIFF CRS for masking
        try:
            tiff_area = transform(transform_from_wgs84_to_tiff, zone_area)
            geoms = prepare_geometry_for_mask(tiff_area)
        except Exception as e:
            print(f"Warning: Error transforming geometry to TIFF CRS: {e}")
            # Fallback: try to make valid first
            try:
                zone_area = make_valid(zone_area)
                tiff_area = transform(transform_from_wgs84_to_tiff, zone_area)
                geoms = prepare_geometry_for_mask(tiff_area)
            except:
                warnings.append(f"Could not process zone {inner_radius_km}-{outer_radius_km}km")
                # Add placeholder result for this zone
                results.append({
                    "vulnerability_threshold": vulnerability,
                    "start_distance": inner_radius_km,
                    "end_distance": outer_radius_km,
                    "zone_population": 0,
                    "estimated_casualties": 0,
                    "error": "Could not process geometry"
                })
                # Skip to next zone
                continue

        try:
            # Read both population (band 1) and country FID (band 2)
            # Dynamic all_touched based on buffer size
            all_touched = outer_radius_km < 2.3  # Enable for small buffers
            
            out_image, out_transform = mask(
                src, 
                geoms, 
                crop=True, 
                nodata=src.nodata, 
                indexes=[1, 2],
                all_touched=all_touched  # Critical for small buffers
            )
            
            population_data = out_image[0]  # Band 1: Population
            country_fid_data = out_image[1]  # Band 2: Country FID
            
            valid_mask = population_data != src.nodata
            
            # Handle FID -1 (no country assigned)
            unassigned_mask = (country_fid_data == -1) & valid_mask
            if np.any(unassigned_mask):
                print(f"Found {np.sum(unassigned_mask)} pixels with FID -1 containing population")
                
                # Get the population in unassigned areas
                unassigned_population = np.sum(population_data[unassigned_mask])
                
                # Find all valid country FIDs in this zone
                valid_fids = np.unique(country_fid_data[(country_fid_data > 0) & valid_mask])
                
                if len(valid_fids) > 0:
                    # Distribute the unassigned population to the nearest valid country
                    # If we're at a coastal area, the most common valid FID is likely the right country
                    fid_counts = np.bincount(country_fid_data[(country_fid_data > 0) & valid_mask].astype(int))
                    most_common_fid = np.argmax(fid_counts) if len(fid_counts) > 0 else -1
                    
                    if most_common_fid > 0:
                        # Replace FID -1 with the most common country FID in the neighborhood
                        print(f"Reassigning FID -1 population ({unassigned_population}) to FID {most_common_fid}")
                        country_fid_data[unassigned_mask] = most_common_fid
            
            zone_population = int(np.sum(population_data[valid_mask]))
            
            # Calculate expected casualties for the zone as a whole
            zone_casualties = int(zone_population * vulnerability)
            total_casualties += zone_casualties

            # Process population by country FID
            unique_fids = np.unique(country_fid_data[valid_mask])
            country_breakdown = []
            
            for fid in unique_fids:
                if fid == src.nodata or fid <= 0:  # Skip nodata and FID -1 or other non-positive FIDs
                    continue
                    
                fid_int = int(fid)
                country_mask = (country_fid_data == fid) & valid_mask
                country_pop = int(np.sum(population_data[country_mask]))
                country_casualties = int(country_pop * vulnerability)
                country_name = COUNTRY_NAMES.get(fid_int, f"Country FID {fid_int}")
                
                country_breakdown.append({
                    "fid": fid_int,
                    "name": country_name,
                    "population": country_pop,
                    "casualties": country_casualties
                })
                
                # Add to overall country statistics using defaultdict
                if not countries_population[fid_int]["name"]: # Initialize name and fid if first time
                    countries_population[fid_int]["fid"] = fid_int
                    countries_population[fid_int]["name"] = country_name

                countries_population[fid_int]["total_population"] += country_pop
                countries_population[fid_int]["total_casualties"] += country_casualties
                countries_population[fid_int]["zone_breakdown"].append({
                    "vulnerability_threshold": vulnerability,
                    "start_distance": inner_radius_km,
                    "end_distance": outer_radius_km,
                    "population": country_pop,
                    "casualties": country_casualties
                })

            results.append({
                "vulnerability_threshold": vulnerability,
                "start_distance": inner_radius_km,
                "end_distance": outer_radius_km,
                "zone_population": zone_population,
                "estimated_casualties": zone_casualties,
                "country_breakdown": country_breakdown
            })

        except ValueError as e:
            results.append({
                "vulnerability_threshold": vulnerability,
                "start_distance": inner_radius_km,
                "end_distance": outer_radius_km,
                "zone_population": 0,
                "estimated_casualties": 0,
                "error": str(e)
            })

        # Update previous buffer for next iteration
        previous_buffer = outer_buffer
    
    # Sort countries by total casualties (descending)
    countries_list = list(countries_population.values())
    countries_list.sort(key=lambda x: x["total_casualties"], reverse=True)
    
    # Verify total matches
    country_total_population = sum(c["total_population"] for c in countries_list)
    country_total_casualties = sum(c["total_casualties"] for c in countries_list)

    zone_total_population = sum(z["zone_population"] for z in results)

    # Handle any discrepancy by adding residual to country with smallest casualties
    if countries_list and (abs(country_total_casualties - total_casualties) > 0 or 
                           abs(country_total_population - zone_total_population) > 0):
        # Find residuals
        casualties_residual = total_casualties - country_total_casualties
        population_residual = zone_total_population - country_total_population
        
        # Find country with smallest casualties count
        min_casualties_country = min(countries_list, key=lambda x: x["total_casualties"])
        
        # Add residuals to this country
        min_casualties_country["total_population"] += population_residual
        min_casualties_country["total_casualties"] += casualties_residual
        
        # Update the country's largest zone to include the residual for internal consistency
        if min_casualties_country["zone_breakdown"]:
            largest_zone = max(min_casualties_country["zone_breakdown"], 
                              key=lambda x: x["population"])
            largest_zone["population"] += population_residual
            largest_zone["casualties"] += casualties_residual
        
        print(f"Added residual of {population_residual} population and {casualties_residual} casualties to {min_casualties_country['name']}")

    # Add warning if there's still a significant discrepancy (should not happen after correction)
    if abs(country_total_casualties - total_casualties) > 1:
        warnings.append(f"Note: Minor discrepancy between total casualties ({total_casualties}) and sum of country casualties ({country_total_casualties})")

    return {
        "zones": results,
        "total_casualties": total_casualties,
        "note": "Population estimates based on a 2.5 arcminute resolution map",
        "warnings": warnings if warnings else None,
        "countries": countries_list
    }

def calculate_population_for_standard_zones(lat, lon, vulnerability_zones, max_distance_km=3000):
    """
    Calculate population in vulnerability zones under 3000km that don't cross the antimeridian
    """
    results = []
    total_casualties = 0
    
    # Sort zones by vulnerability threshold in descending order
    sorted_zones = sorted(vulnerability_zones, key=lambda x: x['threshold'], reverse=True)
    
    # Filter out zones beyond our max distance
    sorted_zones = [z for z in sorted_zones if z['end_distance'] <= max_distance_km]

    if not sorted_zones:
         return {
            "zones": [],
            "total_casualties": 0,
            "note": "No applicable zones within standard calculation limits.",
            "countries": [] # Added for consistency
        }

    # Create common transformers
    wgs84 = CRS.from_epsg(4326)
    
    # Get UTM zone
    zone_number = int((lon + 180) // 6 + 1)
    is_northern = lat >= 0
    utm_epsg = 32600 + zone_number if is_northern else 32700 + zone_number
    
    # Create coordinate transformers
    utm_crs = CRS.from_epsg(utm_epsg)
    transformer_to_utm = Transformer.from_crs(wgs84, utm_crs, always_xy=True)
    
    # Convert WGS84 to UTM
    utm_x, utm_y = transformer_to_utm.transform(lon, lat)
    utm_point = Point(utm_x, utm_y)
    
    # Use the globally loaded RASTER
    if RASTER is None:
        return {
            "zones": [],
            "total_casualties": 0,
            "note": "Error: Population data not available for standard calculation.",
            "countries": [] # Added for consistency
        }
    src = RASTER # Use the global raster object
    
    tiff_crs = src.crs
    
    # Set up transformer from UTM to tiff CRS
    transformer_to_tiff = Transformer.from_crs(utm_crs, tiff_crs, always_xy=True)
    
    # Pre-define transformation function
    def transform_geom_to_tiff_crs(x, y): # Renamed for clarity
        return transformer_to_tiff.transform(x, y)
    
    previous_buffer = None
    
    # Initialize countries_population for this function's scope
    countries_population_standard = defaultdict(lambda: {
        "fid": 0, "name": "", "total_population": 0, "total_casualties": 0, "zone_breakdown": []
    })

    for zone in sorted_zones:
        inner_radius_km = zone['start_distance']
        outer_radius_km = zone['end_distance']
        vulnerability = zone['threshold']

        # Create outer buffer in UTM coordinates
        outer_buffer = utm_point.buffer(outer_radius_km * 1000, resolution=64)
        
        # Calculate actual area for this vulnerability zone
        if previous_buffer:
            # Exclude area covered by higher vulnerability zones
            zone_area = outer_buffer.difference(previous_buffer)
            if not zone_area.is_valid: zone_area = make_valid(zone_area)
            if zone_area.is_empty or not zone_area.is_valid: # Fallback if difference fails
                 zone_area = outer_buffer.buffer(0).difference(previous_buffer.buffer(0))
                 if not zone_area.is_valid: zone_area = make_valid(zone_area)

        else:
            zone_area = outer_buffer

        if zone_area.is_empty: # Skip if zone area becomes empty
            results.append({
                "vulnerability_threshold": vulnerability,
                "start_distance": inner_radius_km,
                "end_distance": outer_radius_km,
                "zone_population": 0,
                "estimated_casualties": 0,
                "country_breakdown": [],
                "error": "Zone area became empty after difference."
            })
            previous_buffer = outer_buffer # Still update previous_buffer
            continue


        # Transform to TIFF's CRS
        try:
            tiff_area = transform(transform_geom_to_tiff_crs, zone_area)
        except Exception as e:
            results.append({
                "vulnerability_threshold": vulnerability,
                "start_distance": inner_radius_km,
                "end_distance": outer_radius_km,
                "zone_population": 0,
                "estimated_casualties": 0,
                "country_breakdown": [],
                "error": f"Failed to transform geometry: {e}"
            })
            previous_buffer = outer_buffer
            continue
        
        # Prepare geometry for masking
        # Using the existing prepare_geometry_for_mask function for consistency
        # Note: prepare_geometry_for_mask expects WGS84 input if it were to be reused directly from the other function.
        # Here, tiff_area is already in TIFF CRS.
        # So, a simplified version or direct construction is needed if tiff_area can be MultiPolygon.
        # Assuming tiff_area from UTM buffer difference is usually a Polygon or clean MultiPolygon.
        
        geoms_for_mask = []
        if isinstance(tiff_area, Polygon):
            poly_dict = {"type": "Polygon", "coordinates": [list(tiff_area.exterior.coords)]}
            if hasattr(tiff_area, 'interiors') and tiff_area.interiors:
                for interior in tiff_area.interiors:
                    poly_dict["coordinates"].append(list(interior.coords))
            geoms_for_mask.append(poly_dict)
        elif isinstance(tiff_area, MultiPolygon):
            for poly in tiff_area.geoms:
                poly_dict = {"type": "Polygon", "coordinates": [list(poly.exterior.coords)]}
                if hasattr(poly, 'interiors') and poly.interiors:
                    for interior in poly.interiors:
                        poly_dict["coordinates"].append(list(interior.coords))
                geoms_for_mask.append(poly_dict)
        else: # Should not happen with buffer operations
             results.append({
                "vulnerability_threshold": vulnerability,
                "start_distance": inner_radius_km,
                "end_distance": outer_radius_km,
                "zone_population": 0,
                "estimated_casualties": 0,
                "country_breakdown": [],
                "error": "Invalid geometry type for masking."
            })
             previous_buffer = outer_buffer
             continue


        # Extract data from raster
        try:
            # Dynamic all_touched based on buffer size
            all_touched = outer_radius_km < 2.3 # Enable for small buffers
            
            out_image, out_transform = mask(
                src,
                geoms_for_mask,
                crop=True,
                nodata=src.nodata,
                indexes=[1, 2],
                all_touched=all_touched  # Critical for small buffers
            )
        except Exception as e: # Catch rasterio mask errors
            results.append({
                "vulnerability_threshold": vulnerability,
                "start_distance": inner_radius_km,
                "end_distance": outer_radius_km,
                "zone_population": 0,
                "estimated_casualties": 0,
                "country_breakdown": [],
                "error": f"Raster masking failed: {e}"
            })
            previous_buffer = outer_buffer
            continue

        population_data = out_image[0]  # Band 1: Population
        country_fid_data = out_image[1]  # Band 2: Country FID
        
        valid_mask = (population_data != src.nodata) & (population_data is not None) # Added None check
        
        # Calculate zone population
        zone_population = int(np.sum(population_data[valid_mask]))
        
        # Calculate casualties
        zone_casualties = int(zone_population * vulnerability)
        total_casualties += zone_casualties
        
        # Process by country
        unique_fids = np.unique(country_fid_data[valid_mask])
        country_breakdown = []
        
        for fid_val in unique_fids:
            if fid_val == src.nodata or fid_val <= 0:  # Skip nodata values and non-positive FIDs
                continue
            
            fid_int = int(fid_val)
            country_mask = (country_fid_data == fid_int) & valid_mask
            country_pop = int(np.sum(population_data[country_mask]))
            country_cas = int(country_pop * vulnerability) # Renamed for clarity
            
            country_name = COUNTRY_NAMES.get(fid_int, f"Country FID {fid_int}")
            
            country_breakdown.append({
                "fid": fid_int,
                "name": country_name,
                "population": country_pop,
                "casualties": country_cas
            })
            
            # Update country statistics
            if not countries_population_standard[fid_int]["name"]:
                countries_population_standard[fid_int]["fid"] = fid_int
                countries_population_standard[fid_int]["name"] = country_name
            
            countries_population_standard[fid_int]["total_population"] += country_pop
            countries_population_standard[fid_int]["total_casualties"] += country_cas
            countries_population_standard[fid_int]["zone_breakdown"].append({
                "vulnerability_threshold": vulnerability,
                "start_distance": inner_radius_km,
                "end_distance": outer_radius_km,
                "population": country_pop,
                "casualties": country_cas
            })
        
        results.append({
            "vulnerability_threshold": vulnerability,
            "start_distance": inner_radius_km,
            "end_distance": outer_radius_km,
            "zone_population": zone_population,
            "estimated_casualties": zone_casualties,
            "country_breakdown": country_breakdown
        })

    countries_list_standard = sorted(
        [data for data in countries_population_standard.values() if data["fid"] != 0], 
        key=lambda x: x["total_casualties"], 
        reverse=True
    )

    return {
        "zones": results,
        "total_casualties": total_casualties,
        "note": "Population estimates based on a 2.5 arcminute resolution map (Simplified Calculation)",
        "countries": countries_list_standard # Added countries list
    }
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
    Main population calculation router that selects optimal processing method.
    
    Analyzes impact location and zone characteristics to determine whether to use
    simplified fast calculation for standard cases or comprehensive processing
    for complex geographical scenarios.
    
    Args:
        lat: Impact latitude in decimal degrees
        lon: Impact longitude in decimal degrees
        vulnerability_zones: List of vulnerability zone definitions
        
    Returns:
        dict: Population and casualty results from appropriate calculation method
    """
    # Analyze scenario complexity to select processing method
    max_radius = max(zone['end_distance'] for zone in vulnerability_zones)
    near_poles = abs(lat) > 84
    near_antimeridian = abs(lon) > 170
    
    # Route to appropriate calculation method based on complexity
    if near_poles or near_antimeridian or max_radius > 3000:
        # Use comprehensive method for complex geographical scenarios
        print(f"Using complex calculation method because: " + 
              (f"near poles ({lat})" if near_poles else "") +
              (f"near antimeridian ({lon})" if near_antimeridian else "") +
              (f"large radius ({max_radius}km)" if max_radius > 3000 else ""))
        
        return calculate_population_in_zones(lat, lon, vulnerability_zones)
    else:
        # Use optimized method for standard cases
        print(f"Using simplified calculation method (lat={lat}, lon={lon}, max_radius={max_radius}km)")
        return calculate_population_for_standard_zones(lat, lon, vulnerability_zones)

def cleanup():
    """
    Clean up resources and close open raster files to prevent memory leaks.
    
    Properly closes any open raster datasets to free system resources
    and prevent file handle accumulation during repeated calculations.
    """
    global RASTER
    
    # Close population raster file if currently open
    if RASTER is not None:
        RASTER.close()
        RASTER = None
        print("Closed raster file")