"""
ENEO Asteroid Impact Simulation Web Application

This Flask-based web application offers an interactive platform for simulating asteroid impact
events and analyzing their potential consequences on Earth. The system meticulously calculates
a range of impact effects, including thermal radiation, seismic activity, airblast waves,
ejecta distribution, and wind damage zones. Furthermore, it provides assessments of the
potential impact on human populations and economic infrastructures.

Key Features:
- Real-time asteroid impact simulation with comprehensive physics modeling
- Population impact assessment using global demographic datasets
- Economic damage calculation based on casualty estimates and GDP data
- Interactive visualization with geographic damage zone mapping
- NASA Sentry API integration for real Near-Earth Object (NEO) data
- Ocean depth querying for tsunami modeling capabilities
- Multi-hazard vulnerability analysis (thermal, seismic, blast, ejecta)

Technical Architecture:
- Flask web framework for RESTful API endpoints
- Geospatial coordinate system handling with antimeridian crossing support
- Spherical geometry calculations for accurate damage zone mapping
- Integration with multiple physics simulation modules
- Population and economic impact modeling

Author: Alexandros Notas
Thesis: Methods of prediction and assessment of asteroid impact consequences
Institution: National Technical University of Athens
Date: June 2025
"""

from flask import Flask, request, jsonify, render_template
from results import run_simulation_full
from population_calculator import calculate_population_in_zones
from shapely.geometry import Point, Polygon
import json
import math
from gdp_calculator import calculate_economic_damage
from utils import get_ocean_depth_from_geotiff, m_to_km
import numpy as np
import random
import requests

    
app = Flask(__name__)

@app.route('/simulate', methods=['POST'])
def simulate():
    """
    Serves as the primary endpoint for conducting asteroid impact simulations.
    
    This endpoint receives asteroid impact parameters, validates them, orchestrates
    the comprehensive simulation analysis, and returns a detailed report. The report
    encompasses calculated damage zones, estimated population impact, and economic
    damage assessments, providing a holistic view of the event's consequences.
    
    Expected JSON Input:
        diameter (float): Diameter of the asteroid in meters (must be ≥ 1.0).
        density (float, optional): Density of the asteroid in kg/m³ (range: 1000-8000, defaults to 3000).
        velocity (float): Entry velocity of the asteroid in km/s (range: 11-72).
        entry_angle (float): Angle of atmospheric entry in degrees (range: 15-90).
        distance (float): Reference distance for analysis in km (must be ≥ 1.0).
        latitude (float, optional): Latitude of the impact point in degrees (defaults to 0).
        longitude (float, optional): Longitude of the impact point in degrees (defaults to 0).
    
    Returns:
        JSON: A structured response containing:
            "results_text" (str): A human-readable summary of the simulation results.
            "results_data" (dict): A dictionary with detailed data, including:
                - vulnerability_analysis: Calculations for various damage zones.
                - population_analysis: Data on the population affected within these zones.
                - economic_analysis: Estimates of economic damage.
                - visualization: Data formatted for geographic visualization.
                - country_visualization: Impact data aggregated by affected countries.
        
    Raises:
        HTTP 400: If input parameters are invalid or fail validation checks.
        HTTP 500: If an internal error occurs during the simulation process.
    """
    try:
        # Parse JSON data from the incoming request
        data = request.get_json()
        
        # Extract and validate simulation parameters, applying defaults where necessary
        diameter = float(data['diameter'])
        density = float(data.get('density', 3000))  # Default: typical rocky asteroid density
        velocity = float(data['velocity'])
        entry_angle = float(data['entry_angle'])
        r_distance = float(data['distance'])
        lat = float(data.get('latitude', 0))    # Default: equator if not provided
        lon = float(data.get('longitude', 0))   # Default: prime meridian if not provided
        
        # Log simulation parameters for debugging and analysis tracking
        print(f"\nSimulation Request - Coordinates: {lat:.6f}°, {lon:.6f}°")
        print(f"Parameters: D={diameter}m, ρ={density}kg/m³, v={velocity}km/s, θ={entry_angle}°")
        
        # Comprehensive input validation based on established physical constraints
        if diameter < 1:
            return jsonify({"error": "Diameter must be 1 meter or greater."}), 400
        if not (1000 <= density <= 8000):
            return jsonify({"error": "Density must be between 1000 and 8000 kg/m³."}), 400
        if not (11 <= velocity <= 72):
            return jsonify({"error": "Entry velocity must be between 11 and 72 km/s."}), 400
        if not (15 <= entry_angle <= 90):
            return jsonify({"error": "Entry angle must be between 15 and 90 degrees."}), 400
        if r_distance < 1:
            return jsonify({"error": "Distance must be 1 km or greater."}), 400
            
    except (KeyError, ValueError):
        # Handle cases where required parameters are missing or have incorrect types
        return jsonify({"error": "Invalid input. Please provide numeric values for all required parameters."}), 400

    # Execute the core impact simulation with validated parameters
    results_text, results_data = run_simulation_full(diameter, density, velocity, entry_angle, r_distance, lat, lon)
    
    # Calculate population impact if geographic coordinates are provided
    if lat != 0 or lon != 0:
        # Extract vulnerability zones from simulation results for population analysis
        vulnerability_zones = results_data.get('vulnerability_analysis', {}).get('zones', [])
        
        # Calculate affected population within each vulnerability zone
        population_data = calculate_population_in_zones(lat, lon, vulnerability_zones)
        results_data['population_analysis'] = population_data

        # Generate country-specific visualization if population data includes country breakdown
        if 'countries' in population_data:
            results_data['country_visualization'] = {
                'countries': population_data['countries']
            }
            
            # Calculate economic impact based on estimated casualties and national GDP data
            economic_data = calculate_economic_damage(population_data['countries'])
            results_data['economic_analysis'] = economic_data

    # Generate visualization data for mapping interface
    visualization_data = generate_visualization_data(
        lat, lon, results_data, 
        results_text, 
        diameter, density, velocity, entry_angle
    )
    results_data['visualization'] = visualization_data

    # Return complete simulation results with text summaries and structured data
    return jsonify({
        "results_text": results_text,
        "results_data": results_data
    })

@app.route('/generate-neo', methods=['GET'])
def generate_neo_endpoint():
    """
    Fetch real Near-Earth Object (NEO) data from NASA's Sentry API for simulation.
    
    This endpoint connects to NASA's Sentry system, which tracks potentially hazardous
    asteroids with non-zero impact probabilities. It retrieves actual asteroid data
    and formats it for use in impact simulations, providing realistic scenarios
    based on known objects in space.
    
    Process:
    1. Query NASA Sentry API for list of tracked objects
    2. Filter for objects larger than 50m diameter (significant impact potential)
    3. Randomly select one object for variety in simulations
    4. Fetch detailed parameters for the selected object
    5. Format data with reasonable assumptions for missing parameters
    
    Returns:
        JSON: Formatted asteroid parameters ready for simulation:
            {
                "name": str,           # Official asteroid designation
                "diameter": float,     # Diameter in meters
                "density": float,      # Density in kg/m³ (assumed value)
                "velocity": float,     # Impact velocity in km/s
                "entry_angle": float   # Entry angle in degrees (random assumption)
            }
    
    Error Responses:
        HTTP 503: Service unavailable (NASA API connection issues)
        HTTP 404: No suitable objects found in Sentry database
        HTTP 500: Data parsing or processing errors
    """
    sentry_api_url = "https://ssd-api.jpl.nasa.gov/sentry.api"
    try:
        # Query NASA Sentry API for complete list of tracked objects
        response_summary = requests.get(sentry_api_url, timeout=15)
        response_summary.raise_for_status()  # Raise exception for HTTP errors
        all_objects = response_summary.json()

        if not all_objects.get("data"):
            return jsonify({"error": "No data returned from Sentry API summary."}), 503

        # Filter for objects with significant impact potential (>50m diameter)
        large_objects = [
            obj for obj in all_objects.get("data", [])
            if float(obj.get("diameter", "0")) * 1000 > 50  # Convert km to m
        ]

        if not large_objects:
            return jsonify({"error": "No suitable objects (>50m) found in Sentry data."}), 404

        # Randomly select an object to provide variety in simulations
        random_object_summary = random.choice(large_objects)
        object_des = random_object_summary["des"]  # Object designation

        # Fetch detailed data for the selected asteroid
        params = {"des": object_des}
        response_object = requests.get(sentry_api_url, params=params, timeout=15)
        response_object.raise_for_status()
        object_details = response_object.json()

        if "error" in object_details or "summary" not in object_details:
            return jsonify({"error": f"Could not retrieve details for object {object_des}."}), 503

        # Extract and process simulation parameters from API response
        summary = object_details.get("summary", {})
        
        # Convert diameter from km to meters, with fallback for invalid data
        try:
            diameter_m = float(summary.get("diameter", "0")) * 1000
        except (ValueError, TypeError):
            diameter_m = 100  # Default to 100m if data is invalid

        # Extract impact velocity with reasonable default
        v_imp = summary.get("v_imp") or 20

        # Use typical rocky asteroid density (Sentry assumption)
        density = 2600
        
        # Generate random entry angle since API doesn't provide this parameter
        assumed_entry_angle = random.uniform(20, 60)

        # Return formatted parameters for simulation use
        return jsonify({
            "name": summary.get("fullname", object_des),
            "diameter": round(diameter_m, 2),
            "density": density,
            "velocity": round(float(v_imp), 2),
            "entry_angle": round(assumed_entry_angle, 2)
        })

    except requests.exceptions.RequestException as e:
        # Handle network connectivity and HTTP errors
        print(f"Error fetching Sentry data: {e}")
        return jsonify({"error": "Could not connect to NASA Sentry API."}), 503
    except (KeyError, IndexError, TypeError, ValueError) as e:
        # Handle data structure and parsing errors
        print(f"Error parsing Sentry API response: {e}")
        return jsonify({"error": "Failed to parse data from NASA Sentry API."}), 500
    except Exception as e:
        # Catch-all for unexpected errors with detailed logging
        print(f"An unexpected error occurred in /generate-neo: {e}")
        return jsonify({"error": "An internal server error occurred."}), 500

@app.route('/')
def index():
    """
    Serve the main application interface.
    
    Renders the primary HTML template that provides the interactive user interface
    for asteroid impact simulations. This includes parameter input forms,
    visualization controls, and results display areas.
    
    Returns:
        HTML: The main application interface template
    """
    return render_template('index.html')

def create_circle_coordinates(center_lat, center_lon, radius_km, points=72):
    """
    Generates geographic coordinates for circular damage zones, carefully considering Earth's spherical geometry.
    
    This function produces accurate circular boundaries on Earth's surface by addressing key complexities:
    - Earth's curvature, managed using great circle calculations.
    - Antimeridian crossing (longitude ±180°), for circles spanning the international date line.
    - Special handling for polar regions to ensure accuracy at high latitudes.
    - Normalization and validation of coordinate systems.
    
    The function uses spherical trigonometry to compute points at a consistent distance from the impact center,
    aiming for a precise representation of damage zones globally.
    
    Args:
        center_lat (float): Latitude of the circle's center, in decimal degrees (-90 to 90).
        center_lon (float): Longitude of the circle's center, in decimal degrees (-180 to 180).
        radius_km (float): Radius of the circle in kilometers.
        points (int, optional): Number of points to generate for the circle's perimeter 
                                (default is 72, providing a point every 5 degrees).
    
    Returns:
        list: A list of coordinates defining the circle. The format depends on the circle's characteristics:
            - For standard circles: A single list of [longitude, latitude] pairs.
            - For circles crossing the antimeridian: A list of two lists, representing multi-polygons.
            - An empty list ([]) if the circle cannot be generated due to invalid parameters.
    
    Mathematical Basis:
        Calculations are based on the spherical law of cosines, fundamental for great circle navigation:
        lat2 = asin(sin(lat1) * cos(d/R) + cos(lat1) * sin(d/R) * cos(bearing))
        lon2 = lon1 + atan2(sin(bearing) * sin(d/R) * cos(lat1), cos(d/R) - sin(lat1) * sin(lat2))
        
        Where:
        - d/R is the angular distance (radius_km / Earth_radius).
        - R is Earth's mean radius (6371 km).
        - bearing is the azimuth angle for each point on the circle.
    """
    from math import sin, cos, pi, asin, radians, degrees, atan2

    # Validate and clamp latitude to valid geographic range
    center_lat = max(-90, min(90, center_lat))
    
    # Use specialized polar coordinate generation for high latitudes
    if abs(center_lat) > 85:
        print(f"Using pole-aware circle generation for latitude {center_lat}")
        return create_polar_circle_coordinates(center_lat, center_lon, radius_km, points)
    
    # Earth's mean radius and angular distance calculation
    R = 6371  # Earth's mean radius in kilometers
    angular_distance = radius_km / R  # Convert linear to angular distance
    
    coordinates = []
    crosses_antimeridian = False
    east_coords = []  # Coordinates east of the antimeridian (positive longitude)
    west_coords = []  # Coordinates west of the antimeridian (negative longitude)
    
    # Generate circle points using spherical trigonometry
    for i in range(points + 1):
        bearing = radians(i * (360 / points))  # Azimuth angle for current point
        lat1 = radians(center_lat)
        lon1 = radians(center_lon)
        
        # Calculate destination point latitude using spherical law of cosines
        lat2 = asin(sin(lat1) * cos(angular_distance) + 
                   cos(lat1) * sin(angular_distance) * cos(bearing))
        
        # Handle longitude calculation with pole protection
        if abs(cos(lat2)) < 1e-10:
            lon2 = lon1  # At poles, longitude is undefined
        else:
            # Calculate destination longitude
            lon2 = lon1 + atan2(sin(bearing) * sin(angular_distance) * cos(lat1),
                                cos(angular_distance) - sin(lat1) * sin(lat2))
        
        # Convert from radians to degrees
        lat2 = degrees(lat2)
        lon2 = degrees(lon2)
        
        # Normalize longitude to [-180, 180] range
        while lon2 > 180:
            lon2 -= 360
        while lon2 < -180:
            lon2 += 360
        
        # Detect antimeridian crossing (longitude jump > 180°)
        if i > 0 and not crosses_antimeridian:
            prev_lon = coordinates[-1][0] if coordinates else None
            if prev_lon is not None and abs(prev_lon - lon2) > 180:
                crosses_antimeridian = True
                
                # Redistribute existing coordinates by hemisphere
                for j in range(len(coordinates)):
                    if coordinates[j][0] >= 0:
                        east_coords.append(coordinates[j])
                    else:
                        west_coords.append(coordinates[j])
                
                # Calculate antimeridian intersection points for clean split
                if prev_lon > 0 and lon2 < 0:
                    # Crossing from east (+180) to west (-180)
                    t = (180 - prev_lon) / ((180 - prev_lon) + (180 + lon2))
                    y_inter = coordinates[-1][1] + t * (lat2 - coordinates[-1][1])
                    
                    east_coords.append([180, y_inter])
                    west_coords.append([-180, y_inter])
                else:
                    # Crossing from west (-180) to east (+180)
                    t = (-180 - prev_lon) / ((-180 - prev_lon) + (180 - lon2))
                    y_inter = coordinates[-1][1] + t * (lat2 - coordinates[-1][1])
                    
                    west_coords.append([-180, y_inter])
                    east_coords.append([180, y_inter])
        
        # Add new point to appropriate coordinate set
        if crosses_antimeridian:
            if lon2 >= 0:
                east_coords.append([lon2, lat2])
            else:
                west_coords.append([lon2, lat2])
        else:
            coordinates.append([lon2, lat2])
    
    # Handle antimeridian-crossing circles with multi-polygon structure
    if crosses_antimeridian:
        # Close polygons by adding first point at the end
        if east_coords and east_coords[0] != east_coords[-1]:
            east_coords.append(east_coords[0].copy())
        if west_coords and west_coords[0] != west_coords[-1]:
            west_coords.append(west_coords[0].copy())
            
        # Ensure minimum polygon validity (4 points including closure)
        if len(east_coords) < 4:
            east_coords = []
        if len(west_coords) < 4:
            west_coords = []
            
        # Return valid coordinate sets as multi-polygon
        result = []
        if east_coords:
            result.append(east_coords)
        if west_coords:
            result.append(west_coords)
            
        return result
    
    # Close single polygon for standard circles
    if coordinates and coordinates[0] != coordinates[-1]:
        coordinates.append(coordinates[0].copy())
        
    return coordinates

def create_polar_circle_coordinates(center_lat, center_lon, radius_km, points=72):
    """
    Generates coordinates for circular damage zones near Earth's poles, effectively managing zones that may cross a pole.
    
    This function is designed to handle the specific geometric challenges of creating circular zones at high latitudes,
    where standard calculations can be less accurate. It uses pole-aware algorithms for a precise depiction of damage zones.
    
    Key features include:
    - Detection and management of pole-crossing scenarios.
    - Coordinate reflection for zones that extend across a pole.
    - Compensation for longitude shifts in regions that have crossed a pole.
    - Validation for extreme polar cases (latitudes > 89.9° or < -89.9°).
    
    Args:
        center_lat (float): Latitude of the circle's center (typically > 85° or < -85°).
        center_lon (float): Longitude of the circle's center.
        radius_km (float): Radius of the circle in kilometers.
        points (int, optional): Number of points for the circle's perimeter (default is 72).
    
    Returns:
        list: Coordinates for the polar circle. The format is:
            - For non-pole-crossing circles: A single list of [longitude, latitude] pairs.
            - For pole-crossing circles: A list of two lists (multi-polygon format).
            - An empty list ([]) for invalid configurations (e.g., radius too large).
    
    Mathematical Approach:
        If a circle crosses a pole, points extending beyond the pole are geometrically reflected:
        - Latitude reflection: lat_reflected = pole_lat - (lat_calculated - pole_lat)
        - Longitude shift: lon_reflected = lon_calculated + 180°
        This method ensures a continuous and accurate representation of the zone across polar regions.
    """
    from math import sin, cos, radians, degrees, pi, asin, atan2
    
    # Earth's mean radius and validation
    R_EARTH = 6371
    
    # Prevent unrealistic global-spanning circles
    if radius_km > 0.95 * R_EARTH:
        print(f"Warning: Polar circle radius {radius_km} km exceeds 95% of Earth's radius - skipping visualization")
        return []
    
    # Determine closest pole and adjust for mathematical stability
    is_north_pole = center_lat > 0
    pole_lat = 90 if is_north_pole else -90
    
    # Adjust extremely polar positions to avoid singularities
    if abs(abs(center_lat) - 90) < 0.1:
        center_lat = 89.9 if is_north_pole else -89.9
    
    # Calculate circle characteristics and pole-crossing potential
    R = 6371  # Earth's radius in kilometers
    angular_distance = radius_km / R  # Angular radius in radians
    
    # Determine if circle extends beyond the pole
    lat_extent = degrees(angular_distance) - abs(90 - abs(center_lat))
    crosses_pole = lat_extent > 0
    
    print(f"Polar circle analysis: center={center_lat}°, radius={radius_km}km, crosses_pole={crosses_pole}")
    
    # Handle pole-crossing circles with hemisphere splitting
    if crosses_pole:
        this_side_coords = []  # Original hemisphere coordinates
        other_side_coords = [] # Opposite hemisphere coordinates after pole crossing
        
        # Generate points using standard spherical calculations
        for i in range(points + 1):
            bearing = radians(i * (360 / points))
            lat1 = radians(center_lat)
            lon1 = radians(center_lon)
            
            # Calculate point at given distance and bearing
            lat2 = asin(sin(lat1) * cos(angular_distance) + 
                       cos(lat1) * sin(angular_distance) * cos(bearing))
            
            # Handle longitude calculation with pole protection
            if abs(cos(lat2)) < 1e-10:
                lon2 = lon1  # Longitude undefined at exact pole
            else:
                lon2 = lon1 + atan2(sin(bearing) * sin(angular_distance) * cos(lat1),
                                    cos(angular_distance) - sin(lat1) * sin(lat2))
            
            # Convert to degrees
            lat2 = degrees(lat2)
            lon2 = degrees(lon2)
            
            # Normalize longitude
            while lon2 > 180:
                lon2 -= 360
            while lon2 < -180:
                lon2 += 360
            
            # Check for pole crossing and reflect coordinates if needed
            if (is_north_pole and lat2 > 90) or (not is_north_pole and lat2 < -90):
                # Reflect point to opposite hemisphere
                excess = lat2 - 90 if is_north_pole else -90 - lat2
                reflected_lat = 90 - excess if is_north_pole else -90 + excess
                
                # Shift longitude by 180° for opposite hemisphere
                other_lon = lon2 + 180
                while other_lon > 180:
                    other_lon -= 360
                    
                other_side_coords.append([other_lon, reflected_lat])
            else:
                # Point remains on original hemisphere
                this_side_coords.append([lon2, lat2])
        
        # Close polygons for valid GeoJSON representation
        if this_side_coords and this_side_coords[0] != this_side_coords[-1]:
            this_side_coords.append(this_side_coords[0].copy())
        if other_side_coords and other_side_coords[0] != other_side_coords[-1]:
            other_side_coords.append(other_side_coords[0].copy())
            
        # Ensure polygon validity (minimum 4 points)
        if len(this_side_coords) < 4:
            this_side_coords = []
        if len(other_side_coords) < 4:
            other_side_coords = []
            
        # Return valid coordinate sets as multi-polygon
        result = []
        if this_side_coords:
            result.append(this_side_coords)
        if other_side_coords:
            result.append(other_side_coords)
        
        print(f"Generated {len(result)} polar coordinate sets")
        return result
    
    # For non-pole-crossing circles, create simple constant-latitude circle
    lat_offset = degrees(angular_distance)
    lat_value = center_lat + lat_offset if is_north_pole else center_lat - lat_offset
    
    # Clamp latitude to valid range
    lat_value = max(-90, min(90, lat_value))
    
    # Generate circle of constant latitude
    coordinates = []
    for i in range(points + 1):
        lon = center_lon + (i * 360 / points)
        while lon > 180:
            lon -= 360
        while lon < -180:
            lon += 360
        coordinates.append([lon, lat_value])
    
    return coordinates

def parse_danger_zones(section_text, lat, lon):
    """
    Extracts and categorizes information about various damage zones from the simulation's textual output.
    
    This function processes structured text to identify different types of hazard zones (e.g., thermal, seismic, airblast)
    and their respective distance ranges. This information is vital for visualizing the spatial extent of these hazards.
    It uses keyword matching and pattern recognition to automatically classify zones based on their descriptive text.
    
    Zone types are classified as follows:
    - Wind zones: Identified by an 'EF' prefix (Enhanced Fujita scale).
    - Tsunami zones: Recognized by amplitude thresholds (e.g., >1km, >100m) and context.
    - Seismic zones: Identified by the keyword 'Richter'.
    - Thermal zones: Detected by keywords like 'burns', 'ignition', or 'Thermal'.
    - Ejecta zones: Characterized by thickness patterns (e.g., ">X m").
    - Airblast zones: Inferred from descriptions of structural damage.
    
    Args:
        section_text (str): The raw text segment from the simulation output describing the zones.
        lat (float): Latitude of the impact point, for generating zone coordinates.
        lon (float): Longitude of the impact point, for generating zone coordinates.
    
    Returns:
        list: A list of parsed zone objects, each containing coordinates, description, distance, and type.
    
    Expected Text Format:
        The input text should contain lines with patterns like: "Description: X-Y km" or "Description: X km".
    """
    zones = []
    lines = section_text.split('\n')

    for line_content in lines:
        line = line_content.strip()
        
        # Skip empty lines and lines without proper format
        if not line or ':' not in line:
            continue

        # Split line into description and distance range components
        description_part = line.split(':', 1)[0].strip()
        range_part_with_km = line.split(':', 1)[-1].strip()

        determined_type_for_line = None

        # Classify zone type using pattern matching and keywords
        
        # Wind zones: Enhanced Fujita scale identification
        if description_part.startswith('EF'):
            determined_type_for_line = 'wind'
            
        # Tsunami zones: Amplitude thresholds in tsunami context
        elif description_part in [">1km", ">100m", ">10m", ">1m"] and "tsunami" in section_text.lower():
             if "km" in range_part_with_km and "None" not in range_part_with_km:
                determined_type_for_line = 'tsunami'
                
        # Seismic zones: Richter scale references
        elif 'Richter' in description_part:
            determined_type_for_line = 'seismic'
        
        # Thermal zones: Heat and fire-related keywords
        elif any(keyword.lower() in description_part.lower() for keyword in 
                ['burns', 'ignition', 'Thermal']):
            determined_type_for_line = 'thermal'
        
        # Ejecta zones: Thickness measurements (e.g., ">10 m")
        elif description_part.startswith(">") and description_part.endswith("m"):
            determined_type_for_line = 'ejecta'
        
        # Airblast zones: Structural damage descriptions
        elif any(keyword.lower() in description_part.lower() for keyword in 
                ["collapse", "distorted", "damage", "shatter", "structural", 
                 "windows", "buildings", "office-type", "wall-bearing", "wood frame"]):
            determined_type_for_line = 'airblast'
        
        # Additional tsunami detection for wave amplitude descriptions
        elif "wave amplitude >" in description_part.lower() or \
             ("tsunami" in section_text.lower() and "zone" in description_part.lower()):
            if "km" in line.split(':', 1)[-1]:
                determined_type_for_line = 'tsunami'

        # Process zones with determined types
        if determined_type_for_line:
            try:
                dist_text = line.split(':', 1)[1].strip()

                # Skip zones without valid distance data
                if 'None' not in dist_text:
                    start_dist_val, end_dist_val = 0.0, 0.0
                    
                    # Parse distance ranges
                    if '-' in dist_text:  # Range format "X-Y km"
                        parts = dist_text.replace('km', '').strip().split('-', 1)
                        start_dist_val = float(parts[0])
                        end_dist_val = float(parts[1] if len(parts) > 1 else parts[0])
                    else:  # Single value format "X km" (implies 0 to X)
                        end_dist_val = float(dist_text.replace('km', '').strip())
                        start_dist_val = 0.0

                    # Validate zone significance and generate coordinates
                    if end_dist_val > start_dist_val or (abs(start_dist_val) < 1e-9 and end_dist_val > 0.01):
                        # Generate geographic coordinates for zone boundary
                        circle_coords = create_circle_coordinates(lat, lon, end_dist_val)
                        
                        # Validate coordinate generation success
                        is_valid_coords = False
                        if circle_coords:
                            if isinstance(circle_coords, list) and len(circle_coords) > 0:
                                if isinstance(circle_coords[0], list):  # Multi-polygon case
                                    is_valid_coords = len(circle_coords[0]) > 0
                                else:  # Single polygon case
                                    if len(circle_coords) > 2 and isinstance(circle_coords[0], list):
                                        is_valid_coords = len(circle_coords[0]) == 2

                        # Add valid zones to results
                        if is_valid_coords:
                            zones.append({
                                'coordinates': circle_coords,
                                'description': description_part,
                                'start_distance': start_dist_val,
                                'end_distance': end_dist_val,
                                'type': determined_type_for_line
                            })
                            
            except ValueError as ve:
                print(f"ValueError parsing zone line '{line}' for type '{determined_type_for_line}': {ve}")
            except Exception as e:
                print(f"Error parsing zone line '{line}' for type '{determined_type_for_line}': {e}")
                
    return zones

def generate_visualization_data(lat, lon, results_data, results_text, diameter, density, velocity, entry_angle):
    """
    Prepares a comprehensive dataset for visualizing impact effects on a map.
    
    This function processes the simulation's outputs to create geographically referenced 
    data for identified damage zones, vulnerability assessments, and other impact 
    characteristics. It combines parsed text results (for basic damage zones)
    with more detailed vulnerability analyses derived from calculated thresholds and distances.
    
    Data is generated for several visualization layers:
    1. Basic damage zones (from textual simulation results).
    2. Vulnerability zones (from physical simulation models).
    3. Combined vulnerability assessments (integrating multiple hazard types).
    4. Type-specific vulnerability zones (e.g., thermal, seismic).
    
    Args:
        lat (float): Latitude of the impact point.
        lon (float): Longitude of the impact point.  
        results_data (dict): Structured data output from the simulation.
        results_text (str): Human-readable textual output from the simulation.
        diameter (float): Asteroid diameter in meters.
        density (float): Asteroid density in kg/m³.
        velocity (float): Asteroid entry velocity in km/s.
        entry_angle (float): Atmospheric entry angle in degrees.
    
    Returns:
        dict: A data structure ready for visualization, organized as:
            {
                'seismic': [zone_objects],           // Seismic damage zones.
                'thermal': [zone_objects],           // Thermal damage zones.
                'airblast': [zone_objects],          // Airblast damage zones.
                'ejecta': [zone_objects],            // Ejecta deposit zones.
                'wind': [zone_objects],              // Wind damage zones.
                'tsunami': [zone_objects],           // Tsunami zones.
                'vulnerability': [zone_objects],      // Combined vulnerability zones.
                'vulnerability_combined': [zone_objects], // Explicitly combined vulnerability zones.
                'vulnerability_thermal': [zone_objects],  // Thermal-specific vulnerability zones.
                'vulnerability_overpressure': [zone_objects], // Overpressure-specific vulnerability zones.
                'vulnerability_seismic': [zone_objects],     // Seismic-specific vulnerability zones.
                'vulnerability_ejecta': [zone_objects],      // Ejecta-specific vulnerability zones.
                'vulnerability_wind': [zone_objects]         // Wind-specific vulnerability zones.
            }
    
    Each zone object has the following structure:
        {
            'coordinates': list,              // Geographic coordinates for the boundary.
            'description': str,               // Human-readable description.
            'start_distance': float,          // Inner radius in km.
            'end_distance': float,           // Outer radius in km.
            'type': str,                     // Category of the zone.
            'threshold': float,              // Vulnerability threshold (for vulnerability zones).
            'specific_vuln_type': str        // Specific type of vulnerability (e.g., 'thermal').
        }
    """
    # Initialize comprehensive visualization data structure
    visualization = {
        'seismic': [], 'thermal': [], 'airblast': [], 'ejecta': [], 'wind': [], 'tsunami': [],
        'vulnerability': [],  # Combined vulnerability zones
        'vulnerability_combined': [],  # Explicitly combined zones
        'vulnerability_thermal': [],   # Type-specific vulnerability zones
        'vulnerability_overpressure': [],
        'vulnerability_seismic': [],
        'vulnerability_ejecta': [],
        'vulnerability_wind': []
    }

    # Parse basic damage zones from structured text results
    split_parts = results_text.split(r'===')  # Split by section delimiters

    # Process each section title-content pair
    for i in range(1, len(split_parts), 2):
        section_title_full = split_parts[i].strip()
        section_content = split_parts[i+1].strip() if (i+1) < len(split_parts) else ""
        
        # Extract zones from section content
        parsed_zones_from_content = parse_danger_zones(section_content, lat, lon)
        
        # Categorize parsed zones into visualization structure
        for zone in parsed_zones_from_content:
            zone_type_from_parser = zone.get('type')
            if zone_type_from_parser and zone_type_from_parser in visualization:
                visualization[zone_type_from_parser].append(zone)
            elif zone_type_from_parser:
                print(f"Debug: Parsed zone type '{zone_type_from_parser}' not in standard categories")

    # Generate advanced vulnerability analysis zones from simulation data
    if 'vulnerability_analysis' in results_data and 'zones' in results_data['vulnerability_analysis']:
        from models import AsteroidImpactSimulation
        
        # Initialize simulation model for detailed vulnerability calculations
        sim = AsteroidImpactSimulation(diameter, velocity, density, entry_angle)
        entry_results = sim.simulate_atmospheric_entry()

        # Determine thermal vulnerability display based on impact velocity
        show_thermal_vulnerability = True
        if entry_results["event_type"] == "ground impact":
            v_check = entry_results["post_breakup_velocity"]
        elif entry_results["event_type"] == "airburst":
            v_check = entry_results["v_breakup"]
        else:
            v_check = entry_results["post_breakup_velocity"]

        # Skip thermal zones for low-velocity impacts
        if m_to_km(v_check) < 15.0:
            show_thermal_vulnerability = False

        # Generate vulnerability threshold range (high to low vulnerability)
        threshold_values = [round(x, 2) for x in np.arange(1.0, 0.04, -0.05)]

        # Process type-specific vulnerability zones
        specific_vuln_types = ['thermal', 'overpressure', 'seismic', 'ejecta', 'wind']

        for vuln_type in specific_vuln_types:
            # Skip thermal zones if velocity condition not met
            if vuln_type == 'thermal' and not show_thermal_vulnerability:
                continue

            previous_max_distance_for_type = 0.0
            visualization_key = f'vulnerability_{vuln_type}'

            # Generate concentric vulnerability zones for each threshold level
            for threshold in threshold_values:
                if threshold < 0.05:  # Skip very low thresholds for clarity
                    continue

                # Calculate maximum distance for current vulnerability level
                current_max_distance = find_specific_vulnerability_distance(
                    sim, entry_results, vuln_type, threshold
                )

                # Skip negligible zones
                if previous_max_distance_for_type == 0.0 and current_max_distance <= 0.01:
                    continue

                # Create zone if distance increased from previous threshold
                if current_max_distance > previous_max_distance_for_type:
                    circle_coords = create_circle_coordinates(lat, lon, current_max_distance, 72)

                    if circle_coords:  # Validate coordinate generation
                        zone_info = {
                            'coordinates': circle_coords,
                            'description': f"{vuln_type.capitalize()}-Induced Vulnerability ≥ {threshold:.2f}",
                            'threshold': threshold,
                            'start_distance': previous_max_distance_for_type,
                            'end_distance': current_max_distance,
                            'specific_vuln_type': vuln_type 
                        }
                        visualization[visualization_key].append(zone_info)
                        previous_max_distance_for_type = current_max_distance
        
        # Process combined vulnerability zones from existing analysis
        vuln_zones_data = results_data['vulnerability_analysis']['zones']
        
        # Sort zones by threshold (descending) for consistent processing
        sorted_vuln_zones = sorted(
            vuln_zones_data, 
            key=lambda x: (-x.get('threshold', 0), -x.get('end_distance', 0))
        )
        
        # Generate combined vulnerability visualization zones
        for zone_data in sorted_vuln_zones:
            if zone_data.get('threshold', 0) >= 0.05 and zone_data.get('end_distance', 0) > 0.01:
                circle_coords = create_circle_coordinates(lat, lon, zone_data['end_distance'], 72)
                
                if circle_coords:  # Validate coordinates
                    zone_info = {
                        'coordinates': circle_coords,
                        'description': f"Combined Vulnerability ≥ {zone_data['threshold']:.2f}",
                        'threshold': zone_data['threshold'],
                        'start_distance': zone_data['start_distance'],
                        'end_distance': zone_data['end_distance']
                    }
                    visualization['vulnerability'].append(zone_info)
                    visualization['vulnerability_combined'].append(zone_info)
    
    return visualization

def find_specific_vulnerability_distance(sim, entry_results, vuln_type, threshold, r_min=0.01, r_max=20000.0, tol=0.01):
    """
    Calculates the maximum distance at which a specific vulnerability type 
    (e.g., thermal, seismic) meets or exceeds a defined threshold.
    
    This function uses a binary search algorithm to efficiently find this boundary 
    distance. Beyond this distance, the vulnerability's intensity drops below the 
    given threshold. It accounts for different impact scenarios (ground impacts vs. airbursts) 
    and calculates relevant vulnerability metrics for each distance evaluated.
    
    The binary search converges quickly, making it suitable for dynamic calculations.
    
    Args:
        sim (AsteroidImpactSimulation): An initialized simulation object.
        entry_results (dict): Results from the atmospheric entry simulation, including
            event_type, post_breakup_velocity, airburst_altitude, v_breakup, and z_star.
        vuln_type (str): The type of vulnerability to assess (e.g., 'thermal', 'overpressure').
        threshold (float): The vulnerability threshold value (0.0 to 1.0).
        r_min (float, optional): Minimum search distance in km (default: 0.01 km).
        r_max (float, optional): Maximum search distance in km (default: 20000.0 km).
        tol (float, optional): Tolerance for search convergence in km (default: 0.01 km).
    
    Returns:
        float: The maximum distance (km) where vulnerability is ≥ threshold. 
               Returns r_min if the threshold isn't met even at close distances.
    
    Mathematical Approach:
        The binary search maintains these conditions:
        - The left boundary (left) of the search interval represents a distance where 
          vulnerability(distance) ≥ threshold.
        - The right boundary (right) represents a distance where vulnerability(distance) < threshold.
        
        The search stops when the interval (right - left) is smaller than the tolerance.
    
    Details on vulnerability calculation:
        - For Ground Impacts: Uses models based on impact energy and crater formation.
        - For Airbursts: Employs models based on burst energy and atmospheric blast wave propagation.
        - Returns 0.0 if a vulnerability type isn't applicable (e.g., ejecta from an airburst).
    """
    from utils import km_to_m, rho_target
    
    def calculate_specific_vulnerability_at_distance(r_km, vuln_type):
        """Calculate vulnerability value at specific distance for given hazard type."""
        D_m = km_to_m(r_km)  # Convert distance to meters
        
        # Ground impact vulnerability calculations
        if entry_results["event_type"] == "ground impact":
            v_surface = entry_results['post_breakup_velocity']
            imp_energy, _ = sim.calculate_impact_energy(v_surface)
            
            if vuln_type == 'thermal':
                phi_ground = sim.calculate_thermal_exposure(imp_energy, r_km)
                return sim.calculate_thermal_vulnerability(phi_ground)
            
            elif vuln_type == 'overpressure':
                p_overpressure = sim.calculate_overpressure_ground_new(D_m, imp_energy)
                return sim.calculate_overpressure_vulnerability(p_overpressure)
            
            elif vuln_type == 'wind':
                p_overpressure = sim.calculate_overpressure_ground_new(D_m, imp_energy)
                wind_velocity = sim.calculate_peak_wind_velocity(p_overpressure)
                return sim.calculate_wind_vulnerability(wind_velocity)
            
            elif vuln_type == 'seismic':
                M = sim.calculate_seismic_magnitude(imp_energy)
                M_eff = sim.calculate_effective_seismic_magnitude(M, r_km)
                return sim.calculate_seismic_vulnerability(M_eff)
            
            elif vuln_type == 'ejecta':
                D_tc = sim.calculate_transient_crater_diameter(v_surface, rho_target, sim.entry_angle_deg)
                t_e = sim.calculate_ejecta_thickness(D_tc, r_km)
                return sim.calculate_ejecta_vulnerability(t_e)
                
        # Airburst vulnerability calculations
        elif entry_results["event_type"] == "airburst":
            z_b = entry_results["airburst_altitude"]
            mass = sim.density * (4.0/3.0) * 3.14159 * ((sim.diameter/2)**3)
            KE_initial = 0.5 * mass * (sim.v0**2)
            KE_post = 0.5 * mass * (entry_results["post_breakup_velocity"]**2)
            KE_internal = KE_initial - KE_post
            airburst_energy = max(KE_post, KE_internal)
            
            if vuln_type == 'thermal':
                phi = sim.calculate_airburst_thermal_flux(airburst_energy, z_b, D_m)
                return sim.calculate_thermal_vulnerability(phi)
            
            elif vuln_type == 'overpressure':
                p_overpressure = sim.calculate_overpressure_airburst(D_m, z_b, airburst_energy, entry_results["z_star"])
                return sim.calculate_overpressure_vulnerability(p_overpressure)
            
            elif vuln_type == 'wind':
                p_overpressure = sim.calculate_overpressure_airburst(D_m, z_b, airburst_energy, entry_results["z_star"])
                wind_velocity = sim.calculate_peak_wind_velocity(p_overpressure)
                return sim.calculate_wind_vulnerability(wind_velocity)
            
            elif vuln_type == 'seismic':
                return 0.0  # No seismic effects for airburst events
            
            elif vuln_type == 'ejecta':
                return 0.0  # No ejecta for airburst events
        
        return 0.0  # Default case

    # Binary search implementation for efficient distance finding
    left, right = r_min, r_max
    best_distance = r_min
    
    # Converge on the boundary where vulnerability meets threshold
    while right - left > tol:
        mid = (left + right) / 2
        if calculate_specific_vulnerability_at_distance(mid, vuln_type) >= threshold:
            best_distance = mid
            left = mid  # Vulnerability still meets threshold, search further out
        else:
            right = mid  # Vulnerability below threshold, search closer in
    
    return best_distance

@app.route('/get_ocean_depth', methods=['POST'])
def get_depth():
    """
    Retrieves the ocean depth for given geographic coordinates.
    
    This endpoint returns bathymetric data (ocean depth) for specified latitude and 
    longitude points. This information is important for tsunami modeling and 
    analyzing impacts over water, as depth significantly influences wave 
    propagation and impact characteristics.
    
    Depth data comes from GeoTIFF bathymetric datasets. This helps distinguish 
    between land and ocean impacts and apply appropriate parameters for
    tsunami generation models.
    
    Expected JSON Input:
        {
            "latitude": float,    // Latitude in decimal degrees (-90 to 90).
            "longitude": float    // Longitude in decimal degrees (-180 to 180).
        }
    
    Returns:
        JSON: A response with:
            {
                "depth": float    // Ocean depth in meters (positive values indicate depth).
                                 // A value of 0 usually means land or very shallow water.
                                 // 'null' can be returned if data is unavailable.
            }
    
    Possible Error Responses:
        HTTP 400: If 'latitude' or 'longitude' are missing or invalid.
        HTTP 404: If coordinates are out of dataset bounds or depth data is unavailable.  
        HTTP 500: If an internal server error occurs.
    
    Usage Example:
        POST /get_ocean_depth
        Content-Type: application/json
        
        {
            "latitude": 35.6762,
            "longitude": 139.6503
        }
        
        Example Response:
        {
            "depth": 15.2 
        }
    """
    try:
        # Parse and validate input coordinates from JSON request
        data = request.get_json()
        lat = float(data['latitude'])
        lon = float(data['longitude'])
        
        # Query bathymetric dataset for ocean depth
        depth_value = get_ocean_depth_from_geotiff(lat, lon)
        
        if depth_value is None:
            # Coordinates out of bounds or dataset access error
            return jsonify({"error": "Could not retrieve depth data (out of bounds or map error)."}), 404
        
        # Return depth value (positive = depth in meters, 0 = land/shallow water)
        return jsonify({"depth": depth_value})
        
    except KeyError:
        # Missing required coordinate parameters
        return jsonify({"error": "Missing latitude or longitude."}), 400
    except ValueError:
        # Invalid coordinate format or values
        return jsonify({"error": "Invalid latitude or longitude."}), 400
    except Exception as e:
        # Unexpected errors with logging for debugging
        print(f"Error in /get_ocean_depth: {e}")
        return jsonify({"error": "An internal error occurred."}), 500

if __name__ == '__main__':
    """
    Application entry point for development server.
    
    Starts the Flask development server with configuration suitable for
    development and testing. The server is configured to:
    - Listen on all network interfaces (0.0.0.0)
    - Use port 5000 as the default
    - Disable debug mode for security
    
    Development Features:
    - Automatic server restart on code changes (if debug=True)
    - Detailed error reporting and interactive debugger
    - Real-time code reloading for rapid development
    
    Production Deployment:
    For production environments, use a robust WSGI server like:
    - Gunicorn: gunicorn -w 4 -b 0.0.0.0:5000 app:app
    - uWSGI: uwsgi --http :5000 --wsgi-file app.py --callable app
    - Apache with mod_wsgi
    
    Security Note:
    Debug mode is disabled (debug=False) to prevent exposure of sensitive
    application details and potential security vulnerabilities.
    """
    app.run(host='0.0.0.0', port=5000, debug=False)