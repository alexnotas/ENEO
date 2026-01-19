"""
ENEO Asteroid Impact Simulation - Map Utility Functions

This module provides geometric and geographic utility functions for the ENEO application.
It handles the generation of coordinates for damage zones, accounting for Earth's
spherical geometry, antimeridian crossings, and polar regions.

Author: Alexandros Notas
Institution: National Technical University of Athens
Date: July 2025
"""

from math import sin, cos, pi, asin, radians, degrees, atan2

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
