"""
ENEO Asteroid Impact Simulation Web Application

This Flask-based web application offers an interactive platform to simulate Near-Earth Object impacts (like asteroids)
and analyze their potential consequences on Earth. The system calculates
blast effects (overpressure and high winds), thermal radiation, seismic activity,
and ejecta distribution. Furthermore, it provides assessments of the
potential impact on human populations and preliminary economic damage.

Technical Architecture:
- Flask web framework for RESTful API endpoints
- Geospatial coordinate system handling (supports antimeridian crossing)
- Accurate damage zone mapping
- Physics simulation modules
- Population and economic impact modeling

Author: Alexandros Notas
Supervisor: Prof. Dimitrios Kaliampakos
Thesis: ENEO - Development of an application for assessing the impacts of a large-scale natural disaster caused by Near-Earth Objects.
Institution: National Technical University of Athens
Date: July 2025
"""

from flask import Flask, request, jsonify, render_template
from src.results import run_simulation_full
from src.population_calculator import calculate_population_in_zones
from shapely.geometry import Point, Polygon
import json
import math
from src.gdp_calculator import calculate_economic_damage
from src.utils import get_ocean_depth_from_geotiff, m_to_km
from src.translation_utils import set_language
from src.visualization_utils import generate_visualization_data
import numpy as np
import random
import requests

    
app = Flask(__name__)

@app.route('/simulate', methods=['POST'])
def simulate():
    """
    Serves as the primary endpoint for conducting asteroid impact simulations.
    
    This endpoint receives the impact parameters, validates them, orchestrates
    the simulation and returns a detailed report to the frontend.
    
    Issues raised:
        HTTP 400: If input parameters are invalid or fail validation checks.
        HTTP 500: If an internal error occurs during the simulation process.
    """
    try:
        # Parse JSON data from the incoming request
        data = request.get_json()
        
        # Extract and validate simulation parameters, applying defaults where necessary
        diameter = float(data['diameter'])
        density = float(data.get('density', 3000)) 
        velocity = float(data['velocity'])
        entry_angle = float(data['entry_angle'])
        r_distance = float(data['distance'])
        lat = float(data.get('latitude', 0))   
        lon = float(data.get('longitude', 0))   
        language = data.get('language', 'en')   
        
        # Log simulation parameters for debugging and analysis tracking
        print(f"\nSimulation Request - Coordinates: {lat:.6f}°, {lon:.6f}°")
        print(f"Parameters: D={diameter}m, ρ={density}kg/m³, v={velocity}km/s, θ={entry_angle}°")
        print(f"Language: {language}")
        
        # Set the language for translations in this request
        set_language(language)
        
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

    # Execute the core simulation
    results_text, results_data = run_simulation_full(diameter, density, velocity, entry_angle, r_distance, lat, lon)
    
    # Calculate population impact
    if lat != 0 or lon != 0:
        # Extract vulnerability zones from simulation results for population analysis
        vulnerability_zones = results_data.get('vulnerability_analysis', {}).get('zones', [])
        
        # Calculate affected population within each vulnerability zone
        population_data = calculate_population_in_zones(lat, lon, vulnerability_zones)
        
        # Add translated country names to the population data
        if 'countries' in population_data:
            from src.translation_utils import get_translation
            for country in population_data['countries']:
                country_name = country.get('name', '')
                # Add translated name as a new field
                country['translated_name'] = get_translation(
                    f'countries.{country_name}', 
                    country_name  # Fallback to original name
                )
        
        results_data['population_analysis'] = population_data

        # Generate country-specific visualization of estimated casualties
        if 'countries' in population_data:
            results_data['country_visualization'] = {
                'countries': population_data['countries']
            }
            
            # Calculate economic impact based on estimated casualties and national GDP data
            economic_data = calculate_economic_damage(population_data['countries'])
            
            # Add translations to economic data as well
            if 'countries' in economic_data:
                for country in economic_data['countries']:
                    country_name = country.get('name', '')
                    country['translated_name'] = get_translation(
                        f'countries.{country_name}', 
                        country_name
                    )
            
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
    Get data for a random Near-Earth Object (NEO) from NASA's Sentry database.
    
    Issues raised:
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
Connection with the main page of the web application.
    """
    return render_template('index.html')

@app.route('/get_ocean_depth', methods=['POST'])
def get_depth():
    """
    Gets the ocean depth at specified geographic coordinates.
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
    """
    app.run(host='0.0.0.0', port=5000, debug=False)
