"""
ENEO Asteroid Impact Simulation - Visualization Utilities

This module handles the generation of visualization data for the frontend map.

Author: Alexandros Notas
Institution: National Technical University of Athens
Date: July 2025
"""

import numpy as np
from src.map_utils import create_circle_coordinates
from src.utils import m_to_km
from src.models import AsteroidImpactSimulation
from src.results import find_specific_vulnerability_distance

def parse_danger_zones(section_text, lat, lon, section_title=None):
    """
    Extracts and categorizes information about various damage zones from the simulation's textual output.
    
    This function processes structured text to identify different types of hazard zones (e.g., thermal, seismic, airblast)
    and their respective distance ranges. This information is vital for visualizing the spatial extent of these hazards.
    It uses keyword matching and pattern recognition to automatically classify zones based on their descriptive text.
    
    """
    zones = []
    lines = section_text.split('\n')
    section_hint = section_title.lower() if section_title else ""
    lower_section_text = section_text.lower()

    section_type_hints = {
        'airblast': ['blast', 'airblast', 'overpressure', 'ζώνες κινδύνου έκρηξης', 'υπερπίεση', 'ωστικό κύμα'],
        'thermal': ['thermal', 'θερμ', 'πυρ', 'φωτιά'],
        'seismic': ['seismic', 'σεισμ', 'richter', 'ρίχτερ'],
        'wind': ['wind', 'ανέμ', 'ef', 'ανεμο'],
        'ejecta': ['ejecta', 'εκτιναγ', 'θραυσ'],
        'tsunami': ['tsunami', 'τσουνάμι']
    }

    thermal_keywords = ['burns', 'ignition', 'thermal', 'εγκαύματα', 'ανάφλεξη', 'θερμ', 'πυρκαγ']
    seismic_keywords = ['σεισμικ', 'σεισμ', 'σεισμική']
    airblast_keywords = [
        "collapse", "distorted", "damage", "shatter", "structural",
        "windows", "buildings", "office-type", "wall-bearing", "wood frame",
        "κατάρρευση", "καταρρεύ", "παραμόρφωση", "παραμορφ", "ζημιά", "διαθραύση", "δομικός",
        "παράθυρα", "κτίρια", "γραφείων", "φέρον τοίχο", "ξύλινο πλαίσιο",
        "γέφυρ", "οχήμ", "εκτοπισ", "ανακατασκευ"
    ]
    tsunami_threshold_labels = [">1km", ">100m", ">10m", ">1m"]

    for line_content in lines:
        line = line_content.strip()
        
        # Skip empty lines and lines without proper format
        if not line or ':' not in line:
            continue

        # Split line into description and distance range components
        description_part = line.split(':', 1)[0].strip()
        range_part_with_km = line.split(':', 1)[-1].strip()
        lower_description = description_part.lower()

        determined_type_for_line = None

        # Classify zone type using pattern matching and keywords
        
        # Wind zones: Enhanced Fujita scale identification
        if description_part.startswith('EF'):
            determined_type_for_line = 'wind'
            
        # Tsunami zones: Amplitude thresholds in tsunami context
        elif description_part in tsunami_threshold_labels and "tsunami" in lower_section_text:
            if "km" in range_part_with_km and "None" not in range_part_with_km:
                determined_type_for_line = 'tsunami'
                
        # Seismic zones: Richter scale references (English and Greek)
        elif 'richter' in lower_description or 'ρίχτερ' in lower_description or \
             any(keyword in lower_description for keyword in seismic_keywords):
            determined_type_for_line = 'seismic'
        
        # Thermal zones: Heat and fire-related keywords (English and Greek)
        elif any(keyword in lower_description for keyword in thermal_keywords):
            determined_type_for_line = 'thermal'
        
        # Ejecta zones: Thickness measurements (e.g., ">10 m")
        elif description_part.startswith(">") and description_part.endswith("m"):
            determined_type_for_line = 'ejecta'
        
        # Airblast zones: Structural damage descriptions (English and Greek)
        elif any(keyword in lower_description for keyword in airblast_keywords):
            determined_type_for_line = 'airblast'
        
        # Additional tsunami detection for wave amplitude descriptions (English and Greek)
        if not determined_type_for_line and (
            "wave amplitude >" in lower_description
            or ("tsunami" in lower_section_text and "zone" in lower_description)
            or ("τσουνάμι" in lower_section_text and ("ζώνη" in lower_description or "zone" in lower_description))
        ):
            if "km" in line.split(':', 1)[-1]:
                determined_type_for_line = 'tsunami'

        # Fall back to section-based inference if keyword detection fails
        if not determined_type_for_line and section_hint:
            for zone_type, hints in section_type_hints.items():
                if any(hint in section_hint for hint in hints):
                    determined_type_for_line = zone_type
                    break

        # Process zones with determined types
        if determined_type_for_line:
            try:
                dist_text = line.split(':', 1)[1].strip()

                # Robust parsing: Skip line if it doesn't clearly end in "km" or contains unparseable values
                if not dist_text.endswith("km") and "km" not in dist_text:
                    # Special check: maybe it's just a value line that isn't a distance (e.g. pressure Pa, m/s etc.)
                    # If this line was categorized by section hint but isn't actually a distance line, skip it silently.
                    continue

                # Skip zones without valid distance data
                if 'None' not in dist_text and dist_text:
                    start_dist_val, end_dist_val = 0.0, 0.0
                    
                    # Parse distance ranges
                    # Remove 'km' first to just get numbers/hyphens
                    # Also handle edge cases where dist_text might have extra trailing chars
                    cleaned_dist = dist_text.replace('km', '').strip()
                    if not cleaned_dist:
                        continue 

                    if '-' in cleaned_dist:  # Range format "X-Y"
                        parts = cleaned_dist.split('-', 1)
                        # Check bits to ensure they are numeric
                        if not parts[0].strip():
                            continue
                        start_dist_val = float(parts[0])
                        # If second part is empty or non-numeric, use same as start? Or fail?
                        # Generally if it parses, assume good.
                        if len(parts) > 1 and parts[1].strip():
                            end_dist_val = float(parts[1])
                        else:
                            end_dist_val = start_dist_val
                    else:  # Single value format "X" (implies 0 to X)
                        if not cleaned_dist.strip(): 
                             continue
                        try:
                            # If cleaned_dist is e.g. "EF5", float() will fail, which is good -> except block.
                            end_dist_val = float(cleaned_dist)
                        except ValueError:
                             # Not a number (e.g. some text residue), skip
                             continue
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
    Prepares a dataset for visualizing impact effects on a map.
    
    This function processes the simulation's outputs to create geographically referenced 
    data for identified damage zones, vulnerability assessments, and other impact 
    characteristics. It combines parsed text results (for basic damage zones)
    with more detailed vulnerability analyses derived from calculated thresholds and distances.
    
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
    split_parts = results_text.split('===')  # Split by section delimiters

    # Process each section title-content pair
    for i in range(1, len(split_parts), 2):
        section_title_full = split_parts[i].strip()
        section_content = split_parts[i+1].strip() if (i+1) < len(split_parts) else ""
        
        # Extract zones from section content
        parsed_zones_from_content = parse_danger_zones(section_content, lat, lon, section_title_full)
        
        # Categorize parsed zones into visualization structure
        for zone in parsed_zones_from_content:
            zone_type_from_parser = zone.get('type')
            if zone_type_from_parser and zone_type_from_parser in visualization:
                visualization[zone_type_from_parser].append(zone)
            elif zone_type_from_parser:
                print(f"Debug: Parsed zone type '{zone_type_from_parser}' not in standard categories")

    # Generate advanced vulnerability analysis zones from simulation data
    if 'vulnerability_analysis' in results_data and 'zones' in results_data['vulnerability_analysis']:
        
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
