"""
GDP Calculator Module for ENEO Asteroid Impact Simulation System

This module handles economic damage calculations for asteroid impact simulations by:
- Loading and processing global GDP per capita data from World Bank datasets
- Providing country name normalization and matching capabilities
- Calculating economic impact based on casualty estimates and productivity loss

The economic damage model estimates financial losses as:
casualties × GDP_per_capita × productivity_years_factor

Key Features:
- Robust country name matching across different naming conventions
- Multiple fallback strategies for country identification (codes, names, fuzzy matching)
- Efficient caching of processed data for repeated lookups
- Comprehensive logging for data quality monitoring

Data Sources:
- World Bank GDP per capita data (NY.GDP.PCAP.CD indicator)
- Country code mapping files for geographic standardization
"""

import pandas as pd
import os
import numpy as np
from pathlib import Path
import logging

# Configure logging for monitoring data loading and processing steps
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
logger = logging.getLogger(__name__)

# Define file paths for GDP data and country code mappings.
# dynamic path resolution
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
MAPS_DIR = os.path.join(BASE_DIR, 'maps') # Assumes 'maps' folder is next to this script

GDP_FILE = os.path.join(MAPS_DIR, "API_NY.GDP.PCAP.CD_DS2_en_csv_v2_85121.csv")
MAPPING_FILE = os.path.join(MAPS_DIR, "excel.xlsx")

# Global variables to store loaded and processed dataframes and lookup dictionaries.
# These are populated by load_gdp_data() to avoid redundant file operations.
merged_df = None
name_variations = {}
code_to_row = {}

# Country name mappings to reconcile differences between datasets.
# This dictionary helps standardize country names from various sources.
COUNTRY_NAME_MAPPINGS = {
    "slovak republic": "slovakia",
    "united states": "united states of america", # Standardizing to a common long form
    "united kingdom": "united kingdom of great britain and northern ireland", # Standardizing
    "russia": "russian federation",
    "ivory coast": "côte d'ivoire", # Handling special characters and common names
    "czechia": "czech republic",
    "north korea": "korea, democratic people's republic of",
    "south korea": "korea, republic of",
    "syria": "syrian arab republic",
    "republic of the congo": "congo", # Distinguishing between Congos
    "democratic republic of the congo": "congo, the democratic republic of the",
    "eswatini": "swaziland", # Reflecting name changes
    "saint barthélemy": "saint barthelemy",
    "saint lucia": "st. lucia", # Abbreviation handling
    "saint kitts and nevis": "st. kitts and nevis",
    "saint vincent and the grenadines": "st. vincent and the grenadines",
    "taiwan": "taiwan, china", # Addressing geopolitical naming conventions
    "laos": "lao pdr",
    "hong kong": "hong kong sar, china",
    "macau": "macao sar, china",
    "gambia": "gambia, the", # Handling "The" in names
    "bahamas": "bahamas, the",
    "uae": "united arab emirates", # Acronym expansion
    "uk": "united kingdom",
    "usa": "united states",
    "brunei": "brunei darussalam",
    "vietnam": "vietnam", # Ensuring consistency for variations
    "viet nam": "vietnam",
    "micronesia": "micronesia, fed. sts.",
    # Additional mappings for countries often missing or named differently:
    "turkey": "turkiye", # Reflecting recent name change
    "palestine": "west bank and gaza", # Common representation in datasets
    "jersey": "channel islands", # Grouping related territories
    "guernsey": "channel islands",
    "åland": "aland islands", # Handling special characters
    "faeroe is.": "faroe islands", # Abbreviation and common name
    "faeroe islands": "faroe islands", 
    "w. sahara": "western sahara", # Abbreviation
    "gibraltar": "gibraltar"
}

def normalize_country_name(name):
    """
    Standardizes a country name to improve matching across different datasets.

    This function converts the name to lowercase, removes leading/trailing whitespace,
    eliminates common prefixes/suffixes (e.g., "Republic of", "The"), and replaces
    certain patterns (e.g., " and " with " & ") to create a more uniform representation.
    It also consults the `COUNTRY_NAME_MAPPINGS` dictionary for specific overrides.

    Args:
        name (str): The raw country name.

    Returns:
        str: The normalized country name. Returns an empty string if the input is invalid.
    """
    if pd.isna(name) or not name:
        return ""
    
    # Convert to lowercase and remove leading/trailing whitespace.
    norm_name = name.lower().strip()
    
    # Define common words and patterns to remove or standardize.
    replacements = {
        "the ": "",
        ", the": "",
        "republic of ": "",
        "the republic of ": "",
        "democratic republic of ": "",
        "united republic of ": "",
        "people's republic of ": "",
        "federation of ": "",
        "federated states of ": "",
        "kingdom of ": "",
        "state of ": "",
        "united kingdom of ": "",
        "union of ": "",
        "commonwealth of ": "",
        "grand duchy of ": "",
        "principality of ": "",
        " and ": " & ", # Standardize conjunctions
        "-": " ",       # Replace hyphens with spaces
        ".": "",        # Remove periods
        ",": ""         # Remove commas
    }
    
    for old, new in replacements.items():
        norm_name = norm_name.replace(old, new)
    
    # Apply specific predefined mappings for known variations.
    if norm_name in COUNTRY_NAME_MAPPINGS:
        return COUNTRY_NAME_MAPPINGS[norm_name]
    
    return norm_name.strip() # Ensure no leading/trailing spaces after replacements

def load_gdp_data():
    """
    Loads and merges GDP per capita data with country mapping information.

    This function reads GDP data from a CSV file and country codes from an Excel file.
    It then attempts to merge these datasets using various keys (ISO_A3, WB_A2, WB_A3,
    and normalized country names) to achieve the best possible coverage. The latest
    available GDP per capita for each country is selected.

    The processed data is stored in global variables (`merged_df`, `name_variations`,
    `code_to_row`) for efficient subsequent lookups.

    Returns:
        pandas.DataFrame: The merged dataframe containing country information and GDP data,
                          or the raw GDP dataframe if merging fails or mapping file is absent.
                          Returns None if a critical error occurs during loading.
    """
    global merged_df, name_variations, code_to_row
    
    # If data is already loaded, return the existing dataframe.
    if merged_df is not None:
        return merged_df
    
    try:
        logger.info("Loading GDP per capita data...")
        # Load GDP per capita data from the World Bank CSV.
        gdp_df = pd.read_csv(GDP_FILE, skiprows=4) # Skip metadata rows
        
        # Filter for the specific "GDP per capita (current US$)" indicator.
        gdp_df = gdp_df[gdp_df["Indicator Code"] == "NY.GDP.PCAP.CD"].copy()
        
        # Define the range of years to search for the latest GDP data.
        years = [str(year) for year in range(1960, 2024)] # Covers 1960 through 2023
        
        # Function to extract the most recent non-null GDP value and its year.
        def get_latest_gdp(row):
            for year in reversed(years): # Iterate from most recent to oldest
                val = row.get(year)
                if pd.notnull(val):
                    return pd.Series({"GDP_value": val, "GDP_year": year})
            return pd.Series({"GDP_value": pd.NA, "GDP_year": pd.NA}) # Return NA if no data found
        
        # Apply the function to find the latest GDP for each country.
        gdp_df[["GDP_value", "GDP_year"]] = gdp_df.apply(get_latest_gdp, axis=1)
        
        # Retain only essential columns from the GDP data.
        gdp_df = gdp_df[["Country Code", "Country Name", "GDP_value", "GDP_year"]].copy()
        
        # Create a normalized name column for improved matching.
        gdp_df["Normalized_Name"] = gdp_df["Country Name"].apply(normalize_country_name)
        
        # Load the country mapping Excel file if it exists.
        if os.path.exists(MAPPING_FILE):
            try:
                logger.info(f"Loading country code mappings from {MAPPING_FILE}...")
                mapping_df = pd.read_excel(MAPPING_FILE)
                
                # Add a normalized name column to the mapping dataframe for consistent merging.
                if "NAME" in mapping_df.columns: # Check for common column names
                    mapping_df["Normalized_Name"] = mapping_df["NAME"].apply(normalize_country_name)
                elif "CNTRY_NAME" in mapping_df.columns:
                    mapping_df["Normalized_Name"] = mapping_df["CNTRY_NAME"].apply(normalize_country_name)
                
                # Initial merge: Attempt to join mapping data with GDP data using ISO_A3 country codes.
                merged_df = pd.merge(
                    mapping_df, 
                    gdp_df,
                    left_on="ISO_A3", 
                    right_on="Country Code", # World Bank uses "Country Code" for ISO codes
                    how="left" # Keep all countries from the mapping file
                )
                
                logger.info(f"Initial merge (ISO_A3) complete. Columns: {list(merged_df.columns)}")
                    
                # Standardize the 'Normalized_Name' column after merging.
                # Pandas appends _x and _y to duplicate column names during merges.
                if "Normalized_Name_x" not in merged_df.columns and "Normalized_Name" in merged_df.columns:
                    merged_df["Normalized_Name_x"] = merged_df["Normalized_Name"]
                elif "Normalized_Name" not in merged_df.columns and "Normalized_Name_x" in merged_df.columns:
                    merged_df["Normalized_Name"] = merged_df["Normalized_Name_x"] # Use _x if _y is from gdp_df
                
                # Iterative merging for countries still missing GDP data:
                # Strategy: Try WB_A2, then WB_A3, then normalized name matching.

                # Attempt 1: Merge using WB_A2 codes for remaining unmatched entries.
                missing_gdp_mask = merged_df["GDP_value"].isna()
                if missing_gdp_mask.any() and "WB_A2" in merged_df.columns:
                    logger.info("Attempting to match remaining countries using WB_A2 codes...")
                    subset_missing_wb_a2 = merged_df[missing_gdp_mask].copy()
                    wb_a2_matches = pd.merge(
                        subset_missing_wb_a2[["WB_A2", "FID", "ISO_A3", "Normalized_Name_x"]], # Select relevant columns
                        gdp_df, # Merge with the original GDP dataframe
                        left_on="WB_A2",
                        right_on="Country Code",
                        how="inner" # Only keep successful matches
                    )
                    # Update the main merged dataframe with these new matches.
                    for _, match_row in wb_a2_matches.iterrows():
                        # Find the original row in merged_df to update
                        idx_to_update = merged_df[
                            (merged_df["ISO_A3"] == match_row["ISO_A3"]) & 
                            (merged_df["GDP_value"].isna()) # Ensure we only update rows that were missing GDP
                        ].index
                        if not idx_to_update.empty:
                            merged_df.loc[idx_to_update, "GDP_value"] = match_row["GDP_value_y"] # _y from gdp_df
                            merged_df.loc[idx_to_update, "GDP_year"] = match_row["GDP_year_y"]
                            merged_df.loc[idx_to_update, "Country Name"] = match_row["Country Name_y"]
                
                # Attempt 2: Merge using WB_A3 codes.
                missing_gdp_mask = merged_df["GDP_value"].isna()
                if missing_gdp_mask.any() and "WB_A3" in merged_df.columns:
                    logger.info("Attempting to match remaining countries using WB_A3 codes...")
                    subset_missing_wb_a3 = merged_df[missing_gdp_mask].copy()
                    wb_a3_matches = pd.merge(
                        subset_missing_wb_a3[["WB_A3", "FID", "ISO_A3", "Normalized_Name_x"]],
                        gdp_df,
                        left_on="WB_A3",
                        right_on="Country Code",
                        how="inner"
                    )
                    for _, match_row in wb_a3_matches.iterrows():
                        idx_to_update = merged_df[
                            (merged_df["ISO_A3"] == match_row["ISO_A3"]) & 
                            (merged_df["GDP_value"].isna())
                        ].index
                        if not idx_to_update.empty:
                            merged_df.loc[idx_to_update, "GDP_value"] = match_row["GDP_value_y"]
                            merged_df.loc[idx_to_update, "GDP_year"] = match_row["GDP_year_y"]
                            merged_df.loc[idx_to_update, "Country Name"] = match_row["Country Name_y"]
                
                # Attempt 3: Merge using normalized country names as a final fallback.
                missing_gdp_mask = merged_df["GDP_value"].isna()
                if missing_gdp_mask.any() and "Normalized_Name_x" in merged_df.columns and "Normalized_Name" in gdp_df.columns:
                    logger.info("Attempting to match remaining countries using normalized names...")
                    # Create a lookup dictionary from normalized names in gdp_df to their GDP data.
                    name_to_gdp_lookup = {
                        row["Normalized_Name"]: {
                            "GDP_value": row["GDP_value"],
                            "GDP_year": row["GDP_year"],
                            "Country_Name": row["Country Name"]
                        }
                        for _, row in gdp_df.iterrows()
                        if pd.notna(row["Normalized_Name"]) and pd.notna(row["GDP_value"])
                    }
                    
                    # Update rows in merged_df that are still missing GDP data.
                    for idx, row_to_update in merged_df[missing_gdp_mask].iterrows():
                        norm_name_from_mapping = row_to_update.get("Normalized_Name_x")
                        if pd.notna(norm_name_from_mapping) and norm_name_from_mapping in name_to_gdp_lookup:
                            gdp_info = name_to_gdp_lookup[norm_name_from_mapping]
                            merged_df.loc[idx, "GDP_value"] = gdp_info["GDP_value"]
                            merged_df.loc[idx, "GDP_year"] = gdp_info["GDP_year"]
                            merged_df.loc[idx, "Country Name"] = gdp_info["Country_Name"]
                
                # Populate lookup dictionaries for faster access in `lookup_gdp`.
                # These dictionaries map various identifiers (FID, ISO codes, names) to country data rows.
                for _, row_data in merged_df.iterrows():
                    if pd.notna(row_data.get("FID")): # Ensure FID exists and is not NaN
                        fid_key = str(int(row_data["FID"])) # Convert FID to string key
                        if pd.notna(row_data.get("GDP_value")): # Store only if GDP data is present
                            code_to_row[fid_key] = row_data
                    
                    # Map various country codes to the row data.
                    for code_key_column in ["ISO_A3", "WB_A2", "WB_A3"]:
                        if pd.notna(row_data.get(code_key_column)):
                            code_to_row[row_data[code_key_column]] = row_data
                    
                    # Map normalized names (from mapping file and GDP file) to the row data.
                    # This allows lookup by different name variations.
                    if pd.notna(row_data.get("NAME")): # Original name from mapping file
                        normalized_map_name = normalize_country_name(row_data["NAME"])
                        name_variations[normalized_map_name] = row_data
                    
                    if pd.notna(row_data.get("Country Name")): # Name from World Bank GDP data
                        normalized_gdp_name = normalize_country_name(row_data["Country Name"])
                        name_variations[normalized_gdp_name] = row_data
                
                logger.info(f"Successfully merged GDP data. Found GDP values for {merged_df['GDP_value'].notna().sum()} out of {len(merged_df)} countries.")
                return merged_df
            
            except Exception as e:
                logger.error(f"Error during data merging process: {e}")
                # Fallback: If merging fails, use only the loaded GDP data.
                merged_df = gdp_df # Ensure merged_df is at least the gdp_df
                return merged_df
        else:
            logger.warning(f"Country mapping file {MAPPING_FILE} not found. Using only GDP data.")
            merged_df = gdp_df # Use only GDP data if mapping file is unavailable.
            return merged_df
            
    except Exception as e:
        logger.error(f"Critical error loading GDP data: {e}")
        return None # Return None if initial GDP file loading fails.

def lookup_gdp(country_code=None, country_name=None, fid=None):
    """
    Looks up GDP per capita and the corresponding year for a country.

    The lookup prioritizes FID, then country codes, then normalized country names.
    It utilizes pre-loaded and processed data for efficiency.

    Args:
        country_code (str, optional): The ISO_A3, WB_A2, or WB_A3 country code.
        country_name (str, optional): The name of the country.
        fid (int or str, optional): The Feature ID (FID) of the country.

    Returns:
        tuple: A tuple containing (gdp_per_capita, gdp_year).
               Returns (None, None) if the country or its GDP data cannot be found.
    """
    # Ensure data is loaded before attempting lookup.
    if merged_df is None:
        load_gdp_data()
        if merged_df is None: # If loading failed critically
            logger.error("GDP data could not be loaded for lookup.")
            return None, None

    # Attempt lookup by FID first, as it's often a reliable unique identifier.
    if fid:
        fid_str_key = str(fid)
        if fid_str_key in code_to_row:
            row = code_to_row[fid_str_key]
            if pd.notna(row.get("GDP_value")):
                return row["GDP_value"], row.get("GDP_year") # Use .get for GDP_year for safety
    
    # Attempt lookup by country code (ISO_A3, WB_A2, WB_A3).
    if country_code and country_code in code_to_row:
        row = code_to_row[country_code]
        if pd.notna(row.get("GDP_value")):
            return row["GDP_value"], row.get("GDP_year")
    
    # Attempt lookup by normalized country name.
    if country_name:
        normalized_name_key = normalize_country_name(country_name)
        if normalized_name_key in name_variations:
            row = name_variations[normalized_name_key]
            if pd.notna(row.get("GDP_value")):
                return row["GDP_value"], row.get("GDP_year")
        
        # Fallback: Direct search in the merged_df using normalized names.
        # This handles cases where name_variations might not be perfectly populated.
        if merged_df is not None:
            # Determine which normalized name column to use from merged_df.
            # It could be 'Normalized_Name_x' (from mapping) or 'Normalized_Name' (from gdp_df).
            norm_col_to_check = None
            if "Normalized_Name_x" in merged_df.columns:
                norm_col_to_check = "Normalized_Name_x"
            elif "Normalized_Name" in merged_df.columns: # Check the one from gdp_df if _x isn't there
                norm_col_to_check = "Normalized_Name"
            
            if norm_col_to_check:
                matches = merged_df[merged_df[norm_col_to_check] == normalized_name_key]
                if not matches.empty:
                    first_match = matches.iloc[0]
                    if pd.notna(first_match.get("GDP_value")):
                        return first_match["GDP_value"], first_match.get("GDP_year")
            
            # Final attempt: Match against the 'Country Name' column from the World Bank data.
            if "Country Name" in merged_df.columns:
                # Normalize the 'Country Name' column on-the-fly for comparison.
                matches_on_gdp_name = merged_df[merged_df["Country Name"].apply(normalize_country_name) == normalized_name_key]
                if not matches_on_gdp_name.empty:
                    first_match_gdp_name = matches_on_gdp_name.iloc[0]
                    if pd.notna(first_match_gdp_name.get("GDP_value")):
                        return first_match_gdp_name["GDP_value"], first_match_gdp_name.get("GDP_year")
                        
    # As a last resort, attempt fuzzy matching (partial string containment).
    # This is less precise and should be used cautiously.
    if country_name and merged_df is not None:
        normalized_name_to_search = normalize_country_name(country_name)
        if len(normalized_name_to_search) > 3: # Avoid overly broad matches for short strings
            
            # Define columns to check for fuzzy matching.
            fuzzy_check_columns = []
            if "Normalized_Name_x" in merged_df.columns: fuzzy_check_columns.append("Normalized_Name_x")
            if "Normalized_Name" in merged_df.columns: fuzzy_check_columns.append("Normalized_Name") # From gdp_df
            if "Country Name" in merged_df.columns: fuzzy_check_columns.append("Country Name") # Raw name from gdp_df
            
            for _, db_row in merged_df.iterrows():
                for col_name in fuzzy_check_columns:
                    db_entry_name = db_row.get(col_name)
                    if pd.notna(db_entry_name) and pd.notna(db_row.get("GDP_value")):
                        # Normalize the database entry name if it's the raw 'Country Name'.
                        current_db_norm_name = normalize_country_name(str(db_entry_name)) if col_name == "Country Name" else str(db_entry_name)
                        
                        # Check for containment in either direction.
                        if (normalized_name_to_search in current_db_norm_name) or \
                           (current_db_norm_name in normalized_name_to_search):
                            logger.info(f"Fuzzy match found for '{country_name}' with '{db_entry_name}'.")
                            return db_row["GDP_value"], db_row.get("GDP_year")
                        
    # If no match is found after all attempts.
    logger.warning(f"GDP data not found for query: code='{country_code}', name='{country_name}', fid='{fid}'.")
    return None, None

def calculate_economic_damage(countries_data, productivity_years=1):
    """
    Estimates economic damage based on casualties and GDP per capita for affected countries.

    The economic damage for each country is calculated as:
    `Number of Casualties * GDP per Capita * Productivity Years Factor`.
    The `productivity_years` factor can be used to model longer-term economic impact.

    Args:
        countries_data (list): A list of dictionaries, where each dictionary represents
                               an affected country and includes at least 'name', 'fid',
                               and 'total_casualties'.
        productivity_years (int or float, optional): A multiplier representing the number
                                                     of years of lost productivity per casualty.
                                                     Defaults to 1.

    Returns:
        dict: A dictionary containing:
            - "countries" (list): Detailed economic damage assessment for each country.
            - "total_economic_damage" (float): Sum of economic damage across all countries.
            - "countries_with_data" (int): Number of countries for which GDP data was found.
            - "countries_without_data" (int): Number of countries lacking GDP data.
            - "note" (str): A brief note about the calculation method.
    """
    # Ensure GDP data is loaded.
    load_gdp_data()
    if merged_df is None:
        logger.error("Cannot calculate economic damage: GDP data failed to load.")
        return { # Return a structure indicating failure
            "countries": [], "total_economic_damage": 0, "countries_with_data": 0,
            "countries_without_data": len(countries_data),
            "note": "Error: GDP data could not be loaded."
        }
    
    damage_results = []
    grand_total_economic_damage = 0.0
    count_countries_with_gdp = 0
    list_countries_without_gdp = []
    
    country_names_for_logging = [c.get('name', 'N/A') for c in countries_data if c.get('name')]
    logger.info(f"Calculating economic damage for {len(countries_data)} countries: {', '.join(country_names_for_logging) if country_names_for_logging else 'No country names provided'}.")
    
    for country_info in countries_data:
        country_name = country_info.get("name", "Unknown Country")
        fid_val = country_info.get("fid")
        casualties_val = country_info.get("total_casualties", 0)
        
        current_country_result = {
            "name": country_name,
            "fid": fid_val,
            "casualties": casualties_val,
            "population": country_info.get("total_population", 0) # Include population for context
        }
        
        # Skip calculation if there are no casualties for this country.
        if casualties_val <= 0:
            current_country_result.update({"gdp_per_capita": None, "gdp_year": None, "economic_damage": 0})
            damage_results.append(current_country_result)
            continue # Move to the next country
        
        # Attempt to retrieve GDP data for the current country.
        gdp_val, gdp_yr = lookup_gdp(
            country_name=country_name, 
            fid=fid_val
            # country_code could also be passed if available in country_info
        )
        
        if gdp_val is not None:
            gdp_per_capita = float(gdp_val)
            current_country_result["gdp_per_capita"] = gdp_per_capita
            current_country_result["gdp_year"] = int(gdp_yr) if pd.notna(gdp_yr) else None
            
            # Calculate economic damage for this country.
            country_economic_damage = casualties_val * gdp_per_capita * productivity_years
            current_country_result["economic_damage"] = country_economic_damage
            
            grand_total_economic_damage += country_economic_damage
            count_countries_with_gdp += 1
            logger.info(f"GDP found for {country_name} (FID: {fid_val}): ${gdp_per_capita:,.2f} (Year: {gdp_yr}). Estimated damage: ${country_economic_damage:,.2f}")
        else:
            # GDP data not found for this country.
            current_country_result.update({"gdp_per_capita": None, "gdp_year": None, "economic_damage": None})
            list_countries_without_gdp.append(country_name)
            logger.warning(f"No GDP data found for {country_name} (FID: {fid_val}). Cannot calculate economic damage.")
        
        damage_results.append(current_country_result)
    
    # Log countries for which GDP data was not found.
    if list_countries_without_gdp:
        countries_to_report = [country if country is not None else "Unknown Country" for country in list_countries_without_gdp]
        if countries_to_report:
            logger.warning(f"Summary: GDP data was not found for the following countries: {', '.join(countries_to_report)}.")
    
    # Sort results by economic damage in descending order for reporting.
    damage_results.sort(key=lambda x: x.get("economic_damage", 0) or 0, reverse=True)
    
    return {
        "countries": damage_results,
        "total_economic_damage": grand_total_economic_damage,
        "countries_with_data": count_countries_with_gdp,
        "countries_without_data": len(list_countries_without_gdp),
        "note": f"Economic damage estimated as: Casualties * GDP per Capita * {productivity_years} (Productivity Years Factor)."
    }