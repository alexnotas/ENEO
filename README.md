<div style="background:#050611; color:#EAF0FF; border:1px solid rgba(255,255,255,0.12); border-radius:18px; padding:28px 22px; margin-bottom:22px; max-width:980px">

<div align="center">

<h1 style="margin:0 0 10px; font-size:44px; letter-spacing:0.08em; color:#FF7A18; text-shadow:0 0 18px rgba(255,122,24,0.25)">ENEO</h1>

<div style="margin:0 0 14px; font-size:16px; letter-spacing:0.22em; opacity:0.9">NEAR-EARTH OBJECT IMPACT SIMULATOR</div>

<div style="height:1px; width:72%; background:linear-gradient(90deg, rgba(255,122,24,0), rgba(255,122,24,0.55), rgba(255,122,24,0)); margin:14px auto 18px"></div>

![ENEO Logo](https://img.shields.io/badge/ENEO-Asteroid%20Impact%20Simulator-blue?style=for-the-badge)
![Python](https://img.shields.io/badge/Python-3.8+-green?style=for-the-badge&logo=python)
![Flask](https://img.shields.io/badge/Flask-2.0+-lightgrey?style=for-the-badge&logo=flask)
![License](https://img.shields.io/badge/License-MIT-yellow?style=for-the-badge)

<p style="margin:18px 0 10px; max-width:860px; line-height:1.7; font-size:15.5px; opacity:0.95">
  <strong>An advanced Near-Earth Object impact simulation platform developed in the School of Mining and Metallurgical Engineering of the National Technical University of Athens</strong>
</p>

</div>

---

## � Overview

**ENEO** is a web-based application for simulating Near-Earth Object (NEO) impact events and analyzing their potential consequences on Earth. Developed as part of a thesis project at the **National Technical University of Athens (NTUA)**, this platform provides comprehensive physics-based modeling of NEO impacts, including atmospheric entry, crater formation, blast effects, seismic impacts, thermal radiation, ejecta and population/economic impact assessments. **ENEO** provides a unified, modular architecture and is the only open-source tool that simulates NEO impact scenarios from atmospheric entry to preliminary socio-economic impact.

### 🎯 Purpose

This simulator allows researchers, students, and space enthusiasts to:
- Model realistic Near-Earth Object (like asteroids) impact scenarios with scientific accuracy
- Visualize damage zones on an interactive global map
- Assess potential population casualties and preliminary economic damage
- Study various impact phenomena
- Integrate real NEO data from NASA's Sentry API
- Use the source code to empower similar projects without needing to code the whole project from the start
- The opportunity to experiment and update physics and vulnerability equations

  
---

## ✨ Features

### Core Capabilities

- **🌍 Interactive Impact Simulation**
  - Physics-based calculations
  - Customizable asteroid parameters (diameter, density, velocity, entry angle)
  - Geographic coordinate input for any location on Earth

- **📊 Comprehensive Impact Analysis**
  - Atmospheric fragmentation and pancake effects
  - Crater formation modeling (transient and final crater dimensions)
  - Airblast overpressure calculations
  - Hurricane-like winds calculations
  - Thermal radiation effects and burn zones
  - Seismic wave propagation and earthquake magnitude
  - Ejecta blanket distribution
  - Preliminary tsunami generation for ocean impacts

- **👥 Population & Economic Impact**
  - Global population effects analysis using gridded population data
  - Vulnerability models for different hazard types
  - Country-by-country casualty estimates
  - Economic damage calculations based on GDP per capita

- **🗺️ Advanced Visualization**
  - Interactive Leaflet-based mapping
  - Color-coded damage zone overlays
  - Antimeridian crossing support for accurate global projections

- **🌐 Multi-language Support**
  - English and Greek translations
  - Dynamic language switching
  - Localized result presentations

- **🛰️ NASA API Integration**
  - Get real asteroid data from NASA Sentry system
  - Pre-populate simulations with known NEO characteristics
  - Stay updated with current asteroid threat assessments

---

## 🏗️ Architecture

### Technology Stack

- **Backend**: Python 3.8+ with Flask web framework
- **Frontend**: HTML5, CSS3, JavaScript (ES6+)
- **Mapping**: Leaflet.js for interactive geographic visualization
- **Data Processing**: NumPy, Pandas, GeoPandas
- **Geospatial**: Shapely, GeoPandas, Rasterio, PyProj
- **Styling**: Bootstrap 5, Material Design Icons

### Project Structure

```
ENEO/
├── app.py                      # Main Flask application & API endpoints
├── app.wsgi                    # WSGI deployment configuration
├── src/                        # Source code modules
│   ├── models.py                   # Core asteroid impact physics models
│   ├── results.py                  # Simulation orchestration and result formatting
│   ├── vulnerability_models.py     # Population vulnerability calculations
│   ├── population_calculator.py    # Population impact assessment
│   ├── gdp_calculator.py           # Economic damage calculations
│   ├── map_utils.py                # Geographic and mapping utilities
│   ├── visualization_utils.py      # Data visualization helpers
│   ├── translation_utils.py        # Multi-language support
│   ├── thresholds.py               # Damage threshold definitions
│   └── utils.py                    # Physical constants and utilities
├── maps/                       # Geographic and demographic data
│   ├── world.shp               # World boundaries shapefile
│   ├── API_SP.POP.TOTL_*.csv   # World Bank population data
│   ├── gdp_data.csv            # GDP per capita data (renamed from API_NY...)
│   ├── country_codes.xlsx      # Country code mapping (renamed from excel.xlsx)
│   └── country_fid_lookup.csv  # Country ID mapping
├── static/                     # Frontend assets
│   ├── css/                    # Stylesheets
│   ├── js/                     # JavaScript modules
│   └── translations/           # Language files (en.json, el.json)
└── templates/                  # HTML templates
    └── index.html              # Main application interface
```

---

## 🚀 Installation

> **⚠️ IMPORTANT: REQUIRED MAP DATA**
>
> The high-resolution map files required for this application are too large for GitHub. You must download them separately.
>
> **Method 1: Automatic Download (Recommended)**
> Run the included script to automatically download and setup the maps:
> ```bash
> python download_data.py
> ```
>
> **Method 2: Manual Download**
> 1. Download the map data from Zenodo: **[https://zenodo.org/records/18302255](https://zenodo.org/records/18302255)**
> 2. Create a folder named `maps` in the root directory of the project:
>    ```bash
>    mkdir maps
>    ```
> 3. Extract/Place all downloaded files (shapefiles, CSVs, etc.) into the `maps/` folder.
>
> *Alternatively, you can download the full repository including all map data from: [https://zenodo.org/records/18326608](https://zenodo.org/records/18326608)*

### Prerequisites

- Python 3.8 or higher
- pip (Python package manager)
- Git (for cloning the repository)

### Step 1: Clone the Repository

```bash
git clone https://github.com/alexnotas/ENEO.git
cd ENEO-main
```

### Step 2: Create a Virtual Environment (Recommended)

```bash
# On Linux/macOS
python3 -m venv venv
source venv/bin/activate

# On Windows
python -m venv venv
venv\Scripts\activate
```

### Step 3: Install Dependencies

```bash
pip install -r requirements.txt
```

**Or install manually (same as `requirements.txt`):**
```bash
pip install flask numpy pandas geopandas shapely rasterio pyproj requests
```

**Required Python Packages:**
- `Flask` - Web framework
- `requests` - HTTP requests (NASA API integration)
- `numpy` - Numerical computations
- `pandas` - Data manipulation
- `shapely` - Geometry operations
- `geopandas` - Geospatial vector data operations (reading `world.shp`, etc.)
- `rasterio` - Raster data reading (population grids, ocean depth GeoTIFF)
- `pyproj` - Coordinate reference system transformations

> Note: you don't import `fiona` directly in this repo. On many systems it is installed automatically as a dependency of `geopandas` (shapefile I/O).

### Step 4: Verify Data Files

Ensure the following data files are present in the `maps/` directory:
- `world.shp` (and associated .dbf, .shx, .prj files) - World country boundaries
- `population_with_country_fid_assigned.tif` - High-resolution gridded population raster
- `ETOPO_2022_v1_60s_N90W180_surface.tif` - Global ocean depth data (ETOPO 2022)
- `gdp_data.csv` - GDP per capita data by country
- `country_codes.xlsx` - Country code mappings
- `country_fid_lookup.csv` - Country FID lookup table

### System Requirements

- **Disk Space**: Approximately 1.5 GB for map data files
- **RAM**: Minimum 4 GB recommended (8 GB for large-scale simulations)
- **Python**: Version 3.8 or higher

---

## 💻 Usage

### Running the Application Locally

1. **Start the Flask development server:**

```bash
python app.py
```

2. **Open your web browser and navigate to:**

```
http://localhost:5000
```
or for all the icons to work:
```
http://127.0.0.1:5000
```
3. **The application should now be running!**

## ✅ Validation & Testing

ENEO separates **unit tests** (fast, automated checks) from **validation scripts** (manual, print-based reproduction of published benchmarks).

### Unit Tests (automated)
Run the unit tests in `tests/unit/` to validate core physics utilities and model invariants:

- `tests/unit/test_utils.py`
- `tests/unit/test_models.py`

### Validation Scripts (manual reproduction)
To reproduce the published/benchmark validation scenarios, run the scripts in `tests/validation/`:

- `tests/validation/validation_physics.py` (Reproduce the physics test results)
- `tests/validation/validation_vulnerability.py` (Reproduce the vulnerability test results)

These validation scripts print values for manual comparison with established benchmarks and are intentionally **not** strict unit tests.

### Using the Simulator

#### Basic Simulation

1. **Enter Asteroid Parameters:**
   - **Diameter**: Size of the asteroid (meters) - e.g., 100m
   - **Density**: Material density (kg/m³) - e.g., 3000 for rocky asteroids
   - **Velocity**: Entry speed (km/s) - range: 11-72 km/s
   - **Entry Angle**: Angle from horizontal (degrees) - range: 15-90°
   - **Distance**: Reference distance for analysis (km)

2. **Set Impact Location:**
   - Enter latitude and longitude coordinates if you run only the Python script
   - Click on the interactive map to select a location if you are running the UI

3. **Run Simulation:**
   - Click the "Run Simulation" button
   - Wait for calculations to complete
  
     
4. **View Results:**
   - **Summary Tab**: Overview of impact effects
   - **Detailed Results**: Comprehensive breakdown of all hazards
   - **Population Impact**: Casualties by country and damage zone
   - **Economic Impact**: Estimated financial losses
   - **Map Visualization**: Geographic representation of damage zones

#### Advanced Features

**Load NASA NEO Data:**
- Click the "NASA Sentry API" button
- Browse known Near-Earth Objects
- Select an asteroid to auto-populate parameters

**Switch Languages:**
- Use the language switcher to toggle between English and Greek

---

## 📚 Documentation

### Physics Models

ENEO implements peer-reviewed scientific models for asteroid impact simulation:

#### 1. **Atmospheric Entry**
- Drag equation modeling with altitude-dependent density
- Pancaking and fragmentation effects
- Airburst detection and altitude calculations

#### 2. **Crater Formation**
- Transient crater dimensions using scaling laws
- Rim collapse and final crater calculations
- Different crater types (simple vs. complex)

#### 3. **Thermal Radiation**
- Fireball luminosity and temperature
- Thermal pulse duration and intensity
- Burn severity zones (1st, 2nd, 3rd degree burns)

#### 4. **Seismic Effects**
- Richter scale magnitude calculation
- Ground motion amplitude
- Structural damage thresholds

#### 5. **Airblast**
- Overpressure decay with distance
- Dynamic pressure calculations
- Wind speed estimates

#### 6. **Ejecta**
- Ballistic fragment distribution
- Ejecta thickness profiles
- Fragment size and velocity

#### 7. **Tsunami(preliminary/under development)** (Ocean Impacts)
- Wave amplitude generation
- Potential for inundation modeling

### Vulnerability Models

Population vulnerability calculations based on research by **C. Rumpf et al. (2017)**:
- Crater proximity fatality rates
- Thermal radiation burn thresholds
- Overpressure injury/fatality curves
- Seismic structural collapse probabilities
- Ejecta fragment impact risk

---

## 📄 License

This project is part of academic research at the **National Technical University of Athens (NTUA) and has a MIT License**. 

---

## 👨‍🎓 Author

**Alexandros Notas**  
School of Mining and Metallurgical Engineering  
National Technical University of Athens (NTUA)

**Thesis:** *ENEO - Development of an application for assessing the impacts of a large-scale natural disaster caused by Near-Earth Objects.*  
**Date:** July 2025

---

## 🙏 Acknowledgments

- **National Technical University of Athens (NTUA)** - School of Mining and Metallurgical Engineering
- **NASA** - For the Sentry API and NEO data
- **Research Community** - For published impact physics models and vulnerability studies

---

## 📞 Contact & Support

- **E-Mail**: alexnotas@metal.ntua.gr
- **Issues**: Please use the GitHub Issues tab for bug reports and feature requests
- **Questions**: For academic inquiries, contact through NTUA channels

---


## 📊 Project Status

**Current Version:** 1.0.0 
**Status:** Active Development  
**Last Updated:** January 2026

---

<div align="center">

</div>

</div>
