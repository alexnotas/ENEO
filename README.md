# ENEO - Evaluating Near Earth Objects Impacts

ENEO is a web-based simulator to evaluate and visualize the impact effects of Near Earth Objects (NEOs). It provides a comprehensive analysis of various impact scenarios, calculating effects like crater size, airblast, thermal radiation, seismic activity, and potential tsunamis. The application also assesses the impact on population and economy based on the selected impact location.

## Features

-   **Interactive Map:** Select an impact location anywhere on Earth using a Leaflet-powered map.
-   **Comprehensive Simulation:** Input asteroid parameters (diameter, density, velocity, entry angle) to run a detailed physical simulation.
-   **Multi-faceted Analysis:** Calculates and displays a wide range of impact effects:
    -   Atmospheric Entry
    -   Crater Dimensions
    -   Airblast (Overpressure)
    -   Thermal Radiation
    -   Seismic Shaking
    -   Ejecta Blanket
    -   Wind Effects
    -   Tsunami potential for oceanic impacts
-   **Vulnerability & Damage Assessment:** Analyzes and visualizes vulnerability zones and damage categories based on established thresholds.
-   **Socio-Economic Impact:** Estimates casualties and economic damage by correlating impact zones with global population and GDP data.
-   **Dynamic Visualization:** Displays hazard zones and impact data dynamically on the map and in detailed result panels.
-   **Responsive UI:** Built with Bootstrap for a seamless experience on both desktop and mobile devices.

## Tech Stack

-   **Backend:** Python, Flask
-   **Frontend:** HTML, CSS, JavaScript
-   **Libraries:**
    -   Leaflet.js for the interactive map
    -   Bootstrap 5 for UI components
    -   NumPy for numerical calculations

## Project Structure

```
.
├── app.py                  # Main Flask application, handles routing and simulation logic
├── models.py               # Core physics simulation models (AsteroidImpactSimulation class)
├── results.py              # Processes and aggregates simulation results
├── vulnerability_models.py # Functions for calculating vulnerability to different effects
├── thresholds.py           # Defines damage thresholds for various impact effects
├── population_calculator.py# Logic for population impact analysis
├── gdp_calculator.py       # Logic for economic impact analysis
├── utils.py                # Utility functions
├── templates/
│   └── index.html          # Main HTML template for the user interface
├── static/
│   ├── css/                # CSS stylesheets
│   └── js/
│       └── script.js       # Frontend JavaScript for interactivity, API calls, and results display
├── maps/                   # Data files (e.g., population, GDP, elevation)
└── app.wsgi                # WSGI entry point for production deployment
```

## Setup and Installation

1.  **Clone the repository:**
    ```sh
    git clone <your-repository-url>
    cd ENEO-main
    ```

2.  **Create and activate a virtual environment:**
    ```sh
    python -m venv venv
    # On Windows
    venv\Scripts\activate
    # On macOS/Linux
    source venv/bin/activate
    ```

3.  **Install dependencies:**
    Create a `requirements.txt` file with the necessary Python packages and install them.
    ````
    # requirements.txt
    Flask
    numpy
    gdal  # Or rasterio, depending on get_ocean_depth_from_geotiff implementation
    # ... other dependencies
    ````
    Then run:
    ```sh
    pip install -r requirements.txt
    ```

4.  **Run the application:**
    ```sh
    flask run
    ```
    The application will be available at `http://127.0.0.1:5000`.

## Usage

1.  Open the web application in your browser.
2.  Use the map to click on a desired impact location. The latitude and longitude will be automatically filled.
3.  Enter the asteroid's parameters in the input form on the left:
    -   Diameter (meters)
    -   Density (kg/m³)
    -   Velocity (km/s)
    -   Entry Angle (degrees)
4.  Click the "Run Simulation" button.
5.  The simulation will run, and the results will be displayed in the panels on the right. Hazard zones will be drawn on the map.
6.  Use the "Zone Display Options" to toggle different hazard visualizations on the map.
7.  Explore the different tabs in the results panel (Overpressure, Thermal, Population, etc.) to see a detailed breakdown of the impact effects.
