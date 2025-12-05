/**
 * ENEO Asteroid Impact Simulation - Client-Side Application Controller
 * 
 * This JavaScript module manages the interactive web interface for asteroid impact simulations.
 * It coordinates user interactions, map visualizations, API communications, and result displays
 * to provide a comprehensive asteroid impact analysis platform.
 * 
 * Key Responsibilities:
 * - Platform-specific UI scaling and responsive design handling
 * - Interactive map management with damage zone visualization
 * - Simulation parameter input and validation
 * - Real-time API communication with backend simulation engine
 * - Dynamic result presentation with tabbed interface
 * - Geographic coordinate handling and zone rendering
 * - Integration with NASA Sentry API for real NEO data
 * 
 * Technical Features:
 * - Cross-platform compatibility with Windows-specific scaling
 * - Leaflet.js integration for advanced mapping capabilities
 * - Dynamic chart generation and effect visualization
 * - Responsive layout adaptation for mobile and desktop
 * - Real-time loading animations and progress tracking
 * 
 * Author: Alexandros Notas
 * Institution: National Technical University of Athens
 * Date: June 2025
 */

// Platform Detection and UI Scaling Management
// Handles Windows-specific scaling issues and responsive design across different screen sizes
document.addEventListener('DOMContentLoaded', function() {
  // Create application wrapper to enable controlled scaling without affecting modals
  const wrapper = document.createElement('div');
  wrapper.id = 'app-wrapper';
  while (document.body.firstChild) {
    wrapper.appendChild(document.body.firstChild);
  }
  document.body.appendChild(wrapper);

  // Detect Windows platform for targeted scaling adjustments
  const isWindows = navigator.userAgent.indexOf('Windows') !== -1;
  
  // Apply initial scaling based on platform and screen size
  applyAppropriateScaling(isWindows, wrapper);
  
  // Monitor window resize events for cross-screen compatibility
  if (isWindows) {
    let lastWidth = window.innerWidth;
    let resizeTimeout;
    
    window.addEventListener('resize', function() {
      // Clear any pending resize handler
      clearTimeout(resizeTimeout);
      
      // Set a timeout to prevent multiple rapid executions
      resizeTimeout = setTimeout(function() {
        // Check if width changed significantly (indicating potential monitor switch)
        const widthChange = Math.abs(window.innerWidth - lastWidth);
        if (widthChange > 200) { // Threshold for significant change
          console.log('Significant width change detected, reapplying scaling');
          applyAppropriateScaling(isWindows, wrapper);
          lastWidth = window.innerWidth;
        }
      }, 300);
    });
  }
  
  // Fix loading overlay to be position fixed
  const loadingOverlay = document.querySelector('.loading-overlay');
  if (loadingOverlay) {
    // Move loading overlay to be direct child of body
    document.body.appendChild(loadingOverlay);
    
    loadingOverlay.style.position = 'fixed';
    loadingOverlay.style.top = '0';
    loadingOverlay.style.left = '0';
    loadingOverlay.style.right = '0';
    loadingOverlay.style.bottom = '0';
    loadingOverlay.style.zIndex = '9999';
  }

  // Move the About modal to be a direct child of the body.
  // This prevents it from being scaled down with the app wrapper on Windows.
  const aboutModal = document.getElementById('aboutModal');
  if (aboutModal) {
    document.body.appendChild(aboutModal);
  }
});

// Function to apply scaling based on current window width
function applyAppropriateScaling(isWindows, wrapper) {
  if (!isWindows) return;
  
  const currentWidth = window.innerWidth;
  const needsScaling = currentWidth >= 992 && currentWidth < 1920;
  
  // Remove any existing scaling from wrapper
  wrapper.style.transform = '';
  wrapper.style.transformOrigin = '';
  wrapper.style.width = '';
  wrapper.style.marginLeft = '';
  wrapper.style.overflowX = '';
  wrapper.style.overflowY = '';
  wrapper.style.height = '';
  
  if (needsScaling) {
    // Apply platform-specific scaling for Windows
    document.documentElement.classList.add('windows-platform');
    
    // Apply transform scale to wrapper instead of body
    wrapper.style.transform = 'scale(0.9)';
    wrapper.style.transformOrigin = 'top center';
    wrapper.style.width = '111.11%'; // 100% ÷ 0.9 to compensate for scaling
    wrapper.style.marginLeft = '-5.55%'; // Center the scaled content
    wrapper.style.overflowX = 'hidden'; // Prevent horizontal scrollbar only
    wrapper.style.overflowY = 'auto'; // Allow vertical scrolling
    wrapper.style.height = 'auto'; // Allow content to determine height
    
    // Set body height to match scaled content
    document.body.style.height = '90vh';
    document.body.style.overflowY = 'auto';
    
    // Force refresh layout to ensure proper scaling
    setTimeout(function() {
      window.dispatchEvent(new Event('resize'));
    }, 100);
  } else {
    // Remove Windows platform class if scaling is no longer needed
    document.documentElement.classList.remove('windows-platform');
    
    // Reset body styles when not scaled
    document.body.style.height = '';
    document.body.style.overflowY = '';
  }
}

// Initialize map with world bounds
const map = L.map('map', {
  minZoom: 2,
  maxZoom: 18,
  maxBounds: [
    [-90, -180],
    [90, 180]
  ],
  maxBoundsViscosity: 1.0
}).setView([0, 0], 2);

// Create custom panes for crater visualization
// Pane for the crater image (non-interactive)
map.createPane('craterImagePane');
map.getPane('craterImagePane').style.zIndex = 450; // Default overlayPane is 400

// Pane for the crater border (interactive)
map.createPane('craterBorderPane');
map.getPane('craterBorderPane').style.zIndex = 451; // Above craterImagePane

// Use ESRI World Map tile layer
L.tileLayer('https://server.arcgisonline.com/ArcGIS/rest/services/World_Street_Map/MapServer/tile/{z}/{y}/{x}', {
  attribution: 'Tiles &copy; Esri',
  noWrap: true,
  bounds: [
    [-90, -180],
    [90, 180]
  ]
}).addTo(map);

let marker;
map.on('click', function(e) {
  // Check if the click is within the map's maxBounds
  if (map.options.maxBounds && !map.options.maxBounds.contains(e.latlng)) {
    alert(t("messages.validPointBoundaries", "Please select a valid point within the map boundaries."));
    return; // Stop further processing if click is out of bounds
  }

  if (marker) {
    map.removeLayer(marker);
  }
  marker = L.marker(e.latlng, {
    zIndexOffset: 1000  // Ensure marker stays on top
  }).addTo(map);
  document.getElementById('latitude').value = e.latlng.lat;
  document.getElementById('longitude').value = e.latlng.lng;

  // Fetch ocean depth but do not display it in a popup
  fetch('/get_ocean_depth', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ latitude: e.latlng.lat, longitude: e.latlng.lng })
  })
  .then(response => {
    if (!response.ok) {
      // Try to parse error json if server sent one
      return response.json().then(errData => {
        throw new Error(errData.error || `HTTP error ${response.status}`);
      }).catch(() => {
        // Fallback if error json parsing fails
        throw new Error(`HTTP error ${response.status}`);
      });
    }
    return response.json();
  })
  .then(data => {
    // Ocean depth is fetched, but we are not creating a popup.
    // You can still use the 'data.depth' if needed elsewhere.
    if (data.depth !== undefined) {
      if (data.depth > 0) {
        console.log(`Ocean Depth: ${data.depth.toFixed(2)} m`);
      } else {
        console.log(`Location: On Land`);
      }
    } else if (data.error) {
      console.log(`Depth: ${data.error}`);
    } else {
      console.log(`Depth data: N/A`);
    }
  })
  .catch(error => {
    console.error('Error fetching ocean depth:', error);
    // Error is logged, but no popup is shown.
  });
});

let zoneLayers = {
  seismic: new L.LayerGroup(),
  thermal: new L.LayerGroup(),
  airblast: new L.LayerGroup(),
  ejecta: new L.LayerGroup(),
  vulnerability: new L.LayerGroup(),
  crater: new L.LayerGroup(),
  wind: new L.LayerGroup(), // Add this line to register the wind layer
  tsunami: new L.LayerGroup()
};

function getZoneColor(type, value) {
  const opacity = 0.6;
  switch (type) {
    case 'vulnerability':
      if (value >= 0.75) return `rgba(255, 0, 0, ${opacity + 0.1})`;       // Darker Red for highest
      if (value >= 0.5) return `rgba(255, 120, 0, ${opacity})`;     // Orange
      if (value >= 0.25) return `rgba(255, 220, 0, ${opacity - 0.1})`;    // Yellow
      return `rgba(0, 255, 0, ${opacity - 0.2})`; // Lighter Green for lowest
    case 'seismic':
      return `rgba(255, 0, 0, ${opacity})`; // Red
    case 'thermal':
      return `rgba(255, 102, 0, ${opacity})`; // Orange
    case 'airblast': // This is for Overpressure
      return `rgba(0, 200, 200, ${opacity})`; // Teal/cyan for Overpressure
    case 'wind':
      return `rgba(0, 102, 255, ${opacity})`; // Blue for wind zones
    case 'ejecta':
      return `rgba(102, 51, 0, ${opacity})`; // Brown
    case 'crater': // Added case for crater circle
      return `rgba(128, 128, 128, ${opacity})`; // Grey for crater
    case 'tsunami':
      return `rgba(0, 100, 255, ${opacity})`; // Blue for tsunami
    default:
      return `rgba(102, 102, 102, ${opacity})`; // Gray
  }
}

function clearAllZones() {
  Object.values(zoneLayers).forEach(layer => {
    if (map.hasLayer(layer)) {
      map.removeLayer(layer);
    }
    layer.clearLayers();
  });
}

function showZones(type, zoneSpecificData) {
  clearAllZones();

  // If Tsunami type is selected for visualization, show nothing and fly to marker.
  // This effectively removes tsunami visualization from the map controls.
  if (type === 'tsunami') {
    if (marker) {
        // Fly to marker to indicate the location is still relevant, but no zones shown.
        map.flyTo(marker.getLatLng(), map.getZoom() < 5 ? 5 : map.getZoom());
    }
    // No zones are drawn for tsunami type.
    return;
  }

  if (type === 'none') {
    if (marker) { 
        map.flyTo(marker.getLatLng(), map.getZoom() < 5 ? 5 : map.getZoom());
    }
    return;
  }
  
  // Determine the actual Leaflet layer to use
  // For specific vulnerability types, use the general 'vulnerability' layer
  const actualZoneLayer = type.startsWith('vulnerability_') ? zoneLayers.vulnerability : zoneLayers[type];
  
  if (!actualZoneLayer) {
      console.warn(`No zone layer configured for type: ${type}`);
      if (marker) map.flyTo(marker.getLatLng(), map.getZoom());
      return;
  }

  // Handle crater type separately
  if (type === 'crater') {
    const fullResultsData = JSON.parse(document.querySelector('#resultsContainer').dataset.lastResults || '{}');
    if (fullResultsData && fullResultsData.crater && 
        fullResultsData.crater.crater_formation && 
        fullResultsData.crater.crater_formation.final_diameter > 0 && marker) {
      
      const craterDiameterMeters = fullResultsData.crater.crater_formation.final_diameter;
      const impactLatLng = marker.getLatLng();
      const craterRadiusMeters = craterDiameterMeters / 2;
      
      zoneLayers.crater.clearLayers();

      const craterCircle = L.circle(impactLatLng, {
        radius: craterRadiusMeters,
        color: '#654321',
        fillColor: '#654321',
        fillOpacity: 1.0,
        weight: 1.5,
        pane: 'craterBorderPane'
      });
      
      const craterLabel = typeof t === 'function' ? t('zones.crater', 'Crater') : 'Crater';
      const craterDiameterLabel = typeof t === 'function' ? t('results.craterDiameter', 'Crater Diameter') : 'Crater Diameter';

      craterCircle.bindPopup(`
        <strong>${craterLabel}</strong><br>
        ${craterDiameterLabel}: ${(craterDiameterMeters / 1000).toFixed(2)} km
      `);
      
      zoneLayers.crater.addLayer(craterCircle);
      map.addLayer(zoneLayers.crater);
      
      const bounds = craterCircle.getBounds();
      map.flyToBounds(bounds, { padding: [30, 30], maxZoom: Math.min(map.getMaxZoom(), 17) });

    } else {
      if (marker) {
          map.flyTo(marker.getLatLng(), map.getZoom());
      }
    }
    return; 
  }
  
  // For other zone types (including specific vulnerabilities)
  if (!zoneSpecificData || (Array.isArray(zoneSpecificData) && zoneSpecificData.length === 0)) {
      if (marker) {
          map.flyTo(marker.getLatLng(), map.getZoom());
      }
      return;
  }

  const initiallyFilteredAndSortedZones = [...zoneSpecificData]
    .filter(zone => {
      if (type === 'seismic') {
        const description = zone.description || '';
        // Match both English "Richter" and Greek "Ρίχτερ"
        const match = description.match(/(\d+)(?:\+|-\d+)\s*(?:Richter|Ρίχτερ)/);
        return match && parseInt(match[1]) >= 4;
      }
      return true;
    })
    .sort((a, b) => b.end_distance - a.end_distance);
  
  if (initiallyFilteredAndSortedZones.length === 0) {
    if (marker) map.flyTo(marker.getLatLng(), map.getZoom());
    return;
  }

  const visualizableZones = initiallyFilteredAndSortedZones.filter(zone => {
      if (type.startsWith('vulnerability_')) {
          if (zone.threshold === 1.0 && zone.end_distance <= 0.01) {
              return false;
          }
      }
      return true;
  });
  
  if (visualizableZones.length === 0) {
    if (marker) map.flyTo(marker.getLatLng(), map.getZoom());
    return;
  }

  const EARTH_RADIUS_KM = 10000;
  const EARTH_DIAMETER_KM = EARTH_RADIUS_KM * 2;
  const maxVisualizationDistance = EARTH_DIAMETER_KM;
  
  const maxDistance = Math.min(visualizableZones[0].end_distance, maxVisualizationDistance);
  
  visualizableZones.forEach(zone => {
    const radiusInMeters = Math.min(zone.end_distance, EARTH_DIAMETER_KM) * 1000;
      
    if (marker) {
      const impactLocation = marker.getLatLng();
      const crossesDateLine = Math.abs(impactLocation.lng) + (radiusInMeters / 111319.9) > 180;
        
      let popupContent = '';
      let zoneColorConfig;

      if (type.startsWith('vulnerability_')) {
        let specificVulnerabilityName;
        // Determine the translated name for the specific vulnerability type
        if (type === 'vulnerability_thermal') {
            specificVulnerabilityName = t('vulnerabilityTypes.thermal', 'Thermal-Induced Vulnerability');
        } else if (type === 'vulnerability_overpressure') {
            specificVulnerabilityName = t('vulnerabilityTypes.overpressure', 'Overpressure-Induced Vulnerability');
        } else if (type === 'vulnerability_wind') {
            specificVulnerabilityName = t('vulnerabilityTypes.wind', 'Wind-Induced Vulnerability');
        } else if (type === 'vulnerability_seismic') {
            specificVulnerabilityName = t('vulnerabilityTypes.seismic', 'Seismic-Induced Vulnerability');
        } else if (type === 'vulnerability_ejecta') {
            specificVulnerabilityName = t('vulnerabilityTypes.ejecta', 'Ejecta-Induced Vulnerability');
        } else { // Default to combined
            specificVulnerabilityName = t('vulnerabilityTypes.combined', 'Combined Vulnerability');
        }

        popupContent = `
          <strong>${specificVulnerabilityName}</strong><br>
          ${t('results.riskLevel', 'Risk Level')}: ${(zone.threshold * 100).toFixed(1)}%<br>
          ${t('results.range', 'Range')}: ${zone.start_distance.toFixed(2)} - ${zone.end_distance.toFixed(2)} km
        `;
        const color = getZoneColor('vulnerability', zone.threshold || 0.5);
        zoneColorConfig = { outline: color, fill: color };
      } else {
        popupContent = `
          <strong>${zone.description}</strong><br>
          ${t('results.range', 'Range')}: ${zone.start_distance.toFixed(2)} - ${zone.end_distance.toFixed(2)} km
        `;
        const fillColor = getZoneColor(type, zone.threshold || 0.5);
        const outlineColor = type === 'thermal' ? '#FF4500' : fillColor;
        zoneColorConfig = { outline: outlineColor, fill: fillColor };
      }
        
      const circle = L.circle(impactLocation, {
        radius: radiusInMeters,
        color: zoneColorConfig.outline,
        fillColor: zoneColorConfig.fill,
        fillOpacity: 0.4,
        weight: 1.5
      });
        
      circle.bindPopup(popupContent);
      actualZoneLayer.addLayer(circle);
        
      if (crossesDateLine) {
        const wrappedLng = impactLocation.lng > 0 ? 
          impactLocation.lng - 360 : 
          impactLocation.lng + 360;
        const wrappedLocation = L.latLng(impactLocation.lat, wrappedLng);
          
        const wrappedCircle = L.circle(wrappedLocation, {
          radius: radiusInMeters,
          color: zoneColorConfig.outline,
          fillColor: zoneColorConfig.fill,
          fillOpacity: 0.4,
          weight: 1.5
        });
          
        wrappedCircle.bindPopup(popupContent);
        actualZoneLayer.addLayer(wrappedCircle);
      }
    }
  });

  map.addLayer(actualZoneLayer);

  if (maxDistance > EARTH_RADIUS_KM) {
    map.setView([0, 0], 1);
    return;
  }

  if (marker) {
    try {
      const center = marker.getLatLng();
      const kmPerLongitudeDegree = 111.32 * Math.cos(center.lat * Math.PI / 180);
      const kmPerLatitudeDegree = 111.32;
      const adjustedMaxDistance = maxDistance * 1.3;
      const latDegrees = adjustedMaxDistance / kmPerLatitudeDegree;
      const lonDegrees = Math.abs(center.lat) > 85 ? 360 : adjustedMaxDistance / kmPerLongitudeDegree;
      
      const bounds = L.latLngBounds(
        [Math.max(-90, center.lat - latDegrees), center.lng - lonDegrees],
        [Math.min(90, center.lat + latDegrees), center.lng + lonDegrees]
      );

      const paddingPixels = Math.max(100, adjustedMaxDistance / 2);
      map.flyToBounds(bounds, {
        padding: [paddingPixels, paddingPixels],
        maxZoom: 10,
        duration: 0.75,
        animate: true,
        easeLinearity: 0.3
      });
      marker.bringToFront();
    } catch (e) {
      console.warn('Error fitting bounds:', e);
      const zoomLevel = Math.min(10, Math.max(1, 14 - Math.log2(maxDistance)));
      map.flyTo(marker.getLatLng(), zoomLevel, {
        duration: 0.75,
        animate: true,
        easeLinearity: 0.3
      });
    }
  }
}

function loadSampleData() {
  // Generate random values within logical ranges
  const diameter = Math.floor(Math.random() * 1451) + 50;    // 50-1500 (changed from 50-500)
  const density = Math.floor(Math.random() * 1001) + 2500;    // 2500-3500
  const velocity = Math.floor(Math.random() * 16) + 17;       // 17-32
  const entry_angle = Math.floor(Math.random() * 31) + 35;    // 35-65
  const distance = Math.floor(Math.random() * 91) + 10;       // 10-100
  
  // Set the form values
  document.getElementById("diameter").value = diameter;
  document.getElementById("density").value = density;
  document.getElementById("velocity").value = velocity;
  document.getElementById("entry_angle").value = entry_angle;
  document.getElementById("distance").value = distance;
}

/**
 * Fetches NEO data from the NASA Sentry API via the backend and populates the form.
 */
async function loadNasaData() {
  const nasaBtn = document.getElementById('nasaDataBtn');
  const originalText = nasaBtn.innerHTML;
  
  // Disable button and show loading state
  nasaBtn.disabled = true;
  nasaBtn.innerHTML = `<span class="spinner-border spinner-border-sm" role="status" aria-hidden="true"></span> ${t("loading.fetching", "Fetching...")}`;

  try {
    const response = await fetch('/generate-neo');
    if (!response.ok) {
      const errorData = await response.json().catch(() => ({ error: 'Failed to fetch data from NASA Sentry API.' }));
      throw new Error(errorData.error || `HTTP Error: ${response.status}`);
    }
    const data = await response.json();

    if (data) {
      // Populate form fields with data from the API
      // Use logical OR (||) to keep existing value if API data is null/undefined
      document.getElementById("diameter").value = data.diameter || document.getElementById("diameter").value;
      document.getElementById("velocity").value = data.velocity || document.getElementById("velocity").value;
      document.getElementById("density").value = data.density || document.getElementById("density").value;
      document.getElementById("entry_angle").value = data.entry_angle || document.getElementById("entry_angle").value;
      
      // Sentry API doesn't provide a reference distance, so we can use a sensible default or keep the existing value
      if (!document.getElementById("distance").value) {
        document.getElementById("distance").value = 100;
      }
      
      // Optionally, display a styled notification with the object's name
      if (data.name) {
        showNasaNotification(data.name);
      }

    } else {
      throw new Error(t("messages.receivedEmptyData", "Received empty data from the server."));
    }
  } catch (error) {
    console.error('Error fetching NASA data:', error);
    alert(`${t("messages.couldNotLoadNasa", "Could not load NASA data:")} ${error.message}`);
  } finally {
    // Restore button to its original state
    nasaBtn.disabled = false;
    nasaBtn.innerHTML = originalText;
  }
}

function parseResultsText(text) {
  const sections = {};
  const splitSections = text.split(/===+ (.*?) ===+/g).slice(1);

  for (let i = 0; i < splitSections.length; i += 2) {
    const sectionName = splitSections[i].trim();
    const content = splitSections[i + 1]
      .trim()
      .split("\n")
      .filter((line) => line.trim())
      .map((line) => {
        const [key, ...values] = line.split(":");
        return {
          key: key?.trim() || "",
          value: values.join(":").trim(),
        };
      });
    sections[sectionName.toLowerCase().replace(/ /g, "_")] = content;
  }
  return sections;
}

function createEffectCard(title, value, distance, subtitle = "") {
  if (
    distance !== undefined && // Check if distance is actually provided
    ["overpressure", "seismic", "thermal", "ejecta", "tsunami"].some((word) => 
      title.toLowerCase().includes(word)
    ) &&
    !title.toLowerCase().includes(` ${t("translations.at", "at")} `)
  ) {
    title = `${title} ${t("translations.at", "at")} ${distance} km`;
  }
  return `
    <div class="col-md-6">
      <div class="effect-card">
        <div class="effect-label">${title}</div>
        <div class="dynamic-value">${value}</div>
        ${subtitle ? `<small>${subtitle}</small>` : ""}
      </div>
    </div>
  `;
}

function createDangerZone(label, range, distance) {
  if (
    ["overpressure", "seismic", "thermal", "ejecta", "tsunami"].some((word) => 
      label.toLowerCase().includes(word)
    ) &&
    !label.toLowerCase().includes(` ${t("translations.at", "at")} `)
  ) {
    label = `${label} ${t("translations.at", "at")} ${distance} km`;
  }
  return `
    <div class="col-12">
      <div class="effect-card">
        <div class="effect-label">${label}</div>
        <div class="dynamic-value">${range}</div>
      </div>
    </div>
  `;
}

function populateComparisonBoxes(comparisons) {
  // Clear all previous comparison content
  document.querySelectorAll('.comparison-box').forEach(box => {
    box.innerHTML = '';
    box.classList.remove('visible');
  });

  if (!comparisons) return;

  // Energy comparison (for Overpressure tab, as it's energy-driven)
  if (comparisons.energy) {
    const energyHtml = `
      <i class="mdi mdi-scale-balance comparison-icon"></i>
      <span class="comparison-text">${comparisons.energy}</span>
    `;
    const overpressureBox = document.getElementById('overpressureComparison');
    if (overpressureBox) {
        overpressureBox.innerHTML = energyHtml;
        overpressureBox.classList.add('visible');
    }
  }

  // Seismic comparison — show the comparison sentence and the nearest historical quake (with magnitude)
  if (comparisons.seismic || comparisons.seismic_event) {
    const seismicBox = document.getElementById('seismicComparison');
    const mainText = comparisons.seismic ? `<div class="comparison-text">${comparisons.seismic}</div>` : '';
    let eventLine = '';
    if (comparisons.seismic_event && comparisons.seismic_event.formatted) {
      eventLine = `<div style="margin-top:6px;font-size:0.95rem;color:var(--text-primary);opacity:0.95"><strong>Compared to:</strong> ${comparisons.seismic_event.formatted}</div>`;
    }
    if (seismicBox) {
      seismicBox.innerHTML = `
        <i class="mdi mdi-pulse comparison-icon"></i>
        <div>
          ${mainText}
          ${eventLine}
        </div>
      `;
      seismicBox.classList.add('visible');
    }
  }

  // Wind comparison — show the comparison sentence and the nearest hurricane with wind speeds
  if (comparisons.wind || comparisons.wind_event) {
    const windBox = document.getElementById('windComparison');
    const mainText = comparisons.wind ? `<div class="comparison-text">${comparisons.wind}</div>` : '';
    let eventLine = '';
    if (comparisons.wind_event && comparisons.wind_event.formatted) {
      eventLine = `<div style="margin-top:6px;font-size:0.95rem;color:var(--text-primary);opacity:0.95"><strong>Comparable storm:</strong> ${comparisons.wind_event.formatted}</div>`;
    }
    if (windBox) {
      windBox.innerHTML = `
        <i class="mdi mdi-weather-windy-variant comparison-icon"></i>
        <div>
          ${mainText}
          ${eventLine}
        </div>
      `;
      windBox.classList.add('visible');
    }
  }
}

// Add this helper function at the top of the script or near handleSubmit
function createImpactOverviewItem(label, value, unit = '') {
  let displayValue = "N/A";
  if (value !== null && typeof value !== 'undefined') {
    if (typeof value === 'number') {
      // Check if the number is very large (likely Joules) or very small for exponential
      if (Math.abs(value) > 1e6 || (Math.abs(value) < 1e-2 && Math.abs(value) !== 0)) {
        displayValue = value.toExponential(2);
      } else if (unit === "km/s" || unit === "km" || unit === "MT") {
        displayValue = value.toFixed(2);
      } else if (unit === "m") {
        displayValue = value.toFixed(0);
      } else { // Default numeric formatting
        displayValue = parseFloat(value.toFixed(3)).toString(); // Handles cases like 0.000
      }
    } else {
      displayValue = value.toString();
    }
    displayValue += ` ${unit}`;
  }

  return `
    <div class="col-lg-4 col-md-6 mb-3">
      <div class="dynamic-value">${displayValue.trim()}</div>
      <div>${label}</div>
    </div>
  `;
}

async function handleSubmit(e, overrideParams = null) {
  e.preventDefault();
  
  // Check if a location has been selected on the map
  const latitude = document.getElementById('latitude').value;
  const longitude = document.getElementById('longitude').value;
  
  if (!latitude || !longitude || !marker) {
    showCustomAlert(t('messages.selectLocationFirst', 'Please select a location on the map before running the simulation.'));
    // Scroll to the map
    const mapElement = document.getElementById('map');
    if (mapElement) {
      mapElement.scrollIntoView({ behavior: 'smooth', block: 'center' });
    }
    return;
  }
  
  const loadingOverlay = document.querySelector(".loading-overlay");
  const resultsContainer = document.getElementById("resultsContainer");
  
  try {
    loadingOverlay.style.display = "flex";
    resultsContainer.classList.add("d-none");
    
    // Array of interesting asteroid impact facts
    const impactFacts = [
      t('facts.chicxulub') || "The Chicxulub impact that killed the dinosaurs released 100 trillion tons of TNT equivalent—over a billion times the energy of the Hiroshima bomb.",
      t('facts.tunguska') || "The Tunguska event in 1908 was caused by a ~50-60 m object that air-burst over Siberia, flattening about 2,150 km² of forest and felling 80 million trees.",
      t('facts.kilometer') || "A kilometer-scale asteroid impacts Earth only about once every 500,000 years on average.",
      t('facts.chelyabinsk') || "The Chelyabinsk meteor in 2013 was only ~17-20 m across but injured ~1,500 people and damaged ~7,200 buildings, mostly from broken glass.",
      t('facts.vredefort') || "South Africa's Vredefort Dome is the largest confirmed terrestrial impact structure, originally ~300 km across, formed 2.023 billion years ago.",
      t('facts.dart') || "NASA's DART mission struck the 170 m-wide moonlet Dimorphos at ~6.6 km/s in 2022, shortening its orbit by 32 minutes—demonstrating kinetic-impactor deflection.",
      t('facts.meteorites') || "Meteorites fall into three broad types: stony (silicate-rich), iron (metallic Ni-Fe), and stony-iron (mixed), each with numerous subgroups.",
      t('facts.kineticEnergy') || "An asteroid's kinetic energy scales as ½mv²—doubling its velocity quadruples its energy—making speed as critical as mass for impact hazards.",
      t('facts.neoCount') || "As of December 2024, over 37,000 Near-Earth Objects have been cataloged, with 2,465 classified as potentially hazardous asteroids.",
      t('facts.mainBelt') || "Although the Main Belt contains millions of asteroids, their average separation is roughly 965,000 km—creating vast voids between objects.",
      t('facts.meteoroidSpeed') || "Meteoroids strike Earth at between 11 and 72 km/s; those above 20 km/s ionize the air ahead, creating bright plasma trails.",
      t('facts.hundredMeter') || "Objects ≥100 m diameter can release tens to hundreds of megatons of TNT equivalent, causing regional devastation.",
      t('facts.earthAccretion') || "Earth accretes ~54 tons of interplanetary material daily, totaling ~19,700 tons per year of dust and micrometeoroids.",
      t('facts.airBurst') || "Bodies <50 m often air-burst (like Tunguska), delivering more widespread damage than crater-forming ground impacts of similar mass.",
      t('facts.barringer') || "The Barringer Crater in Arizona is 1.186 km in diameter and 170 m deep, formed ~50,000 years ago by a ~50 m iron meteorite.",
      t('facts.hoba') || "The Hoba meteorite in Namibia is the largest known intact meteorite at over 60 tonnes, composed almost entirely of iron-nickel alloy.",
      t('facts.neoSurveyor') || "NASA's NEO Surveyor space telescope, planned for launch by 2027, aims to complete the survey of 90% of NEOs ≥140 m.",
      t('facts.phas') || "Of the known near-Earth asteroids, 2,465 are classified as Potentially Hazardous Asteroids—large enough and close enough to pose a significant impact risk.",
      t('facts.torino') || "The Torino scale rates NEO impact threat from 0 (none) to 10 (certain global catastrophe); only asteroid Apophis ever reached level 4.",
      t('facts.dailyMeteoroids') || "Approximately 25 million meteoroids enter Earth's atmosphere each day, ranging from dust-sized particles to occasional boulder-sized objects."
    ];
    
    // Shuffle the array to get a random order
    const shuffledFacts = [...impactFacts].sort(() => 0.5 - Math.random());
    let currentFactIndex = 0;
    
    // Create timeline container with initial fact
    loadingOverlay.innerHTML = `
      <div class="impact-timeline">
        <div class="timeline-progress"></div>
        <div class="timeline-step active" data-step="entry">
          <div class="step-icon"><i class="mdi mdi-meteor"></i></div>
          <div class="step-label">${t('loading.atmosphericEntry') || "Atmospheric Entry"}</div>
          <div class="step-details" id="entryDetails"></div>
        </div>
        <div class="timeline-step" data-step="impact">
          <div class="step-icon"><i class="mdi mdi-explosion"></i></div>
          <div class="step-label">${t('loading.impactEffects') || "Impact Effects"}</div>
          <div class="step-details" id="impactDetails"></div>
        </div>
        <div class="timeline-step" data-step="zones">
          <div class="step-icon"><i class="mdi mdi-map-marker-radius"></i></div>
          <div class="step-label">${t('loading.vulnerabilityZones') || "Vulnerability Zones"}</div>
          <div class="step-details" id="zonesDetails"></div>
        </div>
        <div class="timeline-step" data-step="population">
          <div class="step-icon"><i class="mdi mdi-account-group"></i></div>
          <div class="step-label">${t('loading.populationImpact') || "Population Impact"}</div>
          <div class="step-details" id="populationDetails"></div>
        </div>
      </div>
      <div class="loading-fact" id="factContainer">
        <span class="fact-prefix">${t('loading.didYouKnow') || "Did you know?"}</span>
        <span class="fact-text">${shuffledFacts[currentFactIndex]}</span>
      </div>
    `;
    
    // Set up fact rotation
    const factRotationInterval = setInterval(() => {
      currentFactIndex = (currentFactIndex + 1) % shuffledFacts.length;
      const factContainer = document.getElementById('factContainer');
      if (factContainer) {
        // Add fade-out animation
        factContainer.classList.add('fact-fade');
        
        // Change the text after fade-out and remove the class to fade in
        setTimeout(() => {
          factContainer.querySelector('.fact-text').textContent = shuffledFacts[currentFactIndex];
          factContainer.classList.remove('fact-fade');
        }, 500);
      } else {
        // If container no longer exists (loading complete), clear the interval
        clearInterval(factRotationInterval);
      }
    }, 12000); // Change fact every 12 seconds
    
    // Simulate progress through timeline steps with updated text
    animateTimelineStep('entry', t('loading.calculatingEntry') || 'Calculating entry velocity and atmospheric effects...', 2500);
    
    setTimeout(() => {
      animateTimelineStep('impact', t('loading.simulatingEffects') || 'Simulating seismic, airblast, thermal and ejecta effects...', 3000);
    }, 3000);
    
    setTimeout(() => {
      animateTimelineStep('zones', t('loading.calculatingZones') || 'Calculating vulnerability zones...', 3000);
    }, 6500);
    
    setTimeout(() => {
      animateTimelineStep('population', t('loading.estimatingCasualties') || 'Estimating casualties...', 3000);
    }, 10000);
    
    // Store simulation parameters for potential re-use when language changes
    const simulationParams = overrideParams || {
      diameter: parseFloat(document.getElementById("diameter").value),
      density: parseFloat(document.getElementById("density").value),
      velocity: parseFloat(document.getElementById("velocity").value),
      entry_angle: parseFloat(document.getElementById("entry_angle").value),
      distance: parseFloat(document.getElementById("distance").value),
      latitude: parseFloat(document.getElementById("latitude").value || 0),
      longitude: parseFloat(document.getElementById("longitude").value || 0),
      language: currentLanguage
    };
    
    // Always update language to current language (important when switching languages)
    simulationParams.language = currentLanguage;

    // Make the actual API call
    const response = await fetch("/simulate", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(simulationParams),
    });

    if (!response.ok) throw new Error(`HTTP error! status: ${response.status}`);

    const data = await response.json();
    
    // Store the full results data as a JSON string in a data attribute
    resultsContainer.dataset.lastResults = JSON.stringify(data.results_data);
    
    // Store the simulation parameters for potential re-use when language changes
    resultsContainer.dataset.lastSimulationParams = JSON.stringify(simulationParams);

    // removed: populateComparisonBoxes(data.results_data.comparisons)
    // ---- START: New Impact Overview Population Logic ----
    const impactOverviewContainer = document.getElementById("impactOverviewContent");
    impactOverviewContainer.innerHTML = ''; // Clear previous content

    const eventType = data.results_data.atmospheric_entry.event_type;
    const atmEntry = data.results_data.atmospheric_entry;
    const energyData = data.results_data.energy;
    const craterData = data.results_data.crater?.crater_formation;

    let overviewHTML = '';

    // Common item: Initial Kinetic Energy
    overviewHTML += createImpactOverviewItem(t("results.initialKineticEnergy", "Initial Kinetic Energy (MT TNT)"), energyData.initial_energy_megatons, "MT TNT");

    // Breakup Altitude - show if it's a breakup event (ground impact or airburst)
    if (eventType === "ground impact" || eventType === "airburst") {
        overviewHTML += createImpactOverviewItem(t("results.breakupAltitude", "Breakup Altitude"), atmEntry.breakup_altitude, "m");
    } else if (eventType === "intact") {
        // For intact objects, breakup altitude is not applicable.
        // overviewHTML += createImpactOverviewItem(t("results.breakupAltitude", "Breakup Altitude"), "N/A (Object Intact)"); // Or omit
    }

    if (eventType === "airburst") {
        overviewHTML += createImpactOverviewItem(t("results.airburstAltitude", "Airburst Altitude"), atmEntry.airburst_altitude, "m");
        overviewHTML += createImpactOverviewItem(t("results.postAirburstVelocity", "Post Airburst Velocity"), atmEntry.final_velocity ? atmEntry.final_velocity / 1000 : null, "km/s");
        overviewHTML += createImpactOverviewItem(t("results.airburstEnergy", "Airburst Energy (MT TNT)"), energyData.specific_energy_megatons, "MT TNT");
    } else if (eventType === "ground impact" || eventType === "intact") {
        overviewHTML += createImpactOverviewItem(t("results.impactSpeed", "Impact Speed"), atmEntry.final_velocity ? atmEntry.final_velocity / 1000 : null, "km/s");
        overviewHTML += createImpactOverviewItem(t("results.impactEnergy", "Impact Energy (MT TNT)"), energyData.specific_energy_megatons, "MT TNT");
        overviewHTML += createImpactOverviewItem(t("results.transientCraterDiameter", "Transient Crater Diameter"), craterData?.transient_diameter ? craterData.transient_diameter / 1000 : null, "km");
        overviewHTML += createImpactOverviewItem(t("results.finalCraterDiameter", "Final Crater Diameter"), craterData?.final_diameter ? craterData.final_diameter / 1000 : null, "km");
    } else {
        // Fallback for any other unexpected event types
        overviewHTML += `<div class="col-12"><p class="text-white">Overview for event type: "${eventType}"</p></div>`;
        overviewHTML += createImpactOverviewItem(t("results.finalVelocity", "Final Velocity"), atmEntry.final_velocity ? atmEntry.final_velocity / 1000 : null, "km/s");
    }
    impactOverviewContainer.innerHTML = overviewHTML;
    // ---- END: New Impact Overview Population Logic ----

    // Remove or comment out the old direct DOM manipulations for energyValue, velocityValue, craterValue
    /*
    const energyValueElement = document.getElementById("energyValue");
    const energyLabelElement = document.getElementById("energyLabel"); // Get the label element
    // const eventType = data.results_data.atmospheric_entry.event_type; // Already defined above
    const specificEnergyMt = parseFloat(data.results_data.energy.specific_energy_megatons).toFixed(2);
    const specificEnergyType = data.results_data.energy.specific_energy_type;

    if (specificEnergyType === "Impact" && (eventType === "ground impact" || eventType === "intact")) {
      energyLabelElement.textContent = "Impact Energy (MT)";
      energyValueElement.textContent = `${specificEnergyMt} MT`;
    } else if (specificEnergyType === "Airburst" && eventType === "airburst") {
      energyLabelElement.textContent = "Airburst Energy (MT)";
      energyValueElement.textContent = `${specificEnergyMt} MT`;
    } else {
      energyLabelElement.textContent = "Initial Energy (MT)"; // Fallback
      energyValueElement.textContent = `${parseFloat(data.results_data.energy.initial_energy_megatons).toFixed(2)} MT`;
    }
    
    document.getElementById("velocityValue").textContent = 
      `${(data.results_data.atmospheric_entry.final_velocity / 1000).toFixed(2)} km/s`;
    document.getElementById("craterValue").textContent = 
      data.results_data.crater?.crater_formation?.final_diameter
        ? `${(data.results_data.crater.crater_formation.final_diameter / 1000).toFixed(2)} km`
        : "N/A";
    */

    const sections = parseResultsText(data.results_text);
    populateTab("overpressureEffects", sections.airblast || []); 
    populateTab("thermalEffects", sections.thermal_effects || []);
    populateTab("seismicEffects", sections.seismic_effects || []);
    populateTab("ejectaEffects", sections.ejecta || []);
    populateTab("windEffects", sections.wind_effects || []);
    populateTab("tsunamiEffects", sections.tsunami_effects || []);

    if (sections.vulnerability_models) {
      const vulnContainer = document.getElementById("nav-vulnerability").querySelector('.row.g-4'); // Target the inner row
      const vulnValues = {};
      let totalVuln = 0;
      
      sections.vulnerability_models.forEach(item => {
        if (!item.key.toLowerCase().includes('combined vulnerability')) {
          const value = item.value;
          if (!value.includes("N/A")) {
            const vulnValue = parseFloat(value);
            if (!isNaN(vulnValue)) {
              vulnValues[item.key] = vulnValue;
              totalVuln += vulnValue;
            }
          }
        }
      });

      const vulnerabilityCards = sections.vulnerability_models
        .map(item => {
          if (item.key.toLowerCase().includes('combined vulnerability')) {
            return '';
          }
          const [key, value] = [item.key, item.value];
          let formattedValue = value;
          let vulnValue = 0;
          let relativeContribution = 0;
          
          if (value.includes("N/A")) {
            formattedValue = "Not Applicable";
          } else {
            vulnValue = parseFloat(value);
            if (!isNaN(vulnValue)) {
              relativeContribution = totalVuln > 0 ? (vulnValue / totalVuln) * 100 : 0;
              formattedValue = vulnValue < 0.0001 ? 
                  "< 0.1%" : 
                  `${(vulnValue * 100).toFixed(1)}%`;
            }
          }
          
          let indicatorClass = "bg-success";
          if (!value.includes("N/A") && vulnValue > 0) {
            if (vulnValue >= 0.7) {
              indicatorClass = "bg-danger";
            } else if (vulnValue >= 0.3) {
              indicatorClass = "bg-warning";
            }
          }

          return `
              <div class="col-md-6 mb-3">
                  <div class="effect-card">
                      <div class="effect-label">${key}</div>
                      <div class="vulnerability-value ${value.includes("N/A") ? 'text-muted' : ''}">
                          ${formattedValue}
                          ${!value.includes("N/A") && vulnValue > 0 ? 
                              `<small class="vulnerability-contribution">
                  (${relativeContribution.toFixed(1)}% ${t("vulnerabilityText.ofTotal", "of total vulnerability")})
                              </small>` : 
                              ''}
                      </div>
                      ${!value.includes("N/A") && vulnValue > 0 ? `
                          <div class="vulnerability-scale">
                              <div class="${indicatorClass}" 
                                   style="width: ${relativeContribution}%; height: 100%; border-radius: 2px;">
                              </div>
                          </div>
                      ` : ''}
                  </div>
              </div>
          `;
        })
        .filter(card => card !== '')
        .join("");

      vulnContainer.innerHTML = `
          <div class="row g-4">
              ${vulnerabilityCards}
          </div>
      `;
    }

    if (data.results_data.population_analysis) {
      const populationContainer = document.getElementById("populationResults");
      const totalCasualties = document.getElementById("totalCasualties");
      totalCasualties.textContent = data.results_data.population_analysis.total_casualties.toLocaleString();
      populationContainer.innerHTML = data.results_data.population_analysis.zones
        .map(zone => `
          <div class="col-12">
            <div class="population-card">
              <div class="effect-label">${t('results.vulnerabilityZone') || "Vulnerability Zone"} ${(zone.vulnerability_threshold * 100).toFixed(1)}%</div>
              <div class="population-zone">
                <div>${t('results.distanceRange') || "Distance Range"}: ${zone.start_distance.toFixed(2)} - ${zone.end_distance.toFixed(2)} km</div>
              </div>
            </div>
          </div>
        `)
        .join("");
    }

    // BEGIN: Code from "Εικόνες file" for Affected Countries
    if (data.results_data.population_analysis && data.results_data.population_analysis.countries) {
      const countriesContainer = document.getElementById("countriesResults");
      const countries = data.results_data.population_analysis.countries;
      
      // Sort countries by casualties (descending)
      countries.sort((a, b) => b.total_casualties - a.total_casualties);
      
      // Create HTML for each country
      if (countriesContainer) {
        countriesContainer.innerHTML = countries
          .map(country => {
            // Use translated name from backend if available, otherwise try frontend translation
            const displayName = country.translated_name || 
                               t(`countries.${country.name}`, country.name || `Country FID ${country.fid}`);
            
            return `
              <div class="col-md-6 mb-3">
                <div class="country-card effect-card">
                  <div class="effect-label">${displayName}</div>
                  <div class="country-stats">
                    <div class="row">
                      <div class="col-12">
                        <div class="stat-label">${t('results.totalCasualties', 'Total Casualties')}</div>
                        <div class="stat-value text-danger">${country.total_casualties.toLocaleString()}</div>
                      </div>
                    </div>
                  </div>
                </div>
              </div>
            `;
          })
          .join('');
        
        // If no countries found
        if (countries.length === 0) {
          countriesContainer.innerHTML = `
            <div class="col-12">
              <div class="effect-card text-center">
                <p class="mb-0">${t('messages.noCountryData', 'No country data available')}</p>
              </div>
            </div>
          `;
        }
      }
    }
    // END: Code from "Εικόνες file" for Affected Countries

    // BEGIN: Code from "Εικόνες file" for Economic Impact
    if (data.results_data.economic_analysis) {
      const economicContainer = document.getElementById("economicResults");
      const totalEconomicDamage = document.getElementById("totalEconomicDamage");
      const economicData = data.results_data.economic_analysis;
      
      if (totalEconomicDamage) {
        totalEconomicDamage.textContent = new Intl.NumberFormat('en-US', { 
          style: 'currency', 
          currency: 'USD', 
          notation: 'compact',
          compactDisplay: 'short',
          maximumFractionDigits: 1
        }).format(economicData.total_economic_damage);
      }
      
      if (economicContainer) {
        economicContainer.innerHTML = economicData.countries
          .filter(country => country.economic_damage !== null) 
          .map(country => {
            // Use translated name from backend if available, otherwise try frontend translation
            const displayName = country.translated_name || 
                               t(`countries.${country.name}`, country.name);
            
            const formattedDamage = new Intl.NumberFormat('en-US', { 
              style: 'currency', 
              currency: 'USD', 
              notation: 'compact',
              maximumFractionDigits: 1
            }).format(country.economic_damage);

            const formattedGdpPerCapita = new Intl.NumberFormat('en-US', {
              style: 'currency',
              currency: 'USD',
              maximumFractionDigits: 0
            }).format(country.gdp_per_capita);
            
            return `
              <div class="col-md-6 mb-3">
                <div class="economy-card effect-card">
                  <div class="effect-label">${displayName}</div>
                  <div class="economy-stats">
                    <div class="row">
                      <div class="col-md-6">
                        <div class="stat-label">${t('results.economicDamage', 'Economic Damage')}</div>
                        <div class="stat-value text-danger">${formattedDamage}</div>
                      </div>
                      <div class="col-md-6">
                        <div class="stat-label">${t('results.gdpPerCapita', 'GDP Per Capita')} (${country.gdp_year})</div>
                        <div class="stat-value">${formattedGdpPerCapita}</div>
                      </div>
                    </div>
                  </div>
                </div>
              </div>
            `;
          })
          .join('');
        
        if (economicContainer.innerHTML.trim() === '') {
          economicContainer.innerHTML = `
            <div class="col-12">
              <div class="effect-card text-center">
                <p class="mb-0">${t('messages.noGdpData', 'No GDP data available for affected countries')}</p>
              </div>
            </div>
          `;
        } else {
          economicContainer.innerHTML += `
            <div class="col-12 mt-3">
              <div class="effect-card">
                <div class="effect-label text-white fw-bold">${t('economic.aboutCalculation', 'About Economic Impact Calculation')}</div>
                <p class="mb-2 text-highlight-orange">${t('economic.formula', 'Economic damage is calculated using the formula:')}</p>
                <p class="mb-0"><strong class="text-white">${t('economic.formulaText', 'Casualties × GDP per capita')}</strong></p>
                <p class="mb-0 mt-2"><small class="text-muted text-highlight-orange">${t('economic.explanation', 'This represents the direct economic impact based on casualties and national GDP per capita.')}</small></p>
              </div>
            </div>
          `;
        }
      }
    }
    // END: Code from "Εικόνες file" for Economic Impact

    // The following lines were part of the original active file, ensure they are correctly placed
    // relative to the newly added blocks.
    // This was the old showZones call, it's now handled by the map update logic below.
    // const zoneType = document.querySelector('input[name="zonetype"]:checked').value;
    // if (zoneType !== 'none') {
    //   showZones(zoneType, data.results_data.visualization[zoneType]);
    // }

    if (data.results_data) {
      document.querySelector('#resultsContainer').dataset.lastResults = JSON.stringify(data.results_data);
    }
    
    resultsContainer.classList.remove("d-none");

    // Update map display based on current selection
    const checkedZoneTypeRadio = document.querySelector('input[name="zonetype"]:checked');
    const vulnerabilitySelectorDiv = document.getElementById('vulnerabilitySelectorDiv');
    const vulnerabilityTypeSelect = document.getElementById('vulnerabilityTypeSelect');

    if (checkedZoneTypeRadio && data.results_data && data.results_data.visualization) {
        const mainZoneType = checkedZoneTypeRadio.value;
        if (mainZoneType === 'vulnerability') {
            vulnerabilitySelectorDiv.style.display = 'block';
            const specificVulnerabilityType = vulnerabilityTypeSelect.value;
            if (data.results_data.visualization[specificVulnerabilityType]) {
                showZones(specificVulnerabilityType, data.results_data.visualization[specificVulnerabilityType]);
            } else {
                console.warn(`Initial load: No visualization data for ${specificVulnerabilityType}.`);
                clearAllZones();
            }
        } else if (mainZoneType !== 'none') {
            vulnerabilitySelectorDiv.style.display = 'none';
            if (mainZoneType === 'crater') { // Specifically handle crater
                showZones('crater', null); // Crater data is fetched internally by showZones
            } else if (data.results_data.visualization[mainZoneType]) {
                showZones(mainZoneType, data.results_data.visualization[mainZoneType]);
            } else {
                console.warn(`Initial load: No visualization data for ${mainZoneType}.`);
                clearAllZones();
            }
        } else {
            vulnerabilitySelectorDiv.style.display = 'none';
            showZones('none', null);
        }
    }

    // addSensitivityAnalysisButton(); // Removed call to add sensitivity analysis button

  } catch (error) {
    console.error('Error in simulation:', error);
    alert(`Error: ${error.message}`);
  } finally {
    loadingOverlay.style.display = "none";
    resultsContainer.classList.remove("d-none");
  }
}

function animateTimelineStep(step, detail, duration) {
  const allSteps = document.querySelectorAll('.timeline-step');
  const currentStep = document.querySelector(`.timeline-step[data-step="${step}"]`);
  const detailElement = document.getElementById(`${step}Details`);
  
  // Update active step
  allSteps.forEach(s => s.classList.remove('active'));
  currentStep.classList.add('active');
  
  // Update progress bar
  const stepIndex = Array.from(allSteps).findIndex(s => s.dataset.step === step);
  const progress = ((stepIndex + 1) / allSteps.length) * 100;
  document.querySelector('.timeline-progress').style.width = `${progress}%`;
  
  // Show details with typing effect
  if (detailElement) {
    const text = detail;
    detailElement.textContent = '';
    let i = 0;
    const typeInterval = setInterval(() => {
      if (i < text.length) {
        detailElement.textContent += text.charAt(i);
        i++;
      } else {
        clearInterval(typeInterval);
      }
    }, duration / text.length);
  }
}

function simulateCasualtyCounter() {
  // This function is now empty - we don't need it since the animateTimelineStep
  // already displays the message "Estimating casualties..."
}

function populateTab(tabId, items) {
  const container = document.getElementById(tabId);
  if (!container) return;
  container.innerHTML = ''; // Clear previous content

  // Special handling for the Tsunami tab
  if (tabId === "tsunamiEffects") {
    let tsunamiContentHTML = `
      <div class="col-12">
        <div class="effect-card">
          <h5 class="effect-label">${t("tsunamiDisclaimer.title", "Note on Tsunami Calculations")}</h5>
          <p class="mb-1 mt-2" style="font-size: 0.9rem; line-height: 1.6;">
            ${t("tsunamiDisclaimer.simplified", "The tsunami effects presented here are based on a")} <strong>${t("tsunamiDisclaimer.simplifiedBold", "simplified approximation")}</strong>${t("tsunamiDisclaimer.explanation", ". This model considers initial wave amplitude at the source and basic geometric spreading. The initial wave amplitude calculation is validated against the conditions described in the tsunami algorithm by Rumpf et al. (2016), serving as a basis for future enhancements. More advanced tsunami models were tested but proved too computationally intensive for real-time online purposes; their integration is a potential future improvement.")}
          </p>
          <p class="mb-1 mt-1" style="font-size: 0.9rem; line-height: 1.6;">
            ${t("tsunamiDisclaimer.requirements", "Accurate, operational tsunami modeling is significantly more complex and requires:")}
          </p>
          <ul style="font-size: 0.9rem; line-height: 1.6; margin-top: 0.5rem; margin-bottom: 0.5rem;">
            <li>${t("tsunamiDisclaimer.detailedBathymetry", "Detailed bathymetry (sea floor depth, slope, and specific underwater features)")}</li>
            <li>${t("tsunamiDisclaimer.coastalTopography", "High-resolution coastal topography")}</li>
            <li>${t("tsunamiDisclaimer.waveDynamics", "Consideration of complex wave dynamics (e.g., refraction, diffraction, shoaling, reflection, and run-up)")}</li>
            <li>${t("tsunamiDisclaimer.hydrodynamicEquations", "Advanced hydrodynamic equations and numerical simulations")}</li>
          </ul>
          <p class="mb-0 mt-2" style="font-size: 0.9rem; line-height: 1.6; font-weight: bold;">
            ${t("tsunamiDisclaimer.limitationsBold", "Due to these simplifications, the results from this model are for general informational purposes only and are not visualized on the map.")}
          </p>
        </div>
      </div>
    `;

    // Display textual data if available from 'items'
    let hasTextualData = false;
    if (items && items.length > 0) {
      items.forEach(item => {
        if (item.value && item.value.trim().toLowerCase() !== "none") {
          if (!hasTextualData) {
            tsunamiContentHTML += `
              <div class="col-12 mt-3">
                <div class="effect-card">
                  <h6 class="effect-label">${t("tsunamiDisclaimer.approximateWave", "Approximate Wave Amplitude Distances:")}</h6>
            `;
            hasTextualData = true;
          }
          tsunamiContentHTML += `
            <div style="padding: 0.3rem 0; border-bottom: 1px solid rgba(255,255,255,0.05);">
                <strong>${item.key}:</strong> ${item.value}
            </div>
          `;
        }
      });
      if (hasTextualData) {
        // Close the card if data was added
        tsunamiContentHTML = tsunamiContentHTML.replace(/<\/div>\s*$/, ''); // Remove last closing div if it's just a border
        tsunamiContentHTML += `</div></div>`; // Close the card and col-12
      }
    }

    if (!hasTextualData) {
      // If, after checking items, no textual data was found, add a message to the disclaimer card.
      // This requires finding the disclaimer card and appending to it.
      // A simpler way is to just add another card.
       tsunamiContentHTML += `
         <div class="col-12 mt-3">
           <div class="effect-card">
             <p class="mb-0">${t('tsunamiDisclaimer.noSignificantEffects') || "No significant tsunami effects calculated by the simplified model for the current parameters, or data is not applicable."}</p>
           </div>
         </div>
       `;
    }
    container.innerHTML = tsunamiContentHTML;
    return; // Exit after handling tsunami tab
  }

  // Existing logic for other tabs
  const locationEffects = [];
  const tabSpecificDangerZones = []; 
  
  items.forEach(item => {
    // Skip any entries with "None" value in their value field
    if (item.value.trim().toLowerCase() === "none") {
      return;
    }
    // Also skip if the key suggests it's a section header but has no value,
    // or if the value is effectively an empty range for danger zones.
    if ((item.key.toLowerCase().includes("danger zones") && item.value.trim() === "") ||
        (item.value.includes("-") && item.value.includes("km") && 
         (item.value.match(/0\.0[0-9]-0\.0[0-9]\s*km/) || 
          item.value.match(/0\.00-0\.0[0-9]\s*km/) ||
          item.value.match(/0\.00-0\.1[0-9]\s*km/) || 
          item.value === "0.00-0.01 km"))) {
            // console.log(`Skipping empty/zero range danger zone: ${item.key}: ${item.value}`);
            return; // Skip this item
        }
    
    // Check if this is an airblast wind danger zone entry for the Wind Tab
    if (tabId === "windEffects" && item.key.includes("EF") && item.value.includes("-") && item.value.includes("km")) {
        tabSpecificDangerZones.push(item);
    }
    // Check for Tsunami danger zones
    else if (tabId === "tsunamiEffects" && item.key.toLowerCase().includes("zone") && item.value.includes("km")) {
        tabSpecificDangerZones.push(item);
    }
    else if (item.value.includes("-") && item.value.includes("km") && !item.key.toLowerCase().includes("vulnerability zone")) { // General danger zone format
        tabSpecificDangerZones.push(item);
    } else { // Effect at a specific location
        locationEffects.push(item);
    }
  });

  // First display location-specific effects (applies to all tabs like Overpressure, Thermal, etc.)
  if (locationEffects.length > 0) {
    container.innerHTML += `
      <div class="col-12">
        <h4 class="mt-3 mb-3">${t("results.effectsAtLocation", "Effects at Your Selected Location")}</h4>
      </div>
    `;
    
    const categories = new Map();
    locationEffects.forEach(item => {
      const category = item.key.includes(":") ? 
                     item.key.split(":")[0].trim() : 
                     "General";
      if (!categories.has(category)) {
        categories.set(category, []);
      }
      categories.get(category).push(item);
    });

    categories.forEach((itemsInCategory, category) => {
      if (itemsInCategory.length === 0 || (category === "General" && itemsInCategory.length === 1 && itemsInCategory[0].value.trim().toLowerCase() === "none")) return;
      
      let categoryHasContent = false;
      itemsInCategory.forEach(item => {
        if(item.value.trim().toLowerCase() !== "none") categoryHasContent = true;
      });
      if (!categoryHasContent) return;

      if (category !== "General") {
        container.innerHTML += `
            <div class="col-12">
                <h5 class="mt-2 mb-2">${category}</h5>
            </div>
        `;
      }
      itemsInCategory.forEach(item => {
        if (item.value.trim().toLowerCase() === "none") return; // Skip rendering "None" values here too
        let title = item.key.includes(":") ? 
                    item.key.split(":")[1].trim() : 
                    item.key;

        // Filter out "Sound Intensity" but keep "Sound Pressure Level"
        if (title.toLowerCase() === "sound intensity") {
            return; // Skip this item
        }
        
        // Specifically change "Peak wind velocity" to "Wind velocity" for the wind tab
        if (tabId === "windEffects" && title === "Peak wind velocity") {
          title = t("results.windVelocity", "Wind velocity");
        }
        container.innerHTML += createEffectCard(title, item.value);
      });
    });
  }

  let dangerZoneTitleText = "";
  let dangerZoneDescriptionText = "";

  if (tabId === "overpressureEffects" && tabSpecificDangerZones.length > 0) {
    dangerZoneTitleText = t('effects.overpressureDangerZones') || "Overpressure Danger Zones";
    dangerZoneDescriptionText = t('effects.overpressureDescription') || "Areas affected by different overpressure levels";
  } else if (tabId === "windEffects" && tabSpecificDangerZones.length > 0) {
    dangerZoneTitleText = t('effects.windDangerZones') || "Wind Danger Zones (EF Scale)";
    dangerZoneDescriptionText = t('effects.windDescription') || "Areas affected by different wind speed categories";
  } else if (tabId === "thermalEffects" && tabSpecificDangerZones.length > 0) {
    dangerZoneTitleText = t('effects.thermalDangerZones') || "Thermal Radiation Danger Zones";
    dangerZoneDescriptionText = t('effects.thermalDescription') || "Areas affected by different thermal exposure levels";
  } else if (tabId === "seismicEffects" && tabSpecificDangerZones.length > 0) {
    dangerZoneTitleText = t('effects.seismicHazardZones') || "Seismic Hazard Zones";
    dangerZoneDescriptionText = t('effects.seismicDescription') || "Areas affected by different seismic intensity levels";
  } else if (tabId === "ejectaEffects" && tabSpecificDangerZones.length > 0) {
    dangerZoneTitleText = t('effects.ejectaHazardZones') || "Ejecta Hazard Zones";
    dangerZoneDescriptionText = t('effects.ejectaDescription') || "Areas affected by different ejecta deposit thicknesses";
  } else if (tabId === "tsunamiEffects" && tabSpecificDangerZones.length > 0) {
    dangerZoneTitleText = t('effects.tsunamiHazardZones') || "Tsunami Hazard Zones";
    dangerZoneDescriptionText = t('effects.tsunamiDescription') || "Areas affected by different tsunami wave amplitudes";
  }


  if (dangerZoneTitleText && tabSpecificDangerZones.length > 0) {
    container.innerHTML += `
      <div class="col-12">
        <h4 class="mt-4 pt-3 border-top">${dangerZoneTitleText}</h4>
        <p class="text-muted small mb-3 danger-zone-description">${dangerZoneDescriptionText}</p>
      </div>
    `;
    
    tabSpecificDangerZones.forEach(item => {
      const title = item.key.includes(":") ? 
                  item.key.split(":")[1].trim() : 
                  item.key;
      container.innerHTML += `
        <div class="col-12 mb-2">
          <div class="effect-card danger-zone-card">
            <div class="effect-label">${title}</div>
            <div class="dynamic-value">${item.value}</div>
          </div>
        </div>
      `;
    });
  }

  if (container.innerHTML.trim() === '') {
    container.innerHTML = `
      <div class="col-12">
          <div class="effect-card text-center">
              <p class="mb-0">${t('messages.noDataAvailable') || "No data available for this effect type or at the selected distance."}</p>
          </div>
      </div>
    `;
  }
}

// Add this at the bottom of the file
function handleResponsiveLayout() {
  // Check if we're on mobile
  const isMobile = window.innerWidth < 768;
  const isTablet = window.innerWidth >= 768 && window.innerWidth < 992;
  const isSmallScreen = window.innerWidth < 375;
  const isLandscape = window.innerHeight < 500;
  
  // Adjust map height based on screen size
  const mapElement = document.getElementById('map');
  if (mapElement) {
    if (isLandscape) {
      mapElement.style.height = '400px';
    } else if (isSmallScreen) {
      mapElement.style.height = '300px';
    } else if (isMobile) {
      mapElement.style.height = '400px';
    } else if (isTablet) {
      mapElement.style.height = '500px';
    } else {
      mapElement.style.height = '650px'; // Increased from 400px for desktop
    }
  }
  
  // Ensure the map is properly sized
  if (map) {
    map.invalidateSize();
  }
}

// Call on page load and resize
window.addEventListener('load', handleResponsiveLayout);
window.addEventListener('resize', handleResponsiveLayout);

document.getElementById("simulationForm").addEventListener("submit", handleSubmit);

document.querySelectorAll('input[name="zonetype"]').forEach(radio => {
  radio.addEventListener('change', (e) => {
    const lastResults = document.querySelector('#resultsContainer').dataset.lastResults;
    if (lastResults) {
      try {
        const data = JSON.parse(lastResults);
        const mainZoneType = e.target.value;
        const vulnerabilitySelectorDiv = document.getElementById('vulnerabilitySelectorDiv');
        const vulnerabilityTypeSelect = document.getElementById('vulnerabilityTypeSelect');

        if (mainZoneType === 'vulnerability') {
          vulnerabilitySelectorDiv.style.display = 'block';
          const specificVulnerabilityType = vulnerabilityTypeSelect.value;
          if (data.visualization && data.visualization[specificVulnerabilityType]) {
            showZones(specificVulnerabilityType, data.visualization[specificVulnerabilityType]);
          } else {
            console.warn(`Visualization data for ${specificVulnerabilityType} not found.`);
            clearAllZones();
          }
        } else {
          vulnerabilitySelectorDiv.style.display = 'none';
          if (mainZoneType === 'none') {
            showZones('none', null);
          } else if (mainZoneType === 'crater') { // Specifically handle crater
            showZones('crater', null); // Crater data is fetched internally by showZones
          } else if (data.visualization && data.visualization[mainZoneType]) {
            showZones(mainZoneType, data.visualization[mainZoneType]);
          } else {
            console.warn(`Visualization data for ${mainZoneType} not found.`);
            clearAllZones();
          }
        }
      } catch (error) {
        console.error('Error handling zone type change:', error);
      }
    }
  });
});

// Add event listener for the vulnerability type dropdown
document.addEventListener('DOMContentLoaded', function() {
  const vulnerabilityTypeSelect = document.getElementById('vulnerabilityTypeSelect');
  if (vulnerabilityTypeSelect) {
    vulnerabilityTypeSelect.addEventListener('change', (e) => {
      const lastResults = document.querySelector('#resultsContainer').dataset.lastResults;
      const vulnerabilityRadio = document.querySelector('input[name="zonetype"][value="vulnerability"]');

      if (lastResults && vulnerabilityRadio && vulnerabilityRadio.checked) {
        try {
          const data = JSON.parse(lastResults);
          const specificVulnerabilityType = e.target.value;
          if (data.visualization && data.visualization[specificVulnerabilityType]) {
            showZones(specificVulnerabilityType, data.visualization[specificVulnerabilityType]);
          } else {
            console.warn(`Visualization data for ${specificVulnerabilityType} not found.`);
            clearAllZones();
          }
        } catch (error) {
          console.error('Error handling specific vulnerability type change:', error);
        }
      }
    });
  }
});
HEAD


// Custom Alert Function
function showCustomAlert(message) {
  // Remove any existing alert
  const existingAlert = document.querySelector('.custom-alert-overlay');
  if (existingAlert) {
    existingAlert.remove();
  }

  // Create alert overlay
  const overlay = document.createElement('div');
  overlay.className = 'custom-alert-overlay';
  
  // Create alert box
  const alertBox = document.createElement('div');
  alertBox.className = 'custom-alert-box';
  
  // Create icon
  const icon = document.createElement('div');
  icon.className = 'custom-alert-icon';
  icon.innerHTML = '<i class="mdi mdi-map-marker-alert"></i>';
  
  // Create message
  const messageDiv = document.createElement('div');
  messageDiv.className = 'custom-alert-message';
  messageDiv.textContent = message;
  
  // Create button
  const button = document.createElement('button');
  button.className = 'custom-alert-button';
  button.textContent = 'OK';
  button.onclick = () => {
    overlay.classList.add('fade-out');
    setTimeout(() => overlay.remove(), 300);
  };
  
  // Assemble alert
  alertBox.appendChild(icon);
  alertBox.appendChild(messageDiv);
  alertBox.appendChild(button);
  overlay.appendChild(alertBox);
  
  // Add to body (bypasses app-wrapper scaling on Windows)
  document.body.appendChild(overlay);
  
  // Focus management for accessibility
  button.focus();
  
  // Handle ESC key to close
  const handleEscape = (e) => {
    if (e.key === 'Escape') {
      overlay.classList.add('fade-out');
      setTimeout(() => overlay.remove(), 300);
      document.removeEventListener('keydown', handleEscape);
    }
  };
  document.addEventListener('keydown', handleEscape);
  
  // Trigger animation
  setTimeout(() => overlay.classList.add('show'), 10);
  
  // Auto-dismiss after 5 seconds
  setTimeout(() => {
    if (document.body.contains(overlay)) {
      overlay.classList.add('fade-out');
      setTimeout(() => {
        overlay.remove();
        document.removeEventListener('keydown', handleEscape);
      }, 300);
    }
  }, 5000);
}

// NASA Data Success Notification
function showNasaNotification(neoName) {
  // Remove any existing notification
  const existingNotification = document.querySelector('.nasa-notification-overlay');
  if (existingNotification) {
    existingNotification.remove();
  }

  // Create notification overlay
  const overlay = document.createElement('div');
  overlay.className = 'nasa-notification-overlay';
  
  // Create notification box
  const notificationBox = document.createElement('div');
  notificationBox.className = 'nasa-notification-box';
  
  // Create NASA logo/icon
  const icon = document.createElement('div');
  icon.className = 'nasa-notification-icon';
  icon.innerHTML = '<i class="mdi mdi-satellite-variant"></i>';
  
  // Create success message
  const messageDiv = document.createElement('div');
  messageDiv.className = 'nasa-notification-message';
  
  const label = document.createElement('div');
  label.className = 'nasa-notification-label';
  label.textContent = t('messages.loadedNeoData', 'Loaded data for NEO:');
  
  const neoNameDiv = document.createElement('div');
  neoNameDiv.className = 'nasa-notification-neo-name';
  neoNameDiv.textContent = neoName;
  
  messageDiv.appendChild(label);
  messageDiv.appendChild(neoNameDiv);
  
  // Create close button
  const closeBtn = document.createElement('button');
  closeBtn.className = 'nasa-notification-close';
  closeBtn.innerHTML = '<i class="mdi mdi-close"></i>';
  closeBtn.onclick = () => {
    overlay.classList.add('fade-out');
    setTimeout(() => overlay.remove(), 300);
  };
  
  // Assemble notification
  notificationBox.appendChild(icon);
  notificationBox.appendChild(messageDiv);
  notificationBox.appendChild(closeBtn);
  overlay.appendChild(notificationBox);
  
  // Add to body (bypasses app-wrapper scaling on Windows)
  document.body.appendChild(overlay);
  
  // Focus management for accessibility
  closeBtn.focus();
  
  // Handle ESC key to close
  const handleEscape = (e) => {
    if (e.key === 'Escape') {
      overlay.classList.add('fade-out');
      setTimeout(() => overlay.remove(), 300);
      document.removeEventListener('keydown', handleEscape);
    }
  };
  document.addEventListener('keydown', handleEscape);
  
  // Trigger animation
  setTimeout(() => overlay.classList.add('show'), 10);
  
  // Auto-dismiss after 6 seconds
  setTimeout(() => {
    if (document.body.contains(overlay)) {
      overlay.classList.add('fade-out');
      setTimeout(() => {
        overlay.remove();
        document.removeEventListener('keydown', handleEscape);
      }, 300);
    }
  }, 6000);
}
