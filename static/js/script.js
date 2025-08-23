// Detect Windows platform and apply zoom scaling
document.addEventListener('DOMContentLoaded', function() {
  // Create a wrapper element for the entire application content.
  // This is a workaround to apply scaling transforms for display correction on certain platforms (like Windows)
  // without affecting elements that should remain fixed, like modals or overlays.
  const wrapper = document.createElement('div');
  wrapper.id = 'app-wrapper';
  // Move all existing body content into the new wrapper.
  while (document.body.firstChild) {
    wrapper.appendChild(document.body.firstChild);
  }
  document.body.appendChild(wrapper);

  // --- Platform-Specific Scaling for Windows ---
  // Check if the user agent string indicates a Windows operating system.
  const isWindows = navigator.userAgent.indexOf('Windows') !== -1;
  
  // Apply scaling immediately on load if the platform is Windows.
  applyAppropriateScaling(isWindows, wrapper);
  
  // Add a resize event listener specifically for Windows.
  // This handles cases where the user might move the browser window between monitors
  // with different scaling factors, requiring a re-evaluation of the transform.
  if (isWindows) {
    let lastWidth = window.innerWidth;
    let resizeTimeout;
    
    window.addEventListener('resize', function() {
      // Use a timeout to debounce the resize event, preventing the handler
      // from firing too frequently while the user is resizing the window.
      clearTimeout(resizeTimeout);
      
      resizeTimeout = setTimeout(function() {
        // Check if the window width has changed significantly. A large change
        // is more likely to be a monitor switch than a simple window resize.
        const widthChange = Math.abs(window.innerWidth - lastWidth);
        if (widthChange > 200) { // Threshold for a significant change.
          console.log('Significant width change detected, reapplying scaling');
          applyAppropriateScaling(isWindows, wrapper);
          lastWidth = window.innerWidth; // Update the last known width.
        }
      }, 300); // Wait 300ms after the last resize event to execute.
    });
  }
  
  // --- UI Element Positioning Fixes ---
  // Ensure the loading overlay is a direct child of the body and positioned fixed.
  // This prevents it from being affected by the scaling applied to the app-wrapper.
  const loadingOverlay = document.querySelector('.loading-overlay');
  if (loadingOverlay) {
    document.body.appendChild(loadingOverlay); // Move to be a direct child of <body>.
    
    // Apply styles to make it a full-screen, fixed overlay.
    loadingOverlay.style.position = 'fixed';
    loadingOverlay.style.top = '0';
    loadingOverlay.style.left = '0';
    loadingOverlay.style.right = '0';
    loadingOverlay.style.bottom = '0';
    loadingOverlay.style.zIndex = '9999'; // Ensure it's on top of all other content.
  }

  // Move the "About" modal to be a direct child of the body.
  // This also prevents it from being scaled down along with the main app content on Windows.
  const aboutModal = document.getElementById('aboutModal');
  if (aboutModal) {
    document.body.appendChild(aboutModal);
  }
});

/**
 * Applies or removes CSS scaling to the main app wrapper based on window width.
 * This is primarily intended to correct display issues on Windows systems with
 * certain screen resolutions and scaling settings.
 * @param {boolean} isWindows - Flag indicating if the current OS is Windows.
 * @param {HTMLElement} wrapper - The main application wrapper element.
 */
function applyAppropriateScaling(isWindows, wrapper) {
  if (!isWindows) return; // Only apply scaling on Windows.
  
  const currentWidth = window.innerWidth;
  // Define the resolution range where scaling is needed.
  const needsScaling = currentWidth >= 992 && currentWidth < 1920;
  
  // Reset any existing scaling styles first to ensure a clean state.
  wrapper.style.transform = '';
  wrapper.style.transformOrigin = '';
  wrapper.style.width = '';
  wrapper.style.marginLeft = '';
  wrapper.style.overflowX = '';
  wrapper.style.overflowY = '';
  wrapper.style.height = '';
  
  if (needsScaling) {
    // Add a class to the root element for platform-specific CSS rules.
    document.documentElement.classList.add('windows-platform');
    
    // Apply a 0.9 scale transform to the wrapper.
    wrapper.style.transform = 'scale(0.9)';
    wrapper.style.transformOrigin = 'top center'; // Scale from the top center.
    // Compensate for the scaling by increasing the wrapper's width. 100% / 0.9 = 111.11%
    wrapper.style.width = '111.11%'; 
    // Adjust the margin to re-center the scaled content.
    wrapper.style.marginLeft = '-5.55%'; 
    wrapper.style.overflowX = 'hidden'; // Prevent horizontal scrollbar caused by scaling.
    wrapper.style.overflowY = 'auto';   // Allow vertical scrolling within the wrapper.
    wrapper.style.height = 'auto';
    
    // Adjust body height to ensure the viewport is scrollable.
    document.body.style.height = '90vh';
    document.body.style.overflowY = 'auto';
    
    // Force a layout refresh by dispatching a resize event.
    // This helps ensure all elements correctly adjust to the new scaling.
    setTimeout(function() {
      window.dispatchEvent(new Event('resize'));
    }, 100);
  } else {
    // If scaling is not needed, remove the platform-specific class.
    document.documentElement.classList.remove('windows-platform');
    
    // Reset body styles to their default state.
    document.body.style.height = '';
    document.body.style.overflowY = '';
  }
}

// --- Leaflet Map Initialization ---

// Initialize the Leaflet map within the 'map' div.
const map = L.map('map', {
  minZoom: 2,       // Minimum zoom level (zoomed out).
  maxZoom: 18,      // Maximum zoom level (zoomed in).
  maxBounds: [      // Restrict map panning to the geographical bounds of the world.
    [-90, -180],
    [90, 180]
  ],
  maxBoundsViscosity: 1.0 // Makes the bounds completely solid.
}).setView([0, 0], 2); // Set initial view to the center of the world, zoomed out.

// Create custom panes for rendering crater visuals.
// Panes allow controlling the stacking order of layers on the map.
// Pane for the crater image (non-interactive).
map.createPane('craterImagePane');
map.getPane('craterImagePane').style.zIndex = 450; // Below the border pane.

// Pane for the crater border (interactive, e.g., for popups).
map.createPane('craterBorderPane');
map.getPane('craterBorderPane').style.zIndex = 451; // Above the image pane.

// Use the ESRI World Street Map as the base tile layer.
L.tileLayer('https://server.arcgisonline.com/ArcGIS/rest/services/World_Street_Map/MapServer/tile/{z}/{y}/{x}', {
  attribution: 'Tiles &copy; Esri', // Copyright attribution for the map tiles.
  noWrap: true, // Prevents tiles from repeating horizontally.
  bounds: [ // Ensures tiles are only loaded within the world bounds.
    [-90, -180],
    [90, 180]
  ]
}).addTo(map);

// --- Map Interaction ---

let marker; // Variable to hold the user-placed marker.
// Event listener for clicks on the map.
map.on('click', function(e) {
  // Validate that the click is within the defined map boundaries.
  if (map.options.maxBounds && !map.options.maxBounds.contains(e.latlng)) {
    alert("Please select a valid point within the map boundaries.");
    return; // Stop if the click is outside the valid area.
  }

  // If a marker already exists, remove it before placing a new one.
  if (marker) {
    map.removeLayer(marker);
  }
  // Create a new marker at the clicked location.
  marker = L.marker(e.latlng, {
    zIndexOffset: 1000  // Ensure the marker is always on top of other layers.
  }).addTo(map);
  // Update the hidden form fields with the selected latitude and longitude.
  document.getElementById('latitude').value = e.latlng.lat;
  document.getElementById('longitude').value = e.latlng.lng;

  // Asynchronously fetch the ocean depth for the selected coordinates.
  fetch('/get_ocean_depth', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ latitude: e.latlng.lat, longitude: e.latlng.lng })
  })
  .then(response => {
    if (!response.ok) {
      // If the server returns an error, try to parse the JSON error message.
      return response.json().then(errData => {
        throw new Error(errData.error || `HTTP error ${response.status}`);
      }).catch(() => {
        // Fallback if the error response is not JSON.
        throw new Error(`HTTP error ${response.status}`);
      });
    }
    return response.json();
  })
  .then(data => {
    // The ocean depth is fetched but not displayed in a popup.
    // It's logged to the console for debugging or informational purposes.
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
    // Log any errors during the fetch process to the console.
    console.error('Error fetching ocean depth:', error);
  });
});

// --- Zone Visualization ---

// An object to hold Leaflet LayerGroups for different types of hazard zones.
// This allows toggling them on and off as a group.
let zoneLayers = {
  seismic: new L.LayerGroup(),
  thermal: new L.LayerGroup(),
  airblast: new L.LayerGroup(),
  ejecta: new L.LayerGroup(),
  vulnerability: new L.LayerGroup(),
  crater: new L.LayerGroup(),
  wind: new L.LayerGroup(),
  tsunami: new L.LayerGroup()
};

/**
 * Returns a color for a given hazard zone type and value.
 * @param {string} type - The type of zone (e.g., 'seismic', 'thermal').
 * @param {number} value - A value (e.g., vulnerability threshold) to determine color shade.
 * @returns {string} An RGBA color string.
 */
function getZoneColor(type, value) {
  const opacity = 0.6; // Base opacity for the zones.
  switch (type) {
    case 'vulnerability':
      if (value >= 0.75) return `rgba(255, 0, 0, ${opacity + 0.1})`;       // Darker Red for highest risk
      if (value >= 0.5) return `rgba(255, 120, 0, ${opacity})`;     // Orange for high risk
      if (value >= 0.25) return `rgba(255, 220, 0, ${opacity - 0.1})`;    // Yellow for medium risk
      return `rgba(0, 255, 0, ${opacity - 0.2})`; // Lighter Green for lowest risk
    case 'seismic':
      return `rgba(255, 0, 0, ${opacity})`; // Red
    case 'thermal':
      return `rgba(255, 102, 0, ${opacity})`; // Orange
    case 'airblast': // Overpressure
      return `rgba(0, 200, 200, ${opacity})`; // Teal/cyan
    case 'wind':
      return `rgba(0, 102, 255, ${opacity})`; // Blue
    case 'ejecta':
      return `rgba(102, 51, 0, ${opacity})`; // Brown
    case 'crater':
      return `rgba(128, 128, 128, ${opacity})`; // Grey
    case 'tsunami':
      return `rgba(0, 100, 255, ${opacity})`; // Blue
    default:
      return `rgba(102, 102, 102, ${opacity})`; // Default gray for unknown types.
  }
}

/**
 * Removes all hazard zone layers from the map.
 */
function clearAllZones() {
  Object.values(zoneLayers).forEach(layer => {
    if (map.hasLayer(layer)) {
      map.removeLayer(layer);
    }
    layer.clearLayers(); // Also clear the layers within the LayerGroup.
  });
}

/**
 * Displays hazard zones of a specific type on the map.
 * @param {string} type - The type of zone to display (e.g., 'seismic', 'crater').
 * @param {Array} zoneSpecificData - The data array for the specified zone type.
 */
function showZones(type, zoneSpecificData) {
  clearAllZones(); // Start by clearing any existing zones.

  // Special handling for 'tsunami': do not show any zones, just fly to the marker.
  if (type === 'tsunami') {
    if (marker) {
        map.flyTo(marker.getLatLng(), map.getZoom() < 5 ? 5 : map.getZoom());
    }
    return;
  }

  // If 'none' is selected, just center the map on the marker.
  if (type === 'none') {
    if (marker) { 
        map.flyTo(marker.getLatLng(), map.getZoom() < 5 ? 5 : map.getZoom());
    }
    return;
  }
  
  // Determine the correct Leaflet layer group. Vulnerability subtypes all use the same layer.
  const actualZoneLayer = type.startsWith('vulnerability_') ? zoneLayers.vulnerability : zoneLayers[type];
  
  if (!actualZoneLayer) {
      console.warn(`No zone layer configured for type: ${type}`);
      if (marker) map.flyTo(marker.getLatLng(), map.getZoom());
      return;
  }

  // --- Crater Visualization ---
  // The crater is handled separately as it's a single circle, not a series of zones.
  if (type === 'crater') {
    const fullResultsData = JSON.parse(document.querySelector('#resultsContainer').dataset.lastResults || '{}');
    if (fullResultsData && fullResultsData.crater && 
        fullResultsData.crater.crater_formation && 
        fullResultsData.crater.crater_formation.final_diameter > 0 && marker) {
      
      const craterDiameterMeters = fullResultsData.crater.crater_formation.final_diameter;
      const impactLatLng = marker.getLatLng();
      const craterRadiusMeters = craterDiameterMeters / 2;
      
      zoneLayers.crater.clearLayers();

      // Create a circle representing the crater.
      const craterCircle = L.circle(impactLatLng, {
        radius: craterRadiusMeters,
        color: '#654321',
        fillColor: '#654321',
        fillOpacity: 1.0,
        weight: 1.5,
        pane: 'craterBorderPane' // Render in the custom pane.
      });
      
      // Bind a popup with crater information.
      craterCircle.bindPopup(`
        <strong>Crater</strong><br>
        Diameter: ${(craterDiameterMeters / 1000).toFixed(2)} km
      `);
      
      zoneLayers.crater.addLayer(craterCircle);
      map.addLayer(zoneLayers.crater);
      
      // Zoom the map to fit the crater bounds.
      const bounds = craterCircle.getBounds();
      map.flyToBounds(bounds, { padding: [30, 30], maxZoom: Math.min(map.getMaxZoom(), 17) });

    } else {
      // If no crater data, just fly to the marker.
      if (marker) {
          map.flyTo(marker.getLatLng(), map.getZoom());
      }
    }
    return; // End execution for crater type.
  }
  
  // --- General Zone Visualization ---
  if (!zoneSpecificData || (Array.isArray(zoneSpecificData) && zoneSpecificData.length === 0)) {
      if (marker) {
          map.flyTo(marker.getLatLng(), map.getZoom());
      }
      return;
  }

  // Filter and sort the zones. For seismic, only show zones for Richter 4+.
  // Sort by distance descending to draw larger circles first.
  const initiallyFilteredAndSortedZones = [...zoneSpecificData]
    .filter(zone => {
      if (type === 'seismic') {
        const description = zone.description || '';
        const match = description.match(/(\d+)(?:\+|-\d+)\s*Richter/);
        return match && parseInt(match[1]) >= 4;
      }
      return true;
    })
    .sort((a, b) => b.end_distance - a.end_distance);
  
  if (initiallyFilteredAndSortedZones.length === 0) {
    if (marker) map.flyTo(marker.getLatLng(), map.getZoom());
    return;
  }

  // Further filter out zones that are not visually meaningful (e.g., tiny vulnerability zones).
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

  // Define constants for visualization limits.
  const EARTH_RADIUS_KM = 10000;
  const EARTH_DIAMETER_KM = EARTH_RADIUS_KM * 2;
  const maxVisualizationDistance = EARTH_DIAMETER_KM;
  
  // Get the maximum distance from the sorted zones to use for map bounds calculation.
  const maxDistance = Math.min(visualizableZones[0].end_distance, maxVisualizationDistance);
  
  // Iterate over the visualizable zones and create circles for each.
  visualizableZones.forEach(zone => {
    const radiusInMeters = Math.min(zone.end_distance, EARTH_DIAMETER_KM) * 1000;
      
    if (marker) {
      const impactLocation = marker.getLatLng();
      // Check if the circle will cross the International Date Line.
      const crossesDateLine = Math.abs(impactLocation.lng) + (radiusInMeters / 111319.9) > 180;
        
      let popupContent = '';
      let zoneColorConfig;

      // Create popup content and determine color based on zone type.
      if (type.startsWith('vulnerability_')) {
        let specificVulnerabilityName = "Vulnerability"; // Default name
        if (type === 'vulnerability_thermal') specificVulnerabilityName = "Thermal-Induced Vulnerability";
        else if (type === 'vulnerability_overpressure') specificVulnerabilityName = "Overpressure-Induced Vulnerability";
        else if (type === 'vulnerability_seismic') specificVulnerabilityName = "Seismic-Induced Vulnerability";
        else if (type === 'vulnerability_ejecta') specificVulnerabilityName = "Ejecta-Induced Vulnerability";
        else if (type === 'vulnerability_wind') specificVulnerabilityName = "Wind-Induced Vulnerability";
        else if (type === 'vulnerability_combined') specificVulnerabilityName = "Combined Vulnerability";

        popupContent = `
          <strong>${specificVulnerabilityName} Zone</strong><br>
          Risk Level: ${(zone.threshold * 100).toFixed(1)}%<br>
          Range: ${zone.start_distance.toFixed(2)} - ${zone.end_distance.toFixed(2)} km
        `;
        const color = getZoneColor('vulnerability', zone.threshold || 0.5);
        zoneColorConfig = { outline: color, fill: color };
      } else {
        popupContent = `
          <strong>${zone.description}</strong><br>
          Range: ${zone.start_distance.toFixed(2)} - ${zone.end_distance.toFixed(2)} km
        `;
        const fillColor = getZoneColor(type, zone.threshold || 0.5);
        const outlineColor = type === 'thermal' ? '#FF4500' : fillColor; // Special outline for thermal.
        zoneColorConfig = { outline: outlineColor, fill: fillColor };
      }
        
      // Create the primary circle for the zone.
      const circle = L.circle(impactLocation, {
        radius: radiusInMeters,
        color: zoneColorConfig.outline,
        fillColor: zoneColorConfig.fill,
        fillOpacity: 0.4,
        weight: 1.5
      });
        
      circle.bindPopup(popupContent);
      actualZoneLayer.addLayer(circle);
        
      // If the circle crosses the date line, create a "wrapped" circle on the other side of the map.
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

  // Add the completed layer group to the map.
  map.addLayer(actualZoneLayer);

  // If the effect radius is larger than the Earth's radius, just zoom out to a global view.
  if (maxDistance > EARTH_RADIUS_KM) {
    map.setView([0, 0], 1);
    return;
  }

  // --- Adjust Map View to Fit Zones ---
  if (marker) {
    try {
      const center = marker.getLatLng();
      // Approximate conversion from km to degrees, accounting for latitude.
      const kmPerLongitudeDegree = 111.32 * Math.cos(center.lat * Math.PI / 180);
      const kmPerLatitudeDegree = 111.32;
      const adjustedMaxDistance = maxDistance * 1.3 // Add a 30% buffer.
      const latDegrees = adjustedMaxDistance / kmPerLatitudeDegree;
      // Handle longitude degrees carefully near the poles.
      const lonDegrees = Math.abs(center.lat) > 85 ? 360 : adjustedMaxDistance / kmPerLongitudeDegree;
      
      // Create a bounding box that encompasses the largest zone.
      const bounds = L.latLngBounds(
        [Math.max(-90, center.lat - latDegrees), center.lng - lonDegrees],
        [Math.min(90, center.lat + latDegrees), center.lng + lonDegrees]
      );

      // Use flyToBounds for a smooth animated transition to the new view.
      const paddingPixels = Math.max(100, adjustedMaxDistance / 2);
      map.flyToBounds(bounds, {
        padding: [paddingPixels, paddingPixels],
        maxZoom: 10,
        duration: 0.75,
        animate: true,
        easeLinearity: 0.3
      });
      marker.bringToFront(); // Ensure marker is visible on top of zones.
    } catch (e) {
      // Fallback to a simpler flyTo if bounds calculation fails.
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

/**
 * Populates the input form with random (but plausible) sample data.
 * Useful for demonstrating the application without manual input.
 */
function loadSampleData() {
  // Generate random values within logical ranges for each parameter.
  const diameter = Math.floor(Math.random() * 1451) + 50;    // 50-1500 m
  const density = Math.floor(Math.random() * 1001) + 2500;    // 2500-3500 kg/m^3
  const velocity = Math.floor(Math.random() * 16) + 17;       // 17-32 km/s
  const entry_angle = Math.floor(Math.random() * 31) + 35;    // 35-65 degrees
  const distance = Math.floor(Math.random() * 91) + 10;       // 10-100 km
  
  // Set the form input values.
  document.getElementById("diameter").value = diameter;
  document.getElementById("density").value = density;
  document.getElementById("velocity").value = velocity;
  document.getElementById("entry_angle").value = entry_angle;
  document.getElementById("distance").value = distance;
}

/**
 * Parses a custom-formatted text block into a structured JavaScript object.
 * The text is expected to have sections delineated by "=== Section Name ===".
 * @param {string} text - The raw text to parse.
 * @returns {Object} A structured object representing the parsed sections.
 */
function parseResultsText(text) {
  const sections = {};
  // Split the text by the section headers.
  const splitSections = text.split(/===+ (.*?) ===+/g).slice(1);

  for (let i = 0; i < splitSections.length; i += 2) {
    const sectionName = splitSections[i].trim();
    // Process the content of each section.
    const content = splitSections[i + 1]
      .trim()
      .split("\n") // Split by lines.
      .filter((line) => line.trim()) // Remove empty lines.
      .map((line) => {
        // Split each line into a key-value pair.
        const [key, ...values] = line.split(":");
        return {
          key: key?.trim() || "",
          value: values.join(":").trim(),
        };
      });
    // Store the content in the sections object with a normalized key.
    sections[sectionName.toLowerCase().replace(/ /g, "_")] = content;
  }
  return sections;
}

/**
 * Creates an HTML card for displaying a specific impact effect.
 * @param {string} title - The main title for the card.
 * @param {string} value - The value of the effect.
 * @param {number} distance - The distance at which the effect is measured.
 * @param {string} [subtitle=""] - An optional subtitle.
 * @returns {string} An HTML string for the effect card.
 */
function createEffectCard(title, value, distance, subtitle = "") {
  // Automatically add the distance to the title for certain effect types.
  if (
    distance !== undefined &&
    ["overpressure", "seismic", "thermal", "ejecta", "tsunami"].some((word) => 
      title.toLowerCase().includes(word)
    ) &&
    !title.toLowerCase().includes(" at ")
  ) {
    title = `${title} at ${distance} km`;
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

/**
 * Creates an HTML card for displaying a danger zone's range.
 * @param {string} label - The label for the danger zone.
 * @param {string} range - The range of the zone (e.g., "10-20 km").
 * @param {number} distance - The target distance for context.
 * @returns {string} An HTML string for the danger zone card.
 */
function createDangerZone(label, range, distance) {
  // Similar to createEffectCard, add distance context to the label.
  if (
    ["overpressure", "seismic", "thermal", "ejecta", "tsunami"].some((word) => 
      label.toLowerCase().includes(word)
    ) &&
    !label.toLowerCase().includes(" at ")
  ) {
    label = `${label} at ${distance} km`;
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

/**
 * Creates a formatted HTML item for the main impact overview section.
 * @param {string} label - The label for the data point (e.g., "Impact Energy").
 * @param {string|number} value - The value to display.
 * @param {string} [unit=''] - The unit for the value (e.g., "km", "MT").
 * @returns {string} An HTML string for the overview item.
 */
function createImpactOverviewItem(label, value, unit = '') {
  let displayValue = "N/A";
  if (value !== null && typeof value !== 'undefined') {
    if (typeof value === 'number') {
      // Use exponential notation for very large or very small numbers for readability.
      if (Math.abs(value) > 1e6 || (Math.abs(value) < 1e-2 && Math.abs(value) !== 0)) {
        displayValue = value.toExponential(2);
      } else if (unit === "km/s" || unit === "km" || unit === "MT") {
        displayValue = value.toFixed(2); // Use 2 decimal places for these units.
      } else if (unit === "m") {
        displayValue = value.toFixed(0); // No decimals for meters.
      } else { // Default numeric formatting.
        displayValue = parseFloat(value.toFixed(3)).toString();
      }
    } else {
      displayValue = value.toString();
    }
    displayValue += ` ${unit}`; // Append the unit.
  }

  return `
    <div class="col-lg-4 col-md-6 mb-3">
      <div class="dynamic-value">${displayValue.trim()}</div>
      <div>${label}</div>
    </div>
  `;
}

/**
 * Handles the form submission for the impact simulation.
 * @param {Event} e - The form submission event.
 */
async function handleSubmit(e) {
  e.preventDefault(); // Prevent the default form submission behavior.
  const loadingOverlay = document.querySelector(".loading-overlay");
  const resultsContainer = document.getElementById("resultsContainer");
  
  try {
    // Show the loading overlay and hide the previous results.
    loadingOverlay.style.display = "flex";
    resultsContainer.classList.add("d-none");
    
    // An array of interesting facts about asteroid impacts to display during loading.
    const impactFacts = [
      "The Chicxulub impact that killed the dinosaurs released 100 trillion tons of TNT equivalent—over a billion times the energy of the Hiroshima bomb.",
      "The Tunguska event in 1908 was caused by a ~50-60 m object that air-burst over Siberia, flattening about 2,150 km² of forest and felling 80 million trees.",
      "A kilometer-scale asteroid impacts Earth only about once every 500,000 years on average.",
      "The Chelyabinsk meteor in 2013 was only ~17-20 m across but injured ~1,500 people and damaged ~7,200 buildings, mostly from broken glass.",
      "South Africa's Vredefort Dome is the largest confirmed terrestrial impact structure, originally ~300 km across, formed 2.023 billion years ago.",
      "NASA's DART mission struck the 170 m-wide moonlet Dimorphos at ~6.6 km/s in 2022, shortening its orbit by 32 minutes—demonstrating kinetic-impactor deflection.",
      "Meteorites fall into three broad types: stony (silicate-rich), iron (metallic Ni-Fe), and stony-iron (mixed), each with numerous subgroups.",
      "An asteroid's kinetic energy scales as ½mv²—doubling its velocity quadruples its energy—making speed as critical as mass for impact hazards.",
      "As of December 2024, over 37,000 Near-Earth Objects have been cataloged, with 2,465 classified as potentially hazardous asteroids.",
      "Although the Main Belt contains millions of asteroids, their average separation is roughly 965,000 km—creating vast voids between objects.",
      "Meteoroids strike Earth at between 11 and 72 km/s; those above 20 km/s ionize the air ahead, creating bright plasma trails.",
      "Objects ≥100 m diameter can release tens to hundreds of megatons of TNT equivalent, causing regional devastation.",
      "Earth accretes ~54 tons of interplanetary material daily, totaling ~19,700 tons per year of dust and micrometeoroids.",
      "Bodies <50 m often air-burst (like Tunguska), delivering more widespread damage than crater-forming ground impacts of similar mass.",
      "The Barringer Crater in Arizona is 1.186 km in diameter and 170 m deep, formed ~50,000 years ago by a ~50 m iron meteorite.",
      "The Hoba meteorite in Namibia is the largest known intact meteorite at over 60 tonnes, composed almost entirely of iron-nickel alloy.",
      "NASA's NEO Surveyor space telescope, planned for launch by 2027, aims to complete the survey of 90% of NEOs ≥140 m.",
      "Of the known near-Earth asteroids, 2,465 are classified as Potentially Hazardous Asteroids—large enough and close enough to pose a significant impact risk.",
      "The Torino scale rates NEO impact threat from 0 (none) to 10 (certain global catastrophe); only asteroid Apophis ever reached level 4.",
      "Approximately 25 million meteoroids enter Earth's atmosphere each day, ranging from dust-sized particles to occasional boulder-sized objects."
    ];
    
    // Shuffle the facts array to present them in a random order.
    const shuffledFacts = [...impactFacts].sort(() => 0.5 - Math.random());
    let currentFactIndex = 0;
    
    // Dynamically create the loading timeline UI inside the overlay.
    loadingOverlay.innerHTML = `
      <div class="impact-timeline">
        <div class="timeline-progress"></div>
        <div class="timeline-step active" data-step="entry">
          <div class="step-icon"><i class="mdi mdi-meteor"></i></div>
          <div class="step-label">Atmospheric Entry</div>
          <div class="step-details" id="entryDetails"></div>
        </div>
        <div class="timeline-step" data-step="impact">
          <div class="step-icon"><i class="mdi mdi-explosion"></i></div>
          <div class="step-label">Impact Effects</div>
          <div class="step-details" id="impactDetails"></div>
        </div>
        <div class="timeline-step" data-step="zones">
          <div class="step-icon"><i class="mdi mdi-map-marker-radius"></i></div>
          <div class="step-label">Vulnerability Zones</div>
          <div class="step-details" id="zonesDetails"></div>
        </div>
        <div class="timeline-step" data-step="population">
          <div class="step-icon"><i class="mdi mdi-account-group"></i></div>
          <div class="step-label">Population Impact</div>
          <div class="step-details" id="populationDetails"></div>
        </div>
      </div>
      <div class="loading-fact" id="factContainer">
        <span class="fact-prefix">Did you know?</span>
        <span class="fact-text">${shuffledFacts[currentFactIndex]}</span>
      </div>
    `;
    
    // Set up an interval to rotate the displayed fact every 12 seconds.
    const factRotationInterval = setInterval(() => {
      currentFactIndex = (currentFactIndex + 1) % shuffledFacts.length;
      const factContainer = document.getElementById('factContainer');
      if (factContainer) {
        factContainer.querySelector('.fact-text').textContent = shuffledFacts[currentFactIndex];
      } else {
        // If the container is gone, stop the interval.
        clearInterval(factRotationInterval);
      }
    }, 12000);
    
    // Animate the timeline steps to give the user a sense of progress.
    animateTimelineStep('entry', 'Calculating entry velocity and atmospheric effects...', 2500);
    
    setTimeout(() => {
      animateTimelineStep('impact', 'Modeling thermal, seismic, and airblast waves...', 3000);
    }, 3000);
    
    setTimeout(() => {
      animateTimelineStep('zones', 'Assessing hazard radii and vulnerability...', 3500);
    }, 6500);
    
    setTimeout(() => {
      animateTimelineStep('population', 'Estimating population exposure...', 4000);
    }, 10000);
    
    // Make the actual API call to the backend for simulation.
    const response = await fetch("/simulate", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        // Read values from the form, parsing them as floats.
        diameter: parseFloat(document.getElementById("diameter").value),
        density: parseFloat(document.getElementById("density").value),
        velocity: parseFloat(document.getElementById("velocity").value),
        entry_angle: parseFloat(document.getElementById("entry_angle").value),
        distance: parseFloat(document.getElementById("distance").value),
        latitude: parseFloat(document.getElementById("latitude").value || 0),
        longitude: parseFloat(document.getElementById("longitude").value || 0)
      }),
    });

    // Check if the server response is not OK (e.g., status 400 or 500).
    if (!response.ok) {
      const errorData = await response.json().catch(() => ({ error: 'An unknown server error occurred.' }));
      throw new Error(errorData.error || `Server responded with status ${response.status}`);
    }

    const data = await response.json();
    
    // Store the full JSON results in a data attribute on the results container
    // for later use (e.g., by the showZones function).
    resultsContainer.dataset.lastResults = JSON.stringify(data.results_data);
    
    // ---- START: New Impact Overview Population Logic ----
    // This section populates the main "Impact Overview" panel with key results.
    const impactOverviewContainer = document.getElementById("impactOverviewContent");
    impactOverviewContainer.innerHTML = ''; // Clear previous overview content.
    const eventType = data.results_data.atmospheric_entry.event_type;
    const atmEntry = data.results_data.atmospheric_entry;
    const energyData = data.results_data.energy;
    const craterData = data.results_data.crater?.crater_formation;

    let overviewHTML = '';

    // Common item: Initial Kinetic Energy
    overviewHTML += createImpactOverviewItem("Initial Kinetic Energy (MT TNT)", energyData.initial_energy_megatons, "MT TNT");

    // Breakup Altitude - show if it's a breakup event (ground impact or airburst)
    if (eventType === "ground impact" || eventType === "airburst") {
        overviewHTML += createImpactOverviewItem("Breakup Altitude", atmEntry.breakup_altitude, "m");
    } else if (eventType === "intact") {
        // For intact objects, breakup altitude is not applicable.
        // overviewHTML += createImpactOverviewItem("Breakup Altitude", "N/A (Object Intact)"); // Or omit
    }

    if (eventType === "airburst") {
        overviewHTML += createImpactOverviewItem("Airburst Altitude", atmEntry.airburst_altitude, "m");
        overviewHTML += createImpactOverviewItem("Post Airburst Velocity", atmEntry.final_velocity ? atmEntry.final_velocity / 1000 : null, "km/s");
        overviewHTML += createImpactOverviewItem("Airburst Energy (MT TNT)", energyData.specific_energy_megatons, "MT TNT");
    } else if (eventType === "ground impact" || eventType === "intact") {
        overviewHTML += createImpactOverviewItem("Impact Speed", atmEntry.final_velocity ? atmEntry.final_velocity / 1000 : null, "km/s");
        overviewHTML += createImpactOverviewItem("Impact Energy (MT TNT)", energyData.specific_energy_megatons, "MT TNT");
        overviewHTML += createImpactOverviewItem("Transient Crater Diameter", craterData?.transient_diameter ? craterData.transient_diameter / 1000 : null, "km");
        overviewHTML += createImpactOverviewItem("Final Crater Diameter", craterData?.final_diameter ? craterData.final_diameter / 1000 : null, "km");
    } else {
        // Fallback for any other unexpected event types
        overviewHTML += `<div class="col-12"><p class="text-white">Overview for event type: "${eventType}"</p></div>`;
        overviewHTML += createImpactOverviewItem("Final Velocity", atmEntry.final_velocity ? atmEntry.final_velocity / 1000 : null, "km/s");
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
    populateTab("thermalEffects", sections.thermal_radiation || []);
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
                                  (${relativeContribution.toFixed(1)}% of total vulnerability)
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
              <div class="effect-label">Vulnerability Zone ${(zone.vulnerability_threshold * 100).toFixed(1)}%</div>
              <div class="population-zone">
                <div>Distance Range: ${zone.start_distance.toFixed(2)} - ${zone.end_distance.toFixed(2)} km</div>
              </div>
            </div>
          </div>
        `)
        .join("");
    }

    // BEGIN: Code from "Εικόνες file" for Affected Countries
    if (data.results_data.population_analysis && data.results_data.population_analysis.countries) {
      const countriesContainer = document.getElementById("countriesResults"); // Ensure this ID exists in your HTML
      const countries = data.results_data.population_analysis.countries;
      
      // Sort countries by casualties (descending)
      countries.sort((a, b) => b.total_casualties - a.total_casualties);
      
      // Get total affected population across all countries (variable totalPop is not used in the original snippet's HTML output but kept for consistency)
      const totalPop = countries.reduce((sum, country) => sum + country.total_population, 0); 
      
      // Create HTML for each country
      if (countriesContainer) {
        countriesContainer.innerHTML = countries
          .map(country => {
            return `
              <div class="col-md-6 mb-3">
                <div class="country-card effect-card">
                  <div class="effect-label">${country.name || ('Country FID ' + country.fid)}</div>
                  <div class="country-stats">
                    <div class="row">
                      <div class="col-12">
                        <div class="stat-label">Total Casualties</div>
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
                <p class="mb-0">No country data available</p>
              </div>
            </div>
          `;
        }
      }
    }
    // END: Code from "Εικόνες file" for Affected Countries

    // BEGIN: Code from "Εικόνες file" for Economic Impact
    if (data.results_data.economic_analysis) {
      const economicContainer = document.getElementById("economicResults"); // Ensure this ID exists
      const totalEconomicDamage = document.getElementById("totalEconomicDamage"); // Ensure this ID exists
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
                  <div class="effect-label">${country.name}</div>
                  <div class="economy-stats">
                    <div class="row">
                      <div class="col-md-6">
                        <div class="stat-label">Economic Damage</div>
                        <div class="stat-value text-danger">${formattedDamage}</div>
                      </div>
                      <div class="col-md-6">
                        <div class="stat-label">GDP Per Capita (${country.gdp_year})</div>
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
                <p class="mb-0">No GDP data available for affected countries</p>
              </div>
            </div>
          `;
        } else {
          economicContainer.innerHTML += `
            <div class="col-12 mt-3">
              <div class="effect-card">
                <div class="effect-label text-white fw-bold">About Economic Impact Calculation</div>
                <p class="mb-2 text-highlight-orange">Economic damage is calculated using the formula:</p>
                <p class="mb-0"><strong class="text-white">Casualties × GDP per capita</strong></p>
                <p class="mb-0 mt-2"><small class="text-muted text-highlight-orange">This represents the direct economic impact based on casualties and national GDP per capita.</small></p>
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
          <h5 class="effect-label">Note on Tsunami Calculations</h5>
          <p class="mb-1 mt-2" style="font-size: 0.9rem; line-height: 1.6;">
            The tsunami effects presented here are based on a <strong>simplified approximation</strong>. This model considers initial wave amplitude at the source and basic geometric spreading. The initial wave amplitude calculation is validated against the conditions described in the tsunami algorithm by Rumpf et al. (2016), serving as a basis for future enhancements. More advanced tsunami models were tested but proved too computationally intensive for real-time online purposes; their integration is a potential future improvement.
          </p>
          <p class="mb-1 mt-1" style="font-size: 0.9rem; line-height: 1.6;">
            Accurate, operational tsunami modeling is significantly more complex and requires:
          </p>
          <ul style="font-size: 0.9rem; line-height: 1.6; margin-top: 0.5rem; margin-bottom: 0.5rem;">
            <li>Detailed bathymetry (sea floor depth, slope, and specific underwater features)</li>
            <li>High-resolution coastal topography</li>
            <li>Consideration of complex wave dynamics (e.g., refraction, diffraction, shoaling, reflection, and run-up)</li>
            <li>Advanced hydrodynamic equations and numerical simulations</li>
          </ul>
          <p class="mb-0 mt-2" style="font-size: 0.9rem; line-height: 1.6; font-weight: bold;">
            Due to these simplifications, the results from this model are for general informational purposes only and are not visualized on the map.
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
                  <h6 class="effect-label">Approximate Wave Amplitude Distances:</h6>
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
             <p class="mb-0">No significant tsunami effects calculated by the simplified model for the current parameters, or data is not applicable.</p>
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
        <h4 class="mt-3 mb-3">Effects at Your Selected Location</h4>
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
          title = "Wind velocity";
        }
        container.innerHTML += createEffectCard(title, item.value);
      });
    });
  }

  let dangerZoneTitleText = "";
  let dangerZoneDescriptionText = "";

  if (tabId === "overpressureEffects" && tabSpecificDangerZones.length > 0) {
    dangerZoneTitleText = "Overpressure Danger Zones";
    dangerZoneDescriptionText = "Areas affected by different overpressure levels";
  } else if (tabId === "windEffects" && tabSpecificDangerZones.length > 0) {
    dangerZoneTitleText = "Wind Danger Zones (EF Scale)";
    dangerZoneDescriptionText = "Areas affected by different wind speed categories";
  } else if (tabId === "thermalEffects" && tabSpecificDangerZones.length > 0) {
    dangerZoneTitleText = "Thermal Radiation Danger Zones";
    dangerZoneDescriptionText = "Areas affected by different thermal exposure levels";
  } else if (tabId === "seismicEffects" && tabSpecificDangerZones.length > 0) {
    dangerZoneTitleText = "Seismic Hazard Zones";
    dangerZoneDescriptionText = "Areas affected by different seismic intensity levels";
  } else if (tabId === "ejectaEffects" && tabSpecificDangerZones.length > 0) {
    dangerZoneTitleText = "Ejecta Hazard Zones";
    dangerZoneDescriptionText = "Areas affected by different ejecta deposit thicknesses";
  } else if (tabId === "tsunamiEffects" && tabSpecificDangerZones.length > 0) {
    dangerZoneTitleText = "Tsunami Hazard Zones";
    dangerZoneDescriptionText = "Areas affected by different tsunami wave amplitudes";
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
              <p class="mb-0">No data available for this effect type or at the selected distance.</p>
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