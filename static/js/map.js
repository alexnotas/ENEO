/**
 * ENEO Asteroid Impact Simulation - Map Controller
 * Handles Leaflet map initialization, layers, and zone visualization.
 */

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
  if (typeof map !== 'undefined' && map) {
    map.invalidateSize();
  }
}

// Call on page load and resize
window.addEventListener('load', handleResponsiveLayout);
window.addEventListener('resize', handleResponsiveLayout);
