/**
 * ENEO Asteroid Impact Simulation - Results Display
 * Handles parsing simulation results and generating HTML content for the UI.
 */

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
