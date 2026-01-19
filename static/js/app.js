/**
 * ENEO Asteroid Impact Simulation - Application Controller
 * Entry point for the application, handling form submissions and API interactions.
 */



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
    
    showLoadingScreen(loadingOverlay);
    
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

    // ---- START: Impact Overview Population ----
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

  } catch (error) {
    console.error('Error in simulation:', error);
    alert(`Error: ${error.message}`);
  } finally {
    loadingOverlay.style.display = "none";
    resultsContainer.classList.remove("d-none");
  }
}

// Event Listeners
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
