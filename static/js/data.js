/**
 * ENEO Asteroid Impact Simulation - Data & API
 * Handles data generation and fetching from external APIs (NASA Sentry).
 */

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
