/**
 * ENEO Asteroid Impact Simulation - Loading & Facts
 * Handles the loading screen animations and educational facts.
 */

function showLoadingScreen(loadingOverlay) {
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
}
