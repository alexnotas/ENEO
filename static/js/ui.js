/**
 * ENEO Asteroid Impact Simulation - UI Utilities
 * Handles platform scaling, responsive layout, alerts, and animations.
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
  const isMac = /Macintosh|Mac OS X/.test(navigator.userAgent) && !/iPhone|iPad|iPod/.test(navigator.userAgent);
  const isAndroidTablet = /Android/i.test(navigator.userAgent) && !/Mobile/i.test(navigator.userAgent);

  // Add a dedicated class for macOS-specific map click alignment fixes
  if (isMac) {
    document.documentElement.classList.add('mac-platform');

    // Ensure Leaflet recalculates map dimensions after platform-specific styles are applied
    requestAnimationFrame(() => {
      if (typeof map !== 'undefined' && map) {
        setTimeout(() => map.invalidateSize(true), 0);
      }
    });
  }

  // Add a dedicated class for Android tablet-specific map click alignment fixes
  if (isAndroidTablet) {
    document.documentElement.classList.add('android-tablet-platform');

    // Ensure Leaflet recalculates map dimensions after platform-specific styles are applied
    requestAnimationFrame(() => {
      if (typeof map !== 'undefined' && map) {
        setTimeout(() => map.invalidateSize(true), 0);
      }
    });
  }
  
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
