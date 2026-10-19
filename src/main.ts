/**
 * ============================================================================
 * AegisRadar AESA Signal Lab — Main Interactive Application Controller
 * ============================================================================
 * Orchestrates 3D Three.js WebGL rendering, 2D AESA Beamforming, Pulse Compression,
 * Adaptive CFAR Detection, Doppler MTD filtering, and tactical radar telemetry.
 */

import { Complex } from './dsp/fft_chirp';
import { runCfarDetector } from './dsp/cfar_detector';
import { AESAArray } from './aesa/beamforming_array';
import { computeMtdDopplerSpectrum } from './doppler/mtd_stap';

// Global Three.js declaration from browser bundle / CDN
declare const THREE: any;

// Radar RF & Engineering Parameters
const RADAR_PARAMS = {
  frequencyHz: 9.5e9,        // 9.5 GHz X-Band
  peakPowerWatts: 120e3,     // 120 kW
  antennaGainDbi: 42.5,      // 42.5 dBi
  antennaGainLinear: Math.pow(10, 42.5 / 10),
  pulseWidthSec: 10e-6,      // 10 microseconds
  prfHz: 12500,              // 12.5 kHz (High-PRF)
  noiseFigureDb: 3.2,
  boltzmannK: 1.380649e-23,
  temperatureKelvin: 290,
  bandwidthHz: 20e6,         // 20 MHz chirp bandwidth
  targetDistanceKm: 85.0,    // 85 km current track
};

const SPEED_OF_LIGHT = 299792458; // m/s
const WAVELENGTH = SPEED_OF_LIGHT / RADAR_PARAMS.frequencyHz; // ~0.0316 m (3.16 cm)

// State Variables
let currentAz = 15.0;
let currentEl = 8.0;
let targetRcs = 0.02; // Default: Drone 0.02 m^2
let isAudioMuted = true;
let currentLang = 'en';

// Three.js Scene Variables
let scene: any;
let camera: any;
let renderer: any;
let controls: any;
let aesaGroup: any;
let beamCone: any;
let elementMeshes: any[] = [];

// Instantiated DSP Engines
const aesaArray = new AESAArray({
  numElementsX: 18,
  numElementsY: 18,
  frequencyHz: RADAR_PARAMS.frequencyHz,
  apertureShape: 'Octagon',
  taperType: 'Taylor',
});

// Bilingual Dictionary
const I18N: Record<string, any> = {
  en: {
    subDesc: 'Active Electronically Scanned Array (AESA), Digital Beamforming & Pulse Doppler DSP Lab',
    steeringTitle: 'BEAM STEERING (DIGITAL)',
    lblAzimuth: 'Azimuth Angle (θ_az)',
    lblElevation: 'Elevation Angle (θ_el)',
    ppiTitle: 'PPI RADAR SCOPE (180 km)',
    eqTitle: 'RADAR RANGE EQUATION (Pr & SNR)',
    fftTitle: 'DOPPLER FFT VELOCITY SPECTRUM (MTD)',
  },
  ar: {
    subDesc: 'مختبر مصفوفة المسح الإلكتروني النشط (AESA)، تشكيل الحزم الرقمي ومعالجة نبضات دوبلر',
    steeringTitle: 'توجيه الحزمة الرقمية',
    lblAzimuth: 'زاوية السمت (السمت الأفقي)',
    lblElevation: 'زاوية الارتفاع (الارتفاع الرأسي)',
    ppiTitle: 'شاشة الرادار الدائرية PPI (١٨٠ كم)',
    eqTitle: 'معادلة مدى الرادار (القدرة والـ SNR)',
    fftTitle: 'طيف سرعة دوبلر FFT (فلترة MTD)',
  },
};

// Web Audio API Synthesizer
let audioCtx: AudioContext | null = null;

function initAudio(): void {
  if (!audioCtx) {
    const AudioContextClass = window.AudioContext || (window as any).webkitAudioContext;
    audioCtx = new AudioContextClass();
  }
}

function playRadarPing(): void {
  if (isAudioMuted || !audioCtx) return;
  try {
    const now = audioCtx.currentTime;
    const osc = audioCtx.createOscillator();
    const gain = audioCtx.createGain();

    // High chirp sweeping down slightly
    osc.type = 'sine';
    osc.frequency.setValueAtTime(2400, now);
    osc.frequency.exponentialRampToValueAtTime(800, now + 0.12);

    gain.gain.setValueAtTime(0.08, now);
    gain.gain.exponentialRampToValueAtTime(0.001, now + 0.15);

    osc.connect(gain);
    gain.connect(audioCtx.destination);

    osc.start(now);
    osc.stop(now + 0.16);
  } catch {
    // Ignore audio autoplay restrictions
  }
}

/**
 * Initializes the Three.js 3D WebGL Scene & AESA Array
 */
function initThree(): void {
  const container = document.getElementById('threeContainer');
  if (!container || typeof THREE === 'undefined') return;

  const width = container.clientWidth || 800;
  const height = container.clientHeight || 600;

  scene = new THREE.Scene();
  scene.background = new THREE.Color(0x020617);
  scene.fog = new THREE.FogExp2(0x020617, 0.012);

  camera = new THREE.PerspectiveCamera(45, width / height, 0.1, 1000);
  camera.position.set(0, 18, 28);

  renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
  renderer.setSize(width, height);
  renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
  renderer.shadowMap.enabled = true;
  container.appendChild(renderer.domElement);

  if (typeof THREE.OrbitControls !== 'undefined') {
    controls = new THREE.OrbitControls(camera, renderer.domElement);
    controls.enableDamping = true;
    controls.dampingFactor = 0.05;
    controls.maxPolarAngle = Math.PI / 2 + 0.1;
  }

  // Lighting
  const ambientLight = new THREE.AmbientLight(0x38bdf8, 0.4);
  scene.add(ambientLight);

  const dirLight = new THREE.DirectionalLight(0x10b981, 1.2);
  dirLight.position.set(15, 30, 20);
  scene.add(dirLight);

  const fillLight = new THREE.DirectionalLight(0x06b6d4, 0.6);
  fillLight.position.set(-15, -10, -10);
  scene.add(fillLight);

  // Tactical Coordinate Grid
  const grid = new THREE.GridHelper(50, 50, 0x059669, 0x1e293b);
  grid.position.y = -4;
  scene.add(grid);

  // Build AESA Array Mounting & T/R Element Modules
  buildAesa3DModel();

  // Create Steered Electromagnetic Beam Cone
  buildBeamGeometry();

  window.addEventListener('resize', onWindowResize);
  solveRadarEquation();
  animate();
}

/**
 * Builds the 3D physical AESA Aperture and GaN T/R Modules
 */
function buildAesa3DModel(): void {
  aesaGroup = new THREE.Group();

  // Array backplane chassis
  const chassisGeo = new THREE.CylinderGeometry(9.5, 10.0, 1.2, 8);
  const chassisMat = new THREE.MeshStandardMaterial({
    color: 0x0f172a,
    metalness: 0.85,
    roughness: 0.3,
  });
  const chassis = new THREE.Mesh(chassisGeo, chassisMat);
  chassis.rotation.x = Math.PI / 2;
  aesaGroup.add(chassis);

  // Bezel ring
  const ringGeo = new THREE.TorusGeometry(9.6, 0.25, 8, 8);
  const ringMat = new THREE.MeshStandardMaterial({
    color: 0x10b981,
    emissive: 0x047857,
    emissiveIntensity: 0.4,
  });
  const ring = new THREE.Mesh(ringGeo, ringMat);
  aesaGroup.add(ring);

  // Instanced T/R radiating patch elements
  const elemGeo = new THREE.BoxGeometry(0.55, 0.55, 0.4);
  elementMeshes = [];

  const scale = 2.2; // Visual scaling factor
  for (const el of aesaArray.elements) {
    const elemMat = new THREE.MeshStandardMaterial({
      color: 0x10b981,
      metalness: 0.9,
      roughness: 0.2,
      emissive: 0x065f46,
      emissiveIntensity: 0.3,
    });
    const mesh = new THREE.Mesh(elemGeo, elemMat);
    mesh.position.set(el.x * scale * 100, el.y * scale * 100, 0.6);
    aesaGroup.add(mesh);
    elementMeshes.push(mesh);
  }

  scene.add(aesaGroup);
  updateElementPhaseColors();
}

/**
 * Creates dynamic 3D Steered Electromagnetic Beam Cone
 */
function buildBeamGeometry(): void {
  const coneGeo = new THREE.ConeGeometry(5.0, 24, 32, 1, true);
  coneGeo.translate(0, 12, 0);
  coneGeo.rotateX(Math.PI / 2);

  const coneMat = new THREE.MeshBasicMaterial({
    color: 0x10b981,
    transparent: true,
    opacity: 0.35,
    wireframe: true,
    side: THREE.DoubleSide,
  });

  beamCone = new THREE.Mesh(coneGeo, coneMat);
  scene.add(beamCone);
  updateBeamGeometry();
}

/**
 * Updates 3D Beam vector orientation according to Azimuth & Elevation
 */
function updateBeamGeometry(): void {
  if (!beamCone) return;
  aesaArray.steerBeam(currentAz, currentEl);
  updateElementPhaseColors();

  const azRad = (-currentAz * Math.PI) / 180;
  const elRad = (currentEl * Math.PI) / 180;

  beamCone.rotation.order = 'YXZ';
  beamCone.rotation.y = azRad;
  beamCone.rotation.x = -elRad;
}

/**
 * Shades each AESA T/R element color in real-time according to its steering phase delay
 */
function updateElementPhaseColors(): void {
  for (let i = 0; i < aesaArray.elements.length; i++) {
    const el = aesaArray.elements[i];
    const mesh = elementMeshes[i];
    if (!mesh) continue;

    // Normalize phase [-pi, pi] to [0, 1] for visual hue
    const normPhase = ((el.phaseRad % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI);
    const hue = normPhase / (2 * Math.PI);

    // Dynamic color gradient (Emerald -> Cyan -> Amber)
    mesh.material.color.setHSL(0.35 + hue * 0.25, 0.9, 0.5 * el.amplitude);
    mesh.material.emissive.setHSL(0.35 + hue * 0.25, 0.8, 0.2 * el.amplitude);
  }
}

/**
 * Solves the radar range equation and updates telemetry displays
 */
function solveRadarEquation(): void {
  const pt = RADAR_PARAMS.peakPowerWatts;
  const g = RADAR_PARAMS.antennaGainLinear;
  const lambda = WAVELENGTH;
  const sigma = targetRcs;
  const rMeters = RADAR_PARAMS.targetDistanceKm * 1000;

  // Received Power: Pr = (Pt * G^2 * lambda^2 * sigma) / ((4*pi)^3 * R^4)
  const numerator = pt * g * g * lambda * lambda * sigma;
  const denominator = Math.pow(4 * Math.PI, 3) * Math.pow(rMeters, 4);
  const prWatts = numerator / denominator;
  const prDbm = 10 * Math.log10(Math.max(1e-20, prWatts)) + 30;

  // Thermal Noise Power: N = k * T * B * F
  const noiseWatts =
    RADAR_PARAMS.boltzmannK *
    RADAR_PARAMS.temperatureKelvin *
    RADAR_PARAMS.bandwidthHz *
    Math.pow(10, RADAR_PARAMS.noiseFigureDb / 10);

  // SNR (linear and dB)
  const snrLinear = prWatts / noiseWatts;
  const snrDb = 10 * Math.log10(Math.max(1e-5, snrLinear));

  // Max Detection Range for SNR_min = 13 dB (P_d ~ 0.9, P_fa ~ 1e-6)
  const snrMinLinear = Math.pow(10, 13 / 10);
  const maxRangeMeters = Math.pow(numerator / (denominator / Math.pow(rMeters, 4) * noiseWatts * snrMinLinear), 0.25);
  const maxRangeKm = maxRangeMeters / 1000;

  // Update UI Elements
  const setEl = (id: string, val: string) => {
    const el = document.getElementById(id);
    if (el) el.innerText = val;
  };

  const snrText = `${snrDb >= 0 ? '+' : ''}${snrDb.toFixed(1)} dB (${snrDb > 12 ? 'LOCK' : 'SEARCH'})`;
  setEl('eqPr', `${prDbm.toFixed(1)} dBm`);
  setEl('calcPrVal', `${prDbm.toFixed(1)} dBm`);
  setEl('eqSnr', snrText);
  setEl('calcSnrVal', snrText);
  setEl('eqMaxR', `${maxRangeKm.toFixed(0)} km`);
  setEl('calcRangeVal', `${RADAR_PARAMS.targetDistanceKm.toFixed(1)} km`);
  setEl('calcRcsVal', `${targetRcs.toFixed(3)} m²`);
}

/**
 * Draws the tactical Plan Position Indicator (PPI) Polar Scope
 */
function drawPpi(): void {
  const canvas = document.getElementById('ppiCanvas') as HTMLCanvasElement;
  if (!canvas) return;
  const ctx = canvas.getContext('2d');
  if (!ctx) return;

  const w = canvas.width;
  const h = canvas.height;
  const cx = w / 2;
  const cy = h / 2;
  const maxR = Math.min(cx, cy) - 6;

  // Background clear with slight decay for phosphor persistence
  ctx.fillStyle = 'rgba(2, 6, 23, 0.25)';
  ctx.fillRect(0, 0, w, h);

  // Range rings
  ctx.strokeStyle = 'rgba(16, 185, 129, 0.22)';
  ctx.lineWidth = 1;
  const rings = 4;
  for (let i = 1; i <= rings; i++) {
    ctx.beginPath();
    ctx.arc(cx, cy, (maxR * i) / rings, 0, 2 * Math.PI);
    ctx.stroke();
  }

  // Crosshairs
  ctx.beginPath();
  ctx.moveTo(cx, cy - maxR);
  ctx.lineTo(cx, cy + maxR);
  ctx.moveTo(cx - maxR, cy);
  ctx.lineTo(cx + maxR, cy);
  ctx.stroke();

  // Rotating sweep line
  const time = Date.now() * 0.0018;
  const sweepAngle = time % (2 * Math.PI);

  ctx.strokeStyle = 'rgba(52, 211, 153, 0.85)';
  ctx.lineWidth = 1.8;
  ctx.beginPath();
  ctx.moveTo(cx, cy);
  ctx.lineTo(cx + maxR * Math.cos(sweepAngle), cy + maxR * Math.sin(sweepAngle));
  ctx.stroke();

  // Target Blip (corresponding to targetDistanceKm and currentAz)
  const targetAngleRad = ((currentAz - 90) * Math.PI) / 180;
  const targetDistRatio = Math.min(1.0, RADAR_PARAMS.targetDistanceKm / 180.0);
  const targetX = cx + maxR * targetDistRatio * Math.cos(targetAngleRad);
  const targetY = cy + maxR * targetDistRatio * Math.sin(targetAngleRad);

  // Paint detected target
  ctx.fillStyle = '#34d399';
  ctx.shadowColor = '#10b981';
  ctx.shadowBlur = 10;
  ctx.beginPath();
  ctx.arc(targetX, targetY, 3.5, 0, 2 * Math.PI);
  ctx.fill();
  ctx.shadowBlur = 0;
}

/**
 * Draws the Doppler FFT Velocity Spectrum & CFAR Detection Thresholds
 */
function drawFft(): void {
  const canvas = document.getElementById('fftCanvas') as HTMLCanvasElement;
  if (!canvas) return;
  const ctx = canvas.getContext('2d');
  if (!ctx) return;

  const w = canvas.width;
  const h = canvas.height;

  ctx.fillStyle = 'rgba(2, 6, 23, 0.35)';
  ctx.fillRect(0, 0, w, h);

  // Generate synthetic pulse train with target Doppler return and ground clutter
  const numPulses = 32;
  const targetFd = 38700; // 38.7 kHz Doppler frequency (radial velocity ~ 612 m/s)
  const dt = 1 / RADAR_PARAMS.prfHz;

  const slowTime: Complex[] = [];
  for (let i = 0; i < numPulses; i++) {
    const t = i * dt;
    // Ground clutter at zero Doppler + high speed threat
    const groundClutter = Complex.fromPolar(0.4 + 0.1 * Math.sin(i * 0.3), 0);
    const targetSignal = Complex.fromPolar(1.8, 2 * Math.PI * targetFd * t);
    const noise = Complex.fromPolar(Math.random() * 0.15, Math.random() * 2 * Math.PI);
    slowTime.push(groundClutter.add(targetSignal).add(noise));
  }

  const mtdResult = computeMtdDopplerSpectrum(slowTime, {
    numPulses,
    prfHz: RADAR_PARAMS.prfHz,
    wavelengthMeters: WAVELENGTH,
    windowType: 'Hamming',
  });

  const spectrum = mtdResult.powerSpectrumDb;
  const len = spectrum.length;

  // Run CA-CFAR on power spectrum
  const powerLinear = Array.from(spectrum).map(db => Math.pow(10, db / 10));
  const cfar = runCfarDetector(powerLinear, {
    type: 'CA-CFAR',
    numGuardCells: 2,
    numReferenceCells: 4,
    pfa: 1e-4,
  });

  // Plot CFAR Threshold (Amber line)
  ctx.strokeStyle = '#f59e0b';
  ctx.lineWidth = 1.2;
  ctx.beginPath();
  for (let i = 0; i < len; i++) {
    const x = (i / (len - 1)) * w;
    const threshDb = 10 * Math.log10(Math.max(1e-12, cfar.thresholds[i]));
    const normY = Math.min(1.0, Math.max(0.0, (threshDb + 60) / 60));
    const y = h - normY * (h - 10) - 5;
    if (i === 0) ctx.moveTo(x, y);
    else ctx.lineTo(x, y);
  }
  ctx.stroke();

  // Plot Doppler Spectrum (Emerald curve)
  ctx.strokeStyle = '#10b981';
  ctx.lineWidth = 1.8;
  ctx.beginPath();
  for (let i = 0; i < len; i++) {
    const x = (i / (len - 1)) * w;
    const normY = Math.min(1.0, Math.max(0.0, (spectrum[i] + 60) / 60));
    const y = h - normY * (h - 10) - 5;
    if (i === 0) ctx.moveTo(x, y);
    else ctx.lineTo(x, y);
  }
  ctx.stroke();

  // Highlight CFAR Detections
  ctx.fillStyle = '#ef4444';
  for (const idx of cfar.detectedIndices) {
    const x = (idx / (len - 1)) * w;
    const normY = Math.min(1.0, Math.max(0.0, (spectrum[idx] + 60) / 60));
    const y = h - normY * (h - 10) - 5;
    ctx.beginPath();
    ctx.arc(x, y, 3, 0, 2 * Math.PI);
    ctx.fill();
  }
}

function animate(): void {
  requestAnimationFrame(animate);
  if (controls) controls.update();
  drawPpi();
  drawFft();
  if (renderer && scene && camera) {
    renderer.render(scene, camera);
  }
}

function onWindowResize(): void {
  const container = document.getElementById('threeContainer');
  if (!container || !renderer || !camera) return;
  const width = container.clientWidth;
  const height = container.clientHeight;
  camera.aspect = width / height;
  camera.updateProjectionMatrix();
  renderer.setSize(width, height);
}

/**
 * Sets up tactile DOM Event Listeners
 */
function setupEventListeners(): void {
  // Azimuth Steering Slider
  const sliderAz = document.getElementById('sliderAzimuth') as HTMLInputElement;
  if (sliderAz) {
    sliderAz.addEventListener('input', (e: Event) => {
      currentAz = parseFloat((e.target as HTMLInputElement).value);
      const valAz = document.getElementById('valAzimuth');
      if (valAz) valAz.innerText = `${currentAz >= 0 ? '+' : ''}${currentAz}°`;
      const telSteer = document.getElementById('telSteer');
      if (telSteer) telSteer.innerText = `AZ: ${currentAz >= 0 ? '+' : ''}${currentAz}° | EL: +${currentEl}°`;
      const ppiBearing = document.getElementById('ppiBearing');
      if (ppiBearing) ppiBearing.innerText = `${((currentAz + 360) % 360).toFixed(0).padStart(3, '0')}°`;
      updateBeamGeometry();
    });
  }

  // Elevation Steering Slider
  const sliderEl = document.getElementById('sliderElevation') as HTMLInputElement;
  if (sliderEl) {
    sliderEl.addEventListener('input', (e: Event) => {
      currentEl = parseFloat((e.target as HTMLInputElement).value);
      const valEl = document.getElementById('valElevation');
      if (valEl) valEl.innerText = `${currentEl >= 0 ? '+' : ''}${currentEl}°`;
      const telSteer = document.getElementById('telSteer');
      if (telSteer) telSteer.innerText = `AZ: ${currentAz >= 0 ? '+' : ''}${currentAz}° | EL: +${currentEl}°`;
      updateBeamGeometry();
    });
  }

  // Target RCS Presets
  const presets = [
    { id: 'rcsStealthBtn', val: 0.005 },
    { id: 'rcsDroneBtn', val: 0.02 },
    { id: 'rcsMissileBtn', val: 0.1 },
    { id: 'rcsFighterBtn', val: 3.0 },
  ];

  presets.forEach(p => {
    const btn = document.getElementById(p.id);
    if (!btn) return;
    btn.addEventListener('click', () => {
      targetRcs = p.val;
      presets.forEach(item => {
        const otherBtn = document.getElementById(item.id);
        if (!otherBtn) return;
        if (item.id === p.id) {
          otherBtn.className = 'px-2.5 py-1 text-xs font-semibold rounded bg-emerald-600 text-white shadow-md flex items-center space-x-1';
        } else {
          otherBtn.className = 'px-2.5 py-1 text-xs font-semibold rounded bg-slate-800 text-slate-300 hover:text-emerald-300 flex items-center space-x-1';
        }
      });
      playRadarPing();
      solveRadarEquation();
    });
  });

  // Camera Orbit Reset
  const resetCamBtn = document.getElementById('resetCamBtn');
  if (resetCamBtn) {
    resetCamBtn.addEventListener('click', () => {
      if (camera && controls) {
        camera.position.set(0, 18, 28);
        controls.target.set(0, 0, 0);
      }
    });
  }

  // Audio Mute/Unmute Toggle
  const soundBtn = document.getElementById('soundBtn');
  if (soundBtn) {
    soundBtn.addEventListener('click', () => {
      initAudio();
      isAudioMuted = !isAudioMuted;
      const soundIcon = document.getElementById('soundIcon');
      if (soundIcon) {
        soundIcon.className = isAudioMuted
          ? 'fa-solid fa-volume-xmark'
          : 'fa-solid fa-volume-high text-emerald-400';
      }
    });
  }

  // Language English / Arabic Toggle
  const langToggleBtn = document.getElementById('langToggleBtn');
  if (langToggleBtn) {
    langToggleBtn.addEventListener('click', () => {
      currentLang = currentLang === 'en' ? 'ar' : 'en';
      document.documentElement.dir = currentLang === 'ar' ? 'rtl' : 'ltr';
      langToggleBtn.innerText = currentLang === 'ar' ? 'English' : 'عربي';

      const t = I18N[currentLang];
      const setTxt = (id: string, text: string) => {
        const el = document.getElementById(id);
        if (el) el.innerText = text;
      };

      setTxt('subHeaderDesc', t.subDesc);
      const steeringHeader = document.getElementById('txtSteeringTitle');
      if (steeringHeader) {
        steeringHeader.innerHTML = `<i class="fa-solid fa-compass mr-1"></i> ${t.steeringTitle}`;
      }
      setTxt('lblAzimuth', t.lblAzimuth);
      setTxt('lblElevation', t.lblElevation);
      setTxt('ppiTitle', t.ppiTitle);
      setTxt('eqTitle', t.eqTitle);
      setTxt('fftTitle', t.fftTitle);
    });
  }
}

// Bootstrap on DOMContentLoaded
if (document.readyState === 'loading') {
  document.addEventListener('DOMContentLoaded', () => {
    initThree();
    setupEventListeners();
  });
} else {
  initThree();
  setupEventListeners();
}
