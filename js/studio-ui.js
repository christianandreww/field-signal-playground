/* Shared navigation, learning prompts, and the first simplified tool layout.
   Moves existing controls without replacing their handlers or numerical state. */
(function () {
  'use strict';
  const U = FSP.ui;
  const COURSES = {
    'Electromagnetics': ['EE3101', 'Electromagnetics'],
    'Signals & DSP': ['IE3014', 'Signals & DSP'],
    'Communications': ['IE3017', 'Comms & networks'],
    'Machines': ['EE2005', 'Devices & machines'],
  };
  // [course, navigation title, page title, description, experiment prompt]
  const TOOLS = {
    em: ['EE3101', 'Waves & boundaries', 'Waves at a boundary', 'Explore what happens when a wave meets a different material.', 'Start with Air → Glass, then change the angle. Compare reflected and transmitted power. Try Glass → Air to discover total internal reflection.'],
    fields: ['EE3101', 'Fields & polarization', 'Make invisible fields visible', 'Connect charge distributions, currents, and polarization to the fields they create.', 'Choose a topic, start with its default example, then change one parameter. Compare the numerical result with the closed-form solution.'],
    tline: ['EE3101', 'Transmission lines', 'Follow the impedance', 'Explore reflections and matching on a transmission line.', 'Start with a load equal to the line impedance. Change the load and watch the reflection coefficient move on the Smith chart.'],
    dsp: ['IE3014', 'Biquad lab', 'Shape your signal', 'Connect filter settings and pole–zero positions to what you see and hear.', 'Start with the lowpass filter. Move the cutoff, then edit a pole on the z-plane and compare the response. Audio starts only when you choose to play it.'],
    dsplab: ['IE3014', 'Convolution & spectra', 'See signals combine', 'Explore convolution, windows, and the frequency content of sampled signals.', 'Compare linear and circular convolution, then change a window and look for the trade-off between the main lobe and sidelobes.'],
    dspdesign: ['IE3014', 'Sampling & filter design', 'From samples to systems', 'Explore sampling, z-transforms, FIR filters, and implementation structures.', 'Begin with Sampling & aliasing. Change the sample rate, then use the other topics to connect difference equations, poles, and frequency response.'],
    filters: ['IE3014', 'Butterworth & Chebyshev', 'Design the response', 'Compare classic filters and understand their order, ripple, and transition bands.', 'Keep the cutoff and order the same while switching between Butterworth and Chebyshev. Compare the passband and transition region.'],
    audio: ['IE3014', 'Audio tools', 'Explore sound', 'Tune a note, inspect a spectrum, or listen to a reference tone.', 'Start with the 440 Hz sine tone. Compare its spectrum with a square wave. Microphone and file input remain optional.'],
    comms: ['IE3017', 'Communication channels', 'Connect signal to channel', 'Explore modulation, error rates, capacity, and network addressing.', 'Compare BPSK, QPSK, and 16-QAM at the same Eb/N0. Use the confidence interval to judge how much a finite simulation can tell you.'],
    net: ['IE3017', 'Networks & protocols', 'Make networks make sense', 'Trace coding, delays, access methods, and error detection step by step.', 'Choose a familiar bit sequence in Line coding. Compare its transitions across coding schemes, then explore timing and error detection in the other topics.'],
    magnetics: ['EE2005', 'Magnetics & transformers', 'Connect flux to the circuit', 'Explore magnetic circuits and transformer equivalent circuits.', 'Change the air gap and compare the required magnetomotive force. In the transformer topic, compare efficiency at different loads.'],
    edm2: ['EE2005', 'Transformers & actuators', 'From magnetic field to motion', 'Explore magnetic circuits, force, transformer connections, and energy losses.', 'Start with the default magnetic circuit. Compare the forward calculation and nonlinear solve before exploring force, transformer banks, and all-day efficiency.'],
    rotating: ['EE2005', 'Rotating machines', 'Find the operating point', 'Connect induction and DC machine behavior to the parameters you can change.', 'Begin with the induction motor. Compare its torque–speed curve at different voltages, then explore the load intersection and speed-control options.'],
  };
  const courseNames = Object.fromEntries(Object.values(COURSES));
  const groups = Array.from(document.querySelectorAll('.group-btn'));
  groups.forEach(button => {
    const course = COURSES[button.textContent]; if (!course) return;
    button.replaceChildren(U.el('span', { class: 'group-code', text: course[0] }), U.el('span', { class: 'group-name', text: course[1] }));
  });
  document.querySelectorAll('#tab-picker optgroup').forEach(group => { const course = COURSES[group.label]; if (course) group.label = course[0] + ' · ' + course[1]; });
  const list = document.querySelector('.tabs');
  list.before(U.el('p', { class: 'sidebar-label tools-label', text: 'Explore a tool' }));
  Object.entries(TOOLS).forEach(([id, info]) => {
    const tab = document.getElementById('tab-' + id), panel = document.getElementById('panel-' + id);
    if (!tab || !panel) return;
    tab.textContent = info[1];
    const option = document.querySelector('#tab-picker option[value="tab-' + id + '"]'); if (option) option.textContent = info[1];
    const tip = U.el('details', { class: 'lesson-tip' }, U.el('summary', { text: 'What to try' }), U.el('p', { text: info[4] }));
    panel.prepend(tip);
  });
  const appearance = document.getElementById('appearance-picker');
  appearance.value = FSP.appearance.get();
  appearance.addEventListener('change', () => FSP.appearance.set(appearance.value));
  window.addEventListener('fsp:appearance', () => { appearance.value = FSP.appearance.get(); FSP.redrawActive(); });
  FSP.studio = {
    activate(id) {
      const info = TOOLS[id]; if (!info) return;
      document.documentElement.dataset.subject = info[0];
      document.getElementById('workspace-course').textContent = info[0] + ' / ' + courseNames[info[0]];
      document.getElementById('workspace-title').textContent = info[2];
      document.getElementById('workspace-description').textContent = info[3];
      document.title = info[1] + ' · Field & Signal';
    },
    updateEm(an) {
      const values = [Number.isFinite(an.R) ? (100 * an.R).toFixed(1) + '%' : '—', Number.isFinite(an.T) ? (100 * an.T).toFixed(1) + '%' : '—', an.evanescent ? 'Evanescent' : Number.isFinite(an.thetaT) ? an.thetaT.toFixed(1) + '°' : '—'];
      ['em-reflected', 'em-transmitted', 'em-refracted'].forEach((id, i) => { document.getElementById(id).textContent = values[i]; });
      const text = an.conserv !== 'pass' ? 'In a lossy incident medium, interference also contributes to normal power flow. See the detailed results.' : an.evanescent ? 'The transmitted field is evanescent: it decays away from the boundary.' : an.R > 0.5 ? 'Most of the incident power is reflected at this boundary.' : 'More power enters the second medium than is reflected.';
      document.getElementById('em-takeaway').textContent = text;
    },
  };
  function disclosure(parent, title, elements, className) {
    const details = U.el('details', { class: className || 'control-details' }, U.el('summary', { text: title }));
    parent.appendChild(details); elements.forEach(element => details.appendChild(element)); return details;
  }
  const emPanel = document.getElementById('panel-em'), controls = emPanel.querySelector('.controls'), stage = emPanel.querySelector('.stage');
  const fieldsets = Array.from(controls.children);
  controls.prepend(U.el('p', { class: 'em-controls-label', text: 'Make it yours' }));
  const presets = fieldsets.find(f => f.querySelector('#em-presets'));
  controls.insertBefore(presets, fieldsets[0]);
  presets.querySelector('legend').textContent = 'Start with an example';
  fieldsets[0].querySelector('legend').textContent = 'Wave settings';
  document.querySelector('label[for="em-theta"]').textContent = 'Angle θi';
  document.querySelector('label[for="em-freq"]').textContent = 'Frequency';
  const map = document.getElementById('em-canvas').parentElement;
  const toolbar = U.el('div', { class: 'plot-toolbar' }, U.el('span', { text: 'Wave field · drag to explore' }));
  toolbar.appendChild(document.getElementById('em-play')); map.prepend(toolbar);
  disclosure(controls, 'Material properties', [fieldsets[1], fieldsets[2]]);
  disclosure(controls, 'Display options', [fieldsets[4]]);
  const results = U.el('div', { class: 'em-results', 'aria-label': 'Key results' });
  [['Reflected power', 'em-reflected'], ['Transmitted power', 'em-transmitted'], ['Refracted angle', 'em-refracted']].forEach(([text, id]) => results.appendChild(U.el('div', null, U.el('span', { class: 'em-result-label', text }), U.el('span', { class: 'em-result-value', id, text: '—' }))));
  map.after(results);
  results.after(U.el('p', { class: 'em-takeaway', id: 'em-takeaway' }));
  const detailed = disclosure(stage, 'Detailed results & angle sweep', [document.getElementById('em-hud'), document.getElementById('em-rt').parentElement], 'result-details');
  detailed.addEventListener('toggle', () => { if (detailed.open) FSP.redrawActive(); });
  const deepDive = U.el('div', { class: 'em-deep-dive' });
  deepDive.append(stage.querySelector('details.theory'), detailed);
  emPanel.querySelector('.layout').after(deepDive);
  if (window.matchMedia('(prefers-reduced-motion: reduce)').matches) {
    em.playing = false;
    document.getElementById('em-play').textContent = 'Play';
    document.getElementById('em-play').setAttribute('aria-pressed', 'false');
  }
  // Preserve the existing DSP controls, but do not lead with raw coefficients.
  const dspControls = document.querySelector('#panel-dsp .controls');
  const coefficients = document.getElementById('dsp-coefs').closest('fieldset');
  disclosure(dspControls, 'Filter coefficients & stability', [coefficients]);
  const selected = document.querySelector('.tabs [aria-selected="true"]');
  FSP.studio.activate(selected ? selected.id.slice(4) : 'em');
  FSP.redrawActive();
})();
