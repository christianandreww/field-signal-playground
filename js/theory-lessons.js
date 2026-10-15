/* Small, self-contained lessons. No audio, network access, or automatic animation. */
(function () {
  'use strict';
  const el = FSP.ui.el;
  const style = el('style', { text: '.foundation{margin:16px 0;padding:18px;border:1px solid var(--border);border-radius:16px;background:var(--panel)}.foundation>summary{cursor:pointer;font-weight:650}.foundation p{line-height:1.65;max-width:80ch}.foundation label{display:block;margin:12px 0}.foundation select{max-width:100%}.foundation output{display:block;white-space:pre-wrap;overflow-wrap:anywhere;padding:14px;background:var(--panel2);border-radius:10px;line-height:1.7}.foundation svg{display:block;width:100%;max-width:420px;height:auto;margin:12px auto;color:var(--accent)}.foundation .lesson-step{color:var(--muted);font-size:.85rem;letter-spacing:.04em}.foundation input[type=range]{width:min(240px,100%);vertical-align:middle}' });
  document.head.append(style);
  function lesson(panel, title, theory, worked) {
    const box = el('details', { class: 'foundation' }, el('summary', { text: title }));
    box.append(el('p', { class: 'lesson-step', text: '01 / THEORY' }), el('p', { text: theory }), el('p', { class: 'lesson-step', text: '02 / WORKED EXAMPLE' }), el('p', { text: worked }), el('p', { class: 'lesson-step', text: '03 / EXPLORE' }));
    document.getElementById('panel-' + panel).prepend(box);
    return box;
  }
  function choice(box, name, options, update) {
    const select = el('select', { 'aria-label': name });
    options.forEach(([value, text]) => select.append(el('option', { value, text })));
    box.append(el('label', null, el('span', { text: name + ' ' }), select));
    select.addEventListener('change', update); return select;
  }
  function slider(box, name, min, max, value, step, update) {
    const input = el('input', { type: 'range', min, max, value, step, 'aria-label': name });
    box.append(el('label', null, el('span', { text: name + ' ' }), input));
    input.addEventListener('input', update); return input;
  }
  function output(box) { const out = el('output', { 'aria-live': 'polite', 'aria-atomic': 'true' }); box.append(out); return out; }
  function svg(box) { const node = document.createElementNS('http://www.w3.org/2000/svg', 'svg'); node.setAttribute('viewBox', '0 0 360 260'); node.setAttribute('role', 'img'); box.append(node); return node; }
  function shape(view, tag, attrs) { const node = document.createElementNS(view.namespaceURI, tag); Object.entries(attrs).forEach(([key, value]) => node.setAttribute(key, value)); view.append(node); return node; }

  const fields = lesson('fields', 'Start here · Vector fields, Gauss & Maxwell',
    'The gradient points toward increasing scalar potential; an electrostatic field is E = −∇V. Divergence measures local outward flux, while curl measures local circulation. Gauss’s theorem links a closed-surface flux to volume divergence; Stokes’s theorem links a closed-loop circulation to surface curl.',
    'For the dimensionless field F = (x, y, 0), ∇·F = 2 and ∇×F = 0. Around a circle of radius R, outward flux per unit depth is 2πR². Its circulation is zero. For F = (−y, x, 0), those roles reverse: divergence is zero and the z-component of curl is 2.');
  const fieldKind = choice(fields, 'Field pattern', [['source', 'Outward field (x, y)'], ['rotation', 'Rotating field (−y, x)'], ['uniform', 'Uniform field (1, 0)']], drawField);
  const fieldView = svg(fields), fieldOut = output(fields);
  function drawField() {
    fieldView.replaceChildren(); fieldView.setAttribute('aria-label', fieldKind.selectedOptions[0].text + ' vector arrows');
    shape(fieldView, 'circle', { cx: 180, cy: 130, r: 85, fill: 'none', stroke: 'currentColor', 'stroke-dasharray': '5 5' });
    for (let x = -2; x <= 2; x++) for (let y = -2; y <= 2; y++) {
      const [vx, vy] = fieldKind.value === 'source' ? [x, y] : fieldKind.value === 'rotation' ? [-y, x] : [1, 0];
      if (!vx && !vy) continue;
      const px = 180 + 42*x, py = 130 - 42*y, dx = vx*9, dy = -vy*9, a = Math.atan2(dy, dx);
      shape(fieldView, 'path', { d: `M${px},${py} l${dx},${dy} m${-6*Math.cos(a-.5)},${-6*Math.sin(a-.5)} L${px+dx},${py+dy} l${-6*Math.cos(a+.5)},${-6*Math.sin(a+.5)}`, fill: 'none', stroke: 'currentColor', 'stroke-width': 2 });
    }
    const source = fieldKind.value === 'source', rotation = fieldKind.value === 'rotation';
    fieldOut.textContent = `Dimensionless 2D field, unit-radius circle, unit depth\nDivergence: ${source ? 2 : 0} · Curl z: ${rotation ? 2 : 0}\nOutward flux: ${source ? '2π' : '0'} · Counterclockwise circulation: ${rotation ? '2π' : '0'}\n${rotation ? 'Circulation can exist without outward flux. A static electric field cannot have this pattern.' : source ? 'Positive divergence identifies a distributed source, not a point charge.' : 'Equal inward and outward contributions cancel.'}`;
  }
  drawField();
  fields.append(el('p', { text: 'Gauss’s law in vacuum: ∮ E·dA = Qenclosed/ε₀. For a point charge at the centre of a sphere, symmetry makes E constant on its surface: E = Q/(4πε₀R²). Increasing the radius reduces E, but the total flux stays constant.' }));
  const radius = slider(fields, 'Gaussian sphere radius (m)', .5, 3, 1, .1, gauss), charge = slider(fields, 'Enclosed point charge (nC)', -5, 5, 1, .5, gauss), gaussOut = output(fields);
  function gauss() { const r = +radius.value, q = +charge.value*1e-9, e = q/(4*Math.PI*8.8541878128e-12*r*r); gaussOut.textContent = `R = ${r.toFixed(1)} m · Q = ${charge.value} nC\nSigned radial E = ${e.toFixed(3)} V/m\nClosed-surface flux = ${(q/8.8541878128e-12).toFixed(3)} V·m\nAssumption: isolated central point charge in vacuum. Negative E points inward.`; }
  gauss();
  fields.append(el('p', { text: 'Maxwell connects the pictures: ∇·D = ρfree (electric sources); ∇·B = 0 (no magnetic monopoles); ∇×E = −∂B/∂t (Faraday induction); ∇×H = Jfree + ∂D/∂t (current and displacement current). In vacuum D = ε₀E and B = μ₀H. Integral forms describe flux through closed surfaces and circulation around oriented loops; the local forms describe the same laws at each point. The surface normal and loop direction follow the right-hand rule.' }));

  const fft = lesson('dspdesign', 'Start here · Inside an FFT butterfly',
    'A radix-2 FFT splits a DFT into even- and odd-indexed samples. Each butterfly combines two intermediate values: u + Wv and u − Wv, where W = exp(−j2πk/N). A decimation-in-time implementation starts with bit-reversed input indices.',
    'For x = [1, 0, 0, 0, 0, 0, 0, 0], every DFT bin is 1. With eight constant samples equal to 1, only bin 0 is nonzero and equals 8. These are unnormalised forward DFT values.');
  const signal = choice(fft, 'Eight-sample input', [['impulse', 'Unit impulse'], ['constant', 'Constant = 1'], ['cosine', 'Cosine at bin 1']], drawFFT);
  const stage = slider(fft, 'FFT stage (0 = bit-reversed input)', 0, 3, 0, 1, drawFFT), fftOut = output(fft);
  function drawFFT() {
    const x = Array.from({ length: 8 }, (_, i) => signal.value === 'impulse' ? +(i === 0) : signal.value === 'constant' ? 1 : Math.cos(2*Math.PI*i/8));
    const reverse = [0,4,2,6,1,5,3,7], a = reverse.map(i => [x[i], 0]), steps = [];
    for (let s = 1; s <= +stage.value; s++) {
      const n = 2**s;
      for (let base = 0; base < 8; base += n) for (let k = 0; k < n/2; k++) {
        const angle = -2*Math.PI*k/n, u = a[base+k], b = a[base+k+n/2], v = [Math.cos(angle)*b[0]-Math.sin(angle)*b[1], Math.sin(angle)*b[0]+Math.cos(angle)*b[1]];
        a[base+k] = [u[0]+v[0], u[1]+v[1]]; a[base+k+n/2] = [u[0]-v[0], u[1]-v[1]];
        if (s === +stage.value) steps.push(`${base+k} ↔ ${base+k+n/2}: W = exp(−j2π × ${k}/${n})`);
      }
    }
    const fmt = v => (Math.abs(v) < 1e-9 ? 0 : v).toFixed(3);
    fftOut.textContent = `Input: [${x.map(fmt).join(', ')}]\nBit-reversed indices: ${reverse.join(', ')}\nStage ${stage.value}${+stage.value === 3 ? ' · final bins in normal order' : ' · intermediate values, not frequency bins'}\n${a.map((z, i) => `${i}: ${fmt(z[0])} ${z[1]<-1e-9 ? '−' : '+'} j${fmt(Math.abs(z[1]))}`).join('\n')}\n${steps.join('\n')}`;
  }
  drawFFT();

  const network = lesson('net', 'Start here · A packet’s journey',
    'TCP delivers a byte stream; IP delivers packets between networks; Ethernet delivers frames on a local link. A switch learns source MAC addresses. A router chooses the longest matching destination prefix. For a remote IP destination, the sender resolves the gateway’s MAC address with ARP—not the remote host’s MAC.',
    'Host A is 192.168.1.10/24, its gateway is 192.168.1.1, and server B is 203.0.113.20. A sees that B is outside its subnet, then sends a frame addressed to the gateway. The enclosed IP packet still names B as its destination.');
  const destination = choice(network, 'Destination', [['remote', 'Remote server · 203.0.113.20'], ['local', 'Local host · 192.168.1.20']], drawPacket);
  const hop = slider(network, 'Journey step', 0, 4, 0, 1, drawPacket), packetOut = output(network);
  function drawPacket() {
    const remote = destination.value === 'remote', target = remote ? '203.0.113.20' : '192.168.1.20';
    const steps = [
      `Subnet decision: 192.168.1.10/24 compares destination ${target} with 192.168.1.0/24. ${remote ? 'Different subnet: next hop is gateway 192.168.1.1.' : 'Same subnet: next hop is the destination itself.'}`,
      `Empty ARP cache: broadcast “Who has ${remote ? '192.168.1.1' : target}?” The owner replies with its MAC address. Switches flood the broadcast within the LAN and learn source MACs.`,
      `Encapsulation: Ethernet [source MAC A → ${remote ? 'gateway MAC' : 'host B MAC'}] contains IP [192.168.1.10 → ${target}] containing a TCP segment.`,
      remote ? 'Router removes the incoming Ethernet header. Routes 0.0.0.0/0 and 203.0.113.0/24 both match; /24 wins. It decrements TTL, updates the IPv4 header checksum and creates a frame for the next link.' : 'The switch forwards using the destination MAC table entry, or floods an unknown unicast. Host B accepts the frame; no router is required.',
      'The receiver removes link and IP headers and passes the segment to TCP. TCP uses sequence numbers and acknowledgements to deliver an ordered byte stream; a TCP acknowledgement is not an Ethernet acknowledgement.'
    ];
    packetOut.textContent = `A → switch → ${remote ? 'gateway → routed network → server' : 'local host B'}\nStep ${+hop.value+1} / 5\n${steps[+hop.value]}\n\nModel: IPv4 over Ethernet, no NAT, no VLANs, no loss. ARP is shown only for an initially empty cache. TCP connection setup is outside this example.`;
  }
  drawPacket();

  const machine = lesson('rotating', 'Start here · Synchronous machines',
    'In steady synchronous operation, the rotor turns at the speed of the stator’s rotating magnetic field: ns = 120f/P rpm, with P the number of poles. Unlike an induction motor, a synchronous motor develops steady torque at zero slip. Wound-field rotors use excitation current; permanent-magnet rotors supply flux without rotor excitation current.',
    'A four-pole machine on 50 Hz runs synchronously at 1500 rpm. For a cylindrical-rotor machine with negligible resistance, electromagnetic power is P = 3Vphase Ephase sin(δ)/Xs. Here δ is the signed angle by which internal EMF E leads terminal voltage V: positive power means generation, negative means motoring.');
  const frequency = slider(machine, 'Supply frequency (Hz)', 25, 60, 50, 1, drawMachine);
  const poles = choice(machine, 'Number of poles', [['2', '2'], ['4', '4'], ['6', '6'], ['8', '8']], drawMachine); poles.value = '4';
  const angle = slider(machine, 'Internal EMF lead angle δ (degrees)', -120, 120, -30, 1, drawMachine), machineView = svg(machine), machineOut = output(machine);
  function drawMachine() {
    const ns = 120*(+frequency.value)/(+poles.value), rad = (+angle.value)*Math.PI/180, power = 3*230*230/5*Math.sin(rad);
    machineView.replaceChildren(); machineView.setAttribute('aria-label', `Terminal voltage and internal EMF phasors separated by ${angle.value} degrees`);
    shape(machineView, 'circle', { cx: 150, cy: 130, r: 90, fill: 'none', stroke: 'currentColor', opacity: .25 });
    shape(machineView, 'path', { d: 'M150 130 H240', stroke: 'currentColor', 'stroke-width': 3 });
    shape(machineView, 'path', { d: `M150 130 L${150+90*Math.cos(rad)} ${130-90*Math.sin(rad)}`, stroke: 'currentColor', 'stroke-width': 3, 'stroke-dasharray': '6 3' });
    const label = shape(machineView, 'text', { x: 15, y: 245, fill: 'currentColor', 'font-size': 13 }); label.textContent = 'Solid: V reference · dashed: internal E';
    machineOut.textContent = `Synchronous speed = ${ns.toFixed(0)} rpm · slip = 0\nδ = ${angle.value}° · ${power > 1e-8 ? 'Generating' : power < -1e-8 ? 'Motoring' : 'Zero converted power'}\nElectromagnetic power = ${(power/1000).toFixed(2)} kW\nTorque magnitude = ${Math.abs(power/(ns*2*Math.PI/60)).toFixed(2)} N·m\n${Math.abs(+angle.value)<90 ? 'Positive synchronising stiffness: normal stable branch in this ideal model.' : Math.abs(+angle.value)===90 ? 'Maximum power magnitude: zero synchronising stiffness.' : 'Negative synchronising stiffness: unstable branch in this ideal model.'}\nAssumptions: balanced three-phase, cylindrical rotor, Vphase = Ephase = 230 V RMS, Xs = 5 Ω/phase; resistance, losses and saturation neglected. This power-angle model does not include permanent-magnet saliency or starting transients.`;
  }
  drawMachine();
})();
