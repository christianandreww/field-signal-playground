/* Audio Tools: YIN tuner, scrolling spectrogram, reference tone generator, file player.
   Pure math in FSP.math.audio; the tab owns its own AudioContext. No DOM access at load. */
(function () {
  'use strict';

  /* ================= pure math ================= */
  const NOTE_NAMES = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B'];

  // YIN (de Cheveigne & Kawahara 2002). Returns Hz or null (unvoiced). opts.detail=true -> object with diagnostics.
  function yin(x, fs, opts) {
    opts = opts || {};
    const fmin = opts.fmin || 50, fmax = opts.fmax || 1000, thr = opts.threshold === undefined ? 0.10 : opts.threshold;
    const N = x.length;
    const none = { freq: null, tau: null, dprime: null, threshold: thr };
    const ret = r => (opts.detail ? r : r.freq);
    if (!(fs > 0) || N < 16) return ret(none);
    const tauMin = Math.max(2, Math.floor(fs / fmax));
    let tauMax = Math.ceil(fs / fmin);
    const W = Math.floor(N / 2);                     // integration window
    tauMax = Math.min(tauMax, N - W - 1);
    if (tauMax <= tauMin + 2) return ret(none);
    // difference function d(tau) = sum_j (x[j] - x[j+tau])^2
    const d = new Float64Array(tauMax + 2);
    for (let tau = 1; tau <= tauMax + 1 && tau < N - W + 1; tau++) {
      let s = 0;
      for (let j = 0; j < W; j++) { const e = x[j] - x[j + tau]; s += e * e; }
      d[tau] = s;
    }
    // cumulative mean normalised difference d'(tau)
    const dp = new Float64Array(tauMax + 2); dp[0] = 1;
    let run = 0;
    for (let tau = 1; tau <= tauMax + 1; tau++) {
      run += d[tau];
      dp[tau] = run > 0 ? d[tau] * tau / run : 1;
    }
    // absolute threshold: first tau >= tauMin below thr, then walk down to the local minimum
    let tau = -1;
    for (let t = tauMin; t <= tauMax; t++) {
      if (dp[t] < thr) { while (t + 1 <= tauMax && dp[t + 1] < dp[t]) t++; tau = t; break; }
    }
    if (tau < 0) return ret(none);
    // parabolic interpolation around the minimum
    let tauF = tau;
    if (tau > 1 && tau < tauMax + 1) {
      const a = dp[tau - 1], b = dp[tau], c = dp[tau + 1], den = a - 2 * b + c;
      if (Math.abs(den) > 1e-12) { const off = 0.5 * (a - c) / den; if (Math.abs(off) <= 1) tauF = tau + off; }
    }
    const f = fs / tauF;
    if (!Number.isFinite(f) || f < fmin * 0.98 || f > fmax * 1.02) return ret(none);
    return ret({ freq: f, tau: tauF, dprime: dp[tau], threshold: thr });
  }

  function hzToMidi(f, a4) { a4 = a4 || 440; return f > 0 && Number.isFinite(f) ? 69 + 12 * Math.log2(f / a4) : NaN; }
  function midiToHz(m, a4) { return (a4 || 440) * Math.pow(2, (m - 69) / 12); }
  // cents of f relative to ref (positive = sharp)
  function centsOff(f, ref) { return f > 0 && ref > 0 && Number.isFinite(f) && Number.isFinite(ref) ? 1200 * Math.log2(f / ref) : NaN; }
  // nearest equal-tempered note: {name, octave, midi, ideal (Hz), cents}
  function hzToNote(f, a4) {
    a4 = a4 || 440;
    const m = hzToMidi(f, a4); if (!Number.isFinite(m)) return null;
    const n = Math.round(m), ideal = midiToHz(n, a4);
    return { name: NOTE_NAMES[((n % 12) + 12) % 12], octave: Math.floor(n / 12) - 1, midi: n, ideal, cents: centsOff(f, ideal) };
  }
  // soft-clip transfer curve for the WaveShaper: y = tanh(drive*x)/tanh(drive) over x in [-1,1]
  function softClipCurve(n, drive) {
    n = n || 2048; drive = drive || 1.5; const c = new Float32Array(n), k = Math.tanh(drive);
    for (let i = 0; i < n; i++) { const x = (i / (n - 1)) * 2 - 1; c[i] = Math.tanh(drive * x) / k; }
    return c;
  }
  function dbToLin(db) { return Math.pow(10, db / 20); }
  function linToDb(v) { return v > 1e-12 ? 20 * Math.log10(v) : -Infinity; }
  function rms(x) { let s = 0; for (let i = 0; i < x.length; i++) s += x[i] * x[i]; return x.length ? Math.sqrt(s / x.length) : 0; }
  // spectrogram colour map, t in [0,1] -> [r,g,b]
  function colormap(t) {
    if (!Number.isFinite(t)) t = 0; t = Math.min(1, Math.max(0, t));
    const stops = [[0, 0, 0, 8], [0.25, 40, 10, 110], [0.5, 170, 30, 110], [0.75, 245, 130, 20], [1, 255, 245, 160]];
    for (let i = 1; i < stops.length; i++) if (t <= stops[i][0]) {
      const a = stops[i - 1], b = stops[i], u = (t - a[0]) / (b[0] - a[0]);
      return [a[1] + (b[1] - a[1]) * u, a[2] + (b[2] - a[2]) * u, a[3] + (b[3] - a[3]) * u];
    }
    return [255, 245, 160];
  }
  FSP.math.audio = { NOTE_NAMES, yin, hzToMidi, midiToHz, centsOff, hzToNote, softClipCurve, dbToLin, linToDb, rms, colormap };

  /* ================= tests ================= */
  FSP.registerTests('audio', t => {
    const M = FSP.math.audio, fs = 48000, N = 4096;
    const gen = fn => { const x = new Float64Array(N); for (let i = 0; i < N; i++) x[i] = fn(i / fs); return x; };
    [110, 440, 880].forEach(f0 => {
      const r = M.yin(gen(tt => Math.sin(2 * Math.PI * f0 * tt)), fs);
      t.check('YIN sine ' + f0 + ' Hz', r !== null && t.rel(r, f0, 0.001), 'got ' + r);
    });
    const harm = M.yin(gen(tt => { let s = 0; for (let k = 1; k <= 5; k++) s += Math.sin(2 * Math.PI * 220 * k * tt) / k; return s; }), fs);
    t.check('YIN harmonic series 220 Hz, no octave error', harm !== null && t.rel(harm, 220, 0.001), 'got ' + harm);
    const g = FSP.prng.gaussian(FSP.prng.mulberry32(12345));
    const noisyX = gen(tt => Math.sin(2 * Math.PI * 440 * tt)); for (let i = 0; i < N; i++) noisyX[i] += 0.22 * g();
    const noisy = M.yin(noisyX, fs);
    t.check('YIN 440 Hz + Gaussian noise sigma 0.22 within 2%', noisy !== null && t.rel(noisy, 440, 0.02), 'got ' + noisy);
    const g2 = FSP.prng.gaussian(FSP.prng.mulberry32(777)), wn = new Float64Array(N); for (let i = 0; i < N; i++) wn[i] = g2();
    t.check('YIN white noise -> null', M.yin(wn, fs) === null);
    t.check('YIN silence -> null', M.yin(new Float64Array(N), fs) === null);
    t.check('YIN outside range (30 Hz) -> null', M.yin(gen(tt => Math.sin(2 * Math.PI * 30 * tt)), fs) === null);
    t.check('YIN detail output', (() => { const r = M.yin(gen(tt => Math.sin(2 * Math.PI * 440 * tt)), fs, { detail: true }); return r.freq > 0 && r.dprime < 0.1; })());
    t.check('cents: 440 Hz on A4 = 0', t.near(M.centsOff(440, 440), 0, 1e-9));
    t.check('cents: 466.16 Hz = +100 (+-0.1)', t.near(M.centsOff(466.16, 440), 100, 0.1), String(M.centsOff(466.16, 440)));
    const n = M.hzToNote(440);
    t.check('hzToNote 440 -> A4, 0 cents', n.name === 'A' && n.octave === 4 && t.near(n.cents, 0, 1e-9) && n.midi === 69);
    const n2 = M.hzToNote(261.6256);
    t.check('hzToNote 261.63 -> C4', n2.name === 'C' && n2.octave === 4 && Math.abs(n2.cents) < 0.01);
    const n3 = M.hzToNote(450);
    t.check('hzToNote 450 -> A4 +38.9 cents', n3.name === 'A' && t.near(n3.cents, 38.906, 0.01), String(n3.cents));
    t.check('hzToNote invalid -> null', M.hzToNote(NaN) === null && M.hzToNote(-5) === null);
    const cv = M.softClipCurve(2049, 1.5);
    t.check('soft clip: odd, bounded, endpoints +-1, monotone', t.near(cv[0], -1, 1e-6) && t.near(cv[2048], 1, 1e-6) && t.near(cv[1024], 0, 1e-6) && cv.every((v, i) => i === 0 || v >= cv[i - 1]));
    t.check('dB conversions', t.near(M.dbToLin(-6.0206), 0.5, 1e-4) && t.near(M.linToDb(0.1), -20, 1e-9) && M.linToDb(0) === -Infinity);
  });

  /* ================= UI ================= */
  const S = { src: 'tone', wave: 'sine', freq: 440, level: -12, dbmin: -100, dbmax: -20, log: 1, a4: 440 };
  let ui = null, ac = null, g = null, raf = 0, running = false, active = false, frame = 0;
  let pitch = { f: null, d: null, note: null, hold: 0 };
  let sg = { cv: null, ctx: null };

  function secureNote() {
    if (typeof window === 'undefined') return '';
    if (!window.isSecureContext || !navigator.mediaDevices || !navigator.mediaDevices.getUserMedia)
      return 'Microphone input needs a secure context (https:// or http://localhost). This page cannot use the mic here; the tone generator and file player still work.';
    return '';
  }
  function msg(text, kind) { if (!ui) return; ui.msg.textContent = text || ''; ui.msg.className = 'msg' + (kind ? ' ' + kind : ''); ui.msg.style.display = text ? '' : 'none'; }

  function ensureCtx() {
    if (ac) return true;
    const AC = window.AudioContext || window.webkitAudioContext;
    if (!AC) { msg('Web Audio is not available in this browser.', 'err'); return false; }
    ac = new AC();
    g = {};
    g.bus = ac.createGain();
    g.an = ac.createAnalyser(); g.an.fftSize = 4096; g.an.smoothingTimeConstant = 0.2;
    g.mon = ac.createGain();                        // monitor gain (0 for mic to avoid feedback)
    g.shaper = ac.createWaveShaper(); g.shaper.curve = M().softClipCurve(2048, 1.5); g.shaper.oversample = '2x';
    g.master = ac.createGain(); g.master.gain.value = 0.8;
    g.bus.connect(g.an); g.bus.connect(g.mon); g.mon.connect(g.shaper); g.shaper.connect(g.master); g.master.connect(ac.destination);
    g.td = new Float32Array(g.an.fftSize); g.fd = new Float32Array(g.an.frequencyBinCount);
    return true;
  }
  function M() { return FSP.math.audio; }

  function stopSources() {
    if (!g) return;
    if (g.osc) { try { g.osc.stop(); } catch (e) { /* already stopped */ } try { g.osc.disconnect(); } catch (e) { /* */ } g.osc = null; }
    if (g.toneGain) { try { g.toneGain.disconnect(); } catch (e) { /* */ } g.toneGain = null; }
    if (g.player) { try { g.player.stop(); } catch (e) { /* */ } try { g.player.disconnect(); } catch (e) { /* */ } g.player = null; }
    if (g.micNode) { try { g.micNode.disconnect(); } catch (e) { /* */ } g.micNode = null; }
    if (g.stream) { g.stream.getTracks().forEach(tr => { try { tr.stop(); } catch (e) { /* */ } }); g.stream = null; }
  }
  function stopAll() {
    running = false; stopSources();
    if (raf) { cancelAnimationFrame(raf); raf = 0; }
    pitch = { f: null, d: null, note: null, hold: 0 };
    if (ui) { ui.start.textContent = 'Start'; ui.start.setAttribute('aria-pressed', 'false'); drawTuner(); updateReadouts(); }
  }

  async function start() {
    if (running) { stopAll(); return; }
    msg('');
    if (S.src === 'mic') {
      const sn = secureNote(); if (sn) { msg(sn, 'warn'); return; }
    }
    if (S.src === 'file' && !(g && g.fileBuf)) { msg('Choose an audio file first (decoded in your browser; nothing is uploaded).', 'warn'); return; }
    if (!ensureCtx()) return;
    try { if (ac.state === 'suspended') await ac.resume(); } catch (e) { /* ignore */ }
    stopSources();
    try {
      if (S.src === 'tone') {
        g.osc = ac.createOscillator(); g.osc.type = S.wave; g.osc.frequency.value = S.freq;
        g.toneGain = ac.createGain(); g.toneGain.gain.value = M().dbToLin(S.level);
        g.osc.connect(g.toneGain); g.toneGain.connect(g.bus); g.osc.start();
        g.mon.gain.value = 1;
      } else if (S.src === 'file') {
        g.player = ac.createBufferSource(); g.player.buffer = g.fileBuf; g.player.loop = true;
        g.player.connect(g.bus); g.player.start(); g.mon.gain.value = 1;
      } else {
        let stream;
        try { stream = await navigator.mediaDevices.getUserMedia({ audio: { echoCancellation: false, noiseSuppression: false, autoGainControl: false } }); }
        catch (e) {
          const nm = e && e.name;
          msg(nm === 'NotAllowedError' || nm === 'SecurityError' ? 'Microphone permission was denied. Allow mic access in the browser address bar, or use the tone generator or a file instead.' : nm === 'NotFoundError' ? 'No microphone was found on this device.' : 'Could not open the microphone (' + (e && e.message || nm || 'unknown error') + ').', 'err');
          return;
        }
        g.stream = stream; g.micNode = ac.createMediaStreamSource(stream); g.micNode.connect(g.bus); g.mon.gain.value = 0;
      }
    } catch (e) { msg('Could not start audio: ' + (e && e.message || e), 'err'); stopSources(); return; }
    running = true; ui.start.textContent = 'Stop'; ui.start.setAttribute('aria-pressed', 'true');
    if (active && !raf) raf = requestAnimationFrame(loop);
  }

  async function loadFile(file) {
    if (!file) return;
    if (!ensureCtx()) return;
    try {
      const ab = await file.arrayBuffer();
      const buf = await new Promise((res, rej) => { const p = ac.decodeAudioData(ab, res, rej); if (p && p.catch) p.catch(rej); });
      g.fileBuf = buf;
      ui.fileInfo.textContent = file.name + ' — ' + buf.duration.toFixed(1) + ' s, ' + buf.numberOfChannels + ' ch, ' + buf.sampleRate + ' Hz';
      msg('');
      if (running && S.src === 'file') { stopAll(); start(); }
    } catch (e) { g.fileBuf = null; ui.fileInfo.textContent = 'No file loaded'; msg('Could not decode that file (' + (e && e.message || 'unsupported format') + ').', 'err'); }
  }

  /* ---------- drawing ---------- */
  function cssVar(name, fb) { try { const v = getComputedStyle(document.documentElement).getPropertyValue(name).trim(); return v || fb; } catch (e) { return fb; } }

  function drawTuner() {
    if (!ui) return;
    const { ctx, w, h } = ui.tuner.prep();
    const fg = cssVar('--text', '#e6edf7'), mut = cssVar('--muted', '#8d9bb3'), acc = cssVar('--accent', '#5ab0ff'), good = cssVar('--good', '#4cd08a'), bad = cssVar('--bad', '#ff6b6b');
    ctx.clearRect(0, 0, w, h);
    const cx = w / 2, cy = h - 14, R = Math.max(10, Math.min(w * 0.45, h - 30));
    const ang = c => -Math.PI / 2 + (Math.max(-50, Math.min(50, c)) / 50) * (Math.PI / 3);
    ctx.lineWidth = 2; ctx.strokeStyle = mut; ctx.fillStyle = mut; ctx.font = '11px sans-serif'; ctx.textAlign = 'center';
    ctx.beginPath(); ctx.arc(cx, cy, R, ang(-50), ang(50)); ctx.globalAlpha = 0.5; ctx.stroke(); ctx.globalAlpha = 1;
    [-50, -25, 0, 25, 50].forEach(c => {
      const a = ang(c), r1 = R - (c === 0 ? 14 : 8);
      ctx.beginPath(); ctx.moveTo(cx + Math.cos(a) * r1, cy + Math.sin(a) * r1); ctx.lineTo(cx + Math.cos(a) * R, cy + Math.sin(a) * R); ctx.stroke();
      ctx.fillText((c > 0 ? '+' : '') + c, cx + Math.cos(a) * (R + 12), cy + Math.sin(a) * (R + 12) + 4);
    });
    const n = pitch.note;
    if (n) {
      const a = ang(n.cents), inTune = Math.abs(n.cents) <= 5;
      ctx.strokeStyle = inTune ? good : (Math.abs(n.cents) > 25 ? bad : acc); ctx.lineWidth = 3;
      ctx.beginPath(); ctx.moveTo(cx, cy); ctx.lineTo(cx + Math.cos(a) * (R - 6), cy + Math.sin(a) * (R - 6)); ctx.stroke();
      ctx.fillStyle = fg; ctx.font = 'bold 26px sans-serif'; ctx.fillText(n.name + n.octave, cx, cy - R * 0.35);
    } else {
      ctx.fillStyle = mut; ctx.font = '14px sans-serif'; ctx.fillText(running ? 'no pitch' : 'stopped', cx, cy - R * 0.35);
    }
    ctx.fillStyle = mut; ctx.beginPath(); ctx.arc(cx, cy, 4, 0, 2 * Math.PI); ctx.fill();
  }

  function yAxisFreq(row, H, fs) {            // row 0 = top = highest frequency
    const nyq = fs / 2, fmaxPlot = Math.min(nyq, 20000), u = 1 - row / Math.max(1, H - 1);
    if (S.log) { const f0 = 20; return f0 * Math.pow(fmaxPlot / f0, u); }
    return u * fmaxPlot;
  }

  function sgResizeCheck(pw, ph) {
    if (!sg.cv) sg.cv = document.createElement('canvas');
    if (sg.cv.width !== pw || sg.cv.height !== ph) { sg.cv.width = pw; sg.cv.height = ph; sg.ctx = sg.cv.getContext('2d'); sg.ctx.fillStyle = '#000008'; sg.ctx.fillRect(0, 0, pw, ph); }
  }
  function drawSpectrogram(advance) {
    const { ctx, w, h, dpr } = ui.spec.prep(), cv = ui.spec.cv, pw = cv.width, ph = cv.height;
    sgResizeCheck(pw, ph);
    if (advance && g && g.an) {
      const cols = Math.max(1, Math.round(dpr)), fs = ac.sampleRate, bins = g.fd.length;
      g.an.getFloatFrequencyData(g.fd);
      sg.ctx.setTransform(1, 0, 0, 1, 0, 0);
      sg.ctx.drawImage(sg.cv, -cols, 0);
      const img = sg.ctx.createImageData(cols, ph), range = Math.max(1, S.dbmax - S.dbmin);
      for (let row = 0; row < ph; row++) {
        const f = yAxisFreq(row, ph, fs), bpos = Math.min(bins - 1.001, Math.max(0, f / (fs / 2) * bins)), b0 = Math.floor(bpos), fr = bpos - b0;
        let db = g.fd[b0] * (1 - fr) + g.fd[b0 + 1] * fr; if (!Number.isFinite(db)) db = -200;
        const c = M().colormap((db - S.dbmin) / range);
        for (let k = 0; k < cols; k++) { const o = (row * cols + k) * 4; img.data[o] = c[0]; img.data[o + 1] = c[1]; img.data[o + 2] = c[2]; img.data[o + 3] = 255; }
      }
      sg.ctx.putImageData(img, pw - cols, 0);
    }
    ctx.setTransform(1, 0, 0, 1, 0, 0); ctx.drawImage(sg.cv, 0, 0);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    // frequency axis labels
    const fs = ac ? ac.sampleRate : 48000, ticks = S.log ? [50, 100, 200, 500, 1000, 2000, 5000, 10000] : [0, 2000, 4000, 8000, 12000, 16000, 20000];
    ctx.font = '10px sans-serif'; ctx.textAlign = 'left'; ctx.lineWidth = 1;
    ticks.forEach(f => {
      if (f > Math.min(fs / 2, 20000)) return;
      const fmaxPlot = Math.min(fs / 2, 20000), u = S.log ? Math.log(f / 20) / Math.log(fmaxPlot / 20) : f / fmaxPlot;
      if (u < 0 || u > 1) return;
      const y = Math.round((1 - u) * (h - 1)) + 0.5;
      ctx.strokeStyle = 'rgba(255,255,255,0.18)'; ctx.beginPath(); ctx.moveTo(0, y); ctx.lineTo(w, y); ctx.stroke();
      ctx.fillStyle = 'rgba(0,0,0,0.55)'; ctx.fillRect(2, Math.max(0, y - 11), 34, 11);
      ctx.fillStyle = '#e8eefc'; ctx.fillText(f >= 1000 ? (f / 1000) + 'k' : String(f), 4, Math.max(9, y - 2));
    });
  }

  /* ---------- main loop ---------- */
  function updateReadouts() {
    const n = pitch.note;
    ui.rf.textContent = pitch.f ? pitch.f.toFixed(2) + ' Hz' : '—';
    ui.rn.textContent = n ? n.name + n.octave : '—';
    ui.rc.textContent = n ? (n.cents >= 0 ? '+' : '') + n.cents.toFixed(1) + ' ¢' : '—';
    ui.rt.textContent = n ? n.ideal.toFixed(2) + ' Hz' : '—';
    ui.rd.textContent = pitch.d !== null && Number.isFinite(pitch.d) ? pitch.d.toFixed(3) : '—';
    updateWorking();
  }
  function updateWorking() {
    const n = pitch.note, L = [];
    const fs = ac ? ac.sampleRate : null;
    L.push('Audio context sample rate fs = ' + (fs ? fs + ' Hz' : '— (not started)'));
    L.push('YIN: d(tau) = sum_j (x[j] - x[j+tau])^2,  d\'(tau) = d(tau) * tau / sum_{k<=tau} d(k)');
    L.push('  search tau in [fs/fmax, fs/fmin] = [' + (fs ? Math.floor(fs / 1000) : '—') + ', ' + (fs ? Math.ceil(fs / 50) : '—') + '] samples (fmin 50 Hz, fmax 1000 Hz)');
    L.push('  first tau with d\'(tau) < 0.10, walk down to local minimum, parabolic interpolation, f0 = fs / tau');
    if (pitch.f && n) {
      L.push('  d\' at minimum = ' + (pitch.d !== null ? pitch.d.toFixed(4) : '—') + '  ->  f0 = ' + pitch.f.toFixed(3) + ' Hz');
      L.push('Note: MIDI = 69 + 12 log2(f/' + S.a4 + ') = ' + M().hzToMidi(pitch.f, S.a4).toFixed(3) + '  -> nearest ' + n.name + n.octave + ' (' + n.midi + '), ideal ' + n.ideal.toFixed(3) + ' Hz');
      L.push('Cents = 1200 log2(f / ideal) = ' + n.cents.toFixed(2));
    } else L.push('  no tau below threshold (or signal below gate): unvoiced -> null');
    L.push('Tone: ' + S.wave + ' at ' + S.freq.toFixed(2) + ' Hz, level ' + S.level.toFixed(1) + ' dBFS -> gain ' + M().dbToLin(S.level).toFixed(4));
    L.push('Soft clip (WaveShaper): y = tanh(1.5 x) / tanh(1.5), then master gain 0.8');
    L.push('Spectrogram: AnalyserNode 4096-pt FFT (' + (fs ? (fs / 4096).toFixed(2) : '—') + ' Hz/bin), colour range ' + S.dbmin + ' to ' + S.dbmax + ' dB, ' + (S.log ? 'log (20 Hz..' + (fs ? Math.min(fs / 2, 20000) : '—') + ')' : 'linear') + ' frequency axis');
    ui.work.set(L);
  }

  function loop() {
    raf = 0; if (!active) return;
    frame++;
    if (running && g && g.an) {
      if (frame % 3 === 0) {
        g.an.getFloatTimeDomainData(g.td);
        if (M().rms(g.td) < 0.003) { pitch.hold++; if (pitch.hold > 8) { pitch.f = null; pitch.note = null; pitch.d = null; } }
        else {
          const r = M().yin(g.td, ac.sampleRate, { fmin: 50, fmax: 1000, threshold: 0.10, detail: true });
          if (r.freq) { pitch.f = pitch.f && Math.abs(M().centsOff(r.freq, pitch.f)) < 100 ? pitch.f * 0.5 + r.freq * 0.5 : r.freq; pitch.d = r.dprime; pitch.note = M().hzToNote(pitch.f, S.a4); pitch.hold = 0; }
          else { pitch.hold++; if (pitch.hold > 8) { pitch.f = null; pitch.note = null; pitch.d = null; } }
        }
        drawTuner(); updateReadouts();
      }
      drawSpectrogram(true);
    }
    if (running) raf = requestAnimationFrame(loop);
  }
  function redrawAll() { if (!ui) return; drawTuner(); drawSpectrogram(false); }

  /* ---------- state ---------- */
  function num(v, lo, hi, def) { const x = parseFloat(v); return Number.isFinite(x) ? Math.min(hi, Math.max(lo, x)) : def; }
  function syncUi() {
    if (!ui) return;
    ui.srcSel.value = S.src; ui.waveSel.value = S.wave; ui.fSl.set(S.freq, true); ui.lSl.set(S.level, true);
    ui.minSl.set(S.dbmin, true); ui.maxSl.set(S.dbmax, true); ui.logCb.checked = !!S.log; ui.a4Sl.set(S.a4, true);
    ui.toneBox.style.display = S.src === 'tone' ? '' : 'none'; ui.fileBox.style.display = S.src === 'file' ? '' : 'none';
  }
  FSP.state.bind('audio', {
    get: () => ({ src: S.src, wave: S.wave, freq: +S.freq.toFixed(3), lvl: S.level, dbmin: S.dbmin, dbmax: S.dbmax, log: S.log ? 1 : 0, a4: S.a4 }),
    set: o => {
      o = o || {};
      if (['tone', 'file', 'mic'].indexOf(o.src) >= 0) S.src = o.src;
      if (['sine', 'square', 'sawtooth'].indexOf(o.wave) >= 0) S.wave = o.wave;
      S.freq = num(o.freq, 20, 8000, S.freq); S.level = num(o.lvl, -60, 0, S.level);
      S.dbmin = num(o.dbmin, -140, -40, S.dbmin); S.dbmax = num(o.dbmax, -60, 0, S.dbmax); if (S.dbmax <= S.dbmin + 5) { S.dbmin = -100; S.dbmax = -20; }
      if (o.log !== undefined) S.log = o.log === '0' || o.log === 0 ? 0 : 1;
      S.a4 = num(o.a4, 400, 480, S.a4);
      syncUi();
    },
  });

  /* ---------- tab ---------- */
  function init(panel) {
    const U = FSP.ui, el = U.el;
    ui = {};
    panel.appendChild(el('p', { class: 'note', text: 'Tuner (YIN), scrolling spectrogram, reference tone generator and a looped file player. Uses its own audio context; everything stops when you leave this tab.' }));
    ui.msg = el('div', { class: 'msg', role: 'status' }); ui.msg.style.display = 'none';
    const sn = secureNote(); if (sn) { ui.msg.textContent = sn; ui.msg.className = 'msg warn'; ui.msg.style.display = ''; }

    const src = U.fieldset(panel, 'Source');
    ui.srcSel = U.select(src, 'Input', [['tone', 'Tone generator'], ['file', 'Audio file (looped)'], ['mic', 'Microphone']], S.src, v => { S.src = v; if (running) stopAll(); syncUi(); FSP.state.touch(); });
    const row = el('div', { class: 'row' }); src.appendChild(row);
    ui.start = U.button(row, 'Start', () => { start(); }); ui.start.setAttribute('aria-pressed', 'false');
    ui.toneBox = el('div'); src.appendChild(ui.toneBox);
    ui.waveSel = U.select(ui.toneBox, 'Waveform', [['sine', 'Sine'], ['square', 'Square'], ['sawtooth', 'Sawtooth']], S.wave, v => { S.wave = v; if (g && g.osc) g.osc.type = v; FSP.state.touch(); updateWorking(); });
    ui.fSl = U.slider(ui.toneBox, { label: 'Frequency', min: 20, max: 8000, value: S.freq, log: true, unit: 'Hz', onInput: v => { S.freq = v; if (g && g.osc) g.osc.frequency.setTargetAtTime(v, ac.currentTime, 0.01); FSP.state.touch(); updateWorking(); } });
    ui.lSl = U.slider(ui.toneBox, { label: 'Level', min: -60, max: 0, step: 0.5, value: S.level, unit: 'dBFS', onInput: v => { S.level = v; if (g && g.toneGain) g.toneGain.gain.setTargetAtTime(M().dbToLin(v), ac.currentTime, 0.01); FSP.state.touch(); updateWorking(); } });
    ui.fileBox = el('div'); src.appendChild(ui.fileBox);
    const fi = el('input', { type: 'file', accept: 'audio/*', 'aria-label': 'Audio file' }); fi.addEventListener('change', () => loadFile(fi.files && fi.files[0]));
    ui.fileInfo = el('span', { class: 'note', text: 'No file loaded' });
    ui.fileBox.appendChild(el('div', { class: 'row' }, fi, ui.fileInfo));
    src.appendChild(ui.msg);

    const tun = U.fieldset(panel, 'Tuner (YIN, 50 Hz to 1 kHz)');
    const wrap = el('div', { class: 'canvas-wrap' }); tun.appendChild(wrap);
    ui.tuner = U.canvas(wrap, { height: 150 }); ui.tuner.onResize(() => drawTuner());
    const hud = el('div', { class: 'hud' }); tun.appendChild(hud);
    const mk = (label) => { const v = el('span', { class: 'mono', text: '—' }); hud.appendChild(el('div', null, el('span', { text: label }), v)); return v; };
    ui.rf = mk('Frequency'); ui.rn = mk('Note'); ui.rc = mk('Cents'); ui.rt = mk('Target'); ui.rd = mk('YIN d\'');
    ui.a4Sl = U.slider(tun, { label: 'A4 reference', min: 400, max: 480, step: 1, value: S.a4, unit: 'Hz', onInput: v => { S.a4 = v; FSP.state.touch(); } });

    const sp = U.fieldset(panel, 'Spectrogram');
    const w2 = el('div', { class: 'canvas-wrap' }); sp.appendChild(w2);
    ui.spec = U.canvas(w2, { height: 260 }); ui.spec.onResize(() => drawSpectrogram(false));
    ui.minSl = U.slider(sp, { label: 'Floor', min: -140, max: -40, step: 1, value: S.dbmin, unit: 'dB', onInput: v => { S.dbmin = Math.min(v, S.dbmax - 5); FSP.state.touch(); } });
    ui.maxSl = U.slider(sp, { label: 'Ceiling', min: -60, max: 0, step: 1, value: S.dbmax, unit: 'dB', onInput: v => { S.dbmax = Math.max(v, S.dbmin + 5); FSP.state.touch(); } });
    const lr = el('div', { class: 'row' }); sp.appendChild(lr);
    ui.logCb = el('input', { type: 'checkbox', id: 'audio-log' }); ui.logCb.checked = !!S.log;
    ui.logCb.addEventListener('change', () => { S.log = ui.logCb.checked ? 1 : 0; sg.cv = null; FSP.state.touch(); drawSpectrogram(false); });
    lr.appendChild(ui.logCb); lr.appendChild(el('label', { for: 'audio-log', text: 'Logarithmic frequency axis' }));

    ui.work = U.working(panel);
    syncUi(); updateReadouts();
  }

  FSP.registerTab({
    id: 'audio', title: 'Audio Tools',
    init,
    activate() { active = true; if (ui) { syncUi(); redrawAll(); updateReadouts(); if (running && !raf) raf = requestAnimationFrame(loop); } },
    deactivate() { active = false; stopAll(); if (ac && ac.state === 'running') { try { ac.suspend(); } catch (e) { /* ignore */ } } },
  });
})();
