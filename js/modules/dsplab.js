/* Phase 4: Convolution / DFT / Windows lab. Math in FSP.math.dsplab (pure); UI built lazily in init(). */
(function () {
  'use strict';
  const F = FSP.fft;

  /* ================= pure math ================= */
  const M = {};
  M.WINDOWS = [['rect', 'Rectangular'], ['hann', 'Hann'], ['hamming', 'Hamming'], ['blackman', 'Blackman']];
  // Periodic (DFT-even) windows: w[n] = ... 2πn/N, n = 0..N-1
  M.window = function (type, N) {
    const w = new Float64Array(N), a = 2 * Math.PI / N;
    for (let n = 0; n < N; n++) {
      const c = Math.cos(a * n);
      if (type === 'hann') w[n] = 0.5 - 0.5 * c;
      else if (type === 'hamming') w[n] = 0.54 - 0.46 * c;
      else if (type === 'blackman') w[n] = 0.42 - 0.5 * c + 0.08 * Math.cos(2 * a * n);
      else w[n] = 1;
    }
    return w;
  };
  M.coherentGain = function (w) { let s = 0; for (let i = 0; i < w.length; i++) s += w[i]; return s / w.length; };
  M.enbw = function (w) { let s = 0, s2 = 0; for (let i = 0; i < w.length; i++) { s += w[i]; s2 += w[i] * w[i]; } return w.length * s2 / (s * s); };
  // Zero-padded (length N*pad, pad power of two) magnitude of x (re, optional im) -> {re, im, mag}
  M.paddedSpectrum = function (x, pad) {
    const L = x.length * pad, re = new Float64Array(L), im = new Float64Array(L);
    re.set(x); F.fft(re, im, false);
    const mag = new Float64Array(L); for (let i = 0; i < L; i++) mag[i] = Math.hypot(re[i], im[i]);
    return { re, im, mag };
  };
  const toDb = (m, ref) => { const r = m / ref; return r > 1e-12 ? 20 * Math.log10(r) : -240; };
  // Window magnitude response in dB relative to its peak (DC), bins 0..N*pad/2 at 1/pad-bin spacing.
  M.windowSpectrumDb = function (w, pad) {
    const sp = M.paddedSpectrum(w, pad || 64), ref = sp.mag[0], n = sp.mag.length / 2 + 1, db = new Float64Array(n);
    for (let k = 0; k < n; k++) db[k] = toDb(sp.mag[k], ref);
    return db;
  };
  // Peak sidelobe: highest value after the first null (first local minimum of the magnitude response), dB re peak.
  M.peakSidelobeDb = function (w, pad) {
    pad = pad || 64; const sp = M.paddedSpectrum(w, pad), mag = sp.mag, half = mag.length / 2;
    let k = 1; while (k < half && mag[k + 1] <= mag[k]) k++;
    let best = 0; for (let j = k; j <= half; j++) if (mag[j] > best) best = mag[j];
    return toDb(best, mag[0]);
  };
  M.windowStats = function (type, N) {
    const w = M.window(type, N);
    return { cg: M.coherentGain(w), enbw: M.enbw(w), sidelobeDb: M.peakSidelobeDb(w, 64) };
  };
  // Mainlobe half-width (bins) = distance to first null (first local min)
  M.firstNullBins = function (w, pad) {
    pad = pad || 64; const mag = M.paddedSpectrum(w, pad).mag, half = mag.length / 2;
    let k = 1; while (k < half && mag[k + 1] <= mag[k]) k++; return k / pad;
  };
  // Sinusoid with frequency fb in bins (cycles per N samples): A sin(2π fb n/N + ph)
  M.sinusoid = function (N, fb, A, ph) { const x = new Float64Array(N); for (let n = 0; n < N; n++) x[n] = (A === undefined ? 1 : A) * Math.sin(2 * Math.PI * fb * n / N + (ph || 0)); return x; };
  // Leakage spectrum of a sinusoid under a window, dB re a bin-centred unwindowed sinusoid of same amplitude (0 dB = N*A/2 *CG)
  M.leakageDb = function (N, fb, type, pad) {
    const w = M.window(type, N), x = M.sinusoid(N, fb, 1, 0), y = new Float64Array(N);
    for (let n = 0; n < N; n++) y[n] = x[n] * w[n];
    const sp = M.paddedSpectrum(y, pad || 8), ref = N / 2, half = sp.mag.length / 2 + 1, db = new Float64Array(half);
    for (let k = 0; k < half; k++) db[k] = toDb(sp.mag[k], ref);
    return { db, pad: pad || 8, w };
  };
  // Worst leakage level (dB re peak; binsOnly = only true DFT bins, as an analyzer shows) at bins at least `away` bins from fb
  M.farLeakageDb = function (spec, fb, away, binsOnly) {
    let m = 0, pk = 0; const pad = spec.pad;
    for (let k = 0; k < spec.db.length; k++) { if (binsOnly && k % pad) continue; const v = Math.pow(10, spec.db[k] / 20), b = k / pad; if (v > pk) pk = v; if (Math.abs(b - fb) >= away && v > m) m = v; }
    return m > 0 && pk > 0 ? 20 * Math.log10(m / pk) : -240;
  };
  // General circular convolution of length N (inputs zero-padded or time-aliased modulo N)
  M.circularConv = function (a, b, N) {
    const A = new Float64Array(N), B = new Float64Array(N), y = new Float64Array(N);
    for (let i = 0; i < a.length; i++) A[i % N] += a[i];
    for (let i = 0; i < b.length; i++) B[i % N] += b[i];
    for (let n = 0; n < N; n++) { let s = 0; for (let k = 0; k < N; k++) s += A[k] * B[((n - k) % N + N) % N]; y[n] = s; }
    return y;
  };
  // Parse "1, 2 3;4" -> numbers, or null if junk / empty / too long
  M.parseSeq = function (s, maxLen) {
    if (typeof s !== 'string') return null;
    const parts = s.trim().split(/[\s,;]+/).filter(Boolean); if (!parts.length || parts.length > (maxLen || 16)) return null;
    const out = []; for (const p of parts) { const v = Number(p); if (!Number.isFinite(v)) return null; out.push(v); }
    return out;
  };
  // Terms of y[n] = sum_k x[k] h[(n-k)] for display; circular uses modulo N. Returns [{k, xk, hk, p}] for nonzero overlap.
  M.convTerms = function (x, h, n, circ, N) {
    const out = [];
    for (let k = 0; k < (circ ? N : x.length); k++) {
      const xk = k < x.length ? x[k] : 0; let j = n - k; if (circ) j = ((j % N) + N) % N;
      const hk = (j >= 0 && j < h.length) ? h[j] : 0;
      if (circ ? (k < x.length || false) || true : true) out.push({ k, j, xk, hk, p: xk * hk });
    }
    return out;
  };
  FSP.math.dsplab = M;

  /* ================= tests ================= */
  FSP.registerTests('dsplab', t => {
    const R = FSP.prng.mulberry32(12345), G = FSP.prng.gaussian(R);
    const N = 64, xr = new Float64Array(N), xi = new Float64Array(N); for (let i = 0; i < N; i++) { xr[i] = G(); xi[i] = G(); }
    const d = F.dft(xr, xi), fr = Float64Array.from(xr), fi = Float64Array.from(xi); F.fft(fr, fi, false);
    let e = 0; for (let i = 0; i < N; i++) e = Math.max(e, Math.abs(fr[i] - d.re[i]), Math.abs(fi[i] - d.im[i]));
    t.check('FFT vs direct DFT max error < 1e-9 (N=64, complex random)', e < 1e-9, 'err=' + e.toExponential(2));
    const rr = Float64Array.from(fr), ri = Float64Array.from(fi); F.fft(rr, ri, true);
    e = 0; for (let i = 0; i < N; i++) e = Math.max(e, Math.abs(rr[i] - xr[i]), Math.abs(ri[i] - xi[i]));
    t.check('IFFT(FFT(x)) = x within 1e-12', e < 1e-12, 'err=' + e.toExponential(2));
    let et = 0, ef = 0; for (let i = 0; i < N; i++) { et += xr[i] * xr[i] + xi[i] * xi[i]; ef += fr[i] * fr[i] + fi[i] * fi[i]; }
    t.check('Parseval: sum|x|^2 = (1/N) sum|X|^2 within 1e-9', Math.abs(et - ef / N) < 1e-9, 'diff=' + Math.abs(et - ef / N).toExponential(2));
    const y = F.convolveDirect([1, 2, 3], [0, 1, 0.5]), exp = [0, 1, 2.5, 4, 1.5];
    t.check('convolve [1,2,3]*[0,1,0.5] = [0,1,2.5,4,1.5]', y.length === 5 && exp.every((v, i) => t.near(y[i], v, 1e-12)), Array.from(y).join(','));
    const a = Array.from({ length: 37 }, () => G()), b = Array.from({ length: 21 }, () => G());
    const cd = F.convolveDirect(a, b), cf = F.convolveFFT(a, b); e = 0; for (let i = 0; i < cd.length; i++) e = Math.max(e, Math.abs(cd[i] - cf[i]));
    t.check('FFT linear convolution equals direct within 1e-9', cd.length === cf.length && e < 1e-9, 'err=' + e.toExponential(2));
    // circular conv == IFFT(FFT*FFT) (N=16)
    const c1 = Array.from({ length: 16 }, () => G()), c2 = Array.from({ length: 16 }, () => G());
    const cc = M.circularConv(c1, c2, 16), ar = Float64Array.from(c1), ai = new Float64Array(16), br = Float64Array.from(c2), bi = new Float64Array(16);
    F.fft(ar, ai, false); F.fft(br, bi, false);
    for (let k = 0; k < 16; k++) { const r = ar[k] * br[k] - ai[k] * bi[k], i = ar[k] * bi[k] + ai[k] * br[k]; ar[k] = r; ai[k] = i; }
    F.fft(ar, ai, true); e = 0; for (let i = 0; i < 16; i++) e = Math.max(e, Math.abs(ar[i] - cc[i]));
    t.check('circular conv = IFFT(FFT(a)·FFT(b)) within 1e-9', e < 1e-9, 'err=' + e.toExponential(2));
    // circular N=3 of [1,2,3]*[0,1,0.5] aliases linear [0,1,2.5,4,1.5] -> [4,2.5,4]
    const c3 = M.circularConv([1, 2, 3], [0, 1, 0.5], 3);
    t.check('circular N=3 aliases the linear result to [4,2.5,2.5]', [4, 2.5, 2.5].every((v, i) => t.near(c3[i], v, 1e-12)), Array.from(c3).join(','));
    // windows, N=1024
    const spec = { rect: [-13.26, 1.0, 1.0], hann: [-31.47, 1.5, 0.5], hamming: [-42.67, 1.363, 0.54], blackman: [-58.11, 1.727, 0.42] };
    Object.keys(spec).forEach(k => {
      const s = M.windowStats(k, 1024), x = spec[k];
      t.check(k + ' peak sidelobe ' + x[0] + ' dB (±0.05)', t.near(s.sidelobeDb, x[0], 0.05), 'got ' + s.sidelobeDb.toFixed(3));
      t.check(k + ' ENBW ' + x[1] + ' (±0.005)', t.near(s.enbw, x[1], 0.005), 'got ' + s.enbw.toFixed(4));
      t.check(k + ' coherent gain ' + x[2] + ' (±0.005)', t.near(s.cg, x[2], 0.005), 'got ' + s.cg.toFixed(4));
    });
    // leakage: bin-centred sinusoid under rect has ~no far leakage; off-bin does; Hann reduces it
    const lc = M.leakageDb(256, 20, 'rect', 8), lo = M.leakageDb(256, 20.5, 'rect', 8), lh = M.leakageDb(256, 20.5, 'hann', 8);
    const fc = M.farLeakageDb(lc, 20, 10, true), fo = M.farLeakageDb(lo, 20.5, 10, true), fh = M.farLeakageDb(lh, 20.5, 10, true);
    t.check('on-bin sinusoid (rect) has no leakage at other DFT bins (< -200 dB)', fc < -200, fc.toFixed(1));
    t.check('off-bin sinusoid (rect) leaks into DFT bins >=10 away (-50..-20 dB)', fo > -50 && fo < -20, fo.toFixed(1));
    t.check('Hann cuts off-bin far leakage by > 30 dB vs rect', fo - fh > 30, (fo - fh).toFixed(1));
    t.check('on-bin peak is 0 dB (A·N/2 reference)', t.near(Math.max.apply(null, lc.db), 0, 1e-9));
    t.check('parseSeq accepts "1, 2 3" and rejects junk', M.parseSeq('1, 2 3').length === 3 && M.parseSeq('1,a') === null && M.parseSeq('') === null);
  });

  /* ================= UI ================= */
  const S = { x: '1, 2, 3', h: '0, 1, 0.5', mode: 'linear', nc: 5, step: 0, win: 'hann', wn: 256, span: '8', ln: 256, lf: 20.5, lw: 'hann', lspan: '16' };
  let ui = null, raf = 0, playing = false, lastT = 0;

  function css(name, fb) { try { return getComputedStyle(document.documentElement).getPropertyValue(name).trim() || fb; } catch (e) { return fb; } }
  function colors() { return { bg: css('--panel2', '#0f141f'), grid: css('--border', '#243047'), text: css('--text', '#e6edf3'), muted: css('--muted', '#8b98ab'), a: css('--accent', '#4cc9f0'), b: css('--accent2', '#f9a03f'), c: css('--pink', '#f72585'), ok: css('--ok', '#2ecc71') }; }
  const f3 = v => Number.isFinite(v) ? String(+v.toFixed(3)) : '—';
  const fdb = v => Number.isFinite(v) ? v.toFixed(2) + ' dB' : '—';
  function niceStep(r) { const p = Math.pow(10, Math.floor(Math.log10(r))), m = r / p; return (m < 1.5 ? 1 : m < 3.5 ? 2 : m < 7.5 ? 5 : 10) * p; }

  // Generic plot frame; returns mapping helpers.
  function frame(g, w, h, xr, yr, o) {
    const C = colors(), L = 44, R = 10, T = 8, B = o.xlabel ? 30 : 20, pw = w - L - R, ph = h - T - B;
    const X = v => L + (v - xr[0]) / (xr[1] - xr[0]) * pw, Y = v => T + (1 - (v - yr[0]) / (yr[1] - yr[0])) * ph;
    g.fillStyle = C.bg; g.fillRect(0, 0, w, h); g.font = '11px ' + css('--mono', 'monospace'); g.lineWidth = 1;
    const ys = niceStep((yr[1] - yr[0]) / 5), xs = niceStep((xr[1] - xr[0]) / Math.max(2, Math.floor(pw / 60)));
    g.textAlign = 'right'; g.textBaseline = 'middle';
    for (let v = Math.ceil(yr[0] / ys) * ys; v <= yr[1] + 1e-9; v += ys) { const y = Y(v); g.strokeStyle = C.grid; g.beginPath(); g.moveTo(L, y); g.lineTo(w - R, y); g.stroke(); g.fillStyle = C.muted; g.fillText(f3(v), L - 4, y); }
    g.textAlign = 'center'; g.textBaseline = 'top';
    for (let v = Math.ceil(xr[0] / xs) * xs; v <= xr[1] + 1e-9; v += xs) { const x = X(v); g.strokeStyle = C.grid; g.beginPath(); g.moveTo(x, T); g.lineTo(x, T + ph); g.stroke(); g.fillStyle = C.muted; g.fillText(f3(v), x, T + ph + 3); }
    if (o.xlabel) { g.fillStyle = C.muted; g.fillText(o.xlabel, L + pw / 2, h - 13); }
    g.strokeStyle = C.muted; g.strokeRect(L, T, pw, ph);
    g.save(); g.beginPath(); g.rect(L, T, pw, ph); g.clip();
    return { X, Y, C, end() { g.restore(); } };
  }
  function stems(g, F2, xs, ys, color, off, r) {
    g.strokeStyle = color; g.fillStyle = color; g.lineWidth = 1.5;
    xs.forEach((x, i) => { const px = F2.X(x + (off || 0)), py = F2.Y(ys[i]), p0 = F2.Y(0); g.beginPath(); g.moveTo(px, p0); g.lineTo(px, py); g.stroke(); g.beginPath(); g.arc(px, py, r || 3, 0, 6.2832); g.fill(); });
    g.lineWidth = 1;
  }
  function curve(g, F2, xs, ys, color, lw) {
    g.strokeStyle = color; g.lineWidth = lw || 1.4; g.beginPath(); let first = true;
    for (let i = 0; i < xs.length; i++) { const x = F2.X(xs[i]), y = F2.Y(ys[i]); if (first) { g.moveTo(x, y); first = false; } else g.lineTo(x, y); }
    g.stroke(); g.lineWidth = 1;
  }
  function legend(g, items, x, y) { g.textAlign = 'left'; g.textBaseline = 'middle'; items.forEach(([c, s], i) => { g.fillStyle = c; g.fillRect(x, y + i * 14 - 4, 10, 8); g.fillStyle = colors().text; g.fillText(s, x + 14, y + i * 14); }); }

  /* ---- convolution ---- */
  function convData() {
    const x = M.parseSeq(S.x), h = M.parseSeq(S.h);
    if (!x || !h) return { err: 'Enter 1–16 numbers separated by commas or spaces in both sequences.' };
    const circ = S.mode === 'circular', lin = x.length + h.length - 1;
    const N = circ ? Math.max(1, Math.min(32, Math.round(S.nc))) : lin;
    const y = circ ? M.circularConv(x, h, N) : Array.from(F.convolveDirect(x, h));
    return { x, h, circ, N, lin, y, linY: Array.from(F.convolveDirect(x, h)) };
  }
  function drawConv() {
    if (!ui) return;
    const D = convData();
    if (D.err) { ui.cmsg.textContent = D.err; ui.cmsg.hidden = false; return; }
    ui.cmsg.hidden = true;
    const n = Math.max(0, Math.min(D.N - 1, Math.round(S.step))); ui.cstep.max = D.N - 1; ui.cstep.value = n;
    const terms = M.convTerms(D.x, D.h, n, D.circ, D.N), sum = terms.reduce((a, q) => a + q.p, 0);
    // top plot: x[k] and flipped/shifted h[n-k]
    let c = ui.cv1.prep(); let g = c.ctx;
    const span = D.circ ? D.N : D.lin, lo = D.circ ? -0.7 : -Math.max(1, D.h.length - 1) - 0.7, hi = D.circ ? D.N - 0.3 : span + 0.7;
    const allv = D.x.concat(D.h).map(Math.abs), m = Math.max(1e-9, Math.max.apply(null, allv)) * 1.2;
    let fr = frame(g, c.w, c.h, [lo, hi], [-m, m], { xlabel: 'k' });
    const xs = D.x.map((_, i) => i); stems(g, fr, xs, D.x, fr.C.a, 0, 3.5);
    const hk = [], hv = [];
    for (let k = lo | 0; k <= Math.ceil(hi); k++) { const j = D.circ ? (((n - k) % D.N) + D.N) % D.N : n - k; if (D.circ && (k < 0 || k >= D.N)) continue; if (j >= 0 && j < D.h.length) { hk.push(k); hv.push(D.h[j]); } }
    stems(g, fr, hk, hv, fr.C.b, 0.12, 3);
        g.strokeStyle = fr.C.c; g.beginPath(); g.moveTo(fr.X(n), 0); g.lineTo(fr.X(n), c.h); g.stroke();
    fr.end(); legend(g, [[fr.C.a, 'x[k]'], [fr.C.b, 'h[' + n + '−k] (flipped, slid)']], 52, 16);
    // bottom plot: y[n] with current highlighted
    c = ui.cv2.prep(); g = c.ctx;
    const ym = Math.max(1e-9, Math.max.apply(null, D.y.map(Math.abs))) * 1.2, ymin = Math.min(0, Math.min.apply(null, D.y) * 1.2), ymax = Math.max(0, Math.max.apply(null, D.y) * 1.2);
    fr = frame(g, c.w, c.h, [-0.7, D.N - 0.3], [Math.min(ymin, -ym * 0.05), Math.max(ymax, ym * 0.05)], { xlabel: 'n' });
    const ns = D.y.map((_, i) => i);
    stems(g, fr, ns.filter(i => i > n), D.y.filter((_, i) => i > n).map(() => 0), fr.C.grid, 0, 2);
    stems(g, fr, ns.filter(i => i < n), D.y.filter((_, i) => i < n), fr.C.ok, 0, 3.5);
    stems(g, fr, [n], [D.y[n]], fr.C.c, 0, 5);
    if (D.circ) { const lyl = D.linY; stems(g, fr, lyl.map((_, i) => i).filter(i => i < D.N || true), lyl, fr.C.muted, 0.18, 2); }
    fr.end(); legend(g, D.circ ? [[fr.C.ok, 'y[n] done'], [fr.C.c, 'current'], [fr.C.muted, 'linear result (offset)']] : [[fr.C.ok, 'y[n] done'], [fr.C.c, 'current']], 52, 16);
    // readouts
    const kt = terms.filter(q => q.xk !== 0 || q.hk !== 0).map(q => 'x[' + q.k + ']·h[' + q.j + '] = ' + f3(q.xk) + '·' + f3(q.hk) + ' = ' + f3(q.p));
    ui.cread.textContent = 'n = ' + n + ':  y[' + n + '] = ' + (kt.length ? f3(sum) : '0') + '   (' + (D.circ ? 'circular, N=' + D.N : 'linear, length ' + D.lin) + ')';
    const lines = ['Mode: ' + (D.circ ? 'circular, N = ' + D.N + ' (indices mod N)' : 'linear (length ' + D.x.length + ' + ' + D.h.length + ' − 1 = ' + D.lin + ')'),
      'x = [' + D.x.join(', ') + ']', 'h = [' + D.h.join(', ') + ']', '', 'y[n] = Σ_k x[k]·h[' + (D.circ ? '(n−k) mod N' : 'n−k') + ']', 'Terms for n = ' + n + ':'].concat(kt.map(s => '  ' + s), ['  Σ = ' + f3(sum)], ['', 'y = [' + D.y.map(f3).join(', ') + ']']);
    if (D.circ) {
      lines.push('Linear result = [' + D.linY.map(f3).join(', ') + ']');
      lines.push(D.N < D.lin ? 'N < ' + D.lin + ': the tail wraps around and adds onto the start (time-domain aliasing).' : 'N ≥ ' + D.lin + ': no aliasing, circular result equals linear result (zero padded).');
    }
    lines.push('', 'Check via FFT: ' + (function () { const L = D.circ ? D.N : F.nextPow2(D.lin); return F.isPow2(L) || D.circ ? 'IFFT(FFT(x)·FFT(h)) ' + (D.circ && !F.isPow2(L) ? 'needs power-of-two N here; direct sum shown above' : 'matches the direct sum') : ''; })());
    ui.cwork.set(lines);
  }
  function stopPlay() { playing = false; if (raf) cancelAnimationFrame(raf); raf = 0; if (ui) ui.play.textContent = 'Play'; }
  function tick(ts) {
    if (!playing) return;
    if (ts - lastT > 650) { lastT = ts; const D = convData(); if (D.err) { stopPlay(); return; } S.step = (Math.round(S.step) + 1) % D.N; drawConv(); }
    raf = requestAnimationFrame(tick);
  }

  /* ---- windows ---- */
  function drawWin() {
    if (!ui) return;
    const N = S.wn, w = M.window(S.win, N), pad = 64, db = M.windowSpectrumDb(w, pad), st = M.windowStats(S.win, N);
    let c = ui.wv1.prep(), g = c.ctx;
    let fr = frame(g, c.w, c.h, [0, N - 1], [-0.05, 1.1], { xlabel: 'n' });
    const xs = Array.from(w, (_, i) => i); curve(g, fr, xs, Array.from(w), fr.C.a, 1.8); fr.end();
    c = ui.wv2.prep(); g = c.ctx;
    const maxB = S.span === 'full' ? N / 2 : Math.min(+S.span, N / 2), K = Math.floor(maxB * pad);
    fr = frame(g, c.w, c.h, [0, maxB], [-140, 5], { xlabel: 'frequency (bins, zero-padded ×64)' });
    const bx = [], by = []; const stride = Math.max(1, Math.floor(K / 1500)); for (let k = 0; k <= K; k += stride) { bx.push(k / pad); by.push(Math.max(-140, db[k])); }
    curve(g, fr, bx, by, fr.C.b, 1.4);
    g.strokeStyle = fr.C.c; g.setLineDash([5, 4]); g.beginPath(); g.moveTo(fr.X(0), fr.Y(st.sidelobeDb)); g.lineTo(fr.X(maxB), fr.Y(st.sidelobeDb)); g.stroke(); g.setLineDash([]);
    fr.end(); legend(g, [[fr.C.b, '|W(f)| dB re peak'], [fr.C.c, 'peak sidelobe']], 52, 16);
    ui.wsl.textContent = fdb(st.sidelobeDb); ui.wenbw.textContent = f3(st.enbw) + ' bins'; ui.wcg.textContent = f3(st.cg);
    ui.wnull.textContent = f3(M.firstNullBins(w, pad)) + ' bins';
    ui.wwork.set(['Periodic window, N = ' + N + ':  w[n] = ' + ({ rect: '1', hann: '0.5 − 0.5·cos(2πn/N)', hamming: '0.54 − 0.46·cos(2πn/N)', blackman: '0.42 − 0.5·cos(2πn/N) + 0.08·cos(4πn/N)' })[S.win],
      'Coherent gain  CG = Σw / N = ' + f3(st.cg), 'ENBW = N·Σw² / (Σw)² = ' + f3(st.enbw) + ' bins',
      'Spectrum: zero-pad to ' + N * pad + ' points (×64), |FFT|, dB relative to the DC peak.',
      'Peak sidelobe = highest magnitude beyond the first null (first local minimum), here ' + fdb(st.sidelobeDb) + '.',
      'Mainlobe half-width (first null) = ' + f3(M.firstNullBins(w, pad)) + ' bins.',
      'Reference values at N = 1024: Rect −13.26 dB/1.000/1.000, Hann −31.47/1.500/0.500, Hamming −42.67/1.363/0.540, Blackman −58.11/1.727/0.420.']);
  }

  /* ---- leakage ---- */
  function drawLeak() {
    if (!ui) return;
    const N = S.ln, fb = S.lf, A = M.leakageDb(N, fb, 'rect', 8), B = M.leakageDb(N, fb, S.lw, 8);
    const half = +S.lspan, lo = Math.max(0, fb - half), hi = Math.min(N / 2, fb + half);
    const c = ui.lv.prep(), g = c.ctx, fr = frame(g, c.w, c.h, [lo, hi], [-140, 5], { xlabel: 'frequency (bins)' });
    const pad = A.pad, ka = Math.floor(lo * pad), kb = Math.min(A.db.length - 1, Math.ceil(hi * pad));
    const ix = [], ra = [], rb = []; for (let k = ka; k <= kb; k++) { ix.push(k / pad); ra.push(Math.max(-140, A.db[k])); rb.push(Math.max(-140, B.db[k])); }
    curve(g, fr, ix, ra, fr.C.b, 1.3); curve(g, fr, ix, rb, fr.C.a, 1.6);
    g.strokeStyle = fr.C.c; g.setLineDash([3, 4]); g.beginPath(); g.moveTo(fr.X(fb), 0); g.lineTo(fr.X(fb), c.h); g.stroke(); g.setLineDash([]);
    fr.end(); legend(g, [[fr.C.b, 'rectangular (no window)'], [fr.C.a, M.WINDOWS.find(q => q[0] === S.lw)[1]], [fr.C.c, 'true frequency']], 52, 16);
    const away = Math.min(10, N / 8), fa = M.farLeakageDb(A, fb, away, true), fw = M.farLeakageDb(B, fb, away, true);
    const frac = fb - Math.floor(fb), onbin = Math.abs(frac) < 1e-9;
    const pkA = Math.max.apply(null, A.db), pkB = Math.max.apply(null, B.db);
    ui.lra.textContent = fdb(fa); ui.lrb.textContent = fdb(fw); ui.lpk.textContent = fdb(pkA) + ' / ' + fdb(pkB);
    ui.lbin.textContent = onbin ? 'bin-centred (no leakage with rect)' : 'off-bin by ' + f3(Math.min(frac, 1 - frac)) + ' bin → leakage';
    ui.lwork.set(['x[n] = sin(2π·f·n/N), N = ' + N + ', f = ' + f3(fb) + ' bins (' + (onbin ? 'integer: periodic in the window' : 'non-integer: the record is not an integer number of cycles') + ').',
      'Windowed signal y[n] = w[n]·x[n]; zero-pad ×8, |FFT|, dB re N/2 (the peak of a bin-centred unit sinusoid).',
      'Leakage at DFT bins ≥ ' + away + ' bins from f:  rect ' + fdb(fa) + ',  ' + S.lw + ' ' + fdb(fw) + ' (re each curve\'s own peak).',
      'Window peak level = 20·log10(CG) = ' + fdb(20 * Math.log10(M.coherentGain(B.w))) + ' relative to rect; a window trades mainlobe width for lower sidelobes.']);
  }

  function redrawAll() { drawConv(); drawWin(); drawLeak(); }
  function bindState() {
    FSP.state.bind('dsplab', {
      get: () => ({ x: S.x, h: S.h, mode: S.mode, nc: S.nc, step: S.step, win: S.win, wn: S.wn, span: S.span, ln: S.ln, lf: S.lf, lw: S.lw, lspan: S.lspan }),
      set: o => {
        o = o || {};
        const seq = (k) => { if (typeof o[k] === 'string' && M.parseSeq(o[k])) S[k] = o[k]; };
        seq('x'); seq('h');
        if (o.mode === 'linear' || o.mode === 'circular') S.mode = o.mode;
        const num = (k, lo, hi) => { const v = parseFloat(o[k]); if (Number.isFinite(v)) S[k] = Math.min(hi, Math.max(lo, v)); };
        num('nc', 1, 32); num('step', 0, 64); num('lf', 0, 128);
        const names = M.WINDOWS.map(q => q[0]);
        if (names.indexOf(o.win) >= 0) S.win = o.win; if (names.indexOf(o.lw) >= 0) S.lw = o.lw;
        const pick = (k, list) => { const v = parseInt(o[k], 10); if (list.indexOf(v) >= 0) S[k] = v; };
        pick('wn', [16, 32, 64, 128, 256, 512, 1024]); pick('ln', [64, 128, 256, 512]);
        if (['8', '32', 'full'].indexOf(o.span) >= 0) S.span = o.span;
        if (['4', '16', '64'].indexOf(o.lspan) >= 0) S.lspan = o.lspan;
        S.lf = Math.min(S.lf, S.ln / 2);
        if (ui) syncUI();
      },
    });
  }
  function syncUI() {
    ui.xin.value = S.x; ui.hin.value = S.h; ui.mode.value = S.mode; ui.nc.value = S.nc; ui.ncRow.hidden = S.mode !== 'circular';
    ui.wsel.value = S.win; ui.wnsel.value = String(S.wn); ui.spansel.value = S.span;
    ui.lnsel.value = String(S.ln); ui.lwsel.value = S.lw; ui.lspansel.value = S.lspan; ui.lfs.set(S.lf, true);
    redrawAll();
  }

  const tab = {
    id: 'dsplab', title: 'Convolution & Spectra',
    init(panel) {
      const U = FSP.ui, el = U.el; ui = {};
      const change = () => { FSP.state.touch(); };
      bindState();
      panel.appendChild(el('p', { class: 'note', text: 'Convolution, the DFT and spectral windows: flip-and-slide a sequence, compare window shapes, and watch an off-bin sinusoid leak.' }));
      // convolution
      const fc = U.fieldset(panel, 'Convolution visualizer');
      const grid = el('div', { class: 'grid2' }); fc.appendChild(grid);
      const ctl = el('div'); grid.appendChild(ctl);
      const mk = (label, key) => { const inp = el('input', { type: 'text', 'aria-label': label, size: 18, spellcheck: 'false' }); inp.addEventListener('change', () => { if (M.parseSeq(inp.value)) { S[key] = inp.value.trim(); S.step = 0; change(); drawConv(); } else { ui.cmsg.textContent = label + ': enter 1–16 numbers separated by commas or spaces.'; ui.cmsg.hidden = false; } }); ctl.appendChild(el('div', { class: 'ctl' }, el('label', { text: label }), inp, el('span'))); return inp; };
      ui.xin = mk('x[n]', 'x'); ui.hin = mk('h[n]', 'h');
      ui.mode = U.select(ctl, 'Convolution type', [['linear', 'Linear'], ['circular', 'Circular (mod N)']], S.mode, v => { S.mode = v; S.step = 0; ui.ncRow.hidden = v !== 'circular'; change(); drawConv(); });
      ui.nc = el('input', { type: 'number', min: 1, max: 32, step: 1, 'aria-label': 'Circular length N' });
      ui.nc.addEventListener('change', () => { const v = parseInt(ui.nc.value, 10); if (v >= 1 && v <= 32) { S.nc = v; S.step = 0; change(); drawConv(); } });
      ui.ncRow = el('div', { class: 'ctl' }, el('label', { text: 'Circular length N' }), ui.nc, el('span')); ctl.appendChild(ui.ncRow);
      ui.cstep = el('input', { type: 'range', min: 0, max: 4, step: 1, value: 0, 'aria-label': 'Output index n' });
      ui.cstep.addEventListener('input', () => { stopPlay(); S.step = +ui.cstep.value; change(); drawConv(); });
      ctl.appendChild(el('div', { class: 'ctl' }, el('label', { text: 'Output index n' }), ui.cstep, el('span')));
      const br = el('div', { class: 'row' }); ctl.appendChild(br);
      ui.play = U.button(br, 'Play', () => { if (playing) { stopPlay(); } else { playing = true; lastT = 0; ui.play.textContent = 'Pause'; raf = requestAnimationFrame(tick); } });
      U.button(br, 'Step', () => { stopPlay(); const D = convData(); if (!D.err) { S.step = (Math.round(S.step) + 1) % D.N; change(); drawConv(); } });
      U.button(br, 'Reset', () => { stopPlay(); S.step = 0; change(); drawConv(); });
      ui.cmsg = el('div', { class: 'msg warn', hidden: '' }); ctl.appendChild(ui.cmsg);
      ctl.appendChild(el('div', { class: 'note', text: 'Circular with N shorter than the linear length wraps the tail onto the start (aliasing).' }));
      const plots = el('div'); grid.appendChild(plots);
      const w1 = el('div', { class: 'canvas-wrap' }); plots.appendChild(w1); ui.cv1 = U.canvas(w1, { height: 170 }); ui.cv1.onResize(drawConv);
      const w2 = el('div', { class: 'canvas-wrap' }); w2.style.marginTop = '8px'; plots.appendChild(w2); ui.cv2 = U.canvas(w2, { height: 170 }); ui.cv2.onResize(drawConv);
      ui.cread = el('div', { class: 'readout mono' }); fc.appendChild(ui.cread);
      ui.cwork = U.working(fc);
      // windows
      const fw = U.fieldset(panel, 'Window explorer');
      const g2 = el('div', { class: 'grid2' }); fw.appendChild(g2);
      const c2 = el('div'); g2.appendChild(c2);
      ui.wsel = U.select(c2, 'Window', M.WINDOWS, S.win, v => { S.win = v; change(); drawWin(); });
      ui.wnsel = U.select(c2, 'Length N', [16, 32, 64, 128, 256, 512, 1024].map(String), String(S.wn), v => { S.wn = +v; change(); drawWin(); });
      ui.spansel = U.select(c2, 'Spectrum span', [['8', '0 – 8 bins'], ['32', '0 – 32 bins'], ['full', 'Full (0 – N/2)']], S.span, v => { S.span = v; change(); drawWin(); });
      const hud = el('div', { class: 'hud' }); c2.appendChild(hud);
      const hv = (lab) => { const o = el('span', { class: 'mono', text: '—' }); hud.appendChild(el('div', null, el('span', { text: lab }), o)); return o; };
      ui.wsl = hv('Peak sidelobe'); ui.wenbw = hv('ENBW'); ui.wcg = hv('Coherent gain'); ui.wnull = hv('First null');
      const p2 = el('div'); g2.appendChild(p2);
      const a1 = el('div', { class: 'canvas-wrap' }); p2.appendChild(a1); ui.wv1 = U.canvas(a1, { height: 140 }); ui.wv1.onResize(drawWin);
      const a2 = el('div', { class: 'canvas-wrap' }); a2.style.marginTop = '8px'; p2.appendChild(a2); ui.wv2 = U.canvas(a2, { height: 220 }); ui.wv2.onResize(drawWin);
      ui.wwork = U.working(fw);
      // leakage
      const fl = U.fieldset(panel, 'Off-bin sinusoid leakage');
      const g3 = el('div', { class: 'grid2' }); fl.appendChild(g3);
      const c3 = el('div'); g3.appendChild(c3);
      ui.lnsel = U.select(c3, 'Record length N', ['64', '128', '256', '512'], String(S.ln), v => { S.ln = +v; S.lf = Math.min(S.lf, S.ln / 2); ui.lfs.set(S.lf, true); change(); drawLeak(); });
      ui.lfs = U.slider(c3, { label: 'Frequency (bins)', min: 1, max: 100, step: 0.05, value: S.lf, digits: 5, onInput: v => { S.lf = Math.min(v, S.ln / 2); change(); drawLeak(); } });
      ui.lwsel = U.select(c3, 'Window', M.WINDOWS, S.lw, v => { S.lw = v; change(); drawLeak(); });
      ui.lspansel = U.select(c3, 'View ± bins', ['4', '16', '64'], S.lspan, v => { S.lspan = v; change(); drawLeak(); });
      const bb = el('div', { class: 'row' }); c3.appendChild(bb);
      U.button(bb, 'On-bin', () => { S.lf = Math.round(S.lf); ui.lfs.set(S.lf, true); change(); drawLeak(); });
      U.button(bb, 'Half-bin', () => { S.lf = Math.floor(S.lf) + 0.5; ui.lfs.set(S.lf, true); change(); drawLeak(); });
      const hud2 = el('div', { class: 'hud' }); c3.appendChild(hud2);
      const hv2 = lab => { const o = el('span', { class: 'mono', text: '—' }); hud2.appendChild(el('div', null, el('span', { text: lab }), o)); return o; };
      ui.lbin = hv2('Alignment'); ui.lra = hv2('Far leakage, rect'); ui.lrb = hv2('Far leakage, window'); ui.lpk = hv2('Peak rect / window');
      const a3 = el('div', { class: 'canvas-wrap' }); g3.appendChild(a3); ui.lv = U.canvas(a3, { height: 260 }); ui.lv.onResize(drawLeak);
      ui.lwork = U.working(fl);
      syncUI();
    },
    activate() { if (ui) redrawAll(); },
    deactivate() { stopPlay(); },
  };
  FSP.registerTab(tab);
})();
