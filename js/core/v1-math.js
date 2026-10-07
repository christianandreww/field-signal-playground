/* ===================== CORE (pure math, no DOM) ===================== */
'use strict';
class Complex {
  constructor(re, im) { this.re = re; this.im = im === undefined ? 0 : im; }
  add(b) { return new Complex(this.re + b.re, this.im + b.im); }
  sub(b) { return new Complex(this.re - b.re, this.im - b.im); }
  mul(b) { return new Complex(this.re * b.re - this.im * b.im, this.re * b.im + this.im * b.re); }
  div(b) { const d = b.re * b.re + b.im * b.im; return new Complex((this.re * b.re + this.im * b.im) / d, (this.im * b.re - this.re * b.im) / d); }
  scale(s) { return new Complex(this.re * s, this.im * s); }
  neg() { return new Complex(-this.re, -this.im); }
  conj() { return new Complex(this.re, -this.im); }
  abs() { return Math.hypot(this.re, this.im); }
  abs2() { return this.re * this.re + this.im * this.im; }
  arg() { return Math.atan2(this.im, this.re); }
  sqrt() { // principal branch, Re >= 0, numerically stable
    const r = this.abs();
    if (r === 0) return new Complex(0, 0);
    if (this.re >= 0) { const t = Math.sqrt(0.5 * (r + this.re)); return new Complex(t, this.im / (2 * t)); }
    const t = Math.sqrt(0.5 * (r - this.re));
    const s = this.im < 0 ? -t : t;
    return new Complex(this.im / (2 * s), s);
  }
  exp() { const e = Math.exp(this.re); return new Complex(e * Math.cos(this.im), e * Math.sin(this.im)); }
  isFinite() { return Number.isFinite(this.re) && Number.isFinite(this.im); }
}
const C = (re, im) => new Complex(re, im || 0);
const reMulConj = (a, b) => a.re * b.re + a.im * b.im;          // Re(a b*)
const instVal = (a, phi) => a.re * Math.cos(phi) - a.im * Math.sin(phi); // Re(a e^{jφ})

/* ---------- EM: media and interface ---------- */
const PHYS = { EPS0: 8.8541878128e-12, MU0: 1.25663706212e-6, C0: 299792458 };
PHYS.ETA0 = Math.sqrt(PHYS.MU0 / PHYS.EPS0);

function makeMedium(er, mr, sigma, freq) {
  const omega = 2 * Math.PI * freq;
  const ecRel = C(er, -sigma / (omega * PHYS.EPS0));       // ε_c / ε0 = εr − jσ/(ωε0)
  const nc = C(mr * ecRel.re, mr * ecRel.im).sqrt();        // √(μr ε_c/ε0), Im ≤ 0
  const eta = C(mr, 0).div(ecRel).sqrt().scale(PHYS.ETA0);  // √(μ/ε_c), Re ≥ 0
  const Y = C(1, 0).div(eta);
  return { er, mr, sigma, freq, omega, ecRel, nc, eta, Y, nReal: Math.sqrt(Math.max(0, er * mr)), lossy: sigma > 0 };
}

function solveInterface(thetaDeg, pol, m1, m2) {
  const th = thetaDeg * Math.PI / 180;
  const si = Math.sin(th), ci = Math.cos(th);
  const k0 = 2 * Math.PI;                                    // normalized λ0 = 1
  const k1 = m1.nc.scale(k0), k2 = m2.nc.scale(k0);
  const kx = k1.scale(si), kz1 = k1.scale(ci);
  let kz2 = k2.mul(k2).sub(kx.mul(kx)).sqrt();
  if (kz2.im > 0 || (kz2.im === 0 && kz2.re < 0)) kz2 = kz2.neg(); // decaying / forward branch
  const cosI = C(ci, 0), sinI = C(si, 0);
  const cosT = kz2.div(k2), sinT = kx.div(k2);
  const eta1 = m1.eta, eta2 = m2.eta;
  let r, t, den;
  if (pol === 'TE') {
    den = eta2.mul(cosI).add(eta1.mul(cosT));
    r = eta2.mul(cosI).sub(eta1.mul(cosT)).div(den);
    t = eta2.mul(cosI).scale(2).div(den);
  } else {
    den = eta2.mul(cosT).add(eta1.mul(cosI));
    r = eta2.mul(cosT).sub(eta1.mul(cosI)).div(den);
    t = eta2.mul(cosI).scale(2).div(den);
  }
  if (!(den.abs() > 1e-300) || !r.isFinite() || !t.isFinite()) { r = C(NaN, NaN); t = C(NaN, NaN); }
  const evanescent = !m1.lossy && !m2.lossy && Math.abs(kz2.re) <= 1e-9 * (kz2.abs() + 1e-300) && kz2.im < 0;
  return { thetaDeg, pol, m1, m2, k1, k2, kx, kz1, kz2, cosI, sinI, cosT, sinT, r, t, evanescent };
}

// exp(-j (kx x + sign kz z)) with growth clamped (never NaN/Inf)
function waveFactor(kx, kz, x, z, sign) {
  const pr = kx.re * x + sign * kz.re * z, pi = kx.im * x + sign * kz.im * z;
  const mag = Math.exp(Math.min(40, pi));
  return new Complex(mag * Math.cos(-pr), mag * Math.sin(-pr));
}

// Phasor components at (x, z): TE → {A:Ey, B:Hx, D:Hz}; TM → {A:Hy, B:Ex, D:Ez}
function fieldPhasors(sol, x, z, part) {
  const { pol, kx, kz1, kz2, r, t, m1, m2, cosI, sinI, cosT, sinT } = sol;
  let A = C(0), B = C(0), D = C(0);
  const inc = part === 'inc' || (part === 'total' && z < 0);
  const ref = part === 'ref' || (part === 'total' && z < 0);
  const tr = part === 'tr' || (part === 'total' && z >= 0);
  if (pol === 'TE') {
    const cY1 = cosI.mul(m1.Y), sY1 = sinI.mul(m1.Y), cY2 = cosT.mul(m2.Y), sY2 = sinT.mul(m2.Y);
    if (inc) { const e = waveFactor(kx, kz1, x, z, +1); A = A.add(e); B = B.sub(cY1.mul(e)); D = D.add(sY1.mul(e)); }
    if (ref) { const e = waveFactor(kx, kz1, x, z, -1).mul(r); A = A.add(e); B = B.add(cY1.mul(e)); D = D.add(sY1.mul(e)); }
    if (tr) { const e = waveFactor(kx, kz2, x, z, +1).mul(t); A = A.add(e); B = B.sub(cY2.mul(e)); D = D.add(sY2.mul(e)); }
  } else {
    const H0 = m1.Y, rH = r.neg(), tH = C(1).add(rH);
    const e1c = m1.eta.mul(cosI), e1s = m1.eta.mul(sinI), e2c = m2.eta.mul(cosT), e2s = m2.eta.mul(sinT);
    if (inc) { const h = waveFactor(kx, kz1, x, z, +1).mul(H0); A = A.add(h); B = B.add(e1c.mul(h)); D = D.sub(e1s.mul(h)); }
    if (ref) { const h = waveFactor(kx, kz1, x, z, -1).mul(H0).mul(rH); A = A.add(h); B = B.sub(e1c.mul(h)); D = D.sub(e1s.mul(h)); }
    if (tr) { const h = waveFactor(kx, kz2, x, z, +1).mul(H0).mul(tH); A = A.add(h); B = B.add(e2c.mul(h)); D = D.sub(e2s.mul(h)); }
  }
  return { A, B, D };
}
function poyntingAvg(pol, F) { // ½ Re(E × H*)
  if (pol === 'TE') return { x: 0.5 * reMulConj(F.A, F.D), z: -0.5 * reMulConj(F.A, F.B) };
  return { x: -0.5 * reMulConj(F.D, F.A), z: 0.5 * reMulConj(F.B, F.A) };
}
function poyntingInst(pol, F, phi) { // E(t) × H(t) from real instantaneous fields
  const a = instVal(F.A, phi), b = instVal(F.B, phi), d = instVal(F.D, phi);
  if (pol === 'TE') return { x: a * d, z: -a * b };
  return { x: -d * a, z: b * a };
}

function criticalAngle(m1, m2) {
  if (m1.lossy || m2.lossy || !(m1.nReal > m2.nReal)) return null;
  return Math.asin(m2.nReal / m1.nReal) * 180 / Math.PI;
}
function brewsterAnalytic(pol, m1, m2) {
  if (m1.lossy || m2.lossy) return null;
  const e1 = m1.er, u1 = m1.mr, e2 = m2.er, u2 = m2.mr;
  let num, den;
  if (pol === 'TM') { num = u1 / e1 - u2 / e2; den = u1 / e1 - u1 * e1 / (e2 * e2); }
  else { num = u2 / e2 - u1 / e1; den = u2 / e2 - u1 * u1 / (u2 * e2); }
  if (!(Math.abs(den) > 1e-15)) return null;
  const s2 = num / den;
  if (!(s2 >= 0 && s2 < 1)) return null;
  const th = Math.asin(Math.sqrt(s2)) * 180 / Math.PI;
  const sol = solveInterface(th, pol, m1, m2);            // verify numerically
  if (!sol.r.isFinite() || sol.r.abs() > 1e-6) return null;
  return th;
}
function pseudoBrewster(m1, m2) {
  const f = th => { const s = solveInterface(th, 'TM', m1, m2); return s.r.isFinite() ? s.r.abs() : Infinity; };
  let best = 0, bestV = Infinity;
  for (let th = 0; th <= 89; th += 0.5) { const v = f(th); if (v < bestV) { bestV = v; best = th; } }
  let a = Math.max(0, best - 0.5), b = Math.min(89, best + 0.5);
  const g = (Math.sqrt(5) - 1) / 2;
  let x1 = b - g * (b - a), x2 = a + g * (b - a), f1 = f(x1), f2 = f(x2);
  for (let i = 0; i < 80; i++) {
    if (f1 < f2) { b = x2; x2 = x1; f2 = f1; x1 = b - g * (b - a); f1 = f(x1); }
    else { a = x1; x1 = x2; f1 = f2; x2 = a + g * (b - a); f2 = f(x2); }
  }
  const th = 0.5 * (a + b);
  return { theta: th, rmin: f(th) };
}

function emAnalyze(thetaDeg, pol, m1, m2) {
  const sol = solveInterface(thetaDeg, pol, m1, m2);
  const R = sol.r.isFinite() ? sol.r.abs2() : NaN;
  const inc = fieldPhasors(sol, 0, 0, 'inc'), tr = fieldPhasors(sol, 0, 0, 'tr');
  const SzInc = poyntingAvg(pol, inc).z, SzTr = poyntingAvg(pol, tr).z;
  let T = (Number.isFinite(SzInc) && Math.abs(SzInc) > 1e-300 && Number.isFinite(SzTr)) ? SzTr / SzInc : NaN;
  if (Number.isFinite(T) && Math.abs(T) < 1e-15) T = 0;
  const sum = R + T;
  const lossless = !m1.lossy && !m2.lossy;
  const conserv = m1.lossy ? 'na' : (Number.isFinite(sum) && Math.abs(sum - 1) <= 1e-9 ? 'pass' : 'fail');
  const phaseDeg = sol.r.isFinite() ? argDeg(sol.r) : NaN;
  const delta = sol.kz2.im < -1e-12 ? 1 / (-sol.kz2.im) : Infinity;    // in λ0 units
  const lambda0 = PHYS.C0 / m1.freq;
  const thetaT = sol.evanescent ? null : Math.atan2(sol.kx.re, sol.kz2.re) * 180 / Math.PI;
  let brewster;
  if (lossless) brewster = { kind: 'exact', tm: brewsterAnalytic('TM', m1, m2), te: brewsterAnalytic('TE', m1, m2) };
  else { const p = pseudoBrewster(m1, m2); brewster = { kind: 'pseudo', tm: p.theta, rmin: p.rmin, te: null, label: 'pseudo-Brewster (min |r_TM|)' }; }
  return { sol, r: sol.r, t: sol.t, R, T, sum, conserv, phaseDeg, delta, deltaPhys: delta * lambda0, lambda0, thetaT,
           evanescent: sol.evanescent, thetaC: criticalAngle(m1, m2), brewster, lossless };
}
function argDeg(c) { let d = Math.atan2(c.im, c.re) * 180 / Math.PI; if (d <= -180 + 1e-7) d = 180; return d; }

// Field cache for the heat map: iw×ih internal pixels, `span` wavelengths across, aspect = h/w
function buildFieldCache(sol, iw, ih, span, aspect) {
  const n = iw * ih;
  const Are = new Float32Array(n), Aim = new Float32Array(n), Bre = new Float32Array(n), Bim = new Float32Array(n), Dre = new Float32Array(n), Dim = new Float32Array(n);
  const { kx, kz1, kz2, pol, r, t, m1, m2, cosI, sinI, cosT, sinT } = sol;
  const exr = new Float64Array(iw), exi = new Float64Array(iw);
  for (let i = 0; i < iw; i++) {
    const x = ((i + 0.5) / iw - 0.5) * span;
    const pr = kx.re * x, mg = Math.exp(Math.min(40, kx.im * x));
    exr[i] = mg * Math.cos(-pr); exi[i] = mg * Math.sin(-pr);
  }
  let MI, MR, MT, KB, KD, KB2, KD2;
  if (pol === 'TE') {
    MI = C(1); MR = r; MT = t;
    KB = cosI.mul(m1.Y).neg(); KD = sinI.mul(m1.Y); KB2 = cosT.mul(m2.Y).neg(); KD2 = sinT.mul(m2.Y);
  } else {
    const H0 = m1.Y, rH = r.neg(), tH = C(1).add(rH);
    MI = H0; MR = H0.mul(rH); MT = H0.mul(tH);
    KB = m1.eta.mul(cosI); KD = m1.eta.mul(sinI).neg(); KB2 = m2.eta.mul(cosT); KD2 = m2.eta.mul(sinT).neg();
  }
  let vmaxA = 0, vmaxE = 0;
  for (let j = 0; j < ih; j++) {
    const z = ((j + 0.5) / ih - 0.5) * span * aspect;
    const row = j * iw;
    if (z < 0) {
      const pr = kz1.re * z, pi = kz1.im * z;
      const mI = Math.exp(Math.min(40, pi)), mR = Math.exp(Math.min(40, -pi));
      let fIr = mI * Math.cos(-pr), fIi = mI * Math.sin(-pr), fRr = mR * Math.cos(pr), fRi = mR * Math.sin(pr);
      let a = fIr * MI.re - fIi * MI.im, b = fIr * MI.im + fIi * MI.re; fIr = a; fIi = b;
      a = fRr * MR.re - fRi * MR.im; b = fRr * MR.im + fRi * MR.re; fRr = a; fRi = b;
      for (let i = 0; i < iw; i++) {
        const er = exr[i], ei = exi[i];
        const ur = er * fIr - ei * fIi, ui = er * fIi + ei * fIr;
        const vr = er * fRr - ei * fRi, vi = er * fRi + ei * fRr;
        const sr = ur + vr, si = ui + vi, dr = ur - vr, di = ui - vi;
        const k = row + i;
        Are[k] = sr; Aim[k] = si;
        Bre[k] = KB.re * dr - KB.im * di; Bim[k] = KB.re * di + KB.im * dr;
        Dre[k] = KD.re * sr - KD.im * si; Dim[k] = KD.re * si + KD.im * sr;
      }
    } else {
      const pr = kz2.re * z, mT = Math.exp(Math.min(40, kz2.im * z));
      let fTr = mT * Math.cos(-pr), fTi = mT * Math.sin(-pr);
      const a = fTr * MT.re - fTi * MT.im, b = fTr * MT.im + fTi * MT.re; fTr = a; fTi = b;
      for (let i = 0; i < iw; i++) {
        const er = exr[i], ei = exi[i];
        const hr = er * fTr - ei * fTi, hi = er * fTi + ei * fTr;
        const k = row + i;
        Are[k] = hr; Aim[k] = hi;
        Bre[k] = KB2.re * hr - KB2.im * hi; Bim[k] = KB2.re * hi + KB2.im * hr;
        Dre[k] = KD2.re * hr - KD2.im * hi; Dim[k] = KD2.re * hi + KD2.im * hr;
      }
    }
  }
  for (let k = 0; k < n; k++) {
    if (!Number.isFinite(Are[k]) || !Number.isFinite(Aim[k])) { Are[k] = 0; Aim[k] = 0; }
    if (!Number.isFinite(Bre[k]) || !Number.isFinite(Bim[k])) { Bre[k] = 0; Bim[k] = 0; }
    if (!Number.isFinite(Dre[k]) || !Number.isFinite(Dim[k])) { Dre[k] = 0; Dim[k] = 0; }
    const ma = Math.hypot(Are[k], Aim[k]); if (ma > vmaxA) vmaxA = ma;
    const me = Math.sqrt(Bre[k] * Bre[k] + Bim[k] * Bim[k] + Dre[k] * Dre[k] + Dim[k] * Dim[k]); if (me > vmaxE) vmaxE = me;
  }
  return { iw, ih, span, aspect, pol, Are, Aim, Bre, Bim, Dre, Dim, vmaxA: vmaxA || 1, vmaxE: vmaxE || 1 };
}
/* ---------- DSP: biquad core shared by the UI and the AudioWorklet ---------- */
const BIQUAD_CORE_SRC = `
function df2tStep(xn, c, st) {
  // Direct Form II Transposed; c = [b0, b1, b2, a1, a2], st = [s1, s2]
  const yn = c[0] * xn + st[0];
  st[0] = c[1] * xn - c[3] * yn + st[1];
  st[1] = c[2] * xn - c[4] * yn;
  return yn;
}
function df2tProcess(x, y, n, c, st) {
  for (let i = 0; i < n; i++) y[i] = df2tStep(x[i], c, st);
}`;
const DF2T = new Function(BIQUAD_CORE_SRC + '\nreturn { df2tStep, df2tProcess };')();
const WORKLET_SRC = BIQUAD_CORE_SRC + `
class BiquadDF2TProcessor extends AudioWorkletProcessor {
  constructor() {
    super();
    this.cur = new Float64Array([1, 0, 0, 0, 0]);
    this.target = new Float64Array([1, 0, 0, 0, 0]);
    this.delta = new Float64Array(5);
    this.rampLeft = 0;
    this.states = [];
    this.port.onmessage = (e) => {
      const d = e.data;
      if (d && d.type === 'coefs') {
        const n = Math.max(1, Math.round(sampleRate * 0.005)); // 5 ms linear ramp
        for (let i = 0; i < 5; i++) { this.target[i] = d.c[i]; this.delta[i] = (d.c[i] - this.cur[i]) / n; }
        this.rampLeft = n;
      } else if (d && d.type === 'reset') {
        for (const s of this.states) { s[0] = 0; s[1] = 0; }
      }
    };
  }
  process(inputs, outputs) {
    const inp = inputs[0], out = outputs[0];
    const nOut = out.length, nIn = inp ? inp.length : 0;
    const n = nOut ? out[0].length : 128;
    for (let ch = 0; ch < nOut; ch++) if (!this.states[ch]) this.states[ch] = new Float64Array(2);
    const cur = this.cur;
    for (let i = 0; i < n; i++) {
      if (this.rampLeft > 0) {
        for (let k = 0; k < 5; k++) cur[k] += this.delta[k];
        if (--this.rampLeft === 0) for (let k = 0; k < 5; k++) cur[k] = this.target[k];
      }
      for (let ch = 0; ch < nOut; ch++) {
        const xn = nIn ? inp[Math.min(ch, nIn - 1)][i] : 0;
        out[ch][i] = df2tStep(xn, cur, this.states[ch]);
      }
    }
    for (let ch = 0; ch < nOut; ch++) { const s = this.states[ch]; if (Math.abs(s[0]) < 1e-20) s[0] = 0; if (Math.abs(s[1]) < 1e-20) s[1] = 0; }
    return true;
  }
}
registerProcessor('biquad-df2t', BiquadDF2TProcessor);`;

const FILTER_NAMES = { lp: 'Lowpass', hp: 'Highpass', bp: 'Bandpass (0 dB peak)', notch: 'Notch', peak: 'Peaking EQ' };
function rbjCoefs(type, fc, fs, Q, gainDb) {
  const w0 = 2 * Math.PI * fc / fs, cw = Math.cos(w0), sw = Math.sin(w0);
  const alpha = sw / (2 * Q), A = Math.pow(10, gainDb / 40);
  let b0, b1, b2, a0, a1, a2;
  switch (type) {
    case 'lp': b0 = (1 - cw) / 2; b1 = 1 - cw; b2 = (1 - cw) / 2; a0 = 1 + alpha; a1 = -2 * cw; a2 = 1 - alpha; break;
    case 'hp': b0 = (1 + cw) / 2; b1 = -(1 + cw); b2 = (1 + cw) / 2; a0 = 1 + alpha; a1 = -2 * cw; a2 = 1 - alpha; break;
    case 'bp': b0 = alpha; b1 = 0; b2 = -alpha; a0 = 1 + alpha; a1 = -2 * cw; a2 = 1 - alpha; break;
    case 'notch': b0 = 1; b1 = -2 * cw; b2 = 1; a0 = 1 + alpha; a1 = -2 * cw; a2 = 1 - alpha; break;
    case 'peak': b0 = 1 + alpha * A; b1 = -2 * cw; b2 = 1 - alpha * A; a0 = 1 + alpha / A; a1 = -2 * cw; a2 = 1 - alpha / A; break;
    default: throw new Error('unknown filter type ' + type);
  }
  return { b0: b0 / a0, b1: b1 / a0, b2: b2 / a0, a0: 1, a1: a1 / a0, a2: a2 / a0, a0raw: a0 };
}
function coefArray(c) { return [c.b0, c.b1, c.b2, c.a1, c.a2]; }
function freqResponseAt(c, w) { // H(e^{jw}) as Complex
  const c1 = Math.cos(w), s1 = Math.sin(w), c2 = Math.cos(2 * w), s2 = Math.sin(2 * w);
  const num = C(c.b0 + c.b1 * c1 + c.b2 * c2, -(c.b1 * s1 + c.b2 * s2));
  const den = C(1 + c.a1 * c1 + c.a2 * c2, -(c.a1 * s1 + c.a2 * s2));
  return num.div(den);
}
function magDb(h) { return 20 * Math.log10(Math.max(h.abs(), 1e-6)); } // floor −120 dB
function freqResponse(c, fs, nPts, fmin) {
  nPts = nPts || 1024; fmin = fmin || 20;
  const fmax = fs / 2, freqs = new Float64Array(nPts), mag = new Float64Array(nPts), phase = new Float64Array(nPts);
  let prev = 0, acc = 0;
  for (let i = 0; i < nPts; i++) {
    const f = fmin * Math.pow(fmax / fmin, i / (nPts - 1));
    const h = freqResponseAt(c, 2 * Math.PI * f / fs);
    freqs[i] = f; mag[i] = magDb(h);
    const p = h.abs() < 1e-12 ? prev : h.arg(); // exact zeros (e.g. LP at Nyquist) carry the previous phase
    if (i > 0) { let d = p - prev; while (d > Math.PI) d -= 2 * Math.PI; while (d < -Math.PI) d += 2 * Math.PI; acc += d; } else acc = p;
    prev = p; phase[i] = acc * 180 / Math.PI;
  }
  return { freqs, mag, phase, fmin, fmax };
}
function quadRoots(a, b, c) { // a z² + b z + c = 0
  const tiny = 1e-300;
  if (Math.abs(a) < tiny) {
    if (Math.abs(b) < tiny) return { roots: [], note: 'constant numerator' };
    return { roots: [C(-c / b, 0)], note: 'one finite root (leading coefficient is 0)' };
  }
  const disc = b * b - 4 * a * c, scale = b * b + Math.abs(4 * a * c);
  if (Math.abs(disc) <= 1e-12 * scale) { const x = -b / (2 * a); return { roots: [C(x, 0), C(x, 0)], double: true }; }
  if (disc > 0) {
    const q = -0.5 * (b + Math.sign(b || 1) * Math.sqrt(disc));
    const r1 = q / a, r2 = Math.abs(q) > tiny ? c / q : 0;
    return { roots: [C(r1, 0), C(r2, 0)] };
  }
  const re = -b / (2 * a), im = Math.sqrt(-disc) / (2 * Math.abs(a));
  return { roots: [C(re, im), C(re, -im)], pair: true };
}
function padRoots(q) { // ensure two roots for the editor (a missing zero becomes z = 0: pure delay, |H| unchanged)
  const r = q.roots.slice(); while (r.length < 2) r.push(C(0, 0)); return r;
}
function rootsFromCoefs(c) {
  const qp = quadRoots(1, c.a1, c.a2), qz = quadRoots(c.b0, c.b1, c.b2);
  const poles = padRoots(qp), zeros = padRoots(qz);
  const isReal = arr => arr.every(z => Math.abs(z.im) < 1e-12);
  let K = c.b0;
  if (Math.abs(c.b0) < 1e-300) K = Math.abs(c.b1) > 1e-300 ? c.b1 : c.b2;
  return { poles, zeros, polesReal: isReal(poles), zerosReal: isReal(zeros), K, zeroNote: qz.note || '', poleNote: qp.note || '' };
}
function coefsFromRoots(roots) {
  const [p1, p2] = roots.poles, [z1, z2] = roots.zeros, K = roots.K;
  return { b0: K, b1: -K * (z1.re + z2.re), b2: K * z1.mul(z2).re, a0: 1, a1: -(p1.re + p2.re), a2: p1.mul(p2).re, a0raw: 1 };
}
function poleRadii(c) { return quadRoots(1, c.a1, c.a2).roots.map(p => p.abs()); }
function isStable(c) { const r = poleRadii(c); return [c.b0, c.b1, c.b2, c.a1, c.a2].every(Number.isFinite) && r.length > 0 && r.every(x => x < 1); }
const POLE_RMAX = 0.9999, ZERO_RMAX = 3;
function clampRoot(kind, x, y) {
  const rmax = kind === 'pole' ? POLE_RMAX : ZERO_RMAX, r = Math.hypot(x, y);
  if (!(r <= rmax)) { if (!Number.isFinite(r) || r === 0) return [0, 0]; return [x * rmax / r, y * rmax / r]; }
  return [x, y];
}
// Move root `idx` of kind 'pole'|'zero' to (x, y); enforce conjugate symmetry and real-axis rules. Mutates and returns roots.
function applyRootMove(roots, kind, idx, x, y, snapTol) {
  const arr = kind === 'pole' ? roots.poles : roots.zeros, flag = kind === 'pole' ? 'polesReal' : 'zerosReal';
  let [X, Y] = clampRoot(kind, x, y);
  if (Math.abs(Y) < snapTol) Y = 0;
  const other = 1 - idx;
  if (Y === 0) {
    if (!roots[flag]) { arr[idx] = C(X, 0); arr[other] = C(X, 0); roots[flag] = true; } // pair → double real root
    else arr[idx] = C(X, 0);                                                            // real root slides alone
  } else { arr[idx] = C(X, Y); arr[other] = C(X, -Y); roots[flag] = false; }              // real/pair → conjugate pair
  return roots;
}
function refFreqFor(type, fc, fs) { return type === 'lp' ? 0 : type === 'hp' ? fs / 2 : fc; }
function impulseResponse(c, N) {
  const x = new Float64Array(N), y = new Float64Array(N); x[0] = 1;
  DF2T.df2tProcess(x, y, N, coefArray(c), new Float64Array(2));
  return y;
}
function trimLength(h, eps) { let last = 0; for (let i = 0; i < h.length; i++) if (Math.abs(h[i]) >= eps) last = i; return Math.max(2, last + 1); }

/* ---------- formatting ---------- */
function fmtNum(x, d) { return Number.isFinite(x) ? x.toFixed(d === undefined ? 4 : d) : '—'; }
function fmtDeg(x, d) { return Number.isFinite(x) ? x.toFixed(d === undefined ? 2 : d) + '°' : '—'; }
function fmtSI(v, unit, digits) {
  digits = digits || 3;
  if (!Number.isFinite(v)) return '—';
  if (v === 0) return '0 ' + unit;
  const a = Math.abs(v), P = [[1e12, 'T'], [1e9, 'G'], [1e6, 'M'], [1e3, 'k'], [1, ''], [1e-3, 'm'], [1e-6, 'µ'], [1e-9, 'n'], [1e-12, 'p'], [1e-15, 'f']];
  for (const [s, p] of P) if (a >= s) return (v / s).toPrecision(digits) + ' ' + p + unit;
  return v.toExponential(2) + ' ' + unit;
}
function fmtPolar(c) { return c && c.isFinite() ? c.abs().toFixed(4) + ' ∠ ' + argDeg(c).toFixed(2) + '°' : '—'; }
function fmtComplexN(c) { return c && c.isFinite() ? c.re.toFixed(3) + (c.im < 0 ? ' − j' : ' + j') + Math.abs(c.im).toFixed(3) : '—'; }
/* ---------- acceptance self-tests (console: runSelfTests()) ---------- */
function runV1Tests() {
  const lines = [], results = [];
  const check = (name, ok, detail) => { ok = !!ok; results.push({ name, ok }); const l = (ok ? 'PASS' : 'FAIL') + ' — ' + name + (detail ? '  [' + detail + ']' : ''); lines.push(l); console.log(l); };
  const near = (a, b, tol) => Number.isFinite(a) && Number.isFinite(b) && Math.abs(a - b) <= tol;
  const rnd = (seed => () => { seed |= 0; seed = seed + 0x6D2B79F5 | 0; let t = Math.imul(seed ^ seed >>> 15, 1 | seed); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; })(20240607);
  const F = 1e9;
  // 1
  { const m1 = makeMedium(1, 1, 0, F), m2 = makeMedium(2.25, 1, 0, F);
    const an = emAnalyze(0, 'TE', m1, m2);
    check('EM1 Brewster n1=1,n2=1.5 ≈ 56.3099°', near(an.brewster.tm, 56.3099, 1e-3), fmtDeg(an.brewster.tm, 4));
    check('EM1 θi=0 TE: r=−0.2, t=0.8, R=0.04, T=0.96', near(an.r.re, -0.2, 1e-9) && near(an.r.im, 0, 1e-9) && near(an.t.re, 0.8, 1e-9) && near(an.R, 0.04, 1e-9) && near(an.T, 0.96, 1e-9),
      'r=' + fmtNum(an.r.re, 6) + ' t=' + fmtNum(an.t.re, 6) + ' R=' + fmtNum(an.R, 6) + ' T=' + fmtNum(an.T, 6)); }
  // 2
  { const m1 = makeMedium(2.25, 1, 0, F), m2 = makeMedium(1, 1, 0, F);
    const an = emAnalyze(60, 'TE', m1, m2), an2 = emAnalyze(60, 'TM', m1, m2);
    check('EM2 critical angle n1=1.5→1 ≈ 41.8103°', near(an.thetaC, 41.8103, 1e-3), fmtDeg(an.thetaC, 4));
    check('EM2 TIR at 60°: |r|=1, T≤1e-9, δ finite (both pols)', near(an.r.abs(), 1, 1e-9) && an.T <= 1e-9 && Number.isFinite(an.delta) && an.delta > 0 && near(an2.r.abs(), 1, 1e-9) && an2.T <= 1e-9 && Number.isFinite(an2.delta),
      '|r|=' + fmtNum(an.r.abs(), 10) + ' T=' + an.T.toExponential(2) + ' δ=' + fmtNum(an.delta, 4) + ' λ0'); }
  // 3
  { let worst = 0, nTir = 0, nProp = 0, ok = true;
    for (let i = 0; i < 500; i++) {
      const m1 = makeMedium(1 + 11 * rnd(), 1 + 4 * rnd(), 0, F), m2 = makeMedium(1 + 11 * rnd(), 1 + 4 * rnd(), 0, F);
      const th = 89 * rnd(), pol = rnd() < 0.5 ? 'TE' : 'TM';
      const an = emAnalyze(th, pol, m1, m2);
      if (an.evanescent) nTir++; else nProp++;
      const e = Math.abs(an.sum - 1); if (!(e <= 1e-9)) ok = false; if (e > worst || !Number.isFinite(e)) worst = e;
    }
    check('EM3 R+T=1 within 1e-9 for 500 random lossless cases', ok && nTir > 0 && nProp > 0, 'worst |R+T−1|=' + worst.toExponential(2) + ', TIR cases=' + nTir); }
  // 4
  { let ok = true, bad = '';
    const media = [[1, 1, 0, 2.25, 1, 0], [2.25, 1, 0, 1, 1, 0], [1, 1, 0, 4, 1, 0.05], [1, 1, 0, 1, 1, 5.8e7], [1, 1, 0, 2, 3, 0], [4, 1, 0.1, 1, 1, 0]];
    for (const m of media) for (const th of [0, 89]) for (const pol of ['TE', 'TM']) {
      const an = emAnalyze(th, pol, makeMedium(m[0], m[1], m[2], F), makeMedium(m[3], m[4], m[5], F));
      const vals = [an.r.re, an.r.im, an.t.re, an.t.im, an.R, an.T, an.phaseDeg];
      const txt = [fmtPolar(an.r), fmtPolar(an.t), fmtNum(an.R), fmtNum(an.T), fmtDeg(an.phaseDeg), fmtDeg(an.thetaT), fmtSI(an.deltaPhys, 'm'), fmtDeg(an.brewster.tm), fmtDeg(an.thetaC)].join(' ');
      if (!vals.every(Number.isFinite) || /NaN|Infinity/.test(txt)) { ok = false; bad = JSON.stringify(m) + ' θ=' + th + ' ' + pol; }
    }
    check('EM4 θi=0 and 89° give no NaN/Infinity for any polarization', ok, bad || '6 media × 2 angles × 2 pols'); }
  // 5
  { const m1 = makeMedium(1, 1, 0, F), m2 = makeMedium(4, 1, 0.05, F);
    const a = emAnalyze(30, 'TE', m1, m2), b = emAnalyze(70, 'TM', m1, m2);
    check('EM5 σ2>0: R+T=1 within 1e-9 (TE 30°, TM 70°)', near(a.sum, 1, 1e-9) && near(b.sum, 1, 1e-9), 'TE: ' + a.sum.toFixed(12) + ' TM: ' + b.sum.toFixed(12));
    check('EM5 pseudo-Brewster reported and labelled', b.brewster.kind === 'pseudo' && Number.isFinite(b.brewster.tm) && /pseudo/.test(b.brewster.label), fmtDeg(b.brewster.tm, 3) + ' |r|min=' + b.brewster.rmin.toExponential(2)); }
  // 6
  { const fs = 48000, fc = fs / 10, c = rbjCoefs('lp', fc, fs, 0.7071, 0);
    const atFc = magDb(freqResponseAt(c, 2 * Math.PI * fc / fs)), dc = magDb(freqResponseAt(c, 0)), ny = magDb(freqResponseAt(c, Math.PI));
    check('DSP6 LP Q=0.7071: |H(fc)|=−3.01±0.02 dB, DC=0 dB, Nyquist<−30 dB', near(atFc, -3.01, 0.02) && near(dc, 0, 1e-9) && ny < -30, 'fc:' + atFc.toFixed(3) + ' dc:' + dc.toExponential(1) + ' nyq:' + ny.toFixed(1)); }
  // 7
  { const fs = 48000, fc = 1000;
    const pk = magDb(freqResponseAt(rbjCoefs('peak', fc, fs, 1, 12), 2 * Math.PI * fc / fs));
    const nt = magDb(freqResponseAt(rbjCoefs('notch', fc, fs, 1, 0), 2 * Math.PI * fc / fs));
    check('DSP7 Peaking +12 dB → |H(fc)|=12.00±0.01; Notch |H(fc)|<−80 dB', near(pk, 12, 0.01) && nt < -80, 'peak:' + pk.toFixed(4) + ' notch:' + nt.toFixed(1)); }
  // 8
  { const N = 256, c = rbjCoefs('lp', 4800, 48000, 0.7071, 0), h = impulseResponse(c, N);
    const H = []; for (let k = 0; k < N; k++) H.push(freqResponseAt(c, 2 * Math.PI * k / N));
    let worst = 0;
    for (let n = 0; n < N; n++) { let acc = 0; for (let k = 0; k < N; k++) { const w = 2 * Math.PI * k * n / N; acc += H[k].re * Math.cos(w) - H[k].im * Math.sin(w); } const e = Math.abs(acc / N - h[n]); if (e > worst) worst = e; }
    check('DSP8 DF-II-T impulse response matches IDFT of H within 1e-6 (256 samples)', worst < 1e-6, 'max err=' + worst.toExponential(2)); }
  // 9
  { const cases = [[[C(0.5, 0.3), C(0.5, -0.3)], [C(-1, 0), C(-1, 0)]], [[C(0.6, 0), C(0.6, 0)], [C(0.2, 0.9), C(0.2, -0.9)]], [[C(0.9, 0), C(-0.4, 0)], [C(1, 0), C(-1, 0)]]];
    let worst = 0;
    for (const [p, z] of cases) {
      const c = coefsFromRoots({ poles: p, zeros: z, K: 0.37 }), rr = rootsFromCoefs(c);
      const err = (a, b) => Math.min(a[0].sub(b[0]).abs() + a[1].sub(b[1]).abs(), a[0].sub(b[1]).abs() + a[1].sub(b[0]).abs());
      worst = Math.max(worst, err(p, rr.poles), err(z, rr.zeros), Math.abs(rr.K - 0.37));
    }
    check('DSP9 roots→coefs→roots round trip within 1e-9 (incl. double root)', worst < 1e-9, 'max err=' + worst.toExponential(2)); }
  // 10
  { const roots = rootsFromCoefs(rbjCoefs('lp', 1000, 48000, 0.7071, 0));
    applyRootMove(roots, 'pole', 0, 1.2, 0.3, 1e-9);
    const r1 = poleRadii(coefsFromRoots(roots));
    applyRootMove(roots, 'pole', 1, -5, 0, 1e-9);
    const r2 = poleRadii(coefsFromRoots(roots));
    const ok = r1.every(x => x <= POLE_RMAX + 1e-12) && r2.every(x => x <= POLE_RMAX + 1e-12);
    check('DSP10 dragging a pole to |p|≥1 never yields |p|>0.9999', ok, 'radii: ' + r1.map(x => x.toFixed(6)).join(', ') + ' / ' + r2.map(x => x.toFixed(6)).join(', ')); }
  const passed = results.filter(r => r.ok).length;
  const summary = passed + '/' + results.length + ' checks passed';
  lines.push(summary); console.log(summary);
  return { passed, failed: results.length - passed, results, lines };
}
