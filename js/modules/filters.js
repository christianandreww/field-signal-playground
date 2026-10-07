/* Phase 5: higher-order Butterworth / Chebyshev-I filters (tab `filters`).
   Analog prototype -> frequency transform (LP/HP/BP/BS) -> bilinear with prewarping -> SOS cascade.
   No high-order polynomial is ever formed. Math lives in FSP.math.filters; DOM work only inside init/activate. */
(function () {
  'use strict';

  /* =====================================================================
     PURE MATH
     ===================================================================== */
  const TYPES = { butter: 'Butterworth', cheby1: 'Chebyshev I' };
  const RESPS = { lp: 'Low-pass', hp: 'High-pass', bp: 'Band-pass', bs: 'Band-stop' };
  const NYQ_FRAC = 0.49;   // edges are clamped below 0.49 fs (tan() prewarp blows up at Nyquist)

  const cprod = arr => arr.reduce((a, b) => a.mul(b), C(1, 0));
  const csqrt = z => z.sqrt();
  const cis = (a) => C(Math.cos(a), Math.sin(a));

  // Normalise/validate design parameters (throws Error with a user-readable message).
  function normParams(p) {
    const type = p.type === 'cheby1' ? 'cheby1' : 'butter';
    const resp = RESPS[p.resp] ? p.resp : 'lp';
    const n = Math.round(+p.order);
    if (!(n >= 1 && n <= 10)) throw new Error('Order must be 1–10');
    const fs = +p.fs;
    if (!(fs >= 1000 && fs <= 384000)) throw new Error('Sample rate out of range');
    const rp = type === 'cheby1' ? +p.rp : 0;
    if (type === 'cheby1' && !(rp > 0 && rp <= 12)) throw new Error('Ripple must be in (0, 12] dB');
    const band = resp === 'bp' || resp === 'bs';
    let f1 = +p.f1, f2 = band ? +p.f2 : +p.f1;
    if (!(f1 > 0) || !(f2 > 0)) throw new Error('Frequencies must be > 0');
    const notes = [];
    if (band && f1 > f2) { const t = f1; f1 = f2; f2 = t; notes.push('edges swapped so f1 < f2'); }
    const fmax = NYQ_FRAC * fs;
    if (f1 > fmax) { f1 = fmax; notes.push('f1 clamped to 0.49·fs'); }
    if (f2 > fmax) { f2 = fmax; notes.push('f2 clamped to 0.49·fs'); }
    if (band && !(f2 > f1 * 1.0001)) { f2 = Math.min(fmax, f1 * 1.05); if (!(f2 > f1 * 1.0001)) f1 = f2 / 1.05; notes.push('band widened to f2 ≈ 1.05·f1'); }
    return { type, resp, order: n, fs, f1, f2: band ? f2 : f1, rp, band, notes };
  }

  // Analog low-pass prototype (cutoff 1 rad/s; Chebyshev: ripple edge 1 rad/s). Same conventions as scipy buttap/cheb1ap.
  function prototype(type, n, rp) {
    const poles = []; let k = 1, eps = 0, mu = 0;
    if (type === 'butter') {
      for (let m = -n + 1; m < n; m += 2) poles.push(cis(Math.PI * m / (2 * n)).neg());
    } else {
      eps = Math.sqrt(Math.pow(10, 0.1 * rp) - 1); mu = Math.asinh(1 / eps) / n;
      for (let m = -n + 1; m < n; m += 2) { const th = Math.PI * m / (2 * n); poles.push(C(-Math.sinh(mu) * Math.cos(th), -Math.cosh(mu) * Math.sin(th))); }
      k = cprod(poles.map(p => p.neg())).re; if (n % 2 === 0) k /= Math.sqrt(1 + eps * eps);
    }
    return { poles, zeros: [], k, eps, mu };
  }

  // Analog LP prototype -> LP/HP/BP/BS with warped edges w1,w2 (rad/s). Returns {poles, zeros, k}.
  function transformAnalog(proto, resp, w1, w2) {
    const P = proto.poles, n = P.length; let k = proto.k, poles, zeros = [];
    const wo = resp === 'bp' || resp === 'bs' ? Math.sqrt(w1 * w2) : w1, bw = w2 - w1;
    const split = pl => { const a = [], b = []; pl.forEach(p => { const s = csqrt(p.mul(p).sub(C(wo * wo, 0))); a.push(p.add(s)); b.push(p.sub(s)); }); return a.concat(b); };
    const invK = () => C(1, 0).div(cprod(P.map(p => p.neg()))).re;
    if (resp === 'lp') { poles = P.map(p => p.scale(wo)); k *= Math.pow(wo, n); }
    else if (resp === 'hp') { poles = P.map(p => C(wo, 0).div(p)); for (let i = 0; i < n; i++) zeros.push(C(0, 0)); k *= invK(); }
    else if (resp === 'bp') { poles = split(P.map(p => p.scale(bw / 2))); for (let i = 0; i < n; i++) zeros.push(C(0, 0)); k *= Math.pow(bw, n); }
    else { poles = split(P.map(p => C(bw / 2, 0).div(p))); for (let i = 0; i < n; i++) zeros.push(C(0, wo)); for (let i = 0; i < n; i++) zeros.push(C(0, -wo)); k *= invK(); }
    return { poles, zeros, k };
  }

  function bilinear(an, fs) {
    const f2 = 2 * fs, T = C(f2, 0), map = s => T.add(s).div(T.sub(s));
    const zeros = an.zeros.map(map), poles = an.poles.map(map);
    for (let i = an.zeros.length; i < an.poles.length; i++) zeros.push(C(-1, 0));
    const k = an.k * cprod(an.zeros.map(z => T.sub(z))).div(cprod(an.poles.map(p => T.sub(p)))).re;
    return { zeros, poles, k };
  }

  // Split a root set into conjugate pairs (upper-half representative) and reals. Imag parts of "real" roots are zeroed.
  function splitRoots(arr) {
    const up = [], re = [];
    arr.forEach(z => { if (z.im > 1e-9 * Math.max(1, z.abs())) up.push(z); else if (z.im >= -1e-9 * Math.max(1, z.abs())) re.push(z.re); });
    return { up, re: re.sort((a, b) => a - b) };
  }

  // Digital zpk -> SOS cascade. Poles are paired conjugate-wise (reals together), sections ordered by increasing pole
  // radius (the most resonant section last, as in scipy's zpk2sos); the gain is spread evenly so every section has the
  // same gain magnitude at the reference frequency (DC for LP/BS, Nyquist for HP, centre for BP).
  function zpk2sos(dz, resp, wRef) {
    const pp = splitRoots(dz.poles), zz = splitRoots(dz.zeros);
    if (2 * pp.up.length + pp.re.length !== dz.poles.length || 2 * zz.up.length + zz.re.length !== dz.zeros.length) throw new Error('pole/zero set is not conjugate-symmetric');
    const pSecs = pp.up.map(p => ({ poles: [p, p.conj()], r: p.abs() }));
    const rp = pp.re.slice(), single = [];
    while (rp.length >= 2) { const a = rp.shift(), b = rp.shift(); pSecs.push({ poles: [C(a, 0), C(b, 0)], r: Math.max(Math.abs(a), Math.abs(b)) }); }
    if (rp.length) single.push({ poles: [C(rp[0], 0)], r: Math.abs(rp[0]) });
    pSecs.sort((a, b) => a.r - b.r);
    const secs = single.concat(pSecs);                         // a first-order section (if any) goes first
    const zPairs = zz.up.map(z => [z, z.conj()]), zs = zz.re.slice(), zSingle = [];
    while (zs.length >= 2) { if (zs[0] < -0.5 && zs[zs.length - 1] > 0.5) zPairs.push([C(zs.shift(), 0), C(zs.pop(), 0)]); else zPairs.push([C(zs.shift(), 0), C(zs.shift(), 0)]); }
    if (zs.length) zSingle.push(C(zs[0], 0));
    const n2 = secs.filter(s => s.poles.length === 2).length;
    if (zPairs.length !== n2 || zSingle.length !== single.length) throw new Error('cannot pair zeros with poles');
    let zi = 0;
    secs.forEach(s => { if (s.poles.length === 2) s.zeros = zPairs[zi++]; else s.zeros = zSingle.slice(); });
    // monic coefficients
    const wr = cis(wRef);
    const gains = secs.map(s => {
      const p = s.poles, z = s.zeros;
      s.c = z.length === 2 ? { b0: 1, b1: -(z[0].re + z[1].re), b2: z[0].mul(z[1]).re } : { b0: 1, b1: -z[0].re, b2: 0 };
      s.c.a1 = p.length === 2 ? -(p[0].re + p[1].re) : -p[0].re; s.c.a2 = p.length === 2 ? p[0].mul(p[1]).re : 0;
      return freqResponseAt(s.c, wRef).abs();
    });
    const r = dz.k * gains.reduce((a, b) => a * b, 1), N = secs.length, mag = Math.pow(Math.abs(r), 1 / N);
    secs.forEach((s, i) => { const sc = (mag / gains[i]) * (i === 0 && r < 0 ? -1 : 1); s.c.b0 *= sc; s.c.b1 *= sc; s.c.b2 *= sc; });
    return secs.map(s => ({ b0: s.c.b0, b1: s.c.b1, b2: s.c.b2, a1: s.c.a1, a2: s.c.a2, order: s.poles.length, poles: s.poles, zeros: s.zeros }));
  }

  function design(params) {
    const q = normParams(params), fs = q.fs;
    const warp = f => 2 * fs * Math.tan(Math.PI * f / fs);          // prewarped analog edge, rad/s
    const w1 = warp(q.f1), w2 = warp(q.f2);
    const proto = prototype(q.type, q.order, q.rp);
    const an = transformAnalog(proto, q.resp, w1, w2);
    const dz = bilinear(an, fs);
    const wRef = q.resp === 'hp' ? Math.PI : q.resp === 'bp' ? 2 * Math.atan(Math.sqrt(w1 * w2) / (2 * fs)) : 0;
    const sos = zpk2sos(dz, q.resp, wRef);
    return Object.assign({}, q, { w1, w2, proto, analog: an, digital: dz, sos, wRef });
  }

  /* ---- evaluation ---- */
  function cascadeAt(sos, w) { let h = C(1, 0); for (let i = 0; i < sos.length; i++) h = h.mul(freqResponseAt(sos[i], w)); return h; }
  function zpkAt(dz, w) { const z = cis(w); let h = C(dz.k, 0); dz.zeros.forEach(a => { h = h.mul(z.sub(a)); }); dz.poles.forEach(a => { h = h.div(z.sub(a)); }); return h; }
  const dbOf = (h, floor) => 20 * Math.log10(Math.max(h.abs(), Math.pow(10, (floor === undefined ? -300 : floor) / 20)));
  function groupDelaySec(c, w) {   // samples; NaN on an exact unit-circle zero
    const e1 = cis(-w), e2 = cis(-2 * w);
    const B = C(c.b0, 0).add(e1.scale(c.b1)).add(e2.scale(c.b2)), A = C(1, 0).add(e1.scale(c.a1)).add(e2.scale(c.a2));
    if (B.abs() < 1e-9 * (Math.abs(c.b0) + Math.abs(c.b1) + Math.abs(c.b2) + 1e-300) || A.abs() < 1e-12) return NaN;
    const Sb = e1.scale(c.b1).add(e2.scale(2 * c.b2)), Sa = e1.scale(c.a1).add(e2.scale(2 * c.a2));
    return Sb.div(B).re - Sa.div(A).re;
  }
  function logGrid(fmin, fmax, n) { const f = new Float64Array(n); for (let i = 0; i < n; i++) f[i] = fmin * Math.pow(fmax / fmin, i / (n - 1)); return f; }

  // Full analysis on a frequency grid: total + per-section magnitude (dB), unwrapped phase (deg), group delay (ms).
  function analyze(des, freqs, floorDb) {
    const n = freqs.length, S = des.sos.length, fs = des.fs;
    const mag = new Float64Array(n), phase = new Float64Array(n), gd = new Float64Array(n), secMag = des.sos.map(() => new Float64Array(n));
    const prev = new Float64Array(S), acc = new Float64Array(S);
    for (let i = 0; i < n; i++) {
      const w = 2 * Math.PI * freqs[i] / fs; let tot = C(1, 0), g = 0, ph = 0;
      for (let s = 0; s < S; s++) {
        const h = freqResponseAt(des.sos[s], w); tot = tot.mul(h); secMag[s][i] = dbOf(h, floorDb);
        const a = h.abs() < 1e-12 ? prev[s] : h.arg();
        if (i === 0) acc[s] = a; else { let d = a - prev[s]; while (d > Math.PI) d -= 2 * Math.PI; while (d < -Math.PI) d += 2 * Math.PI; acc[s] += d; }
        prev[s] = a; ph += acc[s]; g += groupDelaySec(des.sos[s], w);
      }
      mag[i] = dbOf(tot, floorDb); phase[i] = ph * 180 / Math.PI; gd[i] = g / fs * 1000;
    }
    return { freqs, mag, phase, gd, secMag };
  }

  // Step response of the cascade (DF-II-T, the same arithmetic as the worklet). N chosen from the slowest pole.
  function stepLength(des) {
    let rmax = 0; des.digital.poles.forEach(p => { rmax = Math.max(rmax, p.abs()); });
    const n = rmax < 1 && rmax > 0 ? Math.log(1e-5) / Math.log(rmax) : 2000;
    return Math.max(64, Math.min(48000, Math.ceil(n * 1.1)));
  }
  function stepResponse(des, N) {
    N = N || stepLength(des); const y = new Float64Array(N), st = des.sos.map(() => new Float64Array(2)), cf = des.sos.map(coefArray);
    for (let i = 0; i < N; i++) { let v = 1; for (let s = 0; s < cf.length; s++) v = DF2T.df2tStep(v, cf[s], st[s]); y[i] = v; }
    return y;
  }

  function maxPoleRadius(des) { let r = 0; des.sos.forEach(s => s.poles.forEach(p => { r = Math.max(r, p.abs()); })); return r; }

  // scipy-compatible SOS rows [b0,b1,b2,1,a1,a2] plus design metadata.
  function toJSON(des) {
    return JSON.stringify({
      design: { type: des.type, response: des.resp, order: des.order, fs: des.fs, f1: des.f1, f2: des.band ? des.f2 : undefined, rippleDb: des.type === 'cheby1' ? des.rp : undefined },
      format: 'sos rows [b0, b1, b2, a0=1, a1, a2], cascade order, scipy.signal.sosfilt compatible',
      sos: des.sos.map(s => [s.b0, s.b1, s.b2, 1, s.a1, s.a2]),
    }, null, 2);
  }

  const f4 = x => Number.isFinite(x) ? x.toPrecision(6) : '—';
  const cstr = z => Number.isFinite(z.re) && Number.isFinite(z.im) ? (z.im === 0 ? f4(z.re) : f4(z.re) + (z.im < 0 ? ' − j' : ' + j') + f4(Math.abs(z.im))) : '—';
  function workingLines(des) {
    const L = [], n = des.order, q = des;
    L.push(TYPES[q.type] + ' ' + RESPS[q.resp] + ', order ' + n + ', fs = ' + q.fs + ' Hz, ' + (q.band ? 'f1 = ' + f4(q.f1) + ' Hz, f2 = ' + f4(q.f2) + ' Hz' : 'fc = ' + f4(q.f1) + ' Hz') + (q.type === 'cheby1' ? ', ripple ' + q.rp + ' dB' : ''));
    L.push('');
    L.push('1. Analog prototype (low-pass, edge = 1 rad/s)');
    if (q.type === 'butter') L.push('   p_m = −exp(jπm/2n), m = −n+1, −n+3, …, n−1   (all |p| = 1)');
    else L.push('   ε = √(10^(Rp/10) − 1) = ' + f4(q.proto.eps) + ',  μ = asinh(1/ε)/n = ' + f4(q.proto.mu) + '\n   p_m = −sinh(μ + jθ_m), θ_m = πm/2n;  DC gain = ' + (n % 2 ? '1' : '1/√(1+ε²) = −' + q.rp + ' dB (even order)'));
    q.proto.poles.forEach((p, i) => L.push('   p' + (i + 1) + ' = ' + cstr(p)));
    L.push('');
    L.push('2. Prewarp + frequency transform (' + RESPS[q.resp] + ')');
    L.push('   Ω = 2·fs·tan(π f / fs):  Ω1 = ' + f4(q.w1) + ' rad/s' + (q.band ? ',  Ω2 = ' + f4(q.w2) + ' rad/s,  Ω0 = √(Ω1Ω2) = ' + f4(Math.sqrt(q.w1 * q.w2)) + ',  BW = ' + f4(q.w2 - q.w1) : ''));
    L.push('   ' + { lp: 's → s/Ω1', hp: 's → Ω1/s', bp: 's → (s² + Ω0²)/(BW·s)   (each prototype pole → 2 poles; order n gives 2n poles)', bs: 's → BW·s/(s² + Ω0²)   (each prototype pole → 2 poles; order n gives 2n poles)' }[q.resp]);
    L.push('   analog poles: ' + q.analog.poles.length + ', finite zeros: ' + q.analog.zeros.length + ', gain k = ' + f4(q.analog.k));
    L.push('');
    L.push('3. Bilinear transform  s = 2fs (z−1)/(z+1)  ⇒  z = (2fs + s)/(2fs − s)');
    L.push('   zeros at s = ∞ map to z = −1;  digital gain k = ' + f4(q.digital.k));
    L.push('');
    L.push('4. SOS cascade (' + q.sos.length + ' sections, poles ordered by increasing radius, gain spread evenly, unity-ish at ' + { lp: 'DC', hp: 'Nyquist', bp: 'centre', bs: 'DC' }[q.resp] + ')');
    q.sos.forEach((s, i) => {
      L.push('   #' + (i + 1) + ' b = [' + [s.b0, s.b1, s.b2].map(f4).join(', ') + ']  a = [1, ' + f4(s.a1) + ', ' + f4(s.a2) + ']');
      L.push('       poles ' + s.poles.map(cstr).join(' ; ') + '   |p| = ' + f4(s.poles[0].abs()));
    });
    L.push('');
    L.push('H(z) = ∏ᵢ (b0ᵢ + b1ᵢ z⁻¹ + b2ᵢ z⁻²)/(1 + a1ᵢ z⁻¹ + a2ᵢ z⁻²). A high-order polynomial is never formed (it would be numerically ill-conditioned for narrow bands or large n).');
    return L;
  }

  FSP.math.filters = { TYPES, RESPS, normParams, prototype, transformAnalog, bilinear, zpk2sos, design, cascadeAt, zpkAt, groupDelaySec, logGrid, analyze, stepLength, stepResponse, maxPoleRadius, toJSON, workingLines, dbOf };

  /* =====================================================================
     TESTS
     ===================================================================== */
  // Reference values from scipy 1.18 (scipy.signal.butter / cheby1, output='sos', sosfreqz): dB at the listed frequencies.
  const REF = [
    { ty: 'butter', r: 'lp', n: 8, rp: 0, fs: 48000, f: [1000, 1500, 4000], fe: [1000], db: [-3.0103, -28.305451, -97.843691], a: [[-1.7578526471777913, 0.7730210883760056], [-1.7887583504227402, 0.804193475715957], [-1.848819839796427, 0.8647732333138347], [-1.9336504795257299, 0.9503358732893509]] },
    { ty: 'cheby1', r: 'hp', n: 7, rp: 0.5, fs: 44100, f: [1000, 2000, 2600, 3000, 4000], fe: [2000], db: [-65.274423, -0.5, -0.016832, -0.446451, -0.095379] },
    { ty: 'butter', r: 'bp', n: 5, rp: 0, fs: 48000, f: [100, 400, 800, 1000, 1200, 1800, 4000, 10000], fe: [800, 1200], db: [-137.507543, -69.84936, -3.0103, 0, -3.0103, -50.188024, -98.25023, -146.2072] },
    { ty: 'cheby1', r: 'bp', n: 4, rp: 2, fs: 16000, f: [100, 150, 300, 390, 1850, 3400, 4000, 5100], fe: [300, 3400], db: [-55.479722, -39.782115, -2.0, -1.982078, -0.039737, -2.0, -18.318792, -39.367669] },
    { ty: 'butter', r: 'bs', n: 9, rp: 0, fs: 48000, f: [100, 450, 900, 1100, 1170, 1650, 4000], fe: [900, 1100], db: [0, 0, -3.0103, -3.0103, -0.000737, 0, 0] },
    { ty: 'cheby1', r: 'bs', n: 10, rp: 1, fs: 48000, f: [100, 500, 1000, 1300, 2500, 4000, 6000, 10000], fe: [1000, 4000], db: [-0.556997, -0.322824, -1.0, -84.272209, -153.905113, -1.0, -0.768158, -0.870208] },
    { ty: 'cheby1', r: 'lp', n: 10, rp: 0.1, fs: 48000, f: [100, 2500, 4000, 5000, 6500, 7500, 10000], fe: [5000], db: [-0.096368, -0.013, -0.089774, -0.1, -46.877846, -66.801605, -104.09695] },
  ];
  const mk = (ty, r, n, fs, f1, f2, rp) => design({ type: ty, resp: r, order: n, fs, f1, f2: f2 || f1, rp: rp || 1 });
  const atDb = (des, f) => dbOf(cascadeAt(des.sos, 2 * Math.PI * f / des.fs), -400);

  FSP.registerTests('filters', t => {
    // --- spec: Butterworth LP, fs 48k, fc 1k
    [[2, -24.48], [4, -48.92], [6, -73.38]].forEach(([n, ref]) => {
      const d = mk('butter', 'lp', n, 48000, 1000), a = atDb(d, 1000), b = atDb(d, 4000);
      t.check('Butterworth LP N=' + n + ': −3.010 dB at fc, ' + ref + ' dB at 4 kHz', t.near(a, -3.010, 0.01) && t.near(b, ref, 0.1), a.toFixed(4) + ' / ' + b.toFixed(3));
    });
    // --- spec: Chebyshev-I order 4, 1 dB
    { const d = mk('cheby1', 'lp', 4, 48000, 1000, 1000, 1), a = atDb(d, 1000), b = atDb(d, 4000), dc = atDb(d, 1e-6);
      t.check('Chebyshev-I N=4 1 dB: −1.000 dB at fc, −60.58 dB at 4 kHz', t.near(a, -1, 0.01) && t.near(b, -60.58, 0.3), a.toFixed(4) + ' / ' + b.toFixed(3));
      t.check('Chebyshev-I even order: DC gain = −ripple (−1 dB)', t.near(dc, -1, 0.005), dc.toFixed(4));
      const d3 = mk('cheby1', 'lp', 5, 48000, 1000, 1000, 1); t.check('Chebyshev-I odd order: DC gain 0 dB', t.near(atDb(d3, 1e-6), 0, 0.005)); }
    // --- spec: stability for every order/type/response, incl. extreme edges
    { let ok = true, worst = 0, cnt = 0, bad = '';
      for (const ty of ['butter', 'cheby1']) for (const r of ['lp', 'hp', 'bp', 'bs']) for (let n = 1; n <= 10; n++) for (const [f1, f2] of [[1000, 3000], [20, 40], [200, 20000], [1000, 1010]]) {
        let d; try { d = mk(ty, r, n, 48000, f1, f2, 3); } catch (e) { ok = false; bad = ty + r + n + ' threw ' + e.message; continue; }
        cnt++; const rm = maxPoleRadius(d); worst = Math.max(worst, rm); if (!(rm < 1) || !d.sos.every(s => [s.b0, s.b1, s.b2, s.a1, s.a2].every(Number.isFinite))) { ok = false; bad = ty + r + n + ' ' + f1 + '/' + f2 + ' r=' + rm; }
      }
      t.check('all poles inside unit circle (' + cnt + ' designs: both types, LP/HP/BP/BS, N=1–10, narrow/wide/extreme edges)', ok, bad || 'max |p| = ' + worst.toFixed(9)); }
    // --- spec: cascade = product of sections; independent check vs direct zpk evaluation
    { let worst = 0, worstProd = 0;
      for (const [ty, r, n, f1, f2] of [['butter', 'lp', 10, 1000], ['cheby1', 'hp', 9, 3000], ['butter', 'bp', 5, 800, 1200], ['cheby1', 'bs', 6, 500, 4000], ['cheby1', 'lp', 4, 1000]]) {
        const d = mk(ty, r, n, 48000, f1, f2, 1);
        for (let i = 0; i < 300; i++) { const f = 20 * Math.pow(1200, i / 299), w = 2 * Math.PI * f / 48000, h = cascadeAt(d.sos, w), z = zpkAt(d.digital, w);
          let pm = 1; d.sos.forEach(s => { pm *= freqResponseAt(s, w).abs(); });
          worst = Math.max(worst, h.sub(z).abs()); worstProd = Math.max(worstProd, Math.abs(h.abs() - pm)); }
      }
      t.check('cascade response equals independent zpk evaluation (abs error < 1e-12)', worst < 1e-12, 'max err ' + worst.toExponential(2));
      t.check('|cascade| equals product of section magnitudes (error < 1e-12)', worstProd < 1e-12, 'max err ' + worstProd.toExponential(2)); }
    // --- scipy references (sosfreqz) including BP/BS and orders 7–10
    REF.forEach(c => {
      const d = mk(c.ty, c.r, c.n, c.fs, c.fe[0], c.fe[1] || c.fe[0], c.rp || 1); let worst = 0, n = 0;
      c.f.forEach((f, i) => { if (c.db[i] > -200) { worst = Math.max(worst, Math.abs(atDb(d, f) - c.db[i])); n++; } });
      t.check('scipy ' + c.ty + ' ' + c.r + ' N=' + c.n + (c.rp ? ' Rp=' + c.rp : '') + ' fs=' + c.fs + ': ' + n + ' dB points within 0.002 dB', worst < 0.002, 'max |Δ| = ' + worst.toExponential(2) + ' dB');
      t.check('scipy ' + c.ty + ' ' + c.r + ' N=' + c.n + ': ' + Math.ceil((c.r === 'lp' || c.r === 'hp' ? c.n : 2 * c.n) / 2) + ' sections', d.sos.length === Math.ceil((c.r === 'lp' || c.r === 'hp' ? c.n : 2 * c.n) / 2));
    });
    { const c = REF[0], d = mk(c.ty, c.r, c.n, c.fs, 1000), mine = d.sos.map(s => [s.a1, s.a2]);
      let worst = 0; mine.forEach((m, i) => { worst = Math.max(worst, Math.abs(m[0] - c.a[i][0]), Math.abs(m[1] - c.a[i][1])); });
      t.check('scipy butter(8,1 kHz) SOS denominators match in order, 1e-12', worst < 1e-12 && mine.length === c.a.length, 'max |Δ| = ' + worst.toExponential(2)); }
    // --- band edges and DC/Nyquist behaviour
    { const bp = mk('butter', 'bp', 4, 48000, 800, 1200), c = Math.sqrt(800 * 1200);
      t.check('Butterworth BP: −3.01 dB at both edges (800/1200 Hz)', t.near(atDb(bp, 800), -3.0103, 0.01) && t.near(atDb(bp, 1200), -3.0103, 0.01));
      const bs = mk('cheby1', 'bs', 5, 48000, 1000, 2000, 0.5);
      t.check('Chebyshev BS: −ripple at edges, DC and Nyquist 0 dB (odd N) , deep notch inside', t.near(atDb(bs, 1000), -0.5, 0.01) && t.near(atDb(bs, 2000), -0.5, 0.01) && t.near(atDb(bs, 1e-6), 0, 0.01) && atDb(bs, 1414.2) < -60, atDb(bs, 1000).toFixed(3) + ' ' + atDb(bs, 1414.2).toFixed(1));
      const hp = mk('butter', 'hp', 6, 48000, 1000);
      t.check('Butterworth HP: −3.01 dB at fc, 0 dB at Nyquist, < −70 dB at 250 Hz', t.near(atDb(hp, 1000), -3.0103, 0.01) && t.near(atDb(hp, 23999.9999), 0, 0.01) && atDb(hp, 250) < -70); }
    // --- time domain / group delay / pairing
    { const d = mk('butter', 'lp', 6, 48000, 1000), y = stepResponse(d);
      t.check('step response of LP settles to DC gain 1', t.near(y[y.length - 1], 1, 1e-4), 'final ' + y[y.length - 1].toFixed(6) + ', N=' + y.length);
      const dh = mk('butter', 'hp', 4, 48000, 1000), yh = stepResponse(dh);
      t.check('step response of HP settles to 0', Math.abs(yh[yh.length - 1]) < 1e-4); }
    { const d = mk('butter', 'lp', 4, 48000, 1000);
      // numeric derivative of the phase of the total cascade vs analytic group delay
      const w = 2 * Math.PI * 500 / 48000, h = 1e-5, dph = (cascadeAt(d.sos, w + h).arg() - cascadeAt(d.sos, w - h).arg()) / (2 * h);
      const tg = d.sos.reduce((a, s) => a + FSP.math.filters.groupDelaySec(s, w), 0);
      t.check('analytic group delay equals −dφ/dω (rel 1e-6)', t.rel(tg, -dph, 1e-6), tg.toFixed(4) + ' samples'); }
    { const a = FSP.math.filters, cases = [['butter', 'lp', 7], ['cheby1', 'bp', 3], ['butter', 'bs', 1], ['cheby1', 'hp', 1]]; let ok = true;
      cases.forEach(([ty, r, n]) => { const d = mk(ty, r, n, 48000, 1000, 2500, 1); const exp = r === 'lp' || r === 'hp' ? Math.ceil(n / 2) : n; if (d.sos.length !== exp) ok = false; });
      t.check('section counts: ceil(N/2) for LP/HP, N for BP/BS', ok);
      const d = mk('butter', 'lp', 5, 48000, 1000), rs = d.sos.map(s => s.poles[0].abs());
      t.check('sections ordered by increasing pole radius; odd order has one 1st-order section', rs.every((r, i) => i === 0 || r >= rs[i - 1] - 1e-12) && d.sos.filter(s => s.order === 1).length === 1 && a.maxPoleRadius(d) < 1); }
    // --- JSON export + validation
    { const d = mk('cheby1', 'bp', 3, 48000, 500, 2000, 1), j = JSON.parse(FSP.math.filters.toJSON(d));
      t.check('SOS JSON round-trips (scipy layout [b0,b1,b2,1,a1,a2])', j.sos.length === d.sos.length && j.sos.every((r, i) => r.length === 6 && r[3] === 1 && r[4] === d.sos[i].a1 && r[0] === d.sos[i].b0) && j.design.rippleDb === 1);
      let threw = 0; [{ order: 0 }, { order: 11 }, { fs: 0 }, { f1: -5 }, { type: 'cheby1', rp: 0 }].forEach(p => { try { design(Object.assign({ type: 'butter', resp: 'lp', order: 4, fs: 48000, f1: 1000, f2: 1000, rp: 1 }, p)); } catch (e) { threw++; } });
      t.check('invalid parameters are rejected', threw === 5, threw + '/5'); }
    { const d = design({ type: 'butter', resp: 'bp', order: 3, fs: 48000, f1: 30000, f2: 100, rp: 1 });
      t.check('edges beyond 0.49·fs are clamped and swapped, design stays stable', d.f1 < d.f2 && d.f2 <= 0.49 * 48000 + 1e-9 && maxPoleRadius(d) < 1 && d.notes.length >= 1, d.f1.toFixed(1) + '/' + d.f2.toFixed(1)); }
  });

  /* =====================================================================
     UI
     ===================================================================== */
  const M = () => FSP.math.filters;
  const S = {
    p: { type: 'butter', resp: 'lp', order: 4, fs: 48000, fc: 1000, f1: 500, f2: 2000, rp: 1, floor: 100 },
    des: null, showSec: [], ui: {}, active: false, built: false, plots: {},
  };
  const A = { ctx: null, mode: '', bus: null, chain: null, master: null, src: null, srcGain: null, srcKind: 'noise', timer: 0, sweepTimer: 0, rebuildTimer: 0, starting: false, fs: 0, des: null };

  function cssVar(name, fb) { try { const v = getComputedStyle(document.documentElement).getPropertyValue(name).trim(); return v || fb; } catch (e) { return fb; } }
  function theme() { return { bg: cssVar('--panel', '#121826'), text: cssVar('--text', '#e6edf3'), muted: cssVar('--muted', '#8b98ab'), grid: cssVar('--border', '#243047'), accent: cssVar('--accent', '#4cc9f0'), a2: cssVar('--accent2', '#f9a03f'), pink: cssVar('--pink', '#f72585'), ok: cssVar('--ok', '#2ecc71'), bad: cssVar('--bad', '#ff5c5c') }; }
  const fmtF = f => !Number.isFinite(f) ? '—' : f >= 1000 ? (+(f / 1000).toPrecision(3)) + 'k' : String(+f.toPrecision(3));
  const fmtDbS = x => Number.isFinite(x) ? (x <= -299 ? '< −300' : x.toFixed(3)) : '—';

  function niceStep(span, target) { const raw = span / Math.max(1, target), p = Math.pow(10, Math.floor(Math.log10(raw))), m = raw / p; return (m < 1.5 ? 1 : m < 3.5 ? 2 : m < 7.5 ? 5 : 10) * p; }
  function ticks(a, b, target) { const st = niceStep(b - a, target), out = []; for (let v = Math.ceil(a / st - 1e-9) * st; v <= b + st * 1e-9; v += st) out.push(Math.abs(v) < st * 1e-9 ? 0 : v); return out; }

  // Draws axes/grid; returns mapping helpers. o: {xmin,xmax,ymin,ymax,logx,xfmt,yfmt,xlabel,ylabel,equal}
  function frame(ctx, w, h, o) {
    const T = theme(), l = 46, r = 10, tp = 8, b = 30, pw = Math.max(10, w - l - r), ph = Math.max(10, h - tp - b);
    ctx.clearRect(0, 0, w, h); ctx.fillStyle = T.bg; ctx.fillRect(0, 0, w, h);
    const lx = o.logx ? Math.log10 : (x => x);
    const X = x => l + (lx(x) - lx(o.xmin)) / (lx(o.xmax) - lx(o.xmin)) * pw, Y = y => tp + (1 - (y - o.ymin) / (o.ymax - o.ymin)) * ph;
    ctx.font = '11px ' + cssVar('--mono', 'monospace'); ctx.lineWidth = 1;
    let xt = [];
    if (o.logx) { for (let d = Math.floor(Math.log10(o.xmin)); d <= Math.ceil(Math.log10(o.xmax)); d++) [1, 2, 5].forEach(m => { const v = m * Math.pow(10, d); if (v >= o.xmin * 0.9999 && v <= o.xmax * 1.0001) xt.push([v, m === 1]); }); }
    else xt = ticks(o.xmin, o.xmax, Math.max(3, Math.floor(pw / 80))).map(v => [v, true]);
    const yt = ticks(o.ymin, o.ymax, Math.max(3, Math.floor(ph / 36)));
    ctx.textBaseline = 'middle'; ctx.textAlign = 'right';
    yt.forEach(v => { const y = Y(v); ctx.strokeStyle = v === 0 ? T.muted : T.grid; ctx.globalAlpha = v === 0 ? 0.8 : 0.6; ctx.beginPath(); ctx.moveTo(l, y); ctx.lineTo(l + pw, y); ctx.stroke(); ctx.globalAlpha = 1; ctx.fillStyle = T.muted; ctx.fillText(o.yfmt ? o.yfmt(v) : String(+v.toPrecision(4)), l - 4, y); });
    ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    xt.forEach(([v, major]) => { const x = X(v); ctx.strokeStyle = T.grid; ctx.globalAlpha = major ? 0.7 : 0.3; ctx.beginPath(); ctx.moveTo(x, tp); ctx.lineTo(x, tp + ph); ctx.stroke(); ctx.globalAlpha = 1; if (major || !o.logx) { ctx.fillStyle = T.muted; ctx.fillText(o.xfmt ? o.xfmt(v) : String(+v.toPrecision(4)), x, tp + ph + 3); } });
    ctx.strokeStyle = T.grid; ctx.strokeRect(l + 0.5, tp + 0.5, pw, ph);
    ctx.fillStyle = T.muted; ctx.textAlign = 'right'; ctx.textBaseline = 'bottom'; if (o.xlabel) ctx.fillText(o.xlabel, l + pw, h - 1);
    if (o.ylabel) { ctx.save(); ctx.translate(11, tp + 2); ctx.rotate(-Math.PI / 2); ctx.textAlign = 'right'; ctx.textBaseline = 'top'; ctx.fillText(o.ylabel, 0, -2); ctx.restore(); }
    return { X, Y, l, t: tp, pw, ph, T, clip() { ctx.save(); ctx.beginPath(); ctx.rect(l, tp, pw, ph); ctx.clip(); } };
  }
  function polyline(ctx, xs, ys, F, color, lw, alpha) {
    ctx.strokeStyle = color; ctx.lineWidth = lw || 1.5; ctx.globalAlpha = alpha === undefined ? 1 : alpha; ctx.beginPath(); let pen = false;
    for (let i = 0; i < xs.length; i++) { const y = ys[i]; if (!Number.isFinite(y)) { pen = false; continue; } const px = F.X(xs[i]), py = F.Y(y); if (pen) ctx.lineTo(px, py); else { ctx.moveTo(px, py); pen = true; } }
    ctx.stroke(); ctx.globalAlpha = 1;
  }
  function vline(ctx, F, x, color) { const px = F.X(x); ctx.strokeStyle = color; ctx.globalAlpha = 0.6; ctx.setLineDash([4, 3]); ctx.beginPath(); ctx.moveTo(px, F.t); ctx.lineTo(px, F.t + F.ph); ctx.stroke(); ctx.setLineDash([]); ctx.globalAlpha = 1; }
  function secColor(i, n) { return 'hsl(' + Math.round(360 * i / Math.max(n, 1) + 20) + ',70%,62%)'; }
  function autoRange(vals, pad, lo0) {
    let lo = Infinity, hi = -Infinity; vals.forEach(v => { if (Number.isFinite(v)) { if (v < lo) lo = v; if (v > hi) hi = v; } });
    if (!Number.isFinite(lo)) { lo = 0; hi = 1; } if (lo0 !== undefined) lo = Math.min(lo, lo0);
    if (hi - lo < 1e-12) { hi = lo + 1; } const d = (hi - lo) * pad; return [lo - d, hi + d];
  }

  function params() { const p = S.p; return { type: p.type, resp: p.resp, order: p.order, fs: p.fs, f1: (p.resp === 'bp' || p.resp === 'bs') ? p.f1 : p.fc, f2: p.f2, rp: p.rp }; }

  function recompute() {
    const msg = S.ui.msg; msg.textContent = ''; msg.className = 'note';
    let d; try { d = M().design(params()); } catch (e) { msg.textContent = e.message; msg.className = 'msg err'; return false; }
    if (d.notes.length) { msg.textContent = 'Adjusted: ' + d.notes.join('; '); msg.className = 'msg warn'; }
    if (S.des && S.des.sos.length !== d.sos.length) S.showSec = d.sos.map(() => false);
    if (S.showSec.length !== d.sos.length) S.showSec = d.sos.map(() => false);
    S.des = d; return true;
  }

  function freqRange(d) { return [10, d.fs / 2]; }

  function drawAll() {
    if (!S.built || !S.des) return;
    ['mag', 'phase', 'gd', 'step', 'splane', 'zplane'].forEach(k => drawPlot(k));
  }
  function drawPlot(k) {
    const pl = S.plots[k], d = S.des; if (!pl || !d || !pl.cv.cv.clientWidth) return;
    const { ctx, w, h } = pl.cv.prep(); pl.draw(ctx, w, h, d);
  }

  function getResp(d) {
    if (S.resp && S.resp.des === d && S.resp.floor === S.p.floor) return S.resp;
    const [a, b] = freqRange(d), fr = M().logGrid(a, b, 1100), r = M().analyze(d, fr, -300);
    return (S.resp = Object.assign(r, { des: d, floor: S.p.floor }));
  }

  function drawMag(ctx, w, h, d) {
    const r = getResp(d), [fa, fb] = freqRange(d), top = 6, bot = -S.p.floor;
    const F = frame(ctx, w, h, { xmin: fa, xmax: fb, ymin: bot, ymax: top, logx: true, xfmt: fmtF, yfmt: v => String(v), xlabel: 'Hz', ylabel: 'dB' });
    F.clip();
    [d.f1, d.band ? d.f2 : null].forEach(f => { if (f) vline(ctx, F, f, F.T.a2); });
    d.sos.forEach((s, i) => { if (S.showSec[i]) polyline(ctx, r.freqs, r.secMag[i], F, secColor(i, d.sos.length), 1, 0.9); });
    polyline(ctx, r.freqs, r.mag, F, F.T.accent, 2.2);
    ctx.restore();
  }
  function drawPhase(ctx, w, h, d) {
    const r = getResp(d), [fa, fb] = freqRange(d), [lo, hi] = autoRange(r.phase, 0.05);
    const F = frame(ctx, w, h, { xmin: fa, xmax: fb, ymin: lo, ymax: hi, logx: true, xfmt: fmtF, yfmt: v => String(Math.round(v)), xlabel: 'Hz', ylabel: 'phase (°)' });
    F.clip(); [d.f1, d.band ? d.f2 : null].forEach(f => { if (f) vline(ctx, F, f, F.T.a2); }); polyline(ctx, r.freqs, r.phase, F, F.T.pink, 2); ctx.restore();
  }
  function drawGd(ctx, w, h, d) {
    const r = getResp(d), [fa, fb] = freqRange(d), fin = Array.from(r.gd).filter(Number.isFinite).sort((a, b) => a - b);
    const med = fin.length ? fin[Math.floor(fin.length / 2)] : 1, cap = Math.max(fin.length ? fin[fin.length - 1] : 1, 0), lim = Math.min(cap, Math.max(10 * Math.abs(med), 0.05));
    const vals = Array.from(r.gd).map(v => Number.isFinite(v) ? Math.min(v, lim * 1.5) : NaN), [lo, hi] = autoRange(vals, 0.08, 0);
    const F = frame(ctx, w, h, { xmin: fa, xmax: fb, ymin: lo, ymax: hi, logx: true, xfmt: fmtF, yfmt: v => String(+v.toPrecision(3)), xlabel: 'Hz', ylabel: 'delay (ms)' });
    F.clip(); [d.f1, d.band ? d.f2 : null].forEach(f => { if (f) vline(ctx, F, f, F.T.a2); }); polyline(ctx, r.freqs, vals, F, F.T.ok, 2); ctx.restore();
  }
  function drawStep(ctx, w, h, d) {
    if (!d._step) d._step = M().stepResponse(d);
    const y = d._step, N = y.length, ms = new Float64Array(N); for (let i = 0; i < N; i++) ms[i] = i / d.fs * 1000;
    const [lo, hi] = autoRange(y, 0.08, 0), xmax = ms[N - 1] || 1;
    const F = frame(ctx, w, h, { xmin: 0, xmax, ymin: Math.min(lo, 0), ymax: Math.max(hi, 1.0), logx: false, xfmt: v => String(+v.toPrecision(3)), yfmt: v => String(+v.toPrecision(3)), xlabel: 'ms', ylabel: 'step' });
    F.clip(); polyline(ctx, ms, y, F, F.T.accent, 1.8); ctx.restore();
  }
  function marker(ctx, x, y, kind, color, r) { ctx.strokeStyle = color; ctx.lineWidth = 1.8; ctx.beginPath(); if (kind === 'x') { ctx.moveTo(x - r, y - r); ctx.lineTo(x + r, y + r); ctx.moveTo(x - r, y + r); ctx.lineTo(x + r, y - r); } else ctx.arc(x, y, r, 0, 2 * Math.PI); ctx.stroke(); }
  function planeDraw(ctx, w, h, poles, zeros, scale, unit, fmt, title) {
    const T = theme(), m = 30, size = Math.min(w, h) - 2 * m, cx = w / 2, cy = h / 2;
    ctx.clearRect(0, 0, w, h); ctx.fillStyle = T.bg; ctx.fillRect(0, 0, w, h);
    const X = x => cx + x / scale * size / 2, Y = y => cy - y / scale * size / 2;
    ctx.font = '11px ' + cssVar('--mono', 'monospace'); ctx.lineWidth = 1; ctx.strokeStyle = T.grid; ctx.fillStyle = T.muted;
    ctx.strokeRect(m - 0.5, cy - size / 2 - 0.5, size + 1, size + 1);
    ctx.beginPath(); ctx.moveTo(m, cy); ctx.lineTo(w - m, cy); ctx.moveTo(cx, cy - size / 2); ctx.lineTo(cx, cy + size / 2); ctx.stroke();
    if (unit) { ctx.strokeStyle = T.muted; ctx.beginPath(); ctx.arc(cx, cy, size / 2 / scale, 0, 2 * Math.PI); ctx.stroke(); }
    ctx.textAlign = 'left'; ctx.textBaseline = 'top'; ctx.fillText(title, m + 3, cy - size / 2 + 3);
    ctx.textAlign = 'center'; ctx.textBaseline = 'top'; [-1, -0.5, 0.5, 1].forEach(f => { ctx.fillText(fmt(f * scale), X(f * scale), cy + 2); });
    ctx.textAlign = 'left'; ctx.textBaseline = 'middle'; [-1, -0.5, 0.5, 1].forEach(f => { if (f !== 0) ctx.fillText(fmt(f * scale), cx + 3, Y(f * scale)); });
    const groups = (arr) => { const g = []; arr.forEach(z => { const e = g.find(q => Math.hypot(q.z.re - z.re, q.z.im - z.im) < 1e-6 * Math.max(1, scale)); if (e) e.n++; else g.push({ z, n: 1 }); }); return g; };
    ctx.save(); ctx.beginPath(); ctx.rect(m, cy - size / 2, size, size); ctx.clip();
    groups(zeros).forEach(g => { const x = X(g.z.re), y = Y(g.z.im); marker(ctx, x, y, 'o', T.ok, 5); if (g.n > 1) { ctx.fillStyle = T.ok; ctx.textAlign = 'left'; ctx.textBaseline = 'bottom'; ctx.fillText('×' + g.n, x + 7, y - 4); } });
    groups(poles).forEach(g => { const x = X(g.z.re), y = Y(g.z.im); marker(ctx, x, y, 'x', T.pink, 5); if (g.n > 1) { ctx.fillStyle = T.pink; ctx.textAlign = 'left'; ctx.textBaseline = 'top'; ctx.fillText('×' + g.n, x + 7, y + 4); } });
    ctx.restore();
  }
  function drawS(ctx, w, h, d) {
    const a = d.analog, all = a.poles.concat(a.zeros), mx = Math.max(1e-9, ...all.map(z => z.abs())) * 1.25 / (2 * Math.PI);
    const scale = [1, 2, 5, 10].map(m => m * Math.pow(10, Math.floor(Math.log10(mx)))).find(v => v >= mx);
    planeDraw(ctx, w, h, a.poles.map(p => p.scale(1 / (2 * Math.PI))), a.zeros.map(p => p.scale(1 / (2 * Math.PI))), scale, false, v => fmtF(v), 's-plane (σ, jω in Hz)');
  }
  function drawZ(ctx, w, h, d) { planeDraw(ctx, w, h, d.digital.poles, d.digital.zeros, 1.25, true, v => String(+v.toPrecision(2)), 'z-plane'); }

  function updateHud() {
    const d = S.des, hud = S.ui.hud; if (!d) return;
    const g = f => atDbF(d, f), rows = [];
    const edge = d.band ? [['|H| at f1', g(d.f1)], ['|H| at f2', g(d.f2)]] : [['|H| at fc', g(d.f1)]];
    rows.push(...edge, ['DC gain (dB)', g(1e-6)], ['Nyquist gain (dB)', g(d.fs / 2 * 0.999999)], [d.band ? 'Centre gain (dB)' : 'Order / sections', d.band ? g(Math.sqrt(d.f1 * d.f2)) : null]);
    hud.textContent = '';
    rows.forEach(([k, v]) => { if (v === null) { hud.appendChild(FSP.ui.el('div', null, FSP.ui.el('span', { text: k }), FSP.ui.el('span', { text: d.order + ' / ' + d.sos.length }))); return; } hud.appendChild(FSP.ui.el('div', null, FSP.ui.el('span', { text: k }), FSP.ui.el('span', { text: fmtDbS(v) + ' dB' }))); });
    hud.appendChild(FSP.ui.el('div', null, FSP.ui.el('span', { text: 'Poles / sections' }), FSP.ui.el('span', { text: d.digital.poles.length + ' / ' + d.sos.length })));
    hud.appendChild(FSP.ui.el('div', null, FSP.ui.el('span', { text: 'Max |pole|' }), FSP.ui.el('span', { text: Number.isFinite(M().maxPoleRadius(d)) ? M().maxPoleRadius(d).toFixed(8) : '—' })));
  }
  function atDbF(d, f) { return M().dbOf(M().cascadeAt(d.sos, 2 * Math.PI * f / d.fs), -300); }

  function updateTable() {
    const d = S.des, tb = S.ui.tbody; if (!d) return; tb.textContent = '';
    d.sos.forEach((s, i) => {
      const ck = FSP.ui.el('input', { type: 'checkbox', 'aria-label': 'Show section ' + (i + 1) + ' magnitude' }); ck.checked = !!S.showSec[i];
      ck.addEventListener('change', () => { S.showSec[i] = ck.checked; drawPlot('mag'); });
      const sw = FSP.ui.el('span', { text: ' ■', style: 'color:' + secColor(i, d.sos.length) });
      const fp = Math.abs(s.poles[0].arg()) * d.fs / (2 * Math.PI), rad = s.poles[0].abs();
      const cells = [s.b0, s.b1, s.b2, s.a1, s.a2, rad, fp].map((v, j) => FSP.ui.el('td', { class: 'mono', style: 'padding:2px 8px 2px 0', text: j === 6 ? (s.order === 2 ? FSP.fmtNum(fp, 5) : s.poles[0].re < 0 ? 'fs/2' : 'DC') : FSP.fmtNum(v, 8) }));
      const tr = FSP.ui.el('tr', null, FSP.ui.el('td', { style: 'white-space:nowrap;padding-right:8px' }, ck, FSP.ui.el('span', { text: ' ' + (i + 1) }), sw), ...cells); tb.appendChild(tr);
    });
  }

  function updateAll(fromUser) {
    if (!S.built) return;
    if (recompute()) {
      S.resp = null; updateHud(); updateTable(); S.ui.working.set(M().workingLines(S.des)); S.ui.json.textContent = M().toJSON(S.des); drawAll(); audioPush();
    }
    if (fromUser) FSP.state.touch();
  }

  function syncVisibility() {
    const p = S.p, band = p.resp === 'bp' || p.resp === 'bs';
    S.ui.fc.el.style.display = band ? 'none' : ''; S.ui.f1.el.style.display = band ? '' : 'none'; S.ui.f2.el.style.display = band ? '' : 'none';
    S.ui.rp.el.style.display = p.type === 'cheby1' ? '' : 'none';
    S.ui.orderNote.textContent = band ? 'Band designs: order n gives 2n poles (n sections), as in scipy.' : 'Low/high-pass: n poles, ⌈n/2⌉ sections.';
  }

  /* ---------- state ---------- */
  function bindState() {
    FSP.state.bind('filters', {
      get() { const p = S.p; return { ty: p.type, rs: p.resp, n: p.order, fs: p.fs, fc: +p.fc.toPrecision(6), f1: +p.f1.toPrecision(6), f2: +p.f2.toPrecision(6), rp: +p.rp.toPrecision(4), fl: p.floor }; },
      set(o) {
        const p = S.p, num = (v, lo, hi) => { const x = FSP.parseSI(v); return Number.isFinite(x) && x >= lo && x <= hi ? x : null; };
        if (o.ty === 'butter' || o.ty === 'cheby1') p.type = o.ty; if (RESPS_OK(o.rs)) p.resp = o.rs;
        let x; if ((x = num(o.n, 1, 10)) !== null) p.order = Math.round(x); if ((x = num(o.fs, 1000, 384000)) !== null) p.fs = x;
        if ((x = num(o.fc, 1, 192000)) !== null) p.fc = x; if ((x = num(o.f1, 1, 192000)) !== null) p.f1 = x; if ((x = num(o.f2, 1, 192000)) !== null) p.f2 = x;
        if ((x = num(o.rp, 0.001, 12)) !== null) p.rp = x; if ((x = num(o.fl, 20, 300)) !== null) p.floor = Math.round(x);
        if (S.built) pushControls();
      },
    });
  }
  function RESPS_OK(v) { return v === 'lp' || v === 'hp' || v === 'bp' || v === 'bs'; }
  function pushControls() {
    const p = S.p, u = S.ui;
    u.type.value = p.type; u.resp.value = p.resp; u.order.set(p.order, true); u.fs.value = String(p.fs); if (!u.fs.value) { u.fs.value = '48000'; p.fs = 48000; }
    u.fc.set(p.fc, true); u.f1.set(p.f1, true); u.f2.set(p.f2, true); u.rp.set(p.rp, true); u.floor.value = String(p.floor); if (!u.floor.value) { u.floor.value = '100'; p.floor = 100; }
    syncVisibility(); updateAll(false);
  }

  /* ---------- audio ---------- */
  function aMsg(text, cls) { const m = S.ui.amsg; m.textContent = text || ''; m.className = text ? 'msg ' + (cls || '') : 'note'; }
  function aState() {
    const b = S.ui.astate, on = A.ctx && A.ctx.state === 'running';
    b.textContent = A.ctx ? A.ctx.state + ' @ ' + A.ctx.sampleRate + ' Hz · ' + A.mode : 'stopped'; b.className = 'badge ' + (on ? 'ok' : A.ctx ? 'warn' : '');
    S.ui.start.disabled = !!A.ctx || A.starting; S.ui.stop.disabled = !A.ctx;
  }
  function audioDesign() {
    const q = params(); q.fs = A.fs; try { return M().design(q); } catch (e) { return null; }
  }
  function makeChain(des) {
    const ctx = A.ctx, out = ctx.createGain(), nodes = []; out.gain.value = 0;
    let prev = A.bus;
    des.sos.forEach(s => {
      let nd;
      if (A.mode === 'AudioWorklet') { nd = new AudioWorkletNode(ctx, 'biquad-df2t', { numberOfInputs: 1, numberOfOutputs: 1, outputChannelCount: [1] }); nd.port.postMessage({ type: 'coefs', c: coefArray(s) }); }
      else nd = ctx.createIIRFilter([s.b0, s.b1, s.b2], [1, s.a1, s.a2]);   // one node per section: 3 coefficients each, far below the 20 limit
      prev.connect(nd); prev = nd; nodes.push(nd);
    });
    prev.connect(out); out.connect(A.master);
    const t = ctx.currentTime; out.gain.setValueAtTime(0, t); out.gain.linearRampToValueAtTime(1, t + 0.04);
    return { nodes, out, n: des.sos.length };
  }
  function dropChain(ch) {
    if (!ch) return; const t = A.ctx.currentTime; ch.out.gain.setValueAtTime(ch.out.gain.value, t); ch.out.gain.linearRampToValueAtTime(0, t + 0.04);
    setTimeout(() => { try { ch.nodes.forEach(n => n.disconnect()); ch.out.disconnect(); } catch (e) { /* already gone */ } }, 150);
  }
  function audioPush() {
    if (!A.ctx) return; const des = audioDesign(); if (!des) return; A.des = des;
    if (A.mode === 'AudioWorklet' && A.chain && A.chain.n === des.sos.length) { des.sos.forEach((s, i) => A.chain.nodes[i].port.postMessage({ type: 'coefs', c: coefArray(s) })); return; }
    clearTimeout(A.rebuildTimer);
    A.rebuildTimer = setTimeout(() => { if (!A.ctx) return; try { const old = A.chain; A.chain = makeChain(A.des); dropChain(old); } catch (e) { aMsg('Could not build filter chain: ' + e.message, 'err'); } }, A.chain && A.mode !== 'AudioWorklet' ? 40 : 0);
  }
  function makeSource() {
    const ctx = A.ctx, kind = S.ui.src.value; stopSource(); A.srcKind = kind;
    if (kind === 'noise') {
      const n = Math.round(ctx.sampleRate * 4), buf = ctx.createBuffer(1, n, ctx.sampleRate), d = buf.getChannelData(0); for (let i = 0; i < n; i++) d[i] = Math.random() * 2 - 1;
      const s = ctx.createBufferSource(); s.buffer = buf; s.loop = true; s.connect(A.srcGain); s.start(); A.src = s;
    } else {
      const o = ctx.createOscillator(); o.type = 'sine'; o.connect(A.srcGain); A.src = o;
      if (kind === 'sine') o.frequency.value = S.ui.sineF.get();
      else {
        const top = Math.min(20000, ctx.sampleRate * 0.45), run = () => { if (!A.ctx) return; const t = ctx.currentTime; o.frequency.cancelScheduledValues(t); o.frequency.setValueAtTime(20, t); o.frequency.exponentialRampToValueAtTime(top, t + 8); };
        run(); A.sweepTimer = setInterval(run, 8200);
      }
      o.start();
    }
    S.ui.sineRow.style.display = kind === 'sine' ? '' : 'none';
  }
  function stopSource() { clearInterval(A.sweepTimer); A.sweepTimer = 0; if (A.src) { try { A.src.stop(); } catch (e) { /* not started */ } try { A.src.disconnect(); } catch (e) { /* ok */ } A.src = null; } }
  async function startAudio() {
    if (A.ctx || A.starting) return; A.starting = true; aState();
    const AC = typeof window !== 'undefined' && (window.AudioContext || window.webkitAudioContext);
    if (!AC) { aMsg('Web Audio API is not available in this browser.', 'err'); A.starting = false; aState(); return; }
    let ctx; try { ctx = new AC({ sampleRate: S.p.fs }); } catch (e) { ctx = new AC(); }
    A.ctx = ctx; A.fs = ctx.sampleRate; A.mode = '';
    ctx.addEventListener('statechange', aState);
    aMsg(ctx.sampleRate !== S.p.fs ? 'AudioContext runs at ' + ctx.sampleRate + ' Hz (design fs ' + S.p.fs + ' Hz); the audio filter is redesigned at the context rate with the same edge frequencies.' : '', 'warn');
    A.master = ctx.createGain(); A.master.gain.value = S.ui.level.get(); A.master.connect(ctx.destination);
    A.srcGain = ctx.createGain(); A.srcGain.gain.value = 1; A.bus = ctx.createGain(); A.srcGain.connect(A.bus);
    try {
      if (!ctx.audioWorklet) throw new Error('AudioWorklet unavailable');
      const url = URL.createObjectURL(new Blob([WORKLET_SRC], { type: 'application/javascript' }));
      try { await ctx.audioWorklet.addModule(url); } finally { URL.revokeObjectURL(url); }
      A.mode = 'AudioWorklet';
    } catch (e) { A.mode = 'IIRFilterNode fallback'; aMsg('AudioWorklet unavailable (' + e.message + '). Using a chain of IIRFilterNodes, rebuilt with a 40 ms crossfade on changes.', 'warn'); }
    if (!A.ctx) { A.starting = false; return; }              // stopped while loading
    A.des = audioDesign();
    if (!A.des) { aMsg('Current parameters cannot be designed at this sample rate.', 'err'); stopAudio(); A.starting = false; return; }
    try { A.chain = makeChain(A.des); makeSource(); await ctx.resume(); } catch (e) { aMsg('Audio start failed: ' + e.message, 'err'); }
    A.starting = false; aState();
  }
  function stopAudio() {
    clearTimeout(A.rebuildTimer); stopSource();
    const ctx = A.ctx; A.ctx = null; A.chain = null; A.bus = null; A.master = null; A.srcGain = null; A.starting = false;
    if (ctx) { try { ctx.close(); } catch (e) { /* ignore */ } }
    if (S.built) aState();
  }

  /* ---------- tab ---------- */
  function copyText(text, done) {
    const fallback = () => { const ta = FSP.ui.el('textarea', { style: 'position:fixed;left:-9999px;top:0', 'aria-hidden': 'true' }); ta.value = text; document.body.appendChild(ta); ta.select(); let ok = false; try { ok = document.execCommand('copy'); } catch (e) { ok = false; } ta.remove(); done(ok); };
    if (typeof navigator !== 'undefined' && navigator.clipboard && navigator.clipboard.writeText) navigator.clipboard.writeText(text).then(() => done(true), fallback); else fallback();
  }

  function init(panel) {
    const U = FSP.ui, u = S.ui, p = S.p;
    const layout = U.el('div', { class: 'layout' }), controls = U.el('aside', { class: 'controls' }), stage = U.el('div', { class: 'stage' });
    layout.appendChild(controls); layout.appendChild(stage); panel.appendChild(layout);

    const fd = U.fieldset(controls, 'Design');
    u.type = U.select(fd, 'Family', [['butter', 'Butterworth'], ['cheby1', 'Chebyshev I']], p.type, v => { p.type = v; syncVisibility(); updateAll(true); });
    u.resp = U.select(fd, 'Response', [['lp', 'Low-pass'], ['hp', 'High-pass'], ['bp', 'Band-pass'], ['bs', 'Band-stop']], p.resp, v => { p.resp = v; syncVisibility(); updateAll(true); });
    u.order = U.slider(fd, { label: 'Order n', min: 1, max: 10, step: 1, value: p.order, digits: 2, onInput: v => { p.order = Math.round(v); updateAll(true); } });
    u.fc = U.slider(fd, { label: 'fc (Hz)', min: 10, max: 48000, value: p.fc, log: true, digits: 5, onInput: v => { p.fc = v; updateAll(true); } });
    u.f1 = U.slider(fd, { label: 'f1 low', min: 10, max: 48000, value: p.f1, log: true, digits: 5, onInput: v => { p.f1 = v; updateAll(true); } });
    u.f2 = U.slider(fd, { label: 'f2 high', min: 10, max: 48000, value: p.f2, log: true, digits: 5, onInput: v => { p.f2 = v; updateAll(true); } });
    u.rp = U.slider(fd, { label: 'Ripple dB', min: 0.01, max: 6, value: p.rp, log: true, digits: 4, onInput: v => { p.rp = v; updateAll(true); } });
    u.fs = U.select(fd, 'fs (Hz)', [['8000', '8000'], ['16000', '16000'], ['44100', '44100'], ['48000', '48000'], ['96000', '96000']], String(p.fs), v => { p.fs = +v; updateAll(true); });
    u.floor = U.select(fd, 'dB floor', [['60', '−60 dB'], ['100', '−100 dB'], ['150', '−150 dB'], ['200', '−200 dB']], String(p.floor), v => { p.floor = +v; S.resp = null; drawAll(); FSP.state.touch(); });
    u.orderNote = U.el('div', { class: 'note' }); fd.appendChild(u.orderNote);
    u.msg = U.el('div', { class: 'note' }); fd.appendChild(u.msg);
    fd.appendChild(U.el('div', { class: 'note', text: 'Edges are the −3.01 dB points (Butterworth) or the ripple-band edges (Chebyshev I), prewarped through the bilinear transform.' }));

    const af = U.fieldset(controls, 'Audio test');
    u.src = U.select(af, 'Source', [['noise', 'White noise'], ['sine', 'Sine'], ['sweep', 'Log sweep 20 Hz→'] ], 'noise', () => { if (A.ctx) makeSource(); });
    u.sineF = U.slider(af, { label: 'Sine Hz', min: 20, max: 20000, value: 1000, log: true, digits: 5, onInput: v => { if (A.ctx && A.src && A.srcKind === 'sine') A.src.frequency.setTargetAtTime(v, A.ctx.currentTime, 0.01); } });
    u.sineRow = u.sineF.el; u.sineRow.style.display = 'none';
    u.level = U.slider(af, { label: 'Level', min: 0, max: 0.5, step: 0.005, value: 0.15, digits: 3, onInput: v => { if (A.master) A.master.gain.setTargetAtTime(v, A.ctx.currentTime, 0.02); } });
    const br = U.el('div', { class: 'row' }); af.appendChild(br);
    u.start = U.button(br, 'Start audio', startAudio); u.start.className = 'btn primary'; u.stop = U.button(br, 'Stop', stopAudio); u.stop.className = 'btn';
    u.astate = U.el('span', { class: 'badge', text: 'stopped' }); br.appendChild(u.astate);
    u.amsg = U.el('div', { class: 'note' }); af.appendChild(u.amsg);
    af.appendChild(U.el('div', { class: 'note', text: 'One worklet biquad per section (5 ms coefficient ramp); IIRFilterNode chain fallback. Keep the level low: noise is broadband.' }));
    u.stop.disabled = true;

    const mkPlot = (parent, key, title, height, draw, square) => {
      const wrap = U.el('div', { class: 'canvas-wrap' }); parent.appendChild(wrap);
      const cv = U.canvas(wrap, { height }); cv.cv.setAttribute('role', 'img'); cv.cv.setAttribute('aria-label', title);
      if (square) { cv.cv.style.height = 'auto'; cv.cv.style.aspectRatio = '1 / 1'; }
      cv.onResize(() => { if (S.active) drawPlot(key); });
      S.plots[key] = { cv, draw }; return wrap;
    };
    const magBox = U.el('div', { class: 'note', text: 'Magnitude (dB): thick line = cascade; tick the sections in the table below to overlay their individual responses. Dashed lines mark the edge frequencies.' });
    stage.appendChild(magBox); mkPlot(stage, 'mag', 'Magnitude response', 280, drawMag);
    const g1 = U.el('div', { class: 'grid2' }); stage.appendChild(g1);
    mkPlot(g1, 'phase', 'Phase response', 220, drawPhase); mkPlot(g1, 'gd', 'Group delay', 220, drawGd);
    mkPlot(stage, 'step', 'Step response', 220, drawStep);
    const g2 = U.el('div', { class: 'grid2' }); stage.appendChild(g2);
    mkPlot(g2, 's' + 'plane', 'Analog s-plane pole-zero plot', 300, drawS, true); mkPlot(g2, 'zplane', 'Digital z-plane pole-zero plot', 300, drawZ, true);
    u.hud = U.el('div', { class: 'hud', 'aria-live': 'off' }); stage.appendChild(u.hud);

    const tf = U.fieldset(stage, 'Second-order sections');
    const tw = U.el('div', { style: 'overflow-x:auto' }); tf.appendChild(tw);
    const tbl = U.el('table', { style: 'border-collapse:collapse;font-size:12px;width:100%' });
    const th = ['#/show', 'b0', 'b1', 'b2', 'a1', 'a2', '|p|', 'f_pole Hz'];
    tbl.appendChild(U.el('thead', null, U.el('tr', null, ...th.map(x => U.el('th', { text: x, style: 'text-align:left;padding:2px 8px 2px 0;color:var(--muted);font-weight:400;white-space:nowrap' })))));
    u.tbody = U.el('tbody'); tbl.appendChild(u.tbody); tw.appendChild(tbl);
    const tr = U.el('div', { class: 'row' }); tf.appendChild(tr);
    U.button(tr, 'Show all', () => { S.showSec = S.showSec.map(() => true); updateTable(); drawPlot('mag'); }).className = 'btn';
    U.button(tr, 'Hide all', () => { S.showSec = S.showSec.map(() => false); updateTable(); drawPlot('mag'); }).className = 'btn';
    const cp = U.button(tr, 'Copy SOS as JSON', () => copyText(S.ui.json.textContent, ok => { u.copyMsg.textContent = ok ? 'Copied to clipboard.' : 'Clipboard blocked: select the JSON below and copy manually.'; })); cp.className = 'btn primary';
    u.copyMsg = U.el('span', { class: 'note' }); tr.appendChild(u.copyMsg);
    tf.appendChild(U.el('div', { class: 'note', text: 'Each row is a Direct Form II Transposed biquad H = (b0 + b1 z⁻¹ + b2 z⁻²)/(1 + a1 z⁻¹ + a2 z⁻²), a0 = 1. JSON uses the scipy layout [b0, b1, b2, 1, a1, a2].' }));
    const jd = U.el('details', { class: 'working-box' }, U.el('summary', { text: 'SOS JSON' })); u.json = U.el('pre', { class: 'working' }); jd.appendChild(u.json); tf.appendChild(jd);

    u.working = U.working(stage);
    S.built = true;
    bindState(); pushControls();
  }

  FSP.registerTab({
    id: 'filters', title: 'Butterworth / Chebyshev',
    init,
    activate() { S.active = true; if (S.built) { syncVisibility(); S.resp = null; drawAll(); } },
    deactivate() { S.active = false; stopAudio(); },
  });
})();
