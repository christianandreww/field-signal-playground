/* FSP module: z-Transform & FIR Design (NTU IE3014, Oppenheim & Schafer "Discrete-Time Signal Processing" 3e).
   Sub-topics: sampling & aliasing, difference equation & z-transform (poles/zeros, ROC, inverse z by partial fractions),
   linear-phase FIR design (window method incl. Kaiser, frequency sampling), filter structures (DF-I, DF-II, TDF-II,
   cascade SOS, parallel, folded linear-phase FIR) and DFT properties (+ quantisation SNR).
   Conventions (O&S): H(z) = Σ b_k z^-k / Σ a_k z^-k with a_0 = 1 after normalisation, i.e.
   y[n] = Σ b_k x[n-k] − Σ_{k≥1} a_k y[n-k]. Frequencies ω in rad/sample, often shown as ω/π.
   Pure math lives in FSP.math.dspdesign (no DOM). All DOM work happens inside init/activate. */
(function () {
  'use strict';
  const PI = Math.PI, TAU = 2 * Math.PI;
  const fin = x => typeof x === 'number' && Number.isFinite(x);

  /* ======================= complex helpers ======================= */
  const C = (re, im) => ({ re: re, im: im || 0 });
  const cadd = (a, b) => C(a.re + b.re, a.im + b.im);
  const csub = (a, b) => C(a.re - b.re, a.im - b.im);
  const cmul = (a, b) => C(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
  const cdiv = (a, b) => { const d = b.re * b.re + b.im * b.im; return C((a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d); };
  const cscale = (a, s) => C(a.re * s, a.im * s);
  const cabs = a => Math.hypot(a.re, a.im);
  const carg = a => Math.atan2(a.im, a.re);
  const cconj = a => C(a.re, -a.im);
  const cpowi = (a, n) => { // integer power (n may be negative)
    if (n === 0) return C(1, 0); const r = Math.pow(cabs(a), n), th = carg(a) * n; return C(r * Math.cos(th), r * Math.sin(th));
  };
  const cx = { C, cadd, csub, cmul, cdiv, cscale, cabs, carg, cconj, cpowi };

  /* ======================= polynomials ======================= */
  // Horner on descending coefficients (real numbers or complex objects) at complex z -> {p, dp}
  function hornerD(c, z) {
    let p = C(0, 0), dp = C(0, 0);
    for (let i = 0; i < c.length; i++) { dp = cadd(cmul(dp, z), p); const ci = typeof c[i] === 'number' ? C(c[i], 0) : c[i]; p = cadd(cmul(p, z), ci); }
    return { p, dp };
  }
  // Multiply complex polynomials (any consistent order)
  function polyMulC(a, b) { const out = []; for (let i = 0; i < a.length + b.length - 1; i++) out.push(C(0, 0)); a.forEach((x, i) => b.forEach((y, j) => { out[i + j] = cadd(out[i + j], cmul(x, y)); })); return out; }
  const toC = arr => arr.map(v => (typeof v === 'number' ? C(v, 0) : v));
  // Expand Π (x − r_i) as descending coefficients (complex) -> monic
  function polyFromRoots(roots) { let p = [C(1, 0)]; roots.forEach(r => { p = polyMulC(p, [C(1, 0), C(-r.re, -r.im)]); }); return p; }

  /* Polynomial roots: Aberth–Ehrlich simultaneous iteration (cubically convergent for simple roots), Newton polish,
     then roots that form a tight cluster are merged into one root of multiplicity m placed at the cluster centroid
     (the centroid of the m computed copies of an m-fold root is accurate to ~machine precision even though each copy
     is only accurate to ~eps^(1/m)). Input: real or complex coefficients, highest power first. Exact zero roots
     (trailing zero coefficients) are extracted exactly. Returns [{r: complex, m: multiplicity}]. */
  function rootsMult(coefDesc, opts) {
    opts = opts || {};
    let c = toC(coefDesc.slice());
    while (c.length && cabs(c[0]) === 0) c.shift();          // leading zeros reduce the degree
    if (c.length === 0) throw new Error('zero polynomial');
    let nZero = 0; while (c.length > 1 && cabs(c[c.length - 1]) === 0) { c.pop(); nZero++; }
    const n = c.length - 1, out = [];
    if (nZero) out.push({ r: C(0, 0), m: nZero });
    if (n === 0) return out;
    const c0 = c[0]; c = c.map(v => cdiv(v, c0));            // monic
    let z;
    if (n === 1) z = [C(-c[1].re, -c[1].im)];
    else {
      const R0 = Math.pow(cabs(c[n]), 1 / n) || 1; z = [];
      for (let k = 0; k < n; k++) { const th = TAU * k / n + 0.4; z.push(C(R0 * Math.cos(th), R0 * Math.sin(th))); }
      for (let it = 0; it < 800; it++) {
        let maxStep = 0;
        for (let i = 0; i < n; i++) {
          const h = hornerD(c, z[i]); if (cabs(h.p) === 0) continue;
          const ratio = cdiv(h.p, h.dp); let s = C(0, 0);
          for (let j = 0; j < n; j++) if (j !== i) { const d = csub(z[i], z[j]); if (cabs(d) > 0) s = cadd(s, cdiv(C(1, 0), d)); }
          const w = cdiv(ratio, csub(C(1, 0), cmul(ratio, s)));
          if (!fin(w.re) || !fin(w.im)) continue;
          z[i] = csub(z[i], w); maxStep = Math.max(maxStep, cabs(w) / Math.max(cabs(z[i]), 1e-300));
        }
        if (maxStep < 1e-15) break;
      }
    }
    // cluster (single linkage) with a relative tolerance; a cluster of m copies is accepted as an m-fold root only if
    // p, p', …, p^(m−1) all (nearly) vanish at the refined centroid (relative to their condition sums); otherwise the
    // roots are kept as distinct simple roots.
    const tol = opts.clusterTol === undefined ? 2e-4 : opts.clusterTol;
    const derivs = [c]; for (let j = 1; j < n; j++) { const p = derivs[j - 1], dpc = []; for (let i = 0; i < p.length - 1; i++) dpc.push(cscale(p[i], p.length - 1 - i)); derivs.push(dpc); }
    const relVal = (p, zz) => { let s = 0; const az = cabs(zz); for (let i = 0; i < p.length; i++) s += cabs(p[i]) * Math.pow(az, p.length - 1 - i); return cabs(hornerD(p, zz).p) / Math.max(s, 1e-300); };
    const newton = (p, z0, it) => { let zz = z0; for (let k = 0; k < it; k++) { const h = hornerD(p, zz); if (cabs(h.dp) === 0) break; const st = cdiv(h.p, h.dp); if (!fin(st.re) || !fin(st.im)) break; zz = csub(zz, st); if (cabs(st) <= 1e-17 * Math.max(1, cabs(zz))) break; } return zz; };
    const used = new Array(n).fill(false);
    for (let i = 0; i < n; i++) {
      if (used[i]) continue; const grp = [i]; used[i] = true;
      for (let g = 0; g < grp.length; g++) for (let j = 0; j < n; j++) if (!used[j]) {
        const sc = Math.max(1, cabs(z[grp[g]])); if (cabs(csub(z[grp[g]], z[j])) < tol * sc) { used[j] = true; grp.push(j); }
      }
      const m = grp.length;
      if (m === 1) { out.push({ r: newton(c, z[i], 4), m: 1 }); continue; }
      let cen = C(0, 0); grp.forEach(k => { cen = cadd(cen, z[k]); }); cen = cscale(cen, 1 / m);
      const ref = newton(derivs[m - 1], cen, 6); // p^(m−1) has a simple root at an m-fold root of p
      let ok = cabs(csub(ref, cen)) < tol * Math.max(1, cabs(cen));
      for (let j = 0; j < m && ok; j++) if (relVal(derivs[j], ref) > 1e-7) ok = false;
      if (ok) out.push({ r: ref, m: m }); else grp.forEach(k => out.push({ r: newton(c, z[k], 4), m: 1 }));
    }
    // clean tiny imaginary parts for real-coefficient input; snap conjugate partners exactly
    const realIn = coefDesc.every(v => typeof v === 'number');
    if (realIn) {
      out.forEach(o => { if (Math.abs(o.r.im) <= 1e-10 * Math.max(1, Math.abs(o.r.re))) o.r = C(o.r.re, 0); });
      out.forEach(o => { if (o.r.im > 0) { let best = null, bd = Infinity; out.forEach(q => { if (q.r.im < 0 && q.m === o.m) { const d = cabs(csub(q.r, cconj(o.r))); if (d < bd) { bd = d; best = q; } } }); if (best && bd < 1e-6 * Math.max(1, cabs(o.r))) { const avg = C((o.r.re + best.r.re) / 2, (o.r.im - best.r.im) / 2); o.r = avg; best.r = cconj(avg); } } });
    }
    return out;
  }
  // Flat list of roots (each repeated m times)
  function polyRoots(coefDesc, opts) { const out = []; rootsMult(coefDesc, opts).forEach(o => { for (let k = 0; k < o.m; k++) out.push(o.r); }); return out; }

  /* ======================= coefficient lists / H(z) ======================= */
  // "1, -0.5 0.25; 1/3" -> numbers (fractions p/q allowed). null on junk/empty/too long.
  function parseList(s, maxLen) {
    if (typeof s !== 'string') return null;
    const parts = s.trim().replace(/^\[|\]$/g, '').split(/[\s,;]+/).filter(Boolean);
    if (!parts.length || parts.length > (maxLen || 64)) return null;
    const out = [];
    for (const p of parts) {
      let v; const m = /^([+-]?(?:\d+\.?\d*|\.\d+)(?:e[+-]?\d+)?)\/((?:\d+\.?\d*|\.\d+)(?:e[+-]?\d+)?)$/i.exec(p);
      if (m) v = Number(m[1]) / Number(m[2]); else v = Number(p);
      if (!fin(v)) return null; out.push(v);
    }
    return out;
  }
  // Normalise (b, a): a_0 -> 1, strip trailing zeros. -> {b, a, M, N, err}
  function normBA(b, a) {
    if (!Array.isArray(b) || !Array.isArray(a) || !b.length || !a.length) return { err: 'Enter at least one b and one a coefficient.' };
    if (!b.concat(a).every(fin)) return { err: 'Coefficients must be finite numbers.' };
    if (a[0] === 0) return { err: 'a₀ must be non-zero (it multiplies y[n]).' };
    const a0 = a[0]; let bb = b.map(v => v / a0), aa = a.map(v => v / a0);
    while (bb.length > 1 && bb[bb.length - 1] === 0) bb.pop();
    while (aa.length > 1 && aa[aa.length - 1] === 0) aa.pop();
    if (bb.every(v => v === 0)) return { err: 'All b coefficients are zero: H(z) = 0.' };
    return { b: bb, a: aa, M: bb.length - 1, N: aa.length - 1, a0: a0 };
  }
  // Poles/zeros in the z-plane (including those at the origin), gain K: H(z) = K z^(N−M) Π(z − z_i)/Π(z − p_i)
  function pzk(b, a) {
    const q = normBA(b, a); if (q.err) throw new Error(q.err);
    const L = Math.max(q.M, q.N), num = q.b.concat(new Array(L - q.M).fill(0)), den = q.a.concat(new Array(L - q.N).fill(0));
    let D = 0; while (num[D] === 0) D++;  // leading zeros of b = pure delay
    const zeros = rootsMult(num), poles = rootsMult(den);
    return { zeros, poles, K: q.b[D], delay: D, b: q.b, a: q.a, M: q.M, N: q.N };
  }
  // Direct recursion y[n] = Σ b_k x[n−k] − Σ_{k≥1} a_k y[n−k] (a normalised).
  function filter(b, a, x) {
    const q = normBA(b, a); if (q.err) throw new Error(q.err); const B = q.b, A = q.a, y = new Float64Array(x.length);
    for (let n = 0; n < x.length; n++) { let s = 0; for (let k = 0; k < B.length && k <= n; k++) s += B[k] * x[n - k]; for (let k = 1; k < A.length && k <= n; k++) s -= A[k] * y[n - k]; y[n] = s; }
    return y;
  }
  const impulse = (b, a, N) => { const x = new Float64Array(N); x[0] = 1; return filter(b, a, x); };
  const stepResp = (b, a, N) => filter(b, a, new Float64Array(N).fill(1));
  // H(e^{jω}) = Σ b_k e^{-jωk} / Σ a_k e^{-jωk}
  function freqAt(b, a, w) {
    let nr = 0, ni = 0, dr = 0, di = 0;
    for (let k = 0; k < b.length; k++) { nr += b[k] * Math.cos(w * k); ni -= b[k] * Math.sin(w * k); }
    for (let k = 0; k < a.length; k++) { dr += a[k] * Math.cos(w * k); di -= a[k] * Math.sin(w * k); }
    return cdiv(C(nr, ni), C(dr, di));
  }
  function freqResp(b, a, npts) {
    const w = new Float64Array(npts), mag = new Float64Array(npts), ph = new Float64Array(npts);
    for (let i = 0; i < npts; i++) { w[i] = PI * i / (npts - 1); const h = freqAt(b, a, w[i]); mag[i] = cabs(h); ph[i] = carg(h); }
    return { w, mag, ph };
  }
  // Stability of the causal system: all poles strictly inside |z| = 1
  function stability(poles) {
    let rmax = 0; poles.forEach(p => { rmax = Math.max(rmax, cabs(p.r)); });
    const onCircle = poles.some(p => Math.abs(cabs(p.r) - 1) < 1e-9);
    return { rmax, verdict: rmax < 1 - 1e-9 ? 'stable' : (onCircle && rmax <= 1 + 1e-9 ? 'marginal' : 'unstable') };
  }

  /* ---------- partial fractions (O&S §3.3.2) ----------
     H(z) = Σ_{r=0}^{M−N} B_r z^−r + Σ_i Σ_{m=1}^{s_i} C_{i,m} / (1 − d_i z^−1)^m
     Work in w = z^−1: B(w) = Q(w)A(w) + R(w); A(w) = a0 Π (1 − d_i w)^{s_i}. For pole d of order s:
     G(w) = (1 − d w)^s R(w)/A(w) = R(w)/D(w), D = a0 Π_{others}; with u = w − 1/d, (1 − d w) = −d u, so
     C_m = g_{s−m} (−d)^{m−s}, where g_j are the Taylor coefficients of G about w = 1/d (O&S eq. 3.47 in series form). */
  function polyDivAsc(num, den) { // ascending-power real arrays; returns {q, r} ascending
    const n = num.slice().reverse(), d = den.slice().reverse(); // descending in w
    const dq = n.length - d.length; if (dq < 0) return { q: [], r: num.slice() };
    const q = new Array(dq + 1).fill(0), r = n.slice();
    for (let i = 0; i <= dq; i++) { const f = r[i] / d[0]; q[i] = f; for (let j = 0; j < d.length; j++) r[i + j] -= f * d[j]; }
    return { q: q.reverse(), r: r.slice(dq + 1).reverse() };
  }
  // Taylor coefficients (ascending) of polynomial p (ascending complex) about point w0, first K terms
  function taylorAt(pAsc, w0, K) {
    let c = pAsc.slice().reverse(); const out = []; // descending
    for (let k = 0; k < K; k++) {
      if (!c.length) { out.push(C(0, 0)); continue; }
      const q = []; let acc = C(0, 0);
      for (let i = 0; i < c.length; i++) { acc = cadd(cmul(acc, w0), c[i]); if (i < c.length - 1) q.push(acc); }
      out.push(acc); c = q;
    }
    return out;
  }
  function partialFractions(b, a) {
    const q = normBA(b, a); if (q.err) throw new Error(q.err);
    if (q.N === 0) return { direct: q.b.slice(), terms: [], clusters: [], b: q.b, a: q.a };
    const dv = polyDivAsc(q.b, q.a), direct = dv.q.length ? dv.q : [];
    let R = dv.q.length ? dv.r : q.b.slice(); while (R.length < q.N) R.push(0); R = R.slice(0, q.N);
    const clusters = rootsMult(q.a); // poles d_i (all non-zero since a_N ≠ 0)
    const terms = [];
    clusters.forEach((cl, i) => {
      const d = cl.r, s = cl.m, w0 = cdiv(C(1, 0), d);
      let Dw = [C(1, 0)]; // ascending in w: Π_{others} (1 − d_j w)^{s_j}
      clusters.forEach((o, j) => { if (j !== i) for (let k = 0; k < o.m; k++) Dw = polyMulC(Dw, [C(1, 0), C(-o.r.re, -o.r.im)]); });
      const Rt = taylorAt(toC(R), w0, s), Dt = taylorAt(Dw, w0, s);
      const g = []; // series division R/D
      for (let k = 0; k < s; k++) { let acc = Rt[k]; for (let j = 1; j <= k; j++) acc = csub(acc, cmul(Dt[j], g[k - j])); g.push(cdiv(acc, Dt[0])); }
      const md = C(-d.re, -d.im);
      for (let m = 1; m <= s; m++) terms.push({ cluster: i, d: d, order: m, mult: s, C: cmul(g[s - m], cpowi(md, m - s)) });
    });
    return { direct, terms, clusters, b: q.b, a: q.a };
  }
  // ROC options: distinct pole radii r_1 < … < r_K → K+1 annuli. idx K = causal (|z| > r_K), 0 = anticausal.
  function rocOptions(pf) {
    const radii = []; pf.clusters.forEach(c => { const r = cabs(c.r); if (!radii.some(x => Math.abs(x - r) < 1e-9 * Math.max(1, r))) radii.push(r); });
    radii.sort((x, y) => x - y); const K = radii.length, out = [];
    for (let i = 0; i <= K; i++) {
      const lo = i === 0 ? 0 : radii[i - 1], hi = i === K ? Infinity : radii[i];
      const kind = i === K ? 'causal (right-sided)' : (i === 0 ? 'anticausal (left-sided)' : 'two-sided');
      out.push({ idx: i, lo, hi, kind, stable: lo < 1 && hi > 1, label: (i === 0 ? '|z| < ' + fmtR(hi) : i === K ? '|z| > ' + fmtR(lo) : fmtR(lo) + ' < |z| < ' + fmtR(hi)) + ' — ' + kind });
    }
    return out;
  }
  const fmtR = r => String(+r.toPrecision(4));
  // binomial-type factor (n+1)(n+2)…(n+m−1)/(m−1)!  (valid for negative n too)
  function binomPoly(n, m) { let v = 1; for (let j = 1; j < m; j++) v *= (n + j) / j; return v; }
  // Closed-form h[n] for the chosen ROC (rocLo = inner radius of the annulus): poles with |d| ≤ rocLo are right-sided.
  function inverseZAt(pf, rocLo, n) {
    let s = C(0, 0);
    if (n >= 0 && n < pf.direct.length) s = C(pf.direct[n], 0);
    pf.terms.forEach(t => {
      const right = cabs(t.d) <= rocLo * (1 + 1e-9) + 1e-300;
      if (right ? n >= 0 : n <= -1) { const v = cscale(cmul(t.C, cpowi(t.d, n)), binomPoly(n, t.order) * (right ? 1 : -1)); s = cadd(s, v); }
    });
    return s.re;
  }
  // Residual of the difference equation for a candidate h[n]: max |Σ a_k h[n−k] − b_n| over n in [n0, n1]
  function diffEqResidual(b, a, hfun, n0, n1) {
    const q = normBA(b, a); let mx = 0;
    for (let n = n0; n <= n1; n++) { let s = 0; for (let k = 0; k < q.a.length; k++) s += q.a[k] * hfun(n - k); const bn = n >= 0 && n < q.b.length ? q.b[n] : 0; mx = Math.max(mx, Math.abs(s - bn)); }
    return mx;
  }

  /* ======================= sampling & aliasing ======================= */
  // Apparent (lowest-alias) frequency of a real tone F sampled at Fs: Fa = |F − k Fs|, k = round(F/Fs), in [0, Fs/2].
  // sign = sign(F − k Fs): when negative the reconstructed cosine has its phase negated.
  function alias(F, Fs) {
    if (!(Fs > 0) || !fin(F) || F < 0) return null;
    const k = Math.round(F / Fs), d = F - k * Fs, Fa = Math.abs(d);
    const nyq = Math.abs(Fa - Fs / 2) < 1e-12 * Fs;
    return { Fa, k, sign: d < 0 ? -1 : 1, atNyquist: nyq, aliased: F > Fs / 2 * (1 + 1e-12) };
  }
  // Rational approximation p/q (lowest terms) of x via continued fractions, |x − p/q| ≤ tol·max(1,|x|), q ≤ maxDen
  function rationalize(x, tol, maxDen) {
    if (!fin(x)) return null; tol = tol === undefined ? 1e-9 : tol; maxDen = maxDen || 1e6;
    const sg = x < 0 ? -1 : 1, ax = Math.abs(x), lim = tol * Math.max(1, ax);
    let h0 = 0, h1 = 1, k0 = 1, k1 = 0, bb = ax;
    for (let i = 0; i < 80; i++) {
      const ai = Math.floor(bb), h2 = ai * h1 + h0, k2 = ai * k1 + k0; if (k2 > maxDen) break;
      h0 = h1; h1 = h2; k0 = k1; k1 = k2;
      if (Math.abs(ax - h1 / k1) <= lim) return { p: sg * h1, q: k1, err: ax - h1 / k1 };
      const fr = bb - ai; if (fr < 1e-15) break; bb = 1 / fr;
    }
    return null;
  }
  function gcd(a, b) { a = Math.abs(a); b = Math.abs(b); while (b) { const t = a % b; a = b; b = t; } return a; }
  const lcm = (a, b) => a / gcd(a, b) * b;
  // Fundamental period of x[n] = cos(2π f n): f = k/N in lowest terms -> N (f = 0 -> N = 1)
  function periodFromRational(p, q) { if (p === 0) return { N: 1, k: 0 }; const g = gcd(p, q); return { N: q / g, k: Math.abs(p) / g }; }
  // Parse ω in rad/sample: "3pi/8", "3*π/8", "0.75pi", "pi/4", "1.2", "-pi". -> {omega, piCoef:{p,q}|null, hasPi}
  function parseOmega(s) {
    if (typeof s !== 'string') return null;
    const m = /^\s*([+-]?(?:\d+\.?\d*|\.\d+)(?:e[+-]?\d+)?)?\s*\*?\s*(pi|π)?\s*(?:\/\s*((?:\d+\.?\d*|\.\d+)))?\s*$/i.exec(s.replace(/\s+/g, ''));
    if (!m || (m[1] === undefined && !m[2])) return null;
    let coef = m[1] === undefined ? 1 : (m[1] === '+' || m[1] === '-' ? NaN : Number(m[1]));
    if (s.trim()[0] === '-' && m[1] === undefined) coef = -1;
    const den = m[3] === undefined ? 1 : Number(m[3]);
    if (!fin(coef) || !fin(den) || den === 0) return null;
    const hasPi = !!m[2], val = coef / den;
    if (!hasPi) return { omega: val, hasPi: false, piCoef: null };
    const r = rationalize(val, 1e-12, 1e7); // decimal input is exactly rational
    return { omega: val * PI, hasPi: true, piCoef: r ? { p: r.p, q: r.q } : null };
  }
  // Periodicity of cos(ω n + φ): returns {periodic, N, k, f:{p,q}, reason}
  function periodicityOmega(po) {
    if (!po) return null;
    if (!po.hasPi) {
      if (po.omega === 0) return { periodic: true, N: 1, k: 0, f: { p: 0, q: 1 }, reason: 'ω = 0: constant sequence.' };
      return { periodic: false, reason: 'ω = ' + po.omega + ' rad is a non-zero rational number of radians, so ω/2π is irrational: no integer N gives ωN = 2πk.' };
    }
    if (!po.piCoef) return { periodic: false, reason: 'ω/π is not a ratio of integers (within 1e-12).' };
    // f = ω/2π = p/(2q)
    const pp = po.piCoef.p, qq = 2 * po.piCoef.q, g = gcd(pp, qq), r = periodFromRational(pp / g, qq / g);
    return { periodic: true, N: r.N, k: r.k, f: { p: pp / g, q: qq / g }, reason: '' };
  }
  // Tone F at rate Fs: f = F/Fs rationalised with tolerance (inputs are finite decimals, so always "rational")
  function periodicityF(F, Fs, tol, maxDen) {
    const f = F / Fs; if (!fin(f)) return null; const r = rationalize(f, tol === undefined ? 1e-9 : tol, maxDen || 1e5);
    if (!r) return { periodic: false, f, reason: 'f = F/Fs = ' + f + ' has no rational form p/q with q ≤ ' + (maxDen || 1e5) + ' within the tolerance: treated as not periodic.' };
    const g = gcd(r.p, r.q), pr = periodFromRational(r.p / g, r.q / g);
    return { periodic: true, N: pr.N, k: pr.k, f, p: r.p / g, q: r.q / g, err: r.err };
  }

  /* ======================= FIR design ======================= */
  // Oppenheim & Schafer 3e Table 7.2 (windows of length M+1, M = order)
  const WIN = {
    rect: { name: 'Rectangular', sidelobe: -13, mlText: '4π/(M+1)', ml: M => 4 * PI / (M + 1), peakErr: -21 },
    hann: { name: 'Hann', sidelobe: -31, mlText: '8π/M', ml: M => 8 * PI / M, peakErr: -44 },
    hamming: { name: 'Hamming', sidelobe: -41, mlText: '8π/M', ml: M => 8 * PI / M, peakErr: -53 },
    blackman: { name: 'Blackman', sidelobe: -57, mlText: '12π/M', ml: M => 12 * PI / M, peakErr: -74 },
    kaiser: { name: 'Kaiser', sidelobe: null, mlText: '(from β)', ml: null, peakErr: null },
  };
  function besselI0(x) { let s = 1, t = 1; const y = x * x / 4; for (let k = 1; k < 300; k++) { t *= y / (k * k); s += t; if (t < 1e-17 * s) break; } return s; }
  // O&S eq. (7.62): β from A = −20 log10 δ
  function kaiserBeta(A) { if (!fin(A)) return NaN; if (A > 50) return 0.1102 * (A - 8.7); if (A >= 21) return 0.5842 * Math.pow(A - 21, 0.4) + 0.07886 * (A - 21); return 0; }
  // O&S eq. (7.63): M = (A − 8)/(2.285 Δω)  (unrounded)
  const kaiserOrderRaw = (A, dw) => (A - 8) / (2.285 * dw);
  // Symmetric window of length M+1, n = 0..M (O&S eqs. 7.47–7.50, 7.59)
  function windowFn(type, M, beta) {
    const L = M + 1, w = new Float64Array(L);
    for (let n = 0; n < L; n++) {
      const c = M > 0 ? Math.cos(TAU * n / M) : 1;
      if (type === 'hann') w[n] = 0.5 - 0.5 * c;
      else if (type === 'hamming') w[n] = 0.54 - 0.46 * c;
      else if (type === 'blackman') w[n] = 0.42 - 0.5 * c + 0.08 * (M > 0 ? Math.cos(2 * TAU * n / M) : 1);
      else if (type === 'kaiser') { const r = M > 0 ? (n - M / 2) / (M / 2) : 0; w[n] = besselI0(beta * Math.sqrt(Math.max(0, 1 - r * r))) / besselI0(beta); }
      else w[n] = 1;
    }
    return w;
  }
  // Band layout from edges e[] (rad/sample, ascending). Returns {pass:[[lo,hi]], stop:[[lo,hi]], cut:[...], dw}
  const RESP = { lp: 'Lowpass', hp: 'Highpass', bp: 'Bandpass', bs: 'Bandstop' };
  function bands(resp, e) {
    const ok = e.every(fin) && e.every((v, i) => v > 0 && v < PI && (i === 0 || v > e[i - 1]));
    if (!ok) return { err: 'Band edges must satisfy 0 < ω₁ < ω₂' + (resp === 'bp' || resp === 'bs' ? ' < ω₃ < ω₄' : '') + ' < π.' };
    if (resp === 'lp') return { pass: [[0, e[0]]], stop: [[e[1], PI]], cut: [(e[0] + e[1]) / 2], dw: e[1] - e[0] };
    if (resp === 'hp') return { stop: [[0, e[0]]], pass: [[e[1], PI]], cut: [(e[0] + e[1]) / 2], dw: e[1] - e[0] };
    if (resp === 'bp') return { stop: [[0, e[0]], [e[3], PI]], pass: [[e[1], e[2]]], cut: [(e[0] + e[1]) / 2, (e[2] + e[3]) / 2], dw: Math.min(e[1] - e[0], e[3] - e[2]) };
    if (resp === 'bs') return { pass: [[0, e[0]], [e[3], PI]], stop: [[e[1], e[2]]], cut: [(e[0] + e[1]) / 2, (e[2] + e[3]) / 2], dw: Math.min(e[1] - e[0], e[3] - e[2]) };
    return { err: 'unknown response' };
  }
  // Ideal impulse response delayed by α = M/2: lowpass sin(ωc(n−α))/(π(n−α)); "allpass" term sin(π(n−α))/(π(n−α)).
  const sincLP = (wc, t) => (Math.abs(t) < 1e-12 ? wc / PI : Math.sin(wc * t) / (PI * t));
  function idealHd(resp, cut, M) {
    const L = M + 1, al = M / 2, h = new Float64Array(L);
    for (let n = 0; n < L; n++) {
      const t = n - al;
      if (resp === 'lp') h[n] = sincLP(cut[0], t);
      else if (resp === 'hp') h[n] = sincLP(PI, t) - sincLP(cut[0], t);
      else if (resp === 'bp') h[n] = sincLP(cut[1], t) - sincLP(cut[0], t);
      else h[n] = sincLP(PI, t) - sincLP(cut[1], t) + sincLP(cut[0], t);
    }
    return h;
  }
  // Order estimate (O&S Table 7.2 main-lobe widths; Kaiser eq. 7.63)
  function estimateOrder(win, dw, A) {
    if (!(dw > 0)) return NaN;
    if (win === 'kaiser') return Math.max(1, Math.ceil(kaiserOrderRaw(A, dw)));
    if (win === 'rect') return Math.max(1, Math.ceil(4 * PI / dw - 1));
    return Math.max(2, Math.ceil((win === 'blackman' ? 12 : 8) * PI / dw));
  }
  // FIR type by symmetry and length parity (O&S §5.7.3)
  function firType(h, tol) {
    tol = tol === undefined ? 1e-9 : tol; const M = h.length - 1; let mx = 0; h.forEach(v => { mx = Math.max(mx, Math.abs(v)); });
    const t = tol * Math.max(mx, 1e-300); let sym = true, anti = true;
    for (let n = 0; n <= M; n++) { if (Math.abs(h[n] - h[M - n]) > t) sym = false; if (Math.abs(h[n] + h[M - n]) > t) anti = false; }
    if (mx === 0) return { type: 0, text: 'all-zero', M };
    const even = M % 2 === 0;
    if (sym) return even ? { type: 1, M, sym: 'symmetric', zeros: 'none forced', bad: [] } : { type: 2, M, sym: 'symmetric', zeros: 'z = −1 (ω = π): H(e^{jπ}) = 0', bad: ['hp', 'bs'] };
    if (anti) return even ? { type: 3, M, sym: 'antisymmetric', zeros: 'z = 1 and z = −1 (ω = 0 and π)', bad: ['lp', 'hp', 'bs'] } : { type: 4, M, sym: 'antisymmetric', zeros: 'z = 1 (ω = 0): H(e^{j0}) = 0', bad: ['lp', 'bs'] };
    return { type: 0, M, sym: 'none', zeros: '—', bad: [], text: 'not linear phase' };
  }
  // Measure |H| against the bands -> {passDev (linear, max ||H|−1|), stopMax (linear), As (dB), Rp (dB peak-to-peak)}
  function measureFIR(h, bd, npts) {
    npts = npts || 2048; let pd = 0, pmax = 0, pmin = Infinity, smax = 0;
    for (let i = 0; i < npts; i++) {
      const w = PI * i / (npts - 1), m = cabs(freqAt(h, [1], w));
      if (bd.pass.some(p => w >= p[0] - 1e-12 && w <= p[1] + 1e-12)) { pd = Math.max(pd, Math.abs(m - 1)); pmax = Math.max(pmax, m); pmin = Math.min(pmin, m); }
      if (bd.stop.some(p => w >= p[0] - 1e-12 && w <= p[1] + 1e-12)) smax = Math.max(smax, m);
    }
    return { passDev: pd, stopMax: smax, As: -20 * Math.log10(Math.max(smax, 1e-300)), Rp: 20 * Math.log10(pmax / Math.max(pmin, 1e-300)) };
  }
  // Window-method design. spec: {resp, edges[rad], win, d1, As(dB), orderMode:'auto'|'manual', M}
  function designWindow(spec) {
    const bd = bands(spec.resp, spec.edges); if (bd.err) return { err: bd.err };
    if (!(spec.d1 > 0 && spec.d1 < 1) || !(spec.As > 0)) return { err: 'Need 0 < δ₁ < 1 and A_s > 0 dB.' };
    const d2 = Math.pow(10, -spec.As / 20), delta = Math.min(spec.d1, d2), A = -20 * Math.log10(delta);
    const notes = []; let M = spec.orderMode === 'manual' ? Math.round(spec.M) : estimateOrder(spec.win, bd.dw, A);
    const Mest = M;
    if (!(M >= 1) || M > 2000) return { err: 'Order M must be between 1 and 2000 (got ' + M + ').' };
    const needEven = spec.resp === 'hp' || spec.resp === 'bs';
    if (needEven && M % 2 === 1) {
      if (spec.orderMode === 'manual') notes.push('M = ' + M + ' is odd → Type II (forced zero at ω = π), which cannot realise a ' + RESP[spec.resp].toLowerCase() + '. Use an even M.');
      else { M += 1; notes.push('M bumped from ' + Mest + ' to ' + M + ' (even) — a Type II filter has a forced zero at ω = π.'); }
    }
    const beta = spec.win === 'kaiser' ? kaiserBeta(A) : null;
    const hd = idealHd(spec.resp, bd.cut, M), w = windowFn(spec.win, M, beta), h = new Float64Array(M + 1);
    for (let n = 0; n <= M; n++) h[n] = hd[n] * w[n];
    if (spec.win !== 'kaiser' && WIN[spec.win].peakErr > -spec.As) notes.push('The ' + WIN[spec.win].name + ' window\'s peak approximation error ≈ ' + WIN[spec.win].peakErr + ' dB cannot reach A_s = ' + spec.As + ' dB at any M (use a stronger window or Kaiser).');
    const meas = measureFIR(h, bd), type = firType(h);
    return { h, hd, w, M, Mest, beta, A, delta, d2, bands: bd, meas, type, notes };
  }
  // Frequency-sampling design (Type I/II linear phase): L = M+1 samples at ω_k = 2πk/L, amplitudes A_k ∈ {1, T, 0},
  // h[n] = (1/L)[A_0 + 2 Σ_{k=1}^{⌊(L−1)/2⌋} A_k cos(2πk(n − M/2)/L)]  (the k = L/2 sample is forced to 0 for L even)
  function designFreqSamp(spec) {
    const bd = bands(spec.resp, spec.edges); if (bd.err) return { err: bd.err };
    const M = Math.round(spec.M); if (!(M >= 2) || M > 2000) return { err: 'Order M must be between 2 and 2000.' };
    const L = M + 1, al = M / 2, T = spec.T === undefined ? 0.4 : spec.T, Ak = [];
    for (let k = 0; k <= Math.floor((L - 1) / 2); k++) {
      const w = TAU * k / L;
      const inP = bd.pass.some(p => w >= p[0] - 1e-12 && w <= p[1] + 1e-12), inS = bd.stop.some(p => w >= p[0] - 1e-12 && w <= p[1] + 1e-12);
      Ak.push({ k, w, A: inP ? 1 : inS ? 0 : T, where: inP ? 'pass' : inS ? 'stop' : 'trans' });
    }
    const h = new Float64Array(L);
    for (let n = 0; n < L; n++) { let s = Ak[0].A; for (let i = 1; i < Ak.length; i++) s += 2 * Ak[i].A * Math.cos(TAU * Ak[i].k * (n - al) / L); h[n] = s / L; }
    const notes = [];
    if (L % 2 === 0 && (spec.resp === 'hp' || spec.resp === 'bs')) notes.push('L = M+1 even → Type II: H(e^{jπ}) is forced to 0, so the passband at π cannot be met. Use an even M.');
    return { h, Ak, M, L, bands: bd, meas: measureFIR(h, bd), type: firType(h), notes };
  }

  /* ======================= filter structures ======================= */
  const nz = c => c.filter(v => v !== 0).length;
  const nontrivial = c => c.filter(v => v !== 0 && v !== 1 && v !== -1).length;
  // Counts for a section with numerator b and denominator a (a[0] = 1)
  function countsBA(b, a, kind) {
    const fb = nontrivial(b), fa = nontrivial(a.slice(1)), terms = nz(b) + nz(a.slice(1));
    const M = b.length - 1, N = a.length - 1;
    return { mult: fb + fa, add: Math.max(0, terms - 1), delay: kind === 'df1' ? M + N : Math.max(M, N) };
  }
  function simDF1(b, a, x) {
    const xs = new Float64Array(b.length), ys = new Float64Array(a.length), y = new Float64Array(x.length);
    for (let n = 0; n < x.length; n++) {
      for (let k = b.length - 1; k > 0; k--) xs[k] = xs[k - 1]; xs[0] = x[n];
      let v = 0; for (let k = 0; k < b.length; k++) v += b[k] * xs[k];     // feed-forward adder chain
      let s = v; for (let k = 1; k < a.length; k++) s -= a[k] * ys[k - 1]; // feedback adder chain
      for (let k = a.length - 1; k > 0; k--) ys[k] = ys[k - 1]; ys[0] = s; y[n] = s;
    }
    return y;
  }
  function simDF2(b, a, x) {
    const L = Math.max(b.length, a.length), w = new Float64Array(L), y = new Float64Array(x.length);
    for (let n = 0; n < x.length; n++) {
      let wn = x[n]; for (let k = 1; k < a.length; k++) wn -= a[k] * w[k]; w[0] = wn;
      let s = 0; for (let k = 0; k < b.length; k++) s += b[k] * w[k]; y[n] = s;
      for (let k = L - 1; k > 0; k--) w[k] = w[k - 1];
    }
    return y;
  }
  function simTDF2(b, a, x) {
    const L = Math.max(b.length, a.length) - 1, s = new Float64Array(L + 1), y = new Float64Array(x.length), bb = k => (k < b.length ? b[k] : 0), aa = k => (k < a.length ? a[k] : 0);
    for (let n = 0; n < x.length; n++) {
      const yn = bb(0) * x[n] + (L > 0 ? s[0] : 0);
      for (let i = 0; i < L; i++) s[i] = s[i + 1] + bb(i + 1) * x[n] - aa(i + 1) * yn;
      y[n] = yn;
    }
    return y;
  }
  // Group roots into real second-order factors: conjugate pairs together, real roots paired (sorted by radius)
  function groupRoots(roots) {
    const cpx = roots.filter(r => r.im > 0), real = roots.filter(r => r.im === 0).sort((p, q) => cabs(p) - cabs(q)), grp = [];
    cpx.forEach(r => grp.push([r, cconj(r)]));
    for (let i = 0; i < real.length; i += 2) grp.push(real.slice(i, i + 2));
    return grp;
  }
  // real coefficients of Π(1 − r_i z^−1) over a group: [1, −Σr, Π r]
  function grpPoly(g) { let p = [C(1, 0)]; g.forEach(r => { p = polyMulC(p, [C(1, 0), C(-r.re, -r.im)]); }); return p.map(c => c.re); }
  // Cascade of second-order sections. Gain K absorbed into section 1's numerator; a pure delay z^−D kept separate.
  function cascadeSOS(b, a) {
    const P = pzk(b, a), D = P.delay;
    const zr = []; P.zeros.forEach(o => { if (cabs(o.r) > 0) for (let k = 0; k < o.m; k++) zr.push(o.r); });
    const pr = []; P.poles.forEach(o => { if (cabs(o.r) > 0) for (let k = 0; k < o.m; k++) pr.push(o.r); });
    const pg = groupRoots(pr).sort((g1, g2) => Math.max(...g1.map(cabs)) - Math.max(...g2.map(cabs))), zg = groupRoots(zr);
    const n = Math.max(pg.length, zg.length, 1), secs = [];
    for (let i = 0; i < n; i++) {
      const pgi = pg[i] || [];
      let bestJ = -1, bestD = Infinity;
      zg.forEach((g, j) => { if (!g) return; const d = pgi.length ? Math.min(...g.map(z => Math.min(...pgi.map(p => cabs(csub(z, p)))))) : 0; if (d < bestD) { bestD = d; bestJ = j; } });
      const zgi = bestJ >= 0 ? zg[bestJ] : []; if (bestJ >= 0) zg[bestJ] = null;
      secs.push({ b: grpPoly(zgi), a: grpPoly(pgi), zeros: zgi, poles: pgi });
    }
    // b-side polynomial is Π(1 − z_i w) times b_D; zeros in w are only the non-zero z-plane zeros; leftover in z-plane
    // zeros/poles at the origin are accounted for by the delay D and the length difference (automatically, in w).
    secs[0].b = secs[0].b.map(v => v * P.K);
    return { secs, K: P.K, delay: D };
  }
  function simCascade(cas, x) {
    let y = Float64Array.from(x);
    if (cas.delay) { const z = new Float64Array(y.length); for (let n = cas.delay; n < y.length; n++) z[n] = y[n - cas.delay]; y = z; }
    cas.secs.forEach(s => { y = simDF2(s.b, s.a, y); });
    return y;
  }
  // Parallel form from partial fractions: direct FIR part + real sections (conjugate pole clusters combined)
  function parallelForm(b, a) {
    const pf = partialFractions(b, a), secs = [];
    pf.clusters.forEach((cl, i) => {
      if (cl.r.im < 0) return; // handled with its conjugate
      const ts = pf.terms.filter(t => t.cluster === i), s = cl.m, d = cl.r;
      // N(w) = Σ_m C_m (1 − d w)^{s−m}, D(w) = (1 − d w)^s (ascending in w)
      let Nw = [C(0, 0)]; const base = [C(1, 0), C(-d.re, -d.im)];
      ts.forEach(t => { let p = [C(1, 0)]; for (let k = 0; k < s - t.order; k++) p = polyMulC(p, base); p = p.map(c => cmul(c, t.C)); while (Nw.length < p.length) Nw.push(C(0, 0)); p.forEach((c, k) => { Nw[k] = cadd(Nw[k], c); }); });
      let Dw = [C(1, 0)]; for (let k = 0; k < s; k++) Dw = polyMulC(Dw, base);
      if (cl.r.im > 0) { // N/D + conj(N)/conj(D) = 2 Re(N·conj(D)) / (D·conj(D))
        const num = polyMulC(Nw, Dw.map(cconj)).map(c => 2 * c.re), den = polyMulC(Dw, Dw.map(cconj)).map(c => c.re);
        secs.push({ b: trimEnd(num), a: den, poles: [d, cconj(d)], mult: s });
      } else secs.push({ b: trimEnd(Nw.map(c => c.re)), a: Dw.map(c => c.re), poles: [d], mult: s });
    });
    return { direct: pf.direct, secs, pf };
  }
  function trimEnd(arr) { const o = arr.slice(); while (o.length > 1 && Math.abs(o[o.length - 1]) < 1e-15) o.pop(); return o; }
  function simParallel(par, x) {
    const y = new Float64Array(x.length);
    par.direct.forEach((c, r) => { for (let n = r; n < x.length; n++) y[n] += c * x[n - r]; });
    par.secs.forEach(s => { const ys = simDF2(s.b, s.a, x); for (let n = 0; n < x.length; n++) y[n] += ys[n]; });
    return y;
  }
  // Linear-phase folded FIR: y[n] = Σ_{k<M/2} h[k](x[n−k] ± x[n−M+k]) (+ h[M/2] x[n−M/2] for even M, symmetric)
  function simFolded(h, x) {
    const ty = firType(h); if (!(ty.type >= 1)) return null;
    const M = h.length - 1, sg = ty.type <= 2 ? 1 : -1, buf = new Float64Array(M + 1), y = new Float64Array(x.length);
    for (let n = 0; n < x.length; n++) {
      for (let k = M; k > 0; k--) buf[k] = buf[k - 1]; buf[0] = x[n];
      let s = 0; for (let k = 0; k < Math.floor((M + 1) / 2); k++) s += h[k] * (buf[k] + sg * buf[M - k]);
      if (M % 2 === 0 && sg > 0) s += h[M / 2] * buf[M / 2];
      y[n] = s;
    }
    return y;
  }
  function foldedCounts(h) {
    const M = h.length - 1, ty = firType(h), half = Math.floor((M + 1) / 2), mid = M % 2 === 0 && ty.type === 1;
    const coefs = h.slice(0, half).concat(mid ? [h[M / 2]] : []);
    return { mult: nontrivial(Array.from(coefs)), add: half + (nz(Array.from(coefs)) - 1), delay: M };
  }
  const maxAbsDiff = (p, q) => { let m = 0; for (let i = 0; i < p.length; i++) m = Math.max(m, Math.abs(p[i] - q[i])); return m; };

  /* ======================= DFT properties ======================= */
  function dft(x, N) { // N-point DFT of x (zero-padded or truncated to N), direct sum, X[k] = Σ x[n] W_N^{kn}, W_N = e^{−j2π/N}
    const re = new Float64Array(N), im = new Float64Array(N);
    for (let k = 0; k < N; k++) { let sr = 0, si = 0; for (let n = 0; n < N && n < x.length; n++) { const a = -TAU * ((k * n) % N) / N; sr += x[n] * Math.cos(a); si += x[n] * Math.sin(a); } re[k] = sr; im[k] = si; }
    return { re, im };
  }
  function idft(re, im) { const N = re.length, xr = new Float64Array(N), xi = new Float64Array(N); for (let n = 0; n < N; n++) { let sr = 0, si = 0; for (let k = 0; k < N; k++) { const a = TAU * ((k * n) % N) / N, c = Math.cos(a), s = Math.sin(a); sr += re[k] * c - im[k] * s; si += re[k] * s + im[k] * c; } xr[n] = sr / N; xi[n] = si / N; } return { re: xr, im: xi }; }
  const padTo = (x, N) => { const o = new Float64Array(N); for (let i = 0; i < Math.min(N, x.length); i++) o[i] = x[i]; return o; };
  const modN = (n, N) => ((n % N) + N) % N;
  function circShift(x, m, N) { const xp = padTo(x, N), o = new Float64Array(N); for (let n = 0; n < N; n++) o[n] = xp[modN(n - m, N)]; return o; }
  function timeRev(x, N) { const xp = padTo(x, N), o = new Float64Array(N); for (let n = 0; n < N; n++) o[n] = xp[modN(-n, N)]; return o; }
  function circConv(x, h, N) { // inputs longer than N are time-aliased modulo N first
    const A = new Float64Array(N), B = new Float64Array(N), y = new Float64Array(N);
    x.forEach((v, i) => { A[i % N] += v; }); h.forEach((v, i) => { B[i % N] += v; });
    for (let n = 0; n < N; n++) { let s = 0; for (let k = 0; k < N; k++) s += A[k] * B[modN(n - k, N)]; y[n] = s; }
    return y;
  }
  function linConv(x, h) { const y = new Float64Array(x.length + h.length - 1); for (let i = 0; i < x.length; i++) for (let j = 0; j < h.length; j++) y[i + j] += x[i] * h[j]; return y; }
  function dtft(x, w) { let r = 0, i = 0; for (let n = 0; n < x.length; n++) { r += x[n] * Math.cos(w * n); i -= x[n] * Math.sin(w * n); } return C(r, i); }
  // B-bit mid-rise uniform quantiser on [−1, 1]: Δ = 2/2^B, q(x) = Δ(⌊x/Δ⌋ + ½) clipped to ±(1 − Δ/2)
  function quantize(x, B) { const D = 2 / Math.pow(2, B), lim = 1 - D / 2; return Math.max(-lim, Math.min(lim, D * (Math.floor(x / D) + 0.5))); }
  // Simulated SNR of a full-scale (A = 1) sinusoid at a frequency incommensurate with the sample rate
  function quantSNR(B, nS, amp) {
    nS = nS || 65536; amp = amp === undefined ? 1 : amp; const f0 = 0.0123456789 * (1 + Math.sqrt(5)) / 2; let ps = 0, pe = 0;
    for (let n = 0; n < nS; n++) { const x = amp * Math.sin(TAU * f0 * n + 0.3), e = quantize(x, B) - x; ps += x * x; pe += e * e; }
    return 10 * Math.log10(ps / pe);
  }
  const snrTheory = B => 6.02 * B + 1.76;

  FSP.math.dspdesign = {
    cx, hornerD, polyMulC, polyFromRoots, rootsMult, polyRoots, parseList, normBA, pzk, filter, impulse, stepResp, freqAt, freqResp, stability,
    polyDivAsc, taylorAt, partialFractions, rocOptions, binomPoly, inverseZAt, diffEqResidual,
    alias, rationalize, gcd, lcm, periodFromRational, parseOmega, periodicityOmega, periodicityF,
    WIN, RESP, besselI0, kaiserBeta, kaiserOrderRaw, windowFn, bands, idealHd, estimateOrder, firType, measureFIR, designWindow, designFreqSamp,
    countsBA, simDF1, simDF2, simTDF2, groupRoots, cascadeSOS, simCascade, parallelForm, simParallel, simFolded, foldedCounts, maxAbsDiff,
    dft, idft, padTo, circShift, timeRev, circConv, linConv, dtft, quantize, quantSNR, snrTheory,
  };

  /* ======================= tests ======================= */

  FSP.math.dspdesign._shared = { C, PI, TAU, cabs, carg, cmul, csub, fin };
})();
