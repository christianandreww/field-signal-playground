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
  FSP.registerTests('dspdesign', t => {
    const D = FSP.math.dspdesign;
    const hasRoot = (rs, re, im, tol) => rs.some(o => Math.abs(o.r.re - re) < tol && Math.abs(o.r.im - im) < tol);
    // roots: z^2 + 1 -> ±j
    { const r = D.rootsMult([1, 0, 1]); t.check('roots of z²+1 are ±j', r.length === 2 && hasRoot(r, 0, 1, 1e-14) && hasRoot(r, 0, -1, 1e-14)); }
    // (z−1)²(z+2) = z³ − 3z + 2: double root 1, simple −2
    { const r = D.rootsMult([1, 0, -3, 2]), d = r.find(o => o.m === 2);
      t.check('z³−3z+2: double root at 1 (multiplicity detected), simple root −2', r.length === 2 && d && Math.abs(d.r.re - 1) < 1e-12 && Math.abs(d.r.im) === 0 && hasRoot(r, -2, 0, 1e-13), JSON.stringify(r.map(o => [o.r.re, o.r.im, o.m]))); }
    // z^10 − 1: ten roots of unity
    { const r = D.polyRoots([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1]); let ok = r.length === 10;
      for (let k = 0; k < 10; k++) { const th = TAU * k / 10; if (!r.some(z => Math.hypot(z.re - Math.cos(th), z.im - Math.sin(th)) < 1e-13)) ok = false; }
      t.check('z¹⁰ − 1: all ten 10th roots of unity to 1e-13', ok); }
    // (z − 0.5)^2 (z² + 0.81): coefficients by hand: (z² − z + 0.25)(z² + 0.81) = z⁴ − z³ + 1.06 z² − 0.81 z + 0.2025
    { const r = D.rootsMult([1, -1, 1.06, -0.81, 0.2025]);
      t.check('(z−0.5)²(z²+0.81): double 0.5 and ±0.9j', r.length === 3 && r.some(o => o.m === 2 && Math.abs(o.r.re - 0.5) < 1e-12) && hasRoot(r, 0, 0.9, 1e-12) && hasRoot(r, 0, -0.9, 1e-12)); }
    // complex-coefficient: (z − (1+2j))(z − (−3+0.5j)) = z² − (−2+2.5j) z + (1+2j)(−3+0.5j) = z² + (2 − 2.5j) z + (−3 − 1 + (0.5 − 6)j) = z² + (2−2.5j)z + (−4 − 5.5j)
    { const r = D.rootsMult([C(1, 0), C(2, -2.5), C(-4, -5.5)]); t.check('complex-coefficient quadratic: roots 1+2j, −3+0.5j', hasRoot(r, 1, 2, 1e-12) && hasRoot(r, -3, 0.5, 1e-12)); }
    // order-10 Wilkinson-type with known roots 0.1..1.0 rebuilt from roots
    { const want = []; for (let k = 1; k <= 8; k++) want.push(C(0.95 * Math.cos(k * 0.7), 0.95 * Math.sin(k * 0.7))); want.push(C(-0.3, 0), C(0.6, 0));
      const p = D.polyFromRoots(want), r = D.polyRoots(p); const ok = want.every(w => r.some(z => cabs(csub(z, w)) < 1e-10));
      t.check('order-10 complex polynomial rebuilt from known roots: all recovered to 1e-10', ok && r.length === 10); }
    // triple root (z+0.5)^3 = z³ + 1.5 z² + 0.75 z + 0.125
    { const r = D.rootsMult([1, 1.5, 0.75, 0.125]); t.check('(z+0.5)³: one root −0.5 of multiplicity 3', r.length === 1 && r[0].m === 3 && Math.abs(r[0].r.re + 0.5) < 1e-10); }
    // impulse response of y[n] = 0.5 y[n−1] + x[n]  (a1 = −0.5)
    { const h = D.impulse([1], [1, -0.5], 30); let e = 0; for (let n = 0; n < 30; n++) e = Math.max(e, Math.abs(h[n] - Math.pow(0.5, n)));
      t.check('h[n] of y[n] = 0.5y[n−1] + x[n] is 0.5^n', e < 1e-15, e.toExponential(2)); }
    // step response of same: Σ 0.5^k = 2(1 − 0.5^{n+1})
    { const s = D.stepResp([1], [1, -0.5], 20); t.check('step response = 2(1 − 0.5^(n+1))', s.every((v, n) => Math.abs(v - 2 * (1 - Math.pow(0.5, n + 1))) < 1e-14)); }
    // partial fractions: 1/((1−0.5z^-1)(1−0.25z^-1)): a = [1, −0.75, 0.125]; A1 = 1/(1 − 0.25/0.5) = 2, A2 = 1/(1 − 0.5/0.25) = −1
    { const pf = D.partialFractions([1], [1, -0.75, 0.125]);
      const c5 = pf.terms.find(x => Math.abs(x.d.re - 0.5) < 1e-12), c25 = pf.terms.find(x => Math.abs(x.d.re - 0.25) < 1e-12);
      t.check('PF of 1/((1−0.5z⁻¹)(1−0.25z⁻¹)): residues 2 and −1', c5 && c25 && Math.abs(c5.C.re - 2) < 1e-12 && Math.abs(c25.C.re + 1) < 1e-12);
      let e = 0; for (let n = 0; n < 25; n++) e = Math.max(e, Math.abs(D.inverseZAt(pf, 0.5, n) - (2 * Math.pow(0.5, n) - Math.pow(0.25, n))));
      t.check('causal inverse = 2(0.5)ⁿ − (0.25)ⁿ', e < 1e-14, e.toExponential(2));
      // anticausal (|z| < 0.25): h[n] = −2(0.5)^n + (0.25)^n for n ≤ −1 ; at n = −1: −4 + 4 = 0; n = −2: −8 + 16 = 8
      t.check('anticausal ROC: h[−1] = 0, h[−2] = 8, h[0] = 0', Math.abs(D.inverseZAt(pf, 0, -1)) < 1e-12 && Math.abs(D.inverseZAt(pf, 0, -2) - 8) < 1e-12 && D.inverseZAt(pf, 0, 0) === 0);
      // two-sided 0.25<|z|<0.5: −2(0.5)^n u[−n−1] − (0.25)^n u[n]; check it satisfies the difference equation everywhere
      const res = D.diffEqResidual([1], [1, -0.75, 0.125], n => D.inverseZAt(pf, 0.25, n), -15, 15);
      t.check('two-sided ROC closed form satisfies the difference equation (residual < 1e-9)', res < 1e-9, res.toExponential(2)); }
    // double pole: 1/(1 − 0.5 z^-1)² = Σ (n+1) 0.5^n
    { const pf = D.partialFractions([1], [1, -1, 0.25]); let e = 0; for (let n = 0; n < 30; n++) e = Math.max(e, Math.abs(D.inverseZAt(pf, 0.5, n) - (n + 1) * Math.pow(0.5, n)));
      t.check('double pole 1/(1−0.5z⁻¹)²: h[n] = (n+1)(0.5)ⁿ', e < 1e-12, e.toExponential(2)); }
    // improper: (1 + 2z^-1 + z^-2)/(1 − 0.5z^-1): long division in w: Q = −2w − 8? check via recursion instead (independent algorithm)
    { const b = [1, 2, 1, 0.5], a = [1, -0.5, 0.06], pf = D.partialFractions(b, a), h = D.impulse(b, a, 40); let e = 0; for (let n = 0; n < 40; n++) e = Math.max(e, Math.abs(D.inverseZAt(pf, 0.5, n) - h[n]));
      t.check('improper H (M ≥ N): direct terms + PF matches the recursion', e < 1e-12 && pf.direct.length === 2, e.toExponential(2)); }
    // Kaiser β: A = 60 -> 0.1102·51.3 = 5.65326; A = 40 -> 0.5842·19^0.4 + 0.07886·19 = 1.89692 + 1.49834 = 3.39526; A = 20 -> 0
    t.check('Kaiser β(A=60) = 5.6533', t.near(D.kaiserBeta(60), 5.65326, 1e-5));
    t.check('Kaiser β(A=40) = 3.3953', t.near(D.kaiserBeta(40), 3.39526, 1e-4));
    t.check('Kaiser β(A=20) = 0', D.kaiserBeta(20) === 0);
    // O&S Example: ωp = 0.4π, ωs = 0.6π, δ = 0.001 -> A = 60, M = 52/(2.285·0.2π) = 36.22 -> 37
    t.check('Kaiser order (O&S example): M = 37', D.estimateOrder('kaiser', 0.2 * PI, 60) === 37);
    // Hamming-windowed lowpass: symmetric, DC gain ≈ 1
    { const r = D.designWindow({ resp: 'lp', edges: [0.3 * PI, 0.5 * PI], win: 'hamming', d1: 0.01, As: 50, orderMode: 'auto' });
      const dc = r.h.reduce((s, v) => s + v, 0); let sym = true; for (let n = 0; n <= r.M; n++) if (Math.abs(r.h[n] - r.h[r.M - n]) > 1e-15) sym = false;
      // Hamming main lobe 8π/M ≤ 0.2π -> M ≥ 40
      t.check('Hamming LP: M = 8π/Δω = 40, symmetric, DC gain ≈ 1 (±0.005)', r.M === 40 && sym && Math.abs(dc - 1) < 0.005, 'M=' + r.M + ' dc=' + dc.toFixed(5));
      t.check('Hamming LP meets ≈ −53 dB class stopband (A_s > 50 dB)', r.meas.As > 50, r.meas.As.toFixed(1)); }
    // highpass auto order is forced even (Type I)
    { const r = D.designWindow({ resp: 'hp', edges: [0.4 * PI, 0.6 * PI], win: 'kaiser', d1: 0.001, As: 60, orderMode: 'auto' });
      t.check('Kaiser HP: M bumped 37 -> 38 (Type I)', r.M === 38 && r.type.type === 1); }
    // FIR type classification
    t.check('type I: [1,2,1]', D.firType([1, 2, 1]).type === 1);
    t.check('type II: [1,1] (zero at ω=π)', D.firType([1, 1]).type === 2 && Math.abs(cabs(D.freqAt([1, 1], [1], PI))) < 1e-15);
    t.check('type III: [1,0,−1] (zeros at 0 and π)', D.firType([1, 0, -1]).type === 3);
    t.check('type IV: [1,−1] (zero at ω=0)', D.firType([1, -1]).type === 4);
    t.check('not linear phase: [1,2,3]', D.firType([1, 2, 3]).type === 0);
    // frequency sampling hits its samples exactly: |H(e^{jω_k})| = A_k
    { const r = D.designFreqSamp({ resp: 'lp', edges: [0.3 * PI, 0.4 * PI], M: 32, T: 0.4 }); let e = 0;
      r.Ak.forEach(s => { e = Math.max(e, Math.abs(cabs(D.freqAt(r.h, [1], s.w)) - s.A)); });
      t.check('frequency sampling: |H(e^{jω_k})| = A_k at all samples', e < 1e-12 && r.type.type === 1, e.toExponential(2)); }
    // structures vs recursion: 4th order with poles 0.8e^{±j0.5}, 0.6e^{±j1.8}, b = [1, 2, 3, 2, 1] and a 5th-order with real poles
    { const mkA = ps => D.polyFromRoots(ps).map(c => c.re);
      const cases = [
        [[1, 2, 3, 2, 1], mkA([C(0.8 * Math.cos(0.5), 0.8 * Math.sin(0.5)), C(0.8 * Math.cos(0.5), -0.8 * Math.sin(0.5)), C(0.6 * Math.cos(1.8), 0.6 * Math.sin(1.8)), C(0.6 * Math.cos(1.8), -0.6 * Math.sin(1.8))])],
        [[0.2, -0.1, 0.4, 0.3, 0, 0.1, 0.05], mkA([C(0.9, 0), C(0.8 * Math.cos(0.7), 0.8 * Math.sin(0.7)), C(0.8 * Math.cos(0.7), -0.8 * Math.sin(0.7)), C(-0.5, 0), C(0.3, 0)])],
        [[0, 0, 1, 0.5], [1, -0.9]]];
      let worst = 0, okAll = true;
      cases.forEach(([b, a]) => {
        const x = new Float64Array(80); x[0] = 1; const ref = D.filter(b, a, x), q = D.normBA(b, a);
        const ys = [D.simDF1(q.b, q.a, x), D.simDF2(q.b, q.a, x), D.simTDF2(q.b, q.a, x), D.simCascade(D.cascadeSOS(b, a), x), D.simParallel(D.parallelForm(b, a), x)];
        ys.forEach(y => { const e = D.maxAbsDiff(y, ref); worst = Math.max(worst, e); if (!(e < 1e-12)) okAll = false; });
      });
      t.check('DF-I, DF-II, TDF-II, cascade SOS, parallel all match the recursion to 1e-12', okAll, worst.toExponential(2)); }
    { const h = [0.1, 0.2, 0.4, 0.2, 0.1], x = new Float64Array(20); x[0] = 1; x[3] = -2; x[7] = 0.5;
      t.check('folded linear-phase FIR equals direct convolution', D.maxAbsDiff(D.simFolded(h, x), D.filter(h, [1], x)) < 1e-15);
      t.check('folded FIR (M=4, type I): 3 multipliers, 4 delays', D.foldedCounts(h).mult === 3 && D.foldedCounts(h).delay === 4); }
    { const cas = D.cascadeSOS([1, 2, 3, 2, 1], [1, -1.2, 0.9, -0.4, 0.1]);
      t.check('cascade: 4th order -> 2 SOS, every a-section has a₀ = 1', cas.secs.length === 2 && cas.secs.every(s => s.a[0] === 1)); }
    // DFT of [1,2,3,4] = [10, −2+2j, −2, −2−2j]
    { const X = D.dft([1, 2, 3, 4], 4), want = [[10, 0], [-2, 2], [-2, 0], [-2, -2]];
      t.check('DFT{1,2,3,4} = [10, −2+2j, −2, −2−2j]', want.every((w, k) => Math.abs(X.re[k] - w[0]) < 1e-12 && Math.abs(X.im[k] - w[1]) < 1e-12)); }
    // circular conv of [1,2,3,4] with [1,1]: y[n] = x[n] + x[(n−1) mod 4] = [5,3,5,7]; N = 5 equals linear [1,3,5,7,4]
    { const y4 = D.circConv([1, 2, 3, 4], [1, 1], 4), y5 = D.circConv([1, 2, 3, 4], [1, 1], 5);
      t.check('4-point circular conv [1,2,3,4]⊛[1,1] = [5,3,5,7]', [5, 3, 5, 7].every((v, i) => Math.abs(y4[i] - v) < 1e-12));
      t.check('N = 5 ≥ L+M−1: circular = linear [1,3,5,7,4]', [1, 3, 5, 7, 4].every((v, i) => Math.abs(y5[i] - v) < 1e-12)); }
    { const x = [1, 2, 3, 4], X = D.dft(x, 4), Xs = D.dft(D.circShift(x, 1, 4), 4); let e = 0;
      for (let k = 0; k < 4; k++) { const w = cmul(C(X.re[k], X.im[k]), C(Math.cos(-TAU * k / 4), Math.sin(-TAU * k / 4))); e = Math.max(e, Math.abs(w.re - Xs.re[k]), Math.abs(w.im - Xs.im[k])); }
      t.check('circular shift by 1 ↔ X[k]·W₄^k', e < 1e-12);
      const r = D.timeRev(x, 4); t.check('time reversal mod 4: [1,4,3,2]', [1, 4, 3, 2].every((v, i) => r[i] === v)); }
    // aliasing 7 kHz at 10 kHz -> 3 kHz (phase flips)
    { const a = D.alias(7000, 10000); t.check('7 kHz sampled at 10 kHz aliases to 3 kHz (phase negated)', a.Fa === 3000 && a.sign === -1 && a.aliased); }
    t.check('12 kHz at Fs = 10 kHz aliases to 2 kHz', D.alias(12000, 10000).Fa === 2000 && D.alias(12000, 10000).sign === 1);
    // periodicity: ω = 3π/8 -> f = 3/16 -> N = 16
    { const p = D.periodicityOmega(D.parseOmega('3pi/8')); t.check('ω = 3π/8 is periodic with N = 16 (k = 3)', p.periodic && p.N === 16 && p.k === 3); }
    t.check('ω = 1 rad is not periodic', D.periodicityOmega(D.parseOmega('1')).periodic === false);
    t.check('ω = π/4 → N = 8; ω = 0.75π → N = 8', D.periodicityOmega(D.parseOmega('π/4')).N === 8 && D.periodicityOmega(D.parseOmega('0.75pi')).N === 8);
    t.check('F = 3 kHz, Fs = 8 kHz -> f = 3/8, N = 8', D.periodicityF(3000, 8000).N === 8);
    t.check('parseOmega rejects junk', D.parseOmega('abc') === null && D.parseOmega('3pi/0') === null);
    // quantisation SNR ≈ 6.02B + 1.76 dB
    t.check('8-bit quantiser SNR ≈ 49.9 dB (±0.3)', t.near(D.quantSNR(8), 6.02 * 8 + 1.76, 0.3), D.quantSNR(8).toFixed(2));
    t.check('12-bit quantiser SNR ≈ 74.0 dB (±0.3)', t.near(D.quantSNR(12), 6.02 * 12 + 1.76, 0.3), D.quantSNR(12).toFixed(2));
    // invalid input
    t.check('normBA rejects a₀ = 0 and empty lists', !!D.normBA([1], [0, 1]).err && !!D.normBA([], [1]).err);
    t.check('parseList: fractions ok, junk -> null', JSON.stringify(D.parseList('1, -1/2 0.25')) === '[1,-0.5,0.25]' && D.parseList('1, x') === null);
    t.check('designWindow rejects unordered edges', !!D.designWindow({ resp: 'lp', edges: [0.5 * PI, 0.3 * PI], win: 'hann', d1: 0.01, As: 40, orderMode: 'auto' }).err);
    t.check('stability: pole at 1.1 unstable, at 0.9 stable', D.stability([{ r: C(1.1, 0), m: 1 }]).verdict === 'unstable' && D.stability([{ r: C(0.9, 0), m: 1 }]).verdict === 'stable');
  });

  /* ======================= UI ======================= */
  const DM = FSP.math.dspdesign;
  const css = n => { try { return getComputedStyle(document.documentElement).getPropertyValue(n).trim() || null; } catch (e) { return null; } };
  function theme() {
    return { bg: css('--panel2') || '#0f141f', text: css('--text') || '#e6edf3', muted: css('--muted') || '#8b98ab', border: css('--border') || '#243047',
      a: css('--accent2') || '#f9a03f', b: css('--ok') || '#2ecc71', c: css('--accent') || '#4cc9f0', pink: css('--pink') || '#f72585',
      warn: css('--warn') || '#ffb347', bad: css('--bad') || '#ff5c5c', mono: css('--mono') || 'monospace' };
  }
  const clean = x => (Math.abs(x) < 1e-12 ? 0 : x);
  const g4 = x => (fin(x) ? String(+clean(x).toPrecision(4)) : '—');
  const g3 = x => (fin(x) ? String(+clean(x).toPrecision(3)) : '—');
  const fmt = (x, d) => (fin(x) ? x.toFixed(d === undefined ? 3 : d) : '—');
  const cstr = c => { const re = clean(c.re), im = clean(c.im); return im === 0 ? g4(re) : g4(re) + (im < 0 ? ' − j' : ' + j') + g4(Math.abs(im)); };
  const SUPD = { 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' };
  const supn = k => String(k).split('').map(d => SUPD[d]).join('');
  const zk = k => (k === 0 ? '' : 'z⁻' + supn(k));
  // Σ c_k z^-k as text
  function polyStr(c) {
    let s = '';
    c.forEach((x, k) => {
      if (x === 0) return; const ax = Math.abs(x), body = k === 0 ? g4(ax) : (ax === 1 ? '' : g4(ax) + ' ') + zk(k);
      s += s === '' ? (x < 0 ? '−' : '') + body : (x < 0 ? ' − ' : ' + ') + body;
    });
    return s === '' ? '0' : s;
  }
  // linear combination text: terms [[coef, 'name'], …]
  function lin(terms) {
    let s = '';
    terms.forEach(t => { const x = t[0]; if (x === 0) return; const ax = Math.abs(x), body = (ax === 1 ? '' : g4(ax) + ' ') + t[1]; s += s === '' ? (x < 0 ? '−' : '') + body : (x < 0 ? ' − ' : ' + ') + body; });
    return s === '' ? '0' : s;
  }
  const idx = k => (k === 0 ? 'n' : 'n−' + k);
  const tickStr = v => (Math.abs(v) < 1e-12 ? '0' : String(+v.toPrecision(4)));
  function niceTicks(lo, hi, n) {
    const span = hi - lo; if (!(span > 0)) return [lo];
    const raw = span / n, mag = Math.pow(10, Math.floor(Math.log10(raw))), f = raw / mag, step = (f < 1.5 ? 1 : f < 3 ? 2 : f < 7 ? 5 : 10) * mag, out = [];
    for (let v = Math.ceil(lo / step - 1e-9) * step; v <= hi + step * 1e-9; v += step) out.push(Math.abs(v) < step * 1e-9 ? 0 : +v.toPrecision(12));
    return out;
  }
  /* generic plot: o {xmin,xmax,ymin,ymax,xlabel,ylabel, series:[{x,y,color,width,dash,stem,base}], rects:[{x0,x1,y0,y1,color}],
     hlines:[{y,x0,x1,color,dash}], vlines:[{x,color}], points:[{x,y,label,color,below}], texts:[{x,y,text,color,align}], legend:[{text,color}]} */
  function plot(c, o) {
    const g = c.prep(), ctx = g.ctx, w = g.w, h = g.h, T = theme();
    const L = w < 400 ? 44 : 54, R = 12, Tp = 10, B = 32, pw = w - L - R, ph = h - Tp - B;
    ctx.clearRect(0, 0, w, h); ctx.fillStyle = T.bg; ctx.fillRect(0, 0, w, h);
    if (pw < 20 || ph < 20) return;
    if (!(o.xmax > o.xmin)) o.xmax = o.xmin + 1; if (!(o.ymax > o.ymin)) o.ymax = o.ymin + 1;
    const X = x => L + (x - o.xmin) / (o.xmax - o.xmin) * pw, Y = y => Tp + (1 - (y - o.ymin) / (o.ymax - o.ymin)) * ph;
    ctx.font = '11px ' + T.mono; ctx.lineWidth = 1; ctx.textBaseline = 'middle'; ctx.textAlign = 'right';
    niceTicks(o.ymin, o.ymax, 5).forEach(v => { const yy = Y(v); ctx.strokeStyle = T.border; ctx.beginPath(); ctx.moveTo(L, yy); ctx.lineTo(L + pw, yy); ctx.stroke(); ctx.fillStyle = T.muted; ctx.fillText(tickStr(v), L - 4, yy); });
    ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    niceTicks(o.xmin, o.xmax, Math.max(3, Math.floor(pw / 70))).forEach(v => { const xx = X(v); ctx.strokeStyle = T.border; ctx.beginPath(); ctx.moveTo(xx, Tp); ctx.lineTo(xx, Tp + ph); ctx.stroke(); ctx.fillStyle = T.muted; ctx.fillText(tickStr(v), xx, Tp + ph + 3); });
    ctx.strokeStyle = T.muted; if (o.ymin < 0 && o.ymax > 0) { ctx.beginPath(); ctx.moveTo(L, Y(0)); ctx.lineTo(L + pw, Y(0)); ctx.stroke(); }
    if (o.xmin < 0 && o.xmax > 0) { ctx.beginPath(); ctx.moveTo(X(0), Tp); ctx.lineTo(X(0), Tp + ph); ctx.stroke(); }
    ctx.fillStyle = T.muted; ctx.textAlign = 'center'; ctx.textBaseline = 'bottom'; ctx.fillText(o.xlabel || '', L + pw / 2, h - 1);
    ctx.save(); ctx.translate(10, Tp + ph / 2); ctx.rotate(-PI / 2); ctx.textBaseline = 'top'; ctx.fillText(o.ylabel || '', 0, 0); ctx.restore();
    ctx.save(); ctx.beginPath(); ctx.rect(L, Tp, pw, ph); ctx.clip();
    (o.rects || []).forEach(r => { ctx.globalAlpha = r.alpha || 0.2; ctx.fillStyle = r.color; const x0 = X(r.x0), x1 = X(r.x1), y1 = Y(r.y1), y0 = Y(r.y0); ctx.fillRect(Math.min(x0, x1), Math.min(y0, y1), Math.abs(x1 - x0), Math.abs(y0 - y1)); ctx.globalAlpha = 1; });
    (o.vlines || []).forEach(v => { ctx.strokeStyle = v.color || T.muted; ctx.setLineDash([4, 4]); ctx.beginPath(); ctx.moveTo(X(v.x), Tp); ctx.lineTo(X(v.x), Tp + ph); ctx.stroke(); ctx.setLineDash([]); });
    (o.hlines || []).forEach(v => { ctx.strokeStyle = v.color || T.muted; ctx.lineWidth = 1.6; ctx.setLineDash(v.dash || [5, 3]); ctx.beginPath(); ctx.moveTo(X(v.x0 === undefined ? o.xmin : v.x0), Y(v.y)); ctx.lineTo(X(v.x1 === undefined ? o.xmax : v.x1), Y(v.y)); ctx.stroke(); ctx.setLineDash([]); ctx.lineWidth = 1; });
    (o.series || []).forEach(s => {
      ctx.strokeStyle = s.color; ctx.fillStyle = s.color; ctx.lineWidth = s.width || 2;
      if (s.stem) {
        const base = Y(s.base || 0), dots = s.x.length <= 90;
        for (let i = 0; i < s.x.length; i++) { if (!fin(s.x[i]) || !fin(s.y[i])) continue; const px = X(s.x[i]), py = Y(s.y[i]); ctx.beginPath(); ctx.moveTo(px, base); ctx.lineTo(px, py); ctx.stroke(); if (dots) { ctx.beginPath(); ctx.arc(px, py, 3, 0, 2 * PI); ctx.fill(); } }
        return;
      }
      ctx.setLineDash(s.dash || []); ctx.beginPath(); let pen = false;
      for (let i = 0; i < s.x.length; i++) { const xv = s.x[i], yv = s.y[i]; if (!fin(xv) || !fin(yv)) { pen = false; continue; } const px = X(xv), py = Y(yv); if (pen) ctx.lineTo(px, py); else { ctx.moveTo(px, py); pen = true; } }
      ctx.stroke(); ctx.setLineDash([]);
    });
    ctx.restore(); ctx.font = '11px ' + T.mono;
    (o.points || []).forEach(p => {
      if (!fin(p.x) || !fin(p.y)) return; const px = X(p.x), py = Y(p.y); if (px < L - 1 || px > L + pw + 1) return;
      ctx.fillStyle = p.color; ctx.strokeStyle = T.bg; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.arc(px, py, 4.5, 0, 2 * PI); ctx.fill(); ctx.stroke();
      if (p.label) { const right = px < L + pw * 0.62; ctx.textAlign = right ? 'left' : 'right'; ctx.textBaseline = p.below ? 'top' : 'bottom'; ctx.fillStyle = T.text; ctx.fillText(p.label, px + (right ? 8 : -8), py + (p.below ? 6 : -6)); }
    });
    (o.texts || []).forEach(t => { ctx.fillStyle = t.color || T.muted; ctx.textAlign = t.align || 'left'; ctx.textBaseline = 'top'; ctx.fillText(t.text, X(t.x), Tp + (t.dy || 3)); });
    (o.legend || []).forEach((lg, i) => { const tw = ctx.measureText(lg.text).width; ctx.globalAlpha = 0.8; ctx.fillStyle = T.bg; ctx.fillRect(L + pw - 8 - tw, Tp + ph - 18 - i * 14, tw + 8, 14); ctx.globalAlpha = 1; ctx.fillStyle = lg.color; ctx.textAlign = 'right'; ctx.textBaseline = 'bottom'; ctx.fillText(lg.text, L + pw - 4, Tp + ph - 4 - i * 14); });
  }
  // pole-zero plot with unit circle. pz: {zeros:[{r,m}], poles:[{r,m}]}
  function pzPlot(c, pz) {
    const g = c.prep(), ctx = g.ctx, w = g.w, h = g.h, T = theme();
    ctx.clearRect(0, 0, w, h); ctx.fillStyle = T.bg; ctx.fillRect(0, 0, w, h);
    const mags = pz.zeros.concat(pz.poles).map(o => cabs(o.r)), R = Math.max(1.5, 1.2 * Math.max.apply(null, mags.concat([0])));
    const rad = Math.min(w, h) / 2 - 20, sc = rad / R, cxp = w / 2, cyp = h / 2; if (rad < 20) return;
    ctx.strokeStyle = T.border; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(cxp - rad, cyp); ctx.lineTo(cxp + rad, cyp); ctx.moveTo(cxp, cyp - rad); ctx.lineTo(cxp, cyp + rad); ctx.stroke();
    ctx.strokeStyle = T.muted; ctx.setLineDash([5, 4]); ctx.beginPath(); ctx.arc(cxp, cyp, sc, 0, 2 * PI); ctx.stroke(); ctx.setLineDash([]);
    ctx.font = '11px ' + T.mono; ctx.fillStyle = T.muted; ctx.textAlign = 'left'; ctx.textBaseline = 'top';
    ctx.fillText('Re', cxp + rad - 16, cyp + 3); ctx.fillText('j Im', cxp + 4, cyp - rad); ctx.fillText('1', cxp + sc + 3, cyp + 3); ctx.fillText('−1', cxp - sc - 16, cyp + 3);
    const P = z => [cxp + z.re * sc, cyp - z.im * sc];
    pz.zeros.forEach(o => { const p = P(o.r); ctx.strokeStyle = T.b; ctx.lineWidth = 2; ctx.beginPath(); ctx.arc(p[0], p[1], 5, 0, 2 * PI); ctx.stroke(); if (o.m > 1) { ctx.fillStyle = T.b; ctx.textAlign = 'left'; ctx.textBaseline = 'bottom'; ctx.fillText('(' + o.m + ')', p[0] + 7, p[1] - 4); } });
    pz.poles.forEach(o => { const p = P(o.r); ctx.strokeStyle = T.pink; ctx.lineWidth = 2; ctx.beginPath(); ctx.moveTo(p[0] - 5, p[1] - 5); ctx.lineTo(p[0] + 5, p[1] + 5); ctx.moveTo(p[0] + 5, p[1] - 5); ctx.lineTo(p[0] - 5, p[1] + 5); ctx.stroke(); if (o.m > 1) { ctx.fillStyle = T.pink; ctx.textAlign = 'left'; ctx.textBaseline = 'top'; ctx.fillText('(' + o.m + ')', p[0] + 7, p[1] + 4); } });
    ctx.fillStyle = T.b; ctx.textAlign = 'left'; ctx.textBaseline = 'top'; ctx.fillText('o zero', 6, 5); ctx.fillStyle = T.pink; ctx.fillText('x pole', 6, 19); ctx.fillStyle = T.muted; ctx.fillText('--- |z| = 1', 6, 33);
  }

  /* ---------- small DOM helpers ---------- */
  const E = function () { return FSP.ui.el.apply(FSP.ui, arguments); };
  function layout(root) { const lay = E('div', { class: 'layout' }), ctl = E('div', { class: 'controls' }), stage = E('div', { class: 'stage' }); lay.appendChild(ctl); lay.appendChild(stage); root.appendChild(lay); return { ctl, stage }; }
  function wrapCanvas(parent, label, height, minW) {
    const wrap = E('div', { class: 'canvas-wrap' }); if (minW) wrap.style.overflowX = 'auto'; parent.appendChild(wrap);
    const c = FSP.ui.canvas(wrap, { height }); c.cv.setAttribute('role', 'img'); c.cv.setAttribute('aria-label', label); c.wrap = wrap; return c;
  }
  const note = t => E('div', { class: 'note', text: t });
  const mkMsg = parent => { const m = E('div', { class: 'msg bad', hidden: '' }); parent.appendChild(m); return m; };
  function showMsg(m, lines, kind) { const t = Array.isArray(lines) ? lines.join('  ') : lines; m.hidden = !t; m.textContent = t || ''; m.className = 'msg ' + (kind || 'bad'); }
  function hudSet(el0, rows) { while (el0.firstChild) el0.removeChild(el0.firstChild); rows.forEach(r => el0.appendChild(E('div', null, E('span', { text: r[0] }), E('span', { text: r[1] })))); }
  function textIn(parent, label, value, onChange) {
    const inp = E('input', { type: 'text', 'aria-label': label, spellcheck: 'false', autocomplete: 'off' }); inp.value = value;
    parent.appendChild(E('div', { class: 'dsp2-tctl' }, E('label', { text: label }), inp));
    inp.addEventListener('input', onChange); return inp;
  }
  const selOpt = (sel, v) => { if (typeof v === 'string' && Array.prototype.some.call(sel.options, o => o.value === v)) sel.value = v; };
  const setNum = (s, v) => { const x = parseFloat(v); if (fin(x)) s.set(x, true); };
  function setStyle() {
    if (document.getElementById('dsp2-style')) return;
    const s = E('style', { id: 'dsp2-style' });
    s.textContent = '.dsp2-tctl{display:grid;grid-template-columns:72px minmax(0,1fr);gap:6px;align-items:center;margin:4px 0}.dsp2-tctl label{color:var(--muted);font-size:12px}.dsp2-tctl input{width:100%;box-sizing:border-box}' +
      '.dsp2-tbl{border-collapse:collapse;font:12px var(--mono);width:100%}.dsp2-tbl th,.dsp2-tbl td{border-bottom:1px solid var(--border);padding:3px 8px;text-align:right;white-space:nowrap}.dsp2-tbl th{color:var(--muted);font-weight:400}.dsp2-tbl tr.sel td{color:var(--accent)}' +
      '.dsp2-scroll{overflow-x:auto;background:var(--panel);border:1px solid var(--border);border-radius:10px;padding:4px 6px}.dsp2-box{font:12px/1.5 var(--mono);white-space:pre-wrap;word-break:break-word;background:var(--panel);border:1px solid var(--border);border-radius:10px;padding:8px 12px}' +
      '#panel-dspdesign select{max-width:100%;min-width:0;width:100%;text-overflow:ellipsis}#panel-dspdesign .grid2{grid-template-columns:repeat(auto-fit,minmax(min(300px,100%),1fr))}.dsp2-ok{color:var(--ok)}.dsp2-no{color:var(--bad)}';
    document.head.appendChild(s);
  }

  /* ======================= 1. sampling & aliasing ======================= */
  function buildSamp(root) {
    const lo = layout(root), ctl = lo.ctl, stage = lo.stage, st = {}, ui = {};
    const chg = () => { upd(); FSP.state.touch(); };
    const fs = FSP.ui.fieldset(ctl, 'Tone & sampling');
    st.F = FSP.ui.slider(fs, { label: 'Tone F', min: 1, max: 100000, value: 7000, unit: 'Hz', log: true, digits: 6, onInput: chg });
    st.Fs = FSP.ui.slider(fs, { label: 'Fs', min: 1, max: 100000, value: 10000, unit: 'Hz', log: true, digits: 6, onInput: chg });
    st.ph = FSP.ui.slider(fs, { label: 'Phase φ', min: -180, max: 180, step: 1, value: 40, unit: '°', onInput: chg });
    fs.appendChild(note('x(t) = cos(2πFt + φ), sampled at t = nTs, Ts = 1/Fs. Values accept SI suffixes (7k = 7000).'));
    const fo = FSP.ui.fieldset(ctl, 'Periodicity of cos(ωn)');
    ui.om = textIn(fo, 'ω (rad)', '3pi/8', chg);
    fo.appendChild(note('Write 3pi/8, 0.75pi, pi/4 or a plain number of radians such as 1.2. The sampled tone above has ω = 2πF/Fs.'));
    ui.msg = mkMsg(ctl);
    ui.c1 = wrapCanvas(stage, 'Continuous tone, its samples and the reconstructed alias', 250); ui.c1.onResize(upd);
    ui.c2 = wrapCanvas(stage, 'Spectrum of the sampled tone with replicas', 220); ui.c2.onResize(upd);
    ui.hud = E('div', { class: 'hud' }); stage.appendChild(ui.hud); ui.work = FSP.ui.working(stage);
    function upd() {
      if (!ui.hud) return;
      const T = theme(), F = st.F.get(), Fs = st.Fs.get(), ph = st.ph.get(), phr = ph * PI / 180, al = DM.alias(F, Fs), f = F / Fs, om = TAU * f;
      // --- time plot in units of sample index
      const cyc = f * 16, cx0 = [], cy0 = [], series = [], texts = [];
      if (cyc <= 250) { const np = Math.min(8000, Math.max(600, Math.ceil(cyc * 24))); for (let i = 0; i <= np; i++) { const tau = 16 * i / np; cx0.push(tau); cy0.push(Math.cos(TAU * f * tau + phr)); } series.push({ x: cx0, y: cy0, color: T.c, width: 1.5 }); }
      else texts.push({ x: 0.2, y: 0, text: 'tone is >15 cycles per sample: continuous curve omitted', color: T.warn });
      if (al && al.aliased) { const ax = [], ay = []; for (let i = 0; i <= 800; i++) { const tau = 16 * i / 800; ax.push(tau); ay.push(Math.cos(TAU * (al.Fa / Fs) * tau + al.sign * phr)); } series.push({ x: ax, y: ay, color: T.pink, width: 1.6, dash: [6, 4] }); }
      const sx = [], sy = []; for (let n = 0; n <= 16; n++) { sx.push(n); sy.push(Math.cos(TAU * f * n + phr)); }
      series.push({ x: sx, y: sy, color: T.a, stem: true, width: 1.6 });
      plot(ui.c1, { xmin: 0, xmax: 16, ymin: -1.9, ymax: 1.3, xlabel: 't·Fs  (= sample index n)', ylabel: 'amplitude', series, texts,
        legend: [{ text: 'x[n] samples', color: T.a }, { text: 'x(t) = cos(2πFt+φ)', color: T.c }].concat(al && al.aliased ? [{ text: 'alias cos(2πFa t ± φ)', color: T.pink }] : []) });
      // --- spectrum: lines at ±F + m·Fs, amplitude 1/2 each
      const xr = Math.min(Math.max(1.5 * Fs, 1.15 * F), 6 * Fs), orig = { x: [], y: [], color: T.c, stem: true, width: 2.4 }, base = { x: [], y: [], color: T.pink, stem: true, width: 2 }, other = { x: [], y: [], color: T.muted, stem: true, width: 1.4 };
      [1, -1].forEach(s => { const f0 = s * F, m0 = Math.ceil((-xr - f0) / Fs), m1 = Math.floor((xr - f0) / Fs); for (let m = Math.max(m0, m1 - 400); m <= m1; m++) { const fr = f0 + m * Fs, tgt = m === 0 ? orig : (Math.abs(fr) <= Fs / 2 + 1e-9 * Fs ? base : other); tgt.x.push(fr); tgt.y.push(0.5); } });
      plot(ui.c2, { xmin: -xr, xmax: xr, ymin: 0, ymax: 0.72, xlabel: 'f (Hz)', ylabel: 'amplitude', rects: [{ x0: -Fs / 2, x1: Fs / 2, y0: 0, y1: 0.72, color: T.c, alpha: 0.1 }],
        series: [other, base, orig], texts: [{ x: 0, y: 0, text: 'baseband −Fs/2…Fs/2', color: T.muted, align: 'center' }],
        legend: [{ text: 'original ±F', color: T.c }, { text: 'alias in baseband', color: T.pink }, { text: 'replicas ±F + m·Fs', color: T.muted }] });
      // --- results
      const pF = DM.periodicityF(F, Fs), po = DM.parseOmega(ui.om.value), pw = po ? DM.periodicityOmega(po) : null, msgs = [];
      if (!po) msgs.push('ω: could not read "' + ui.om.value + '" (try 3pi/8, 0.4pi or 1.2).');
      showMsg(ui.msg, msgs);
      const nyq = al.atNyquist ? ' (exactly Fs/2: the samples depend on φ; amplitude seen = |cos φ|)' : '';
      hudSet(ui.hud, [
        ['Nyquist rate 2F', g4(2 * F) + ' Hz'], ['Fs', g4(Fs) + ' Hz  → ' + (Fs > 2 * F ? 'OK, Fs > 2F' : 'ALIASING, Fs ≤ 2F')],
        ['Nyquist freq Fs/2', g4(Fs / 2) + ' Hz'], ['apparent freq Fa', g4(al.Fa) + ' Hz' + (al.aliased ? '  (alias of ' + g4(F) + ' Hz)' : '')],
        ['digital f = F/Fs', g4(f) + ' cycles/sample'], ['ω = 2πF/Fs', g4(om) + ' rad/sample = ' + g4(om / PI) + 'π'],
        ['reconstructed phase', (al.sign < 0 && al.aliased ? '−φ = ' + g4(-ph) : g4(ph)) + '°'],
        ['sampled period N', pF && pF.periodic ? pF.N + ' samples (' + pF.k + ' cycle' + (pF.k === 1 ? '' : 's') + ')' : 'not periodic'],
        ['ω entry: ' + ui.om.value, pw ? (pw.periodic ? 'periodic, N = ' + pw.N + ' (ω/2π = ' + pw.f.p + '/' + pw.f.q + ')' : 'NOT periodic') : '—'],
      ]);
      const L = [];
      L.push('1. Sampling theorem: no aliasing needs Fs > 2F.  2F = 2 × ' + g4(F) + ' = ' + g4(2 * F) + ' Hz,  Fs = ' + g4(Fs) + ' Hz  →  ' + (Fs > 2 * F ? 'satisfied.' : 'violated, the tone folds.'));
      L.push('   Nyquist frequency Fs/2 = ' + g4(Fs) + ' / 2 = ' + g4(Fs / 2) + ' Hz.' + nyq);
      L.push('2. Digital frequency: f = F/Fs = ' + g4(F) + ' / ' + g4(Fs) + ' = ' + g4(f) + ' cycles/sample;  ω = 2πf = ' + g4(om) + ' rad/sample (' + g4(om / PI) + 'π).');
      L.push('3. Alias: the sample sequence cannot distinguish F from F − k·Fs.  k = round(F/Fs) = round(' + g4(F / Fs) + ') = ' + al.k + '.');
      L.push('   Fa = |F − k·Fs| = |' + g4(F) + ' − ' + al.k + ' × ' + g4(Fs) + '| = ' + g4(al.Fa) + ' Hz (always in 0…Fs/2).');
      if (al.aliased) L.push('   F − k·Fs = ' + g4(F - al.k * Fs) + ' Hz ' + (al.sign < 0 ? '< 0, so cos(2π(−Fa)nTs + φ) = cos(2πFa nTs − φ): the reconstructed tone has phase −φ = ' + g4(-ph) + '°.' : '> 0, phase unchanged.'));
      else L.push('   F ≤ Fs/2, so Fa = F: ideal reconstruction returns the original tone.');
      L.push('4. Ideal D/A reconstruction (low-pass at Fs/2) outputs cos(2π·' + g4(al.Fa) + '·t + (' + g4(al.aliased ? al.sign * ph : ph) + '°)).');
      L.push('5. Spectrum of the sampled signal repeats every Fs: lines at ±F + m·Fs = ±' + g4(F) + ' + m × ' + g4(Fs) + ' Hz (amplitude ½ each for a real cosine); the one inside −Fs/2…Fs/2 is the alias at ±' + g4(al.Fa) + ' Hz.');
      if (pF && pF.periodic) L.push('6. Periodicity of x[n] = cos(2πf n + φ): periodic iff f = p/q rational. f = ' + g4(f) + ' = ' + pF.p + '/' + pF.q + ' in lowest terms → N = q = ' + pF.N + ' samples (k = p = ' + pF.k + ' cycles per period): N·f = ' + pF.N + ' × ' + pF.p + '/' + pF.q + ' = ' + pF.k + ' an integer.');
      else L.push('6. Periodicity: f = F/Fs = ' + g4(f) + ' is not a ratio of small integers → treated as not periodic.');
      if (po && pw) { L.push('7. Entry ω = ' + ui.om.value + ':  ' + (pw.periodic ? 'ω/2π = ' + pw.f.p + '/' + pw.f.q + ' rational → periodic with N = ' + pw.N + ' (smallest N with ωN = 2πk, k = ' + pw.k + ').' : pw.reason)); }
      ui.work.set(L);
    }
    return {
      update: upd,
      get: () => ({ F: st.F.get(), Fs: st.Fs.get(), ph: st.ph.get(), om: ui.om.value }),
      set(o) { setNum(st.F, o.F); setNum(st.Fs, o.Fs); setNum(st.ph, o.ph); if (typeof o.om === 'string' && o.om.length < 40) ui.om.value = o.om; },
    };
  }

  /* ======================= 2. difference equation & z-transform ======================= */
  const PRESETS = [
    ['custom', 'Custom'], ['lp1', 'First-order IIR: y[n] = 0.9y[n−1] + x[n]'], ['res', 'Resonator: poles 0.9e^{±jπ/4}'], ['dp', 'Double pole 1/(1−0.5z⁻¹)²'],
    ['two', 'Two real poles 1/((1−0.5z⁻¹)(1−0.25z⁻¹))'], ['ma', 'FIR moving average (4 taps)'], ['bad', 'Unstable: pole at 1.2'], ['imp', 'Improper: (1+2z⁻¹+z⁻²+0.5z⁻³)/(1−0.5z⁻¹+0.06z⁻²)'],
  ];
  const PRESET_VAL = { lp1: ['1', '1, -0.9'], res: ['1', '1, -1.2728, 0.81'], dp: ['1', '1, -1, 0.25'], two: ['1', '1, -0.75, 0.125'], ma: ['0.25, 0.25, 0.25, 0.25', '1'], bad: ['1', '1, -1.2'], imp: ['1, 2, 1, 0.5', '1, -0.5, 0.06'] };
  // closed-form text of h[n] for the chosen ROC
  function closedForm(pf, roc) {
    const parts = []; pf.direct.forEach((c, r) => { if (c !== 0) parts.push(g4(c) + ' δ[n' + (r ? '−' + r : '') + ']'); });
    const done = new Set();
    pf.terms.forEach((t, i) => {
      if (done.has(i) || cabs(t.C) < 1e-10) return; const right = cabs(t.d) <= roc.lo * (1 + 1e-9) + 1e-300, step = right ? 'u[n]' : 'u[−n−1]', sg = right ? 1 : -1;
      let fac = ''; if (t.order > 1) { const f = []; for (let j = 1; j < t.order; j++) f.push('(n+' + j + ')'); fac = f.join('') + (t.order > 2 ? '/' + [1, 1, 2, 6, 24, 120, 720][t.order - 1] : '') + ' '; }
      if (t.mult === 1 && Math.abs(t.d.im) > 1e-9) {
        if (t.d.im < 0) return; const j = pf.terms.findIndex((u, q) => q !== i && u.mult === 1 && Math.abs(u.d.re - t.d.re) < 1e-9 && Math.abs(u.d.im + t.d.im) < 1e-9); if (j >= 0) done.add(j);
        parts.push((sg < 0 ? '−' : '') + g4(2 * cabs(t.C)) + ' (' + g4(cabs(t.d)) + ')ⁿ cos(' + g4(Math.atan2(t.d.im, t.d.re)) + 'n ' + (carg(t.C) < 0 ? '− ' : '+ ') + g4(Math.abs(carg(t.C))) + ') ' + step);
      } else {
        const cs = Math.abs(t.C.im) > 1e-9 ? '(' + cstr(t.C) + ')' : g4((sg < 0 ? -1 : 1) * t.C.re), ds = Math.abs(t.d.im) > 1e-9 ? '(' + cstr(t.d) + ')' : '(' + g4(t.d.re) + ')';
        parts.push((Math.abs(t.C.im) > 1e-9 && sg < 0 ? '−' : '') + cs + ' ' + fac + ds + 'ⁿ ' + step);
      }
    });
    return parts.length ? parts.join('\n       + ').replace(/\+ -/g, '− ') : '0';
  }
  function buildZ(root) {
    const lo = layout(root), ctl = lo.ctl, stage = lo.stage, st = {}, ui = {}; let rocIdx = -1, rocSig = '';
    const chg = () => { upd(); FSP.state.touch(); };
    const fs = FSP.ui.fieldset(ctl, 'System');
    ui.preset = FSP.ui.select(fs, 'Preset', PRESETS, 'lp1', v => { if (PRESET_VAL[v]) { ui.b.value = PRESET_VAL[v][0]; ui.a.value = PRESET_VAL[v][1]; rocIdx = -1; } chg(); });
    ui.b = textIn(fs, 'b (num)', '1', () => { ui.preset.value = 'custom'; rocIdx = -1; chg(); });
    ui.a = textIn(fs, 'a (den)', '1, -0.9', () => { ui.preset.value = 'custom'; rocIdx = -1; chg(); });
    fs.appendChild(note('H(z) = Σ b_k z⁻ᵏ / Σ a_k z⁻ᵏ, i.e. a₀y[n] = Σ b_k x[n−k] − Σ_{k≥1} a_k y[n−k]. Separate by commas or spaces; fractions like 1/3 are fine.'));
    ui.roc = FSP.ui.select(fs, 'ROC', [['0', '—']], '0', v => { rocIdx = parseInt(v, 10); chg(); });
    st.n = FSP.ui.slider(fs, { label: 'Samples', min: 8, max: 64, step: 1, value: 24, unit: 'n', onInput: chg });
    ui.scale = FSP.ui.select(fs, '|H| axis', [['db', 'dB'], ['lin', 'linear']], 'db', chg);
    ui.msg = mkMsg(ctl); ui.warn = E('div', { class: 'msg warn', hidden: '' }); ctl.appendChild(ui.warn);
    const g2 = E('div', { class: 'grid2' }); stage.appendChild(g2);
    ui.cpz = wrapCanvas(g2, 'Pole-zero plot with unit circle', 300); ui.cpz.onResize(upd);
    ui.ch = wrapCanvas(g2, 'Impulse response h[n] for the chosen ROC', 300); ui.ch.onResize(upd);
    ui.cs = wrapCanvas(stage, 'Step response of the causal system', 200); ui.cs.onResize(upd);
    const g3_ = E('div', { class: 'grid2' }); stage.appendChild(g3_);
    ui.cm = wrapCanvas(g3_, 'Magnitude response', 230); ui.cm.onResize(upd);
    ui.cp = wrapCanvas(g3_, 'Phase response', 230); ui.cp.onResize(upd);
    ui.hud = E('div', { class: 'hud' }); stage.appendChild(ui.hud);
    ui.closed = E('div', { class: 'dsp2-box' }); stage.appendChild(ui.closed);
    ui.work = FSP.ui.working(stage);
    function upd() {
      if (!ui.hud) return;
      const T = theme(), bl = DM.parseList(ui.b.value), al = DM.parseList(ui.a.value);
      if (!bl || !al) { showMsg(ui.msg, 'Coefficients must be numbers separated by commas/spaces (fractions like 1/3 allowed, at most 64 values).'); return; }
      const q = DM.normBA(bl, al); if (q.err) { showMsg(ui.msg, q.err); return; }
      let pz, pf; try { pz = DM.pzk(q.b, q.a); pf = DM.partialFractions(q.b, q.a); } catch (e) { showMsg(ui.msg, 'Could not factorise this system: ' + e.message); return; }
      showMsg(ui.msg, '');
      const roc = DM.rocOptions(pf), K = roc.length - 1; if (rocIdx < 0 || rocIdx > K) rocIdx = K;
      const sig = roc.map(r => r.label).join('|'); if (sig !== rocSig) { rocSig = sig; while (ui.roc.firstChild) ui.roc.removeChild(ui.roc.firstChild); roc.forEach(r => ui.roc.appendChild(E('option', { value: String(r.idx), text: q.N === 0 ? 'all z except 0 — causal FIR' : r.label }))); }
      ui.roc.value = String(rocIdx); const rc = roc[rocIdx], Nn = Math.round(st.n.get());
      const stab = DM.stability(pz.poles), causal = rocIdx === K;
      // h[n] for chosen ROC
      const n0 = causal ? 0 : rocIdx === 0 ? -(Nn - 1) : -Math.floor(Nn / 2), n1 = causal ? Nn - 1 : rocIdx === 0 ? 1 : Math.ceil(Nn / 2), hx = [], hy = [];
      for (let n = n0; n <= n1; n++) { hx.push(n); hy.push(DM.inverseZAt(pf, rc.lo, n)); }
      const hm = Math.max.apply(null, hy.map(Math.abs).concat([1e-12])), hmin = Math.min.apply(null, hy.concat([0])), hmax = Math.max.apply(null, hy.concat([0]));
      const sres = DM.diffEqResidual(q.b, q.a, n => DM.inverseZAt(pf, rc.lo, n), -15, 15), tol = 1e-8 * Math.max(1, Math.min(hm, 1e12));
      let cerr = NaN; if (causal) { const hr = DM.impulse(q.b, q.a, Nn); cerr = 0; for (let n = 0; n < Nn; n++) cerr = Math.max(cerr, Math.abs(hr[n] - hy[n])); }
      pzPlot(ui.cpz, pz);
      plot(ui.ch, { xmin: n0 - 0.7, xmax: n1 + 0.7, ymin: hmin - 0.1 * (hmax - hmin || 1), ymax: hmax + 0.12 * (hmax - hmin || 1), xlabel: 'n', ylabel: 'h[n]', series: [{ x: hx, y: hy, color: T.c, stem: true, width: 1.8 }] });
      const sr = DM.stepResp(q.b, q.a, Nn), sx = [], sy = Array.from(sr); for (let n = 0; n < Nn; n++) sx.push(n);
      const smx = Math.max.apply(null, sy.concat([0])), smn = Math.min.apply(null, sy.concat([0]));
      plot(ui.cs, { xmin: -0.7, xmax: Nn - 0.3, ymin: smn - 0.1 * (smx - smn || 1), ymax: smx + 0.12 * (smx - smn || 1), xlabel: 'n  (causal system, zero initial state)', ylabel: 's[n]', series: [{ x: sx, y: sy, color: T.a, stem: true, width: 1.6 }] });
      const fr = DM.freqResp(q.b, q.a, 512), wx = Array.from(fr.w).map(v => v / PI);
      const finMag = Array.from(fr.mag).filter(fin), db = Array.from(fr.mag).map(m => (fin(m) ? 20 * Math.log10(Math.max(m, 1e-9)) : NaN)), dbf = db.filter(fin);
      const isDb = ui.scale.value === 'db';
      if (isDb) { const hi = Math.max.apply(null, dbf.concat([-60])) + 5, lo2 = Math.max(Math.min.apply(null, dbf.concat([0])) - 3, hi - 90); plot(ui.cm, { xmin: 0, xmax: 1, ymin: lo2, ymax: hi, xlabel: 'ω/π', ylabel: '|H(e^{jω})| (dB)', series: [{ x: wx, y: db, color: T.c, width: 2 }] }); }
      else plot(ui.cm, { xmin: 0, xmax: 1, ymin: 0, ymax: (Math.max.apply(null, finMag.concat([1e-6])) * 1.1), xlabel: 'ω/π', ylabel: '|H(e^{jω})|', series: [{ x: wx, y: Array.from(fr.mag), color: T.c, width: 2 }] });
      const ph = Array.from(fr.ph); for (let i = 1; i < ph.length; i++) { if (!fin(ph[i]) || !fin(ph[i - 1])) continue; while (ph[i] - ph[i - 1] > PI) { for (let j = i; j < ph.length; j++) ph[j] -= TAU; } while (ph[i] - ph[i - 1] < -PI) { for (let j = i; j < ph.length; j++) ph[j] += TAU; } }
      const phn = ph.map(v => v / PI), pf0 = phn.filter(fin);
      plot(ui.cp, { xmin: 0, xmax: 1, ymin: Math.min.apply(null, pf0.concat([0])) - 0.1, ymax: Math.max.apply(null, pf0.concat([0])) + 0.1, xlabel: 'ω/π', ylabel: '∠H (rad/π, unwrapped)', series: [{ x: wx, y: phn, color: T.pink, width: 2 }] });
      // warnings
      const w = []; if (stab.verdict !== 'stable') w.push('Causal system is ' + stab.verdict.toUpperCase() + ': max |pole| = ' + g4(stab.rmax) + (stab.verdict === 'marginal' ? ' (on the unit circle, h[n] does not decay, H(e^{jω}) has a singularity)' : ' ≥ 1, h[n] grows without bound.'));
      if (!rc.stable) w.push('The chosen ROC does not contain |z| = 1, so this h[n] is not absolutely summable and H(e^{jω}) (plotted from the formula) is not its DTFT.');
      showMsg(ui.warn, w, 'warn');
      ui.closed.textContent = 'ROC: ' + (q.N === 0 ? 'all z except z = 0 (FIR)' : rc.label) + (rc.stable ? '   [contains the unit circle → stable]' : '   [excludes the unit circle → unstable]') + '\nh[n] = ' + closedForm(pf, rc) + '\ncheck: difference-equation residual over −15…15 = ' + (sres < tol ? sres.toExponential(1) + '  OK' : sres.toExponential(1) + '  (large: numerical growth)') + (causal ? ';  max|closed form − recursion| = ' + cerr.toExponential(1) : '');
      hudSet(ui.hud, [
        ['order (M, N)', q.M + ', ' + q.N], ['stability (causal)', stab.verdict + ' (max|p| = ' + g4(stab.rmax) + ')'],
        ['zeros', pz.zeros.map(o => cstr(o.r) + (o.m > 1 ? ' (×' + o.m + ')' : '')).join(',  ') || 'none'], ['poles', pz.poles.map(o => cstr(o.r) + (o.m > 1 ? ' (×' + o.m + ')' : '')).join(',  ') || 'none'],
        ['gain K', g4(pz.K) + (pz.delay ? '  with z⁻' + supn(pz.delay) : '')], ['DC gain H(1)', g4(DM.freqAt(q.b, q.a, 0).re)],
        ['|H(e^{jπ})|', g4(cabs(DM.freqAt(q.b, q.a, PI)))], ['ROC chosen', q.N === 0 ? 'all z ≠ 0' : rc.label],
      ]);
      // working
      const L = [];
      L.push('Difference equation (a₀ = 1): ' + lin(q.a.map((c, k) => [c, 'y[' + idx(k) + ']'])) + ' = ' + lin(q.b.map((c, k) => [c, 'x[' + idx(k) + ']'])));
      L.push('   → y[n] = ' + lin(q.b.map((c, k) => [c, 'x[' + idx(k) + ']']).concat(q.a.map((c, k) => [-c, 'y[' + idx(k) + ']']).slice(1))));
      L.push('z-transform (shift property x[n−k] ↔ z⁻ᵏX(z)):  H(z) = Y/X = (' + polyStr(q.b) + ') / (' + polyStr(q.a) + ')');
      const fac = o => (cabs(o.r) === 0 ? 'z' : '(z − ' + (Math.abs(o.r.im) > 1e-12 ? '(' + cstr(o.r) + ')' : g4(o.r.re)) + ')') + (o.m > 1 ? supn(o.m) : '');
      L.push('Multiply numerator and denominator by z^' + Math.max(q.M, q.N) + ' (max(M, N)) to get polynomials in z, then factor:  H(z) = K ∏(z − zᵢ) / ∏(z − pᵢ) = ' + g4(pz.K) + (pz.delay ? ' z⁻' + supn(pz.delay) : '') + ' · ' + (pz.zeros.map(fac).join(' ') || '1') + ' / (' + (pz.poles.map(fac).join(' ') || '1') + ')');
      L.push('Poles (roots of the denominator):  ' + (pz.poles.map(o => 'p = ' + cstr(o.r) + ', |p| = ' + g4(cabs(o.r)) + (o.m > 1 ? ' (double+ ×' + o.m + ')' : '')).join(';  ') || 'none (FIR: only the trivial poles at z = 0 for M > 0)'));
      L.push('Zeros:  ' + (pz.zeros.map(o => 'z = ' + cstr(o.r) + (o.m > 1 ? ' (×' + o.m + ')' : '')).join(';  ') || 'none'));
      L.push('Stability of the causal system: BIBO stable iff all poles lie inside |z| = 1.  max|p| = ' + g4(stab.rmax) + ' → ' + stab.verdict + '.');
      L.push('');
      if (q.N === 0) L.push('FIR: h[n] = b_n, n = 0…' + q.M + ' = [' + q.b.map(g4).join(', ') + '];  ROC is everything except z = 0.');
      else {
        L.push('Inverse z-transform by partial fractions (O&S 3.3.2):  H(z) = ' + (pf.direct.length ? 'Σ B_r z⁻ʳ + ' : '') + 'Σ_i Σ_m C_{i,m} / (1 − d_i z⁻¹)^m');
        if (pf.direct.length) L.push('  M ≥ N: polynomial part Σ B_r z⁻ʳ = ' + polyStr(pf.direct) + '   (δ-terms at n = 0…' + (pf.direct.length - 1) + ')');
        pf.terms.forEach(t => { const simple = t.mult === 1; L.push('  pole d = ' + cstr(t.d) + (t.mult > 1 ? ' (multiplicity ' + t.mult + ', term m = ' + t.order + ')' : '') + ':  C = ' + cstr(t.C) + (simple ? '   [C = (1 − d z⁻¹) H(z) at z = d]' : '   [repeated pole: O&S eq. 3.47, derivatives of (1−d w)^s H in w = z⁻¹]')); });
        L.push('ROC choices come from the pole radii: ' + roc.map(r => r.label).join('  |  '));
        L.push('Chosen ROC: ' + rc.label + '.  A pole of radius ≤ ' + g4(rc.lo) + ' (inside the ROC\'s inner edge) gives a right-sided term;  outside it, a left-sided term:');
        L.push('  right-sided: dⁿu[n] ↔ 1/(1 − d z⁻¹), |z| > |d|;     left-sided: −dⁿu[−n−1] ↔ 1/(1 − d z⁻¹), |z| < |d|');
        L.push('  repeated: (n+1)…(n+m−1)/(m−1)! · dⁿ u[n] ↔ 1/(1 − d z⁻¹)^m');
        L.push('h[n] = ' + closedForm(pf, rc));
        const ex = []; for (let n = Math.max(n0, -3); n <= Math.min(n1, 4); n++) ex.push('h[' + n + '] = ' + g4(DM.inverseZAt(pf, rc.lo, n))); L.push('Values: ' + ex.join(', '));
        L.push('Verification: substituting this h[n] into Σ a_k h[n−k] = b_n for n = −15…15 leaves a worst residual of ' + sres.toExponential(1) + '.');
      }
      L.push('');
      L.push('Frequency response (only exists when the ROC contains the unit circle): H(e^{jω}) = Σ b_k e^{−jωk} / Σ a_k e^{−jωk}.  H(1) = ' + g4(DM.freqAt(q.b, q.a, 0).re) + ',  |H(e^{jπ})| = ' + g4(cabs(DM.freqAt(q.b, q.a, PI))) + '.');
      L.push('Step response s[n] = Σ_{k=0}^{n} h[k] (causal); final value if stable = H(1) = ' + g4(DM.freqAt(q.b, q.a, 0).re) + '.');
      ui.work.set(L);
    }
    return {
      update: upd,
      get: () => ({ b: ui.b.value, a: ui.a.value, roc: rocIdx, n: st.n.get(), sc: ui.scale.value, pre: ui.preset.value }),
      set(o) {
        if (typeof o.b === 'string' && o.b.length < 400) ui.b.value = o.b; if (typeof o.a === 'string' && o.a.length < 400) ui.a.value = o.a;
        const r = parseInt(o.roc, 10); rocIdx = fin(r) && r >= 0 && r < 50 ? r : -1; setNum(st.n, o.n); selOpt(ui.scale, o.sc); selOpt(ui.preset, o.pre);
      },
    };
  }

  /* ======================= 3. FIR design ======================= */
  const RDEF = { lp: ['0.3, 0.4', 'ωp, ωs (×π)'], hp: ['0.5, 0.6', 'ωs, ωp (×π)'], bp: ['0.2, 0.3, 0.6, 0.7', 'ωs1 ωp1 ωp2 ωs2 (×π)'], bs: ['0.3, 0.4, 0.6, 0.7', 'ωp1 ωs1 ωs2 ωp2 (×π)'] };
  const WINEQ = {
    rect: 'w[n] = 1, 0 ≤ n ≤ M', hann: 'w[n] = 0.5 − 0.5 cos(2πn/M)', hamming: 'w[n] = 0.54 − 0.46 cos(2πn/M)',
    blackman: 'w[n] = 0.42 − 0.5 cos(2πn/M) + 0.08 cos(4πn/M)', kaiser: 'w[n] = I₀(β√(1 − ((n−M/2)/(M/2))²)) / I₀(β)',
  };
  function buildFIR(root) {
    const lo = layout(root), ctl = lo.ctl, stage = lo.stage, st = {}, ui = {};
    const chg = () => { upd(); FSP.state.touch(); };
    const f0 = FSP.ui.fieldset(ctl, 'Method & response');
    ui.mode = FSP.ui.select(f0, 'Method', [['win', 'Window method'], ['fs', 'Frequency sampling']], 'win', chg);
    ui.resp = FSP.ui.select(f0, 'Response', [['lp', 'Lowpass'], ['hp', 'Highpass'], ['bp', 'Bandpass'], ['bs', 'Bandstop']], 'lp', v => { ui.edges.value = RDEF[v][0]; chg(); });
    ui.edges = textIn(f0, RDEF.lp[1], RDEF.lp[0], chg);
    f0.appendChild(note('Band edges in ω/π, ascending. LP: passband edge then stopband edge. HP: stopband then passband. BP: stop, pass, pass, stop. BS: pass, stop, stop, pass.'));
    ui.win = FSP.ui.select(f0, 'Window', [['rect', 'Rectangular'], ['hann', 'Hann'], ['hamming', 'Hamming'], ['blackman', 'Blackman'], ['kaiser', 'Kaiser']], 'hamming', chg);
    const f1 = FSP.ui.fieldset(ctl, 'Tolerances & order');
    st.d1 = FSP.ui.slider(f1, { label: 'Passband δ₁', min: 0.0001, max: 0.3, value: 0.01, log: true, digits: 4, onInput: chg });
    st.As = FSP.ui.slider(f1, { label: 'Stop As', min: 20, max: 120, step: 1, value: 50, unit: 'dB', onInput: chg });
    ui.om = FSP.ui.select(f1, 'Order', [['auto', 'Auto (from spec)'], ['manual', 'Manual M']], 'auto', chg);
    st.M = FSP.ui.slider(f1, { label: 'Order M', min: 2, max: 200, step: 1, value: 40, onInput: chg });
    st.T = FSP.ui.slider(f1, { label: 'Transition T', min: 0, max: 1, step: 0.01, value: 0.4, onInput: chg });
    f1.appendChild(note('Filter length L = M + 1 taps, delay α = M/2. Frequency sampling uses the manual M and the transition sample T.'));
    ui.msg = mkMsg(ctl); ui.warn = E('div', { class: 'msg warn', hidden: '' }); ctl.appendChild(ui.warn);
    ui.ch = wrapCanvas(stage, 'FIR impulse response h[n]', 220); ui.ch.onResize(upd);
    ui.cw = wrapCanvas(stage, 'Window or frequency samples', 200); ui.cw.onResize(upd);
    ui.cm = wrapCanvas(stage, 'Magnitude response in dB with specification', 290); ui.cm.onResize(upd);
    ui.hud = E('div', { class: 'hud' }); stage.appendChild(ui.hud); ui.work = FSP.ui.working(stage);
    function upd() {
      if (!ui.hud) return;
      const T = theme(), isWin = ui.mode.value === 'win', resp = ui.resp.value, fe = DM.parseList(ui.edges.value, 8), need = resp === 'lp' || resp === 'hp' ? 2 : 4;
      ui.edges.previousSibling.textContent = RDEF[resp][1]; ui.edges.setAttribute('aria-label', RDEF[resp][1]); ui.win.parentNode.hidden = !isWin; ui.om.parentNode.hidden = !isWin; st.d1.el.hidden = !isWin; st.As.el.hidden = !isWin; st.T.el.hidden = isWin;
      st.M.el.hidden = isWin && ui.om.value === 'auto';
      if (!fe || fe.length !== need) { showMsg(ui.msg, 'Enter ' + need + ' band edges (ω/π, between 0 and 1) for ' + DM.RESP[resp].toLowerCase() + ', e.g. ' + RDEF[resp][0] + '.'); return; }
      const edges = fe.map(v => v * PI);
      const spec = isWin ? { resp, edges, win: ui.win.value, d1: st.d1.get(), As: st.As.get(), orderMode: ui.om.value, M: st.M.get() } : { resp, edges, M: st.M.get(), T: st.T.get() };
      const r = isWin ? DM.designWindow(spec) : DM.designFreqSamp(spec);
      if (r.err) { showMsg(ui.msg, r.err); return; }
      showMsg(ui.msg, '');
      if (isWin && spec.orderMode === 'auto' && r.M <= 200) st.M.set(r.M, true);
      const w2 = r.notes.slice(); if (r.type.bad && r.type.bad.indexOf(resp) >= 0) w2.push('FIR Type ' + roman(r.type.type) + ' cannot realise a ' + DM.RESP[resp].toLowerCase() + ' response (forced zeros: ' + r.type.zeros + ').');
      showMsg(ui.warn, w2, 'warn');
      const bd = r.bands, M = r.M, hx = [], hy = Array.from(r.h); for (let n = 0; n <= M; n++) hx.push(n);
      const hm = Math.max.apply(null, hy.map(Math.abs).concat([1e-9])), series = [];
      if (isWin) series.push({ x: hx, y: Array.from(r.hd), color: T.muted, dash: [4, 3], width: 1.2 });
      series.push({ x: hx, y: hy, color: T.c, stem: true, width: M > 60 ? 1.1 : 1.8 });
      plot(ui.ch, { xmin: -0.7, xmax: M + 0.7, ymin: -hm * 0.45 - 0.02 * hm, ymax: hm * 1.12, xlabel: 'n  (center of symmetry α = ' + g4(M / 2) + ')', ylabel: 'h[n]', series, legend: [{ text: 'h[n]', color: T.c }].concat(isWin ? [{ text: 'ideal hd[n] (delayed)', color: T.muted }] : []) });
      if (isWin) plot(ui.cw, { xmin: -0.7, xmax: M + 0.7, ymin: -0.05, ymax: 1.1, xlabel: 'n', ylabel: 'w[n]', series: [{ x: hx, y: Array.from(r.w), color: T.a, stem: M <= 60, width: 1.6 }], legend: [{ text: DM.WIN[spec.win].name + ' window', color: T.a }] });
      else {
        const np = 400, ax = [], ay = []; for (let i = 0; i <= np; i++) { const w = PI * i / np; ax.push(w / PI); ay.push(cabs(DM.freqAt(r.h, [1], w))); }
        plot(ui.cw, { xmin: 0, xmax: 1, ymin: -0.05, ymax: Math.max(1.15, Math.max.apply(null, ay) * 1.05), xlabel: 'ω/π', ylabel: '|H| (linear)', series: [{ x: ax, y: ay, color: T.c, width: 1.8 }],
          points: r.Ak.map(s => ({ x: s.w / PI, y: s.A, color: s.where === 'trans' ? T.warn : s.where === 'pass' ? T.b : T.pink })), legend: [{ text: 'samples A_k at ω_k = 2πk/L', color: T.b }] });
      }
      const np = 1024, mx = [], dbv = []; for (let i = 0; i < np; i++) { const w = PI * i / (np - 1); mx.push(w / PI); dbv.push(20 * Math.log10(Math.max(cabs(DM.freqAt(r.h, [1], w)), 1e-9))); }
      const floorDb = -Math.min(140, Math.max(90, (isWin ? spec.As : 60) + 45)), rects = [], hl = [];
      if (isWin) {
        const up = 20 * Math.log10(1 + spec.d1), dn = 20 * Math.log10(1 - spec.d1);
        bd.stop.forEach(b => { rects.push({ x0: b[0] / PI, x1: b[1] / PI, y0: -spec.As, y1: 10, color: T.bad, alpha: 0.18 }); hl.push({ y: -spec.As, x0: b[0] / PI, x1: b[1] / PI, color: T.bad }); });
        bd.pass.forEach(b => { rects.push({ x0: b[0] / PI, x1: b[1] / PI, y0: up, y1: 10, color: T.bad, alpha: 0.18 }, { x0: b[0] / PI, x1: b[1] / PI, y0: floorDb, y1: dn, color: T.bad, alpha: 0.18 }); });
      }
      plot(ui.cm, { xmin: 0, xmax: 1, ymin: floorDb, ymax: 6, xlabel: 'ω/π', ylabel: '|H(e^{jω})| (dB)', series: [{ x: mx, y: dbv, color: T.c, width: 2 }], rects, hlines: hl, vlines: edges.map(e => ({ x: e / PI })),
        legend: isWin ? [{ text: '|H|', color: T.c }, { text: 'forbidden region (spec mask)', color: T.bad }] : [{ text: '|H|', color: T.c }] });
      // readouts
      const m = r.meas, ty = r.type, rows = [['order M / taps L', M + ' / ' + (M + 1)], ['FIR type', ty.type ? 'Type ' + roman(ty.type) + ' (' + ty.sym + ', M ' + (M % 2 ? 'odd' : 'even') + ')' : 'not linear phase'],
        ['forced zeros', ty.zeros || '—'], ['group delay', g4(M / 2) + ' samples (constant)'], ['measured As', g4(m.As) + ' dB'], ['passband max |ΔH|', g4(m.passDev) + '  (' + g4(20 * Math.log10(1 + m.passDev)) + ' dB)'], ['passband ripple Rp', g4(m.Rp) + ' dB p-p']];
      if (isWin) {
        const W = DM.WIN[spec.win], ok1 = m.As >= spec.As - 0.05, ok2 = m.passDev <= spec.d1 * 1.02;
        rows.push(['As spec ' + g4(spec.As) + ' dB', ok1 ? 'MET' : 'NOT MET'], ['δ₁ spec ' + g4(spec.d1), ok2 ? 'MET' : 'NOT MET'], ['A = −20log δ', g4(r.A) + ' dB'], ['transition Δω', g4(bd.dw / PI) + 'π rad'],
          ['window', W.name + (W.sidelobe ? ', sidelobe ' + W.sidelobe + ' dB, main lobe ' + W.mlText : ', β = ' + g4(r.beta))]);
      } else rows.push(['frequency samples', r.Ak.length + ' (k = 0…' + (r.Ak.length - 1) + ')'], ['transition value T', g4(spec.T)]);
      hudSet(ui.hud, rows);
      // working
      const L = [], e = edges.map(v => g4(v / PI) + 'π'), pi = x => g4(x / PI) + 'π', cuts = bd.cut;
      L.push('Specification: ' + DM.RESP[resp] + ', band edges ' + e.join(', ') + ' (' + RDEF[resp][1].replace(' (×π)', '') + ').');
      if (isWin) {
        const W = DM.WIN[spec.win], dl = Math.min(spec.d1, r.d2);
        L.push('1. Tolerances: δ₂ = 10^(−As/20) = 10^(−' + g4(spec.As) + '/20) = ' + g4(r.d2) + ';  δ₁ = ' + g4(spec.d1) + ';  the window method gives equal ripple in both bands, so δ = min(δ₁, δ₂) = ' + g4(dl) + '.');
        L.push('   A = −20 log₁₀ δ = ' + g4(r.A) + ' dB.');
        L.push('2. Transition width Δω = ' + (resp === 'bp' || resp === 'bs' ? 'min of the two transition bands = ' : '') + g4(bd.dw) + ' rad = ' + pi(bd.dw) + '.');
        if (spec.win === 'kaiser') {
          L.push('3. Kaiser (O&S eqs. 7.62–7.63): β = ' + (r.A > 50 ? '0.1102 (A − 8.7)' : r.A >= 21 ? '0.5842 (A − 21)^0.4 + 0.07886 (A − 21)' : '0 (A < 21)') + ' = ' + g4(r.beta) + ';  M = (A − 8) / (2.285 Δω) = (' + g4(r.A) + ' − 8) / (2.285 × ' + g4(bd.dw) + ') = ' + g4(DM.kaiserOrderRaw(r.A, bd.dw)) + ' → M = ' + r.Mest + ' (round up).');
        } else {
          L.push('3. ' + W.name + ' window: peak approximation error ≈ ' + W.peakErr + ' dB (' + (W.peakErr <= -r.A ? 'enough for A = ' + g4(r.A) + ' dB' : 'NOT enough for A = ' + g4(r.A) + ' dB') + '), main-lobe width ' + W.mlText + ' = Δω → M ≥ ' + (spec.win === 'rect' ? '4π/Δω − 1' : (spec.win === 'blackman' ? 12 : 8) + 'π/Δω') + ' = ' + g4(spec.win === 'rect' ? 4 * PI / bd.dw - 1 : (spec.win === 'blackman' ? 12 : 8) * PI / bd.dw) + ' → M = ' + r.Mest + (spec.orderMode === 'manual' ? ' (estimate; you chose M = ' + M + ')' : '') + '.');
        }
        if (r.M !== r.Mest) L.push('   Parity: ' + (resp === 'hp' || resp === 'bs' ? 'high-pass/band-stop need H(e^{jπ}) ≠ 0 so M must be even (Type I) → M = ' + r.M + '.' : 'M changed to ' + r.M + '.'));
      } else {
        L.push('1. Frequency sampling (O&S §7.5): L = M + 1 = ' + r.L + ' samples of the desired amplitude on ω_k = 2πk/L = ' + g4(2 / r.L) + 'π·k; only k = 0…' + Math.floor((r.L - 1) / 2) + ' are needed (conjugate symmetry A_{L−k} = A_k).');
        L.push('   Samples: ' + r.Ak.slice(0, 24).map(s => 'k=' + s.k + ' (' + pi(s.w) + ') → ' + g4(s.A) + (s.where === 'trans' ? ' [T]' : '')).join(';  ') + (r.Ak.length > 24 ? ' …' : ''));
        L.push('   Every sample with ω_k strictly between the stopband and passband edges takes the value T = ' + g4(spec.T) + ' (none ⇒ plain 1/0 samples).');
      }
      L.push((isWin ? '4' : '2') + '. Type and delay: M = ' + M + ' is ' + (M % 2 ? 'odd' : 'even') + ', h symmetric ⇒ Type ' + roman(ty.type || 0) + (ty.zeros && ty.zeros !== 'none forced' ? ' (forced zero: ' + ty.zeros + ')' : '') + '. H(e^{jω}) = e^{−jωM/2} A(ω), group delay α = M/2 = ' + g4(M / 2) + ' samples.');
      if (isWin) {
        L.push('5. Cut-off' + (cuts.length > 1 ? 's' : '') + ' (centre of each transition band): ' + cuts.map((c, i) => 'ωc' + (cuts.length > 1 ? i + 1 : '') + ' = ' + pi(c)).join(', ') + '.');
        const sc = (n, c) => 'sin(' + n + '·(n−α))/(π(n−α))';
        const ideal = resp === 'lp' ? 'hd[n] = sin(ωc(n − α)) / (π(n − α)),   hd[α] = ωc/π = ' + g4(cuts[0] / PI)
          : resp === 'hp' ? 'hd[n] = δ[n − α] − sin(ωc(n − α)) / (π(n − α)),   hd[α] = 1 − ωc/π = ' + g4(1 - cuts[0] / PI)
            : resp === 'bp' ? 'hd[n] = [sin(ωc2(n − α)) − sin(ωc1(n − α))] / (π(n − α)),   hd[α] = (ωc2 − ωc1)/π = ' + g4((cuts[1] - cuts[0]) / PI)
              : 'hd[n] = δ[n − α] − [sin(ωc2(n − α)) − sin(ωc1(n − α))] / (π(n − α)),   hd[α] = 1 − (ωc2 − ωc1)/π = ' + g4(1 - (cuts[1] - cuts[0]) / PI);
        L.push('6. Ideal delayed response (α = ' + g4(M / 2) + '):  ' + ideal + (M % 2 ? '  (α is not an integer: no δ[n−α] sample, use the sinc form)' : ''));
        L.push('7. Window: ' + WINEQ[spec.win] + ',  0 ≤ n ≤ M = ' + M + (spec.win === 'kaiser' ? ',  β = ' + g4(r.beta) : '') + ';   h[n] = hd[n]·w[n].');
        const nn = Math.min(M, 4), ex = []; for (let n = 0; n <= nn; n++) ex.push('h[' + n + '] = ' + g4(r.hd[n]) + ' × ' + g4(r.w[n]) + ' = ' + g4(r.h[n])); L.push('   e.g. ' + ex.join(';  ') + (M > nn ? ' …' : ''));
        L.push('   DC gain ΣH = ' + g4(hy.reduce((s, v) => s + v, 0)) + ' (should be ' + (resp === 'lp' || resp === 'bs' ? '≈ 1' : '≈ 0') + ').');
        L.push('8. Check: measured stopband attenuation ' + g4(m.As) + ' dB vs required ' + g4(spec.As) + ' dB;  worst passband deviation ' + g4(m.passDev) + ' vs δ₁ = ' + g4(spec.d1) + '  →  ' + (m.As >= spec.As - 0.05 && m.passDev <= spec.d1 * 1.02 ? 'specification met.' : 'NOT met, increase M or choose a stronger window (Blackman/Kaiser).'));
      } else {
        L.push('3. Impulse response: h[n] = (1/L) [A₀ + 2 Σ_{k=1}^{⌊(L−1)/2⌋} A_k cos(2πk(n − M/2)/L)],  n = 0…M  (L = ' + r.L + ').');
        const nn = Math.min(M, 4), ex = []; for (let n = 0; n <= nn; n++) ex.push('h[' + n + '] = ' + g4(r.h[n])); L.push('   ' + ex.join(', ') + (M > nn ? ' …' : '') + ';   |H(e^{jω_k})| reproduces A_k exactly at the samples (see the dots).');
        L.push('4. Check: stopband attenuation ' + g4(m.As) + ' dB, passband deviation ' + g4(m.passDev) + '.  A larger transition sample T (try 0.3–0.6) trades transition width for stopband attenuation.');
      }
      L.push('Zeros forced by the FIR type: Type I none; Type II z = −1; Type III z = ±1; Type IV z = +1.  Allowed: I any; II LP/BP; III BP; IV HP/BP.');
      ui.work.set(L);
    }
    return {
      update: upd,
      get: () => ({ mode: ui.mode.value, resp: ui.resp.value, ed: ui.edges.value, win: ui.win.value, d1: st.d1.get(), As: st.As.get(), om: ui.om.value, M: st.M.get(), T: st.T.get() }),
      set(o) {
        selOpt(ui.mode, o.mode); selOpt(ui.resp, o.resp); selOpt(ui.win, o.win); selOpt(ui.om, o.om); if (typeof o.ed === 'string' && o.ed.length < 80) ui.edges.value = o.ed;
        setNum(st.d1, o.d1); setNum(st.As, o.As); setNum(st.M, o.M); setNum(st.T, o.T);
      },
    };
  }
  const roman = n => ['—', 'I', 'II', 'III', 'IV'][n] || '—';

  /* ======================= 4. filter structures ======================= */
  // tiny block-diagram toolkit on a 2D context
  function Dg(ctx, T) {
    const d = {
      line(pts, arrow, color, dash) {
        ctx.strokeStyle = color || T.text; ctx.fillStyle = color || T.text; ctx.lineWidth = 1.4; ctx.setLineDash(dash || []); ctx.beginPath(); ctx.moveTo(pts[0][0], pts[0][1]); for (let i = 1; i < pts.length; i++) ctx.lineTo(pts[i][0], pts[i][1]); ctx.stroke(); ctx.setLineDash([]);
        if (arrow) { const p = pts[pts.length - 1], q = pts[pts.length - 2], a = Math.atan2(p[1] - q[1], p[0] - q[0]); ctx.beginPath(); ctx.moveTo(p[0], p[1]); ctx.lineTo(p[0] - 7 * Math.cos(a - 0.4), p[1] - 7 * Math.sin(a - 0.4)); ctx.lineTo(p[0] - 7 * Math.cos(a + 0.4), p[1] - 7 * Math.sin(a + 0.4)); ctx.closePath(); ctx.fill(); }
      },
      dot(x, y) { ctx.fillStyle = T.text; ctx.beginPath(); ctx.arc(x, y, 3, 0, 2 * PI); ctx.fill(); },
      adder(x, y, sign) { ctx.fillStyle = T.bg; ctx.strokeStyle = T.b; ctx.lineWidth = 1.8; ctx.beginPath(); ctx.arc(x, y, 9, 0, 2 * PI); ctx.fill(); ctx.stroke(); ctx.strokeStyle = T.b; ctx.lineWidth = 1.6; ctx.beginPath(); ctx.moveTo(x - 4, y); ctx.lineTo(x + 4, y); if (sign !== '-') { ctx.moveTo(x, y - 4); ctx.lineTo(x, y + 4); } ctx.stroke(); },
      box(x, y, w, h, label, color, lines) {
        ctx.fillStyle = T.bg; ctx.strokeStyle = color || T.c; ctx.lineWidth = 1.8; ctx.fillRect(x - w / 2, y - h / 2, w, h); ctx.strokeRect(x - w / 2, y - h / 2, w, h);
        ctx.fillStyle = color || T.c; ctx.font = '12px ' + T.mono; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
        if (lines) { const n = lines.length; lines.forEach((ln, i) => { ctx.fillStyle = i === 0 ? (color || T.c) : T.text; ctx.font = (i === 0 ? 'bold 12px ' : '11px ') + T.mono; ctx.fillText(ln, x, y - h / 2 + (i + 0.5) * h / n); }); } else ctx.fillText(label, x, y);
      },
      delay(x, y) { d.box(x, y, 34, 20, 'z⁻¹', T.c); },
      gain(x, y, v) { const s = g3(v), w = Math.max(32, 7 * s.length + 10); d.box(x, y, w, 20, s, T.a); },
      txt(x, y, s, color, align) { ctx.fillStyle = color || T.muted; ctx.font = '11px ' + T.mono; ctx.textAlign = align || 'left'; ctx.textBaseline = 'middle'; ctx.fillText(s, x, y); },
    };
    return d;
  }
  const PITCH = 54;
  const rowsH = L => 46 + L * PITCH + 44;
  function drawDF2(d, b, a, oy) {
    const M = b.length - 1, N = a.length - 1, L = Math.max(M, N), y0 = oy + 46, yk = k => y0 + k * PITCH, xl = 120, xgl = 205, xc = 300, xgr = 395, xr = 480, xin = 14, xout = 575;
    if (L > 0) d.line([[xc, y0], [xc, yk(L)]], true);
    if (N >= 1) {
      for (let k = 1; k <= N; k++) d.line([[xc, yk(k)], [xl, yk(k)]], true);
      d.line([[xl, yk(N)], [xl, y0]], true); d.line([[xin, y0], [xl, y0]], true); d.line([[xl, y0], [xc, y0]], true);
    } else d.line([[xin, y0], [xc, y0]], true);
    for (let k = 0; k <= M; k++) d.line([[xc, yk(k)], [xr, yk(k)]], true);
    if (M >= 1) d.line([[xr, yk(M)], [xr, y0]], true);
    d.line([[xr, y0], [xout, y0]], true);
    for (let k = 0; k <= L; k++) d.dot(xc, yk(k));
    for (let k = 0; k < L; k++) d.delay(xc, yk(k) + PITCH / 2);
    for (let k = 1; k <= N; k++) d.gain(xgl, yk(k), -a[k]);
    for (let k = 0; k < N; k++) d.adder(xl, yk(k));
    for (let k = 0; k <= M; k++) d.gain(xgr, yk(k), b[k]);
    for (let k = 0; k < M; k++) d.adder(xr, yk(k));
    d.txt(xin, y0 - 16, 'x[n]'); d.txt(xout - 30, y0 - 16, 'y[n]'); d.txt(xc + 8, y0 - 14, 'w[n]'); d.txt(xc + 8, yk(1) - 12 + (L > 0 ? 0 : 99), 'w[n−1]');
    d.txt(xgl, oy + 18, 'feedback: −a_k', undefined, 'center'); d.txt(xgr, oy + 18, 'feed-forward: b_k', undefined, 'center');
  }
  function drawDF1(d, b, a, oy) {
    const M = b.length - 1, N = a.length - 1, R = Math.max(M, N), y0 = oy + 46, yk = k => y0 + k * PITCH, xin = 14, xX = 70, xBg = 140, xC = 205, xD = 320, xAg = 395, xY = 470, xout = 575;
    d.line([[xin, y0], [xX, y0]], true); if (M >= 1) d.line([[xX, y0], [xX, yk(M)]], true);
    for (let k = 0; k <= M; k++) d.line([[xX, yk(k)], [xC, yk(k)]], true);
    if (M >= 1) d.line([[xC, yk(M)], [xC, y0]], true);
    if (N >= 1) {
      d.line([[xC, y0], [xD, y0]], true); d.line([[xD, yk(N)], [xD, y0]], true); d.line([[xD, y0], [xY, y0]], true); d.line([[xY, y0], [xY, yk(N)]], true);
      for (let k = 1; k <= N; k++) d.line([[xY, yk(k)], [xD, yk(k)]], true);
      d.line([[xY, y0], [xout, y0]], true);
    } else d.line([[xC, y0], [xout, y0]], true);
    for (let k = 0; k <= M; k++) d.dot(xX, yk(k)); if (N >= 1) for (let k = 0; k <= N; k++) d.dot(xY, yk(k));
    for (let k = 0; k < M; k++) d.delay(xX, yk(k) + PITCH / 2);
    for (let k = 0; k < N; k++) d.delay(xY, yk(k) + PITCH / 2);
    for (let k = 0; k <= M; k++) d.gain(xBg, yk(k), b[k]); for (let k = 0; k < M; k++) d.adder(xC, yk(k));
    for (let k = 1; k <= N; k++) d.gain(xAg, yk(k), -a[k]); for (let k = 0; k < N; k++) d.adder(xD, yk(k));
    d.txt(xin, y0 - 16, 'x[n]'); d.txt(xout - 30, y0 - 16, 'y[n]'); d.txt(xC + 6, y0 - 14, 'v[n]'); d.txt(xBg, oy + 18, 'b_k', undefined, 'center'); d.txt(xAg, oy + 18, '−a_k', undefined, 'center');
    if (R === 0) d.txt(xin, y0 + 20, '(pure gain)');
  }
  function drawTDF2(d, b, a, oy) {
    const M = b.length - 1, N = a.length - 1, L = Math.max(M, N), y0 = oy + 46, yk = k => y0 + k * PITCH, xin = 14, xl = 70, xgl = 170, xc = 300, xgr = 430, xr = 520, xout = 585;
    d.line([[xin, y0], [xl, y0]], true); if (M >= 1) d.line([[xl, y0], [xl, yk(M)]], true);
    for (let k = 0; k <= M; k++) d.line([[xl, yk(k)], [xc, yk(k)]], true);
    if (L > 0) d.line([[xc, yk(L)], [xc, y0]], true);
    d.line([[xc, y0], [xout, y0]], true); if (N >= 1) d.line([[xr, y0], [xr, yk(N)]], true);
    for (let k = 1; k <= N; k++) d.line([[xr, yk(k)], [xc, yk(k)]], true);
    for (let k = 0; k <= M; k++) d.dot(xl, yk(k)); if (N >= 1) for (let k = 0; k <= N; k++) d.dot(xr, yk(k));
    for (let k = 0; k < L; k++) d.delay(xc, yk(k) + PITCH / 2);
    for (let k = 0; k <= M; k++) d.gain(xgl, yk(k), b[k]); for (let k = 1; k <= N; k++) d.gain(xgr, yk(k), -a[k]);
    for (let k = 0; k <= L; k++) { const ins = (k <= M ? 1 : 0) + (k >= 1 && k <= N ? 1 : 0) + (k < L ? 1 : 0); if (ins >= 2) d.adder(xc, yk(k)); }
    d.txt(xin, y0 - 16, 'x[n]'); d.txt(xout - 30, y0 - 16, 'y[n]'); d.txt(xgl, oy + 18, 'b_k', undefined, 'center'); d.txt(xgr, oy + 18, '−a_k', undefined, 'center'); d.txt(xc + 24, y0 + PITCH / 2, 's₁[n−1]');
  }
  function secLabel(s) { return ['b: ' + s.b.map(g3).join(' '), 'a: ' + s.a.map(g3).join(' ')]; }
  function drawCascade(d, cas, oy) {
    const y = oy + 50, n = cas.secs.length; let x = 60; d.txt(8, y - 22, 'x[n]'); d.line([[8, y], [x, y]], true);
    if (cas.delay) { d.box(x + 28, y, 56, 30, 'z⁻' + supn(cas.delay), T0.c); d.line([[x + 56, y], [x + 80, y]], true); x += 80; }
    cas.secs.forEach((s, i) => { d.box(x + 70, y, 140, 56, '', T0.c, ['H' + (i + 1) + '(z)'].concat(secLabel(s))); d.line([[x + 140, y], [x + 175, y]], true); x += 175; });
    d.txt(x + 4, y - 22, 'y[n]'); d.txt(60, oy + 98, n + ' SOS in cascade, K = ' + g3(cas.K) + ' in H1. Section 1 as DF-II:');
    drawDF2(d, cas.secs[0].b, cas.secs[0].a, oy + 110);
  }
  function drawParallel(d, par, oy) {
    const rows = []; if (par.direct.length) rows.push({ t: 'Σ B_r z⁻ʳ', l: ['direct FIR', 'B: ' + par.direct.map(g3).join(' ')] });
    par.secs.forEach((s, i) => rows.push({ t: 'H' + (i + 1), l: ['H' + (i + 1) + '(z)'].concat(secLabel(s)) }));
    const P = 78, ym = oy + 20 + (rows.length - 1) * P / 2 + 30, xs = 40, xb = 250, xsum = 480; d.txt(8, ym - 22, 'x[n]'); d.line([[8, ym], [xs, ym]]); d.dot(xs, ym);
    d.line([[xs, oy + 50], [xs, oy + 50 + (rows.length - 1) * P]]);
    rows.forEach((r, i) => { const y = oy + 50 + i * P; d.line([[xs, y], [xb - 80, y]], true); d.box(xb, y, 160, 62, '', T0.c, r.l); d.line([[xb + 80, y], [xsum, y]]); });
    d.line([[xsum, oy + 50], [xsum, oy + 50 + (rows.length - 1) * P]]); d.line([[xsum, ym], [570, ym]], true); d.adder(xsum, ym); d.txt(528, ym - 20, 'y[n]');
  }
  function drawFold(d, h, oy) {
    const M = h.length - 1, ty = DM.firType(h), anti = ty.type >= 3, K = Math.floor((M + 1) / 2), dx = 46, x0 = 40, xk = k => x0 + k * dx, yt = oy + 62, step = 36, ybase = yt + 52, yl = k => ybase + (K - 1 - k) * step, midGain = M % 2 === 0 && ty.type === 1, ybus = ybase + (K - 1) * step + 90;
    d.line([[8, yt], [xk(M) + 6, yt]], true); d.txt(8, yt - 18, 'x[n]');
    for (let k = 0; k <= M; k++) { d.dot(xk(k), yt); d.txt(xk(k), yt - 16, k === 0 ? '' : '−' + k, undefined, 'center'); if (k < M) d.delay(xk(k) + dx / 2, yt); }
    const joins = [];
    for (let k = 0; k < K; k++) {
      d.line([[xk(k), yt], [xk(k), yl(k) - 9]], true); d.line([[xk(M - k), yt], [xk(M - k), yl(k)], [xk(k) + 9, yl(k)]], true);
      d.line([[xk(k), yl(k) + 9], [xk(k), yl(k) + 20]], true); joins.push({ x: xk(k), y: yl(k) + 30, v: h[k] });
    }
    if (midGain) { const xm = xk(M / 2); d.line([[xm, yt], [xm, ybus - 60]], true); joins.push({ x: xm, y: ybus - 50, v: h[M / 2], mid: true }); }
    joins.sort((p, q) => p.x - q.x);
    joins.forEach((j, i) => { d.line([[j.x, j.y + 10], [j.x, ybus - 9 * (i > 0 ? 1 : 0)]], i > 0); if (i === 0) { /* corner */ } });
    d.line([[joins[0].x, ybus], [xk(M) + 50, ybus]], true);
    for (let k = 0; k < K; k++) d.adder(xk(k), yl(k), anti ? '-' : '+');
    joins.forEach((j, i) => { d.gain(j.x, j.y, j.v); if (i > 0) d.adder(j.x, ybus); });
    d.txt(xk(M) + 10, ybus - 18, 'y[n]'); d.txt(8, oy + 12, anti ? 'Antisymmetric: first-adder input from the late tap is subtracted (−), so each adder forms x[n−k] − x[n−M+k]' : 'Symmetric: each adder forms x[n−k] + x[n−M+k], then one multiplier h[k] (half the multipliers)');
  }
  const T0 = { c: '#4cc9f0' }; // overwritten with the live theme before each draw
  function buildStruct(root) {
    const lo = layout(root), ctl = lo.ctl, stage = lo.stage, ui = {}; let strucV = 'df2';
    const chg = () => { upd(); FSP.state.touch(); };
    const fs = FSP.ui.fieldset(ctl, 'System');
    ui.b = textIn(fs, 'b (num)', '1, 2, 3, 2, 1', chg);
    ui.a = textIn(fs, 'a (den)', '1, -1.2, 0.9, -0.4, 0.1', chg);
    fs.appendChild(note('H(z) = Σ b_k z⁻ᵏ / Σ a_k z⁻ᵏ (a₀ is normalised to 1). For an FIR use a = 1.'));
    ui.kind = FSP.ui.select(fs, 'Structure', [['df1', 'Direct form I'], ['df2', 'Direct form II (canonical)'], ['tdf2', 'Transposed DF-II'], ['casc', 'Cascade of 2nd-order sections'], ['par', 'Parallel (partial fractions)'], ['lp', 'Linear-phase FIR (folded)']], 'df2', v => { strucV = v; chg(); });
    fs.appendChild(note('A multiplication by 0 or ±1 is not counted. Each structure is simulated on 120 random samples and compared with the direct recursion.'));
    ui.msg = mkMsg(ctl); ui.warn = E('div', { class: 'msg warn', hidden: '' }); ctl.appendChild(ui.warn);
    ui.cd = wrapCanvas(stage, 'Block diagram of the chosen structure', 300, true); ui.cd.onResize(upd);
    ui.tblw = E('div', { class: 'dsp2-scroll' }); stage.appendChild(ui.tblw);
    ui.work = FSP.ui.working(stage);
    let lastW = 0, lastH = 0;
    function upd() {
      if (!ui.tblw) return;
      const T = theme(), bl = DM.parseList(ui.b.value), al = DM.parseList(ui.a.value);
      if (!bl || !al) { showMsg(ui.msg, 'Coefficients must be numbers separated by commas/spaces (fractions like 1/3 allowed).'); return; }
      const q = DM.normBA(bl, al); if (q.err) { showMsg(ui.msg, q.err); return; }
      if (Math.max(q.M, q.N) > 12) { showMsg(ui.msg, 'Order above 12 is not drawn (M = ' + q.M + ', N = ' + q.N + '). Reduce the order; the cascade form is the practical choice for high orders.'); return; }
      showMsg(ui.msg, '');
      const b = q.b, a = q.a, M = q.M, N = q.N, L = Math.max(M, N), warn = [];
      let x = new Float64Array(120), s = 12345; for (let i = 0; i < x.length; i++) { s = (s * 1103515245 + 12345) % 2147483648; x[i] = s / 2147483648 * 2 - 1; }
      const ref = DM.filter(b, a, x), sc = Math.max(1, Math.max.apply(null, Array.from(ref).map(Math.abs))), R = {};
      const run = (k, fn, cnt) => { try { const y = fn(); R[k] = { err: DM.maxAbsDiff(y, ref) / sc, c: cnt() }; } catch (e) { R[k] = { na: e.message }; } };
      let cas = null, par = null;
      const lpOk = N === 0 && DM.firType(b).type >= 1;
      run('df1', () => DM.simDF1(b, a, x), () => DM.countsBA(b, a, 'df1')); run('df2', () => DM.simDF2(b, a, x), () => DM.countsBA(b, a, 'df2')); run('tdf2', () => DM.simTDF2(b, a, x), () => DM.countsBA(b, a, 'df2'));
      run('casc', () => { cas = DM.cascadeSOS(b, a); return DM.simCascade(cas, x); }, () => { const c = { mult: 0, add: 0, delay: cas.delay }; cas.secs.forEach(z => { const k = DM.countsBA(z.b, z.a, 'df2'); c.mult += k.mult; c.add += k.add; c.delay += k.delay; }); return c; });
      run('par', () => { par = DM.parallelForm(b, a); return DM.simParallel(par, x); }, () => { const c = { mult: 0, add: 0, delay: Math.max(0, par.direct.length - 1) }; par.secs.forEach(z => { const k = DM.countsBA(z.b, z.a, 'df2'); c.mult += k.mult; c.add += k.add; c.delay += k.delay; }); c.mult += par.direct.filter(v => v !== 0 && v !== 1 && v !== -1).length; c.add += par.secs.length + (par.direct.length ? 1 : 0) - 1 + Math.max(0, par.direct.length - 1); return c; });
      if (lpOk) run('lp', () => DM.simFolded(b, x), () => DM.foldedCounts(b)); else R.lp = { na: 'needs an FIR with a symmetric/antisymmetric h[n]' };
      const names = { df1: 'Direct form I', df2: 'Direct form II', tdf2: 'Transposed DF-II', casc: 'Cascade of SOS', par: 'Parallel', lp: 'Linear-phase (folded)' };
      while (ui.tblw.firstChild) ui.tblw.removeChild(ui.tblw.firstChild);
      const tb = E('table', { class: 'dsp2-tbl' }), hd = E('tr'); ['structure', 'multipliers', 'adders', 'delays', 'max |y − y_ref| / scale'].forEach(t => hd.appendChild(E('th', { text: t }))); tb.appendChild(hd);
      Object.keys(names).forEach(k => { const r = R[k], tr = E('tr', { class: k === strucV ? 'sel' : '' }); tr.appendChild(E('td', { text: names[k], style: 'text-align:left' }));
        if (r.na) { const td = E('td', { text: r.na }); td.colSpan = 4; tr.appendChild(td); } else { [r.c.mult, r.c.add, r.c.delay].forEach(v => tr.appendChild(E('td', { text: String(v) }))); tr.appendChild(E('td', { text: r.err.toExponential(1), class: r.err < 1e-9 ? 'dsp2-ok' : 'dsp2-no' })); }
        tb.appendChild(tr); });
      ui.tblw.appendChild(tb);
      // diagram
      const k = strucV; Object.assign(T0, T);
      let W = 620, H = rowsH(L);
      if (k === 'casc' && cas) { W = Math.max(620, 80 + 175 * cas.secs.length + 60); H = 110 + rowsH(2); }
      else if (k === 'par' && par) { H = 50 + Math.max(1, par.secs.length + (par.direct.length ? 1 : 0)) * 78 + 20; }
      else if (k === 'lp') { W = 40 + q.M * 46 + 110; const K2 = Math.floor((q.M + 1) / 2); H = 62 + 52 + (K2 - 1) * 36 + 130; if (!lpOk) { W = 620; H = 120; } }
      if (lastW !== W || lastH !== H) { ui.cd.cv.style.minWidth = W + 'px'; ui.cd.cv.style.height = H + 'px'; lastW = W; lastH = H; }
      const g = ui.cd.prep(), ctx = g.ctx; ctx.clearRect(0, 0, g.w, g.h); ctx.fillStyle = T.bg; ctx.fillRect(0, 0, g.w, g.h);
      const d = Dg(ctx, T);
      const msgs = [];
      try {
        if (k === 'df1') drawDF1(d, b, a, 0); else if (k === 'df2') drawDF2(d, b, a, 0); else if (k === 'tdf2') drawTDF2(d, b, a, 0);
        else if (k === 'casc') { if (cas) drawCascade(d, cas, 0); else msgs.push('Cascade failed: ' + R.casc.na); }
        else if (k === 'par') { if (par) drawParallel(d, par, 0); else msgs.push('Parallel failed: ' + R.par.na); }
        else if (lpOk) drawFold(d, b, 0); else msgs.push('Linear-phase folding needs an FIR (a = 1) whose h[n] is symmetric or antisymmetric.');
      } catch (e) { msgs.push('Could not draw: ' + e.message); }
      if (msgs.length) { ctx.fillStyle = T.warn; ctx.font = '12px ' + T.mono; ctx.textAlign = 'left'; ctx.textBaseline = 'top'; ctx.fillText(msgs[0], 12, 14); }
      if (N > 0 && DM.stability(DM.pzk(b, a).poles).verdict !== 'stable') warn.push('The recursive system is not BIBO stable: errors are normalised to the output scale, but a real fixed-point implementation would overflow.');
      showMsg(ui.warn, warn, 'warn');
      // working
      const W2 = [], cnt = R[k] && R[k].c;
      W2.push('H(z) = (' + polyStr(b) + ') / (' + polyStr(a) + '),   M = ' + M + ', N = ' + N + '.');
      if (k === 'df1') { W2.push('Direct form I: two separate delay lines.'); W2.push('  v[n] = ' + lin(b.map((c, i) => [c, 'x[' + idx(i) + ']']))); W2.push('  y[n] = v[n] ' + (N ? '+ ' + lin(a.map((c, i) => [-c, 'y[' + idx(i) + ']']).slice(1)) : '')); W2.push('  Cost: (M+1) + N = ' + (M + 1 + N) + ' multiplier slots, M + N = ' + (M + N) + ' adders, M + N = ' + (M + N) + ' delays (counting 0/±1 as free: ' + (cnt ? cnt.mult : '—') + ' mult).'); }
      else if (k === 'df2') { W2.push('Direct form II (canonical): one shared delay line w[n] (O&S §6.3.2).'); W2.push('  w[n] = x[n] ' + (N ? '+ ' + lin(a.map((c, i) => [-c, 'w[' + idx(i) + ']']).slice(1)) : '')); W2.push('  y[n] = ' + lin(b.map((c, i) => [c, 'w[' + idx(i) + ']']))); W2.push('  Cost: M + N + 1 = ' + (M + N + 1) + ' multiplier slots, M + N = ' + (M + N) + ' adders, max(M, N) = ' + L + ' delays (' + (cnt ? cnt.mult : '—') + ' non-trivial mult).'); }
      else if (k === 'tdf2') { W2.push('Transposed DF-II: reverse every arrow of DF-II, swap input and output (same H(z), same cost; states are partial sums).'); W2.push('  y[n] = b₀ x[n] + s₁[n−1];   s_k[n] = b_k x[n] − a_k y[n] + s_{k+1}[n−1],  k = 1…' + L + '  (s_{' + (L + 1) + '} = 0)'); W2.push('  Cost: ' + (cnt ? cnt.mult : '—') + ' non-trivial multipliers, ' + (cnt ? cnt.add : '—') + ' adders, ' + L + ' delays.'); }
      else if (k === 'casc' && cas) { W2.push('Cascade: factor numerator and denominator into real second-order factors (conjugate pairs together), H(z) = K z⁻ᴰ Π H_i(z), H_i = (1 + b₁z⁻¹ + b₂z⁻²)/(1 + a₁z⁻¹ + a₂z⁻²). Poles sorted by radius, each zero pair matched to the nearest pole pair to keep intermediate gains moderate.'); W2.push('  K = ' + g4(cas.K) + (cas.delay ? ',  D = ' + cas.delay + ' pure delays' : '')); cas.secs.forEach((z, i) => W2.push('  H' + (i + 1) + '(z) = (' + polyStr(z.b) + ') / (' + polyStr(z.a) + ')')); W2.push('  Each section is a DF-II block (2 delays). Less sensitive to coefficient quantisation than the direct forms.'); }
      else if (k === 'par' && par) { W2.push('Parallel: partial fractions in z⁻¹, conjugate poles combined into real second-order sections, H(z) = Σ B_r z⁻ʳ + Σ H_i(z).'); if (par.direct.length) W2.push('  direct part: ' + polyStr(par.direct)); par.secs.forEach((z, i) => W2.push('  H' + (i + 1) + '(z) = (' + polyStr(z.b) + ') / (' + polyStr(z.a) + ')   poles ' + z.poles.map(cstr).join(', '))); W2.push('  Output = sum of all branches. Errors in one section do not propagate to the others.'); }
      else if (lpOk) { const ty = DM.firType(b); W2.push('Linear-phase FIR, Type ' + roman(ty.type) + ' (' + ty.sym + ', M = ' + M + '): h[n] = ' + (ty.type <= 2 ? '+' : '−') + 'h[M−n].'); W2.push('  y[n] = Σ_{k<' + Math.floor((M + 1) / 2) + '} h[k] (x[n−k] ' + (ty.type <= 2 ? '+' : '−') + ' x[n−' + M + '+k])' + (M % 2 === 0 && ty.type === 1 ? ' + h[' + M / 2 + '] x[n−' + M / 2 + ']' : '')); W2.push('  Folded cost: ' + DM.foldedCounts(b).mult + ' multipliers (direct form: ' + (M + 1) + '), ' + DM.foldedCounts(b).add + ' adders, ' + M + ' delays.'); }
      W2.push('');
      W2.push('Verification: every structure is run on 120 random samples; the table shows max|y − y_ref| / max(1, max|y_ref|) relative to the direct recursion y[n] = Σb_k x[n−k] − Σa_k y[n−k]. A value ≲ 1e-12 confirms the diagram realises the same H(z).');
      ui.work.set(W2);
    }
    return {
      update: upd,
      get: () => ({ b: ui.b.value, a: ui.a.value, k: strucV }),
      set(o) { if (typeof o.b === 'string' && o.b.length < 400) ui.b.value = o.b; if (typeof o.a === 'string' && o.a.length < 400) ui.a.value = o.a; selOpt(ui.kind, o.k); strucV = ui.kind.value; },
    };
  }

  /* ======================= 5. DFT properties ======================= */
  let snrTable = null; // simulated full-scale SNR for B = 2…16 (computed once, lazily)
  function buildDFT(root) {
    const lo = layout(root), ctl = lo.ctl, stage = lo.stage, st = {}, ui = {};
    const chg = () => { upd(); FSP.state.touch(); };
    const fs = FSP.ui.fieldset(ctl, 'Sequences & DFT size');
    ui.x = textIn(fs, 'x[n]', '1, 2, 3, 4', chg); ui.h = textIn(fs, 'h[n]', '1, 1', chg);
    st.N = FSP.ui.slider(fs, { label: 'DFT size N', min: 1, max: 64, step: 1, value: 4, onInput: chg });
    st.m = FSP.ui.slider(fs, { label: 'Shift m', min: -16, max: 16, step: 1, value: 1, onInput: chg });
    fs.appendChild(note('x[n] and h[n] start at n = 0. N ≥ their lengths; N larger than len(x) means zero padding, which samples the same DTFT more finely.'));
    const fq = FSP.ui.fieldset(ctl, 'Quantiser (mid-rise, ±1 full scale)');
    st.B = FSP.ui.slider(fq, { label: 'Bits B', min: 2, max: 16, step: 1, value: 8, onInput: chg });
    st.amp = FSP.ui.slider(fq, { label: 'Amplitude', min: 0.01, max: 1, value: 1, log: true, digits: 3, onInput: chg });
    ui.msg = mkMsg(ctl);
    ui.c1 = wrapCanvas(stage, 'DFT magnitude samples over the DTFT', 240); ui.c1.onResize(upd);
    ui.c2 = wrapCanvas(stage, 'Circular shift', 190); ui.c2.onResize(upd);
    ui.c3 = wrapCanvas(stage, 'Circular versus linear convolution', 200); ui.c3.onResize(upd);
    ui.tw = E('div', { class: 'dsp2-scroll' }); stage.appendChild(ui.tw);
    ui.hud = E('div', { class: 'hud' }); stage.appendChild(ui.hud);
    ui.c4 = wrapCanvas(stage, 'Quantisation SNR versus bits', 220); ui.c4.onResize(upd);
    ui.work = FSP.ui.working(stage);
    function upd() {
      if (!ui.hud) return;
      const T = theme(), xl = DM.parseList(ui.x.value, 64), hl = DM.parseList(ui.h.value, 64);
      if (!xl || !hl) { showMsg(ui.msg, 'x[n] and h[n] must be lists of numbers (commas/spaces, at most 64 values, fractions allowed).'); return; }
      const N = Math.round(st.N.get()), need = Math.max(xl.length, hl.length);
      if (N < need) { showMsg(ui.msg, 'N = ' + N + ' is smaller than the longest sequence (' + need + '). Increase N (or shorten x, h); otherwise the DFT would only see a truncated sequence.'); return; }
      showMsg(ui.msg, '');
      const m = Math.round(st.m.get()), X = DM.dft(xl, N), H = DM.dft(hl, N), xp = DM.padTo(xl, N), hp = DM.padTo(hl, N);
      // spectrum
      const cxs = [], cys = [], npt = 512; for (let i = 0; i <= npt; i++) { const kk = N * i / npt; cxs.push(kk); cys.push(cabs(DM.dtft(xl, TAU * kk / N))); }
      const kx = [], ky = []; for (let k = 0; k < N; k++) { kx.push(k); ky.push(Math.hypot(X.re[k], X.im[k])); }
      const ymx = Math.max.apply(null, cys.concat(ky).concat([1e-9])) * 1.12;
      plot(ui.c1, { xmin: -0.5, xmax: N + 0.5, ymin: 0, ymax: ymx, xlabel: 'k  (ω = 2πk/N; k = N ↔ ω = 2π)', ylabel: '|X|', series: [{ x: cxs, y: cys, color: T.muted, width: 1.4 }, { x: kx, y: ky, color: T.c, stem: true, width: 2 }], legend: [{ text: '|X[k]| (N-point DFT)', color: T.c }, { text: '|X(e^{jω})| (DTFT)', color: T.muted }] });
      // circular shift
      const xs = DM.circShift(xl, m, N), nn = []; for (let n = 0; n < N; n++) nn.push(n);
      const off = a => a.map(v => v);
      const nA = nn.map(v => v - 0.1), nB = nn.map(v => v + 0.1), xmax = Math.max.apply(null, Array.from(xp).concat(Array.from(xs)).map(Math.abs).concat([1e-9])) * 1.2;
      plot(ui.c2, { xmin: -0.7, xmax: N - 0.3, ymin: -xmax, ymax: xmax, xlabel: 'n', ylabel: 'value', series: [{ x: nA, y: off(Array.from(xp)), color: T.c, stem: true, width: 2 }, { x: nB, y: Array.from(xs), color: T.pink, stem: true, width: 2 }], legend: [{ text: 'x[n] (zero-padded to N)', color: T.c }, { text: 'x[(n − ' + m + ')ₙ]', color: T.pink }] });
      // convolution
      const yc = DM.circConv(xl, hl, N), yl = DM.linConv(xl, hl), Lmin = xl.length + hl.length - 1, cxn = [], lxn = [];
      for (let n = 0; n < N; n++) cxn.push(n + 0.12); for (let n = 0; n < yl.length; n++) lxn.push(n - 0.12);
      const ym = Math.max.apply(null, Array.from(yc).concat(Array.from(yl)).map(Math.abs).concat([1e-9])) * 1.2;
      plot(ui.c3, { xmin: -0.7, xmax: Math.max(N, yl.length) - 0.3, ymin: Array.from(yc).concat(Array.from(yl)).some(v => v < 0) ? -ym : -0.05 * ym, ymax: ym, xlabel: 'n', ylabel: 'y[n]', series: [{ x: lxn, y: Array.from(yl), color: T.a, stem: true, width: 2 }, { x: cxn, y: Array.from(yc), color: T.b, stem: true, width: 2 }], legend: [{ text: 'linear x*h (length ' + Lmin + ')', color: T.a }, { text: N + '-point circular x⊛h', color: T.b }] });
      // table
      while (ui.tw.firstChild) ui.tw.removeChild(ui.tw.firstChild);
      const tb = E('table', { class: 'dsp2-tbl' }), hd = E('tr'); ['k', 'x[k]', 'Re X[k]', 'Im X[k]', '|X[k]|', '∠X[k] (°)'].forEach(t => hd.appendChild(E('th', { text: t }))); tb.appendChild(hd);
      for (let k = 0; k < N; k++) { const tr = E('tr'), mag = Math.hypot(X.re[k], X.im[k]); [String(k), g4(xp[k]), g4(X.re[k]), g4(X.im[k]), g4(mag), mag < 1e-9 ? '—' : g4(Math.atan2(X.im[k], X.re[k]) * 180 / PI)].forEach(t => tr.appendChild(E('td', { text: t }))); tb.appendChild(tr); }
      ui.tw.appendChild(tb);
      // properties (numerical verification)
      let eShift = 0, eConv = 0, eRev = 0, eInv = 0;
      const Xs = DM.dft(xs, N); for (let k = 0; k < N; k++) { const a = -TAU * ((((k * m) % N) + N) % N) / N, wr = Math.cos(a), wi = Math.sin(a), pr = X.re[k] * wr - X.im[k] * wi, pi2 = X.re[k] * wi + X.im[k] * wr; eShift = Math.max(eShift, Math.abs(pr - Xs.re[k]), Math.abs(pi2 - Xs.im[k])); }
      const Yr = new Float64Array(N), Yi = new Float64Array(N); for (let k = 0; k < N; k++) { Yr[k] = X.re[k] * H.re[k] - X.im[k] * H.im[k]; Yi[k] = X.re[k] * H.im[k] + X.im[k] * H.re[k]; }
      const yi = DM.idft(Yr, Yi); for (let n = 0; n < N; n++) eConv = Math.max(eConv, Math.abs(yi.re[n] - yc[n]), Math.abs(yi.im[n]));
      const xr = DM.timeRev(xl, N), Xr = DM.dft(xr, N); for (let k = 0; k < N; k++) { const kk = (N - k) % N; eRev = Math.max(eRev, Math.abs(Xr.re[k] - X.re[kk]), Math.abs(Xr.im[k] - X.im[kk])); }
      const xi = DM.idft(X.re, X.im); for (let n = 0; n < N; n++) eInv = Math.max(eInv, Math.abs(xi.re[n] - xp[n]), Math.abs(xi.im[n]));
      let Ex = 0, EX = 0; for (let n = 0; n < N; n++) { Ex += xp[n] * xp[n]; EX += X.re[n] * X.re[n] + X.im[n] * X.im[n]; } EX /= N;
      // quantiser
      const B = Math.round(st.B.get()), amp = st.amp.get(), delta = 2 / Math.pow(2, B), snrSim = DM.quantSNR(B, 65536, amp), snrTh = DM.snrTheory(B) + 20 * Math.log10(amp);
      if (!snrTable) { snrTable = []; for (let b = 2; b <= 16; b++) snrTable.push(DM.quantSNR(b, 65536, 1)); }
      const bx = [], th = []; for (let b = 2; b <= 16; b++) { bx.push(b); th.push(DM.snrTheory(b)); }
      plot(ui.c4, { xmin: 1.5, xmax: 16.5, ymin: 0, ymax: 110, xlabel: 'bits B', ylabel: 'SNR (dB), full-scale sine', series: [{ x: bx, y: th, color: T.muted, width: 1.4, dash: [5, 3] }, { x: bx, y: snrTable, color: T.c, stem: true, width: 1.8 }], points: [{ x: B, y: snrSim, color: T.pink, label: g4(snrSim) + ' dB (A=' + g3(amp) + ')', below: true }], legend: [{ text: 'simulated (A = 1)', color: T.c }, { text: '6.02B + 1.76', color: T.muted }] });
      const okc = v => (v < 1e-9 ? 'OK' : 'ERROR');
      hudSet(ui.hud, [
        ['length x / h / N', xl.length + ' / ' + hl.length + ' / ' + N], ['zero padding', N > xl.length ? (N - xl.length) + ' zeros appended to x' : 'none'],
        ['Parseval Σ|x|² / (1/N)Σ|X|²', g4(Ex) + ' / ' + g4(EX)], ['IDFT(DFT(x)) = x', 'max err ' + eInv.toExponential(1) + ' ' + okc(eInv)],
        ['shift X·W^{km}', 'max err ' + eShift.toExponential(1) + ' ' + okc(eShift)], ['reversal X[(−k)ₙ]', 'max err ' + eRev.toExponential(1) + ' ' + okc(eRev)],
        ['IDFT(X·H) = x⊛h', 'max err ' + eConv.toExponential(1) + ' ' + okc(eConv)], ['circular = linear?', N >= Lmin ? 'yes (N ≥ ' + Lmin + ')' : 'NO, time aliasing (need N ≥ ' + Lmin + ')'],
        ['Δ = 2/2^B', g4(delta)], ['noise power Δ²/12', g4(delta * delta / 12)], ['SNR theory', g4(snrTh) + ' dB'], ['SNR simulated', g4(snrSim) + ' dB'],
      ]);
      const L = [], W = 'e^{−j2π/' + N + '}';
      L.push('DFT definition: X[k] = Σ_{n=0}^{N−1} x[n] W_N^{kn},  W_N = e^{−j2π/N},  N = ' + N + (N > xl.length ? ' (x zero-padded with ' + (N - xl.length) + ' zero' + (N - xl.length > 1 ? 's' : '') + ')' : '') + '.  Samples of the DTFT: X[k] = X(e^{jω}) at ω = 2πk/N.');
      L.push('x = [' + Array.from(xp).map(g4).join(', ') + ']');
      const kshow = Math.min(N, 4);
      for (let k = 0; k < kshow; k++) { const terms = []; for (let n = 0; n < Math.min(N, 4); n++) terms.push(g4(xp[n]) + '·W^' + ((k * n) % N)); L.push('X[' + k + '] = ' + terms.join(' + ') + (N > 4 ? ' + …' : '') + ' = ' + cstr(C(X.re[k], X.im[k])) + '   (|X| = ' + g4(Math.hypot(X.re[k], X.im[k])) + ')'); }
      if (N > kshow) L.push('… full table above (k = 0…' + (N - 1) + '). W_N^{kn} depends only on (kn mod N): W^0 = 1' + (N % 2 === 0 ? ', W^' + N / 2 + ' = −1' : '') + (N % 4 === 0 ? ', W^' + N / 4 + ' = −j' : '') + '.');
      L.push('X[0] = Σ x[n] = ' + g4(Ex === 0 ? 0 : xp.reduce((s, v) => s + v, 0)) + ' (DC);  Parseval: Σ|x[n]|² = ' + g4(Ex) + ' = (1/N)Σ|X[k]|² = ' + g4(EX) + '.');
      L.push('');
      L.push('Circular shift: DFT{x[(n − m)ₙ]} = W_N^{km} X[k].  m = ' + m + ': x[(n−' + m + ')ₙ] = [' + Array.from(xs).map(g4).join(', ') + ']; check against X[k]·e^{−j2πk·' + m + '/' + N + '}: max error ' + eShift.toExponential(1) + '.');
      L.push('Time reversal x[(−n)ₙ] = [' + Array.from(xr).map(g4).join(', ') + '] ↔ X[(−k)ₙ] (error ' + eRev.toExponential(1) + ').');
      L.push('Circular convolution: y[n] = Σ_{k=0}^{N−1} x[k] h[(n−k)ₙ] ↔ Y[k] = X[k]H[k].  y = [' + Array.from(yc).map(g4).join(', ') + ']  (via IDFT(X·H): error ' + eConv.toExponential(1) + ').');
      L.push('Linear convolution x*h has length Lx + Lh − 1 = ' + xl.length + ' + ' + hl.length + ' − 1 = ' + Lmin + ': [' + Array.from(yl).map(g4).join(', ') + ']. ' + (N >= Lmin ? 'N ≥ ' + Lmin + ', so no wrap-around: circular = linear.' : 'N = ' + N + ' < ' + Lmin + ': the last ' + (Lmin - N) + ' outputs wrap onto n = 0…' + (Lmin - N - 1) + ' (time aliasing).'));
      L.push('');
      L.push('Quantisation (B = ' + B + ' bits over ±1, mid-rise): Δ = 2/2^B = 2/' + Math.pow(2, B) + ' = ' + g4(delta) + ';  error uniform in ±Δ/2 → noise power σe² = Δ²/12 = ' + g4(delta * delta / 12) + '.');
      L.push('  Sine of amplitude A = ' + g3(amp) + ': signal power A²/2 = ' + g4(amp * amp / 2) + ';  SNR = 10 log₁₀((A²/2)/(Δ²/12)) = 6.02B + 1.76 + 20 log₁₀A = 6.02×' + B + ' + 1.76 + (' + g4(20 * Math.log10(amp)) + ') = ' + g4(snrTh) + ' dB.  Simulated (65536 samples, incommensurate frequency): ' + g4(snrSim) + ' dB.  Each extra bit adds ≈ 6.02 dB.');
      ui.work.set(L);
    }
    return {
      update: upd,
      get: () => ({ x: ui.x.value, h: ui.h.value, N: st.N.get(), m: st.m.get(), B: st.B.get(), amp: st.amp.get() }),
      set(o) { if (typeof o.x === 'string' && o.x.length < 400) ui.x.value = o.x; if (typeof o.h === 'string' && o.h.length < 400) ui.h.value = o.h; setNum(st.N, o.N); setNum(st.m, o.m); setNum(st.B, o.B); setNum(st.amp, o.amp); },
    };
  }

  /* ======================= tab ======================= */
  const host = { active: false, sub: 'samp' };
  let subs = null; const panels = {}, btns = {};
  function showSub(id, silent) {
    host.sub = id;
    Object.keys(panels).forEach(k => { panels[k].hidden = k !== id; btns[k].setAttribute('aria-pressed', String(k === id)); btns[k].classList.toggle('active', k === id); });
    if (host.active) subs[id].update();
    if (!silent) FSP.state.touch();
  }
  FSP.registerTab({
    id: 'dspdesign', title: 'z-Transform & FIR Design',
    init(panel) {
      setStyle();
      const bar = E('div', { class: 'seg row', role: 'group', 'aria-label': 'z-transform and FIR design topic' }); panel.appendChild(bar);
      const defs = [['samp', 'Sampling & aliasing'], ['z', 'Difference eq. & z-transform'], ['fir', 'FIR design'], ['struct', 'Filter structures'], ['dft', 'DFT properties']];
      defs.forEach(d => { btns[d[0]] = E('button', { type: 'button', class: 'seg-btn', 'aria-pressed': 'false', text: d[1], onclick: () => showSub(d[0]) }); bar.appendChild(btns[d[0]]); });
      defs.forEach(d => { panels[d[0]] = E('div', { hidden: '' }); panel.appendChild(panels[d[0]]); });
      subs = { samp: buildSamp(panels.samp), z: buildZ(panels.z), fir: buildFIR(panels.fir), struct: buildStruct(panels.struct), dft: buildDFT(panels.dft) };
      const pre = { samp: 'sp_', z: 'zt_', fir: 'fd_', struct: 'sx_', dft: 'df_' };
      FSP.state.bind('dspdesign', {
        get() { const o = { sub: host.sub }; Object.keys(subs).forEach(k => { const g = subs[k].get(); Object.keys(g).forEach(q => { o[pre[k] + q] = g[q]; }); }); return o; },
        set(o) {
          o = o || {};
          Object.keys(subs).forEach(k => { const r = {}; Object.keys(o).forEach(q => { if (q.indexOf(pre[k]) === 0) r[q.slice(pre[k].length)] = o[q]; }); try { subs[k].set(r); } catch (e) { /* ignore bad state */ } });
          showSub(o.sub && panels[o.sub] ? o.sub : 'samp', true);
        },
      });
      showSub(host.sub, true);
    },
    activate() { host.active = true; subs[host.sub].update(); },
    deactivate() { host.active = false; },
  });
})();
