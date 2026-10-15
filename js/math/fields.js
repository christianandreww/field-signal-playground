/* FSP module: Fields & Polarization (EE3101 Engineering Electromagnetics).
   Sub-topics: Coulomb charge distributions, Biot–Savart / Ampère, Faraday & displacement current,
   plane waves & polarization, normal incidence & standing waves.
   Notation follows Hayt: phasors with e^{jωt}, waves e^{-γz}, γ = α + jβ, η intrinsic impedance.
   Pure math lives in FSP.math.fields (no DOM). All DOM work happens inside init/activate. */
(function () {
  'use strict';
  const PI = Math.PI;
  const EPS0 = 8.8541878128e-12, MU0 = 1.25663706212e-6, C0 = 299792458;
  const KE = 1 / (4 * PI * EPS0);            // 1/(4πε0) ≈ 8.98755e9 m/F
  const ETA0 = Math.sqrt(MU0 / EPS0);        // ≈ 376.73 Ω
  const fin = Number.isFinite;
  const D2R = PI / 180, R2D = 180 / PI;

  /* ======================= complex helpers ======================= */
  const cx = (re, im) => ({ re: re, im: im || 0 });
  const cadd = (a, b) => cx(a.re + b.re, a.im + b.im);
  const csub = (a, b) => cx(a.re - b.re, a.im - b.im);
  const cmul = (a, b) => cx(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
  const cdiv = (a, b) => { const d = b.re * b.re + b.im * b.im; return cx((a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d); };
  const cscale = (a, s) => cx(a.re * s, a.im * s);
  const cabs = a => Math.hypot(a.re, a.im);
  const carg = a => Math.atan2(a.im, a.re);
  const cexp = a => { const e = Math.exp(a.re); return cx(e * Math.cos(a.im), e * Math.sin(a.im)); };
  function csqrt(a) { // principal branch (Re >= 0)
    const r = cabs(a); if (r === 0) return cx(0, 0);
    if (a.re >= 0) { const t = Math.sqrt(0.5 * (r + a.re)); return cx(t, a.im / (2 * t)); }
    const t = Math.sqrt(0.5 * (r - a.re)), s = a.im < 0 ? -t : t; return cx(a.im / (2 * s), s);
  }
  function ctanh(a) {
    if (a.re > 20) return cx(1, 0); if (a.re < -20) return cx(-1, 0);
    const d = Math.cosh(2 * a.re) + Math.cos(2 * a.im); return cx(Math.sinh(2 * a.re) / d, Math.sin(2 * a.im) / d);
  }
  const ccosh = a => cx(Math.cosh(a.re) * Math.cos(a.im), Math.sinh(a.re) * Math.sin(a.im));
  const csinh = a => cx(Math.sinh(a.re) * Math.cos(a.im), Math.cosh(a.re) * Math.sin(a.im));
  const cfin = a => fin(a.re) && fin(a.im);

  /* ======================= numerical integration ======================= */
  // Vector adaptive Simpson: f(x) -> array of m numbers. Returns {v, evals, capped}.
  function asimpV(f, a, b, opt) {
    opt = opt || {}; const tol = opt.tol || 1e-9, maxDepth = opt.depth || 40, cap = opt.cap || 400000;
    let evals = 0, capped = false; const bud = opt.budget || null;   // optional shared budget {left} across nested calls
    const F = x => { evals++; if (bud) bud.left--; return f(x); };
    const m = (a + b) / 2, fa = F(a), fm = F(m), fb = F(b), n = fa.length;
    const simp = (h, f0, f1, f2) => { const r = new Array(n); for (let i = 0; i < n; i++) r[i] = h / 6 * (f0[i] + 4 * f1[i] + f2[i]); return r; };
    const whole = simp(b - a, fa, fm, fb);
    // absolute tolerance from a coarse 5-point estimate (scale of the answer)
    const q1 = F(a + (b - a) / 4), q3 = F(a + 3 * (b - a) / 4);
    let scale = 0; const coarse = new Array(n);
    for (let i = 0; i < n; i++) { coarse[i] = (b - a) / 12 * (fa[i] + 4 * q1[i] + 2 * fm[i] + 4 * q3[i] + fb[i]); scale = Math.max(scale, Math.abs(coarse[i]), Math.abs(whole[i])); }
    const eps0 = tol * Math.max(scale, opt.floor || 1e-300);
    function rec(a, b, fa, fm, fb, wh, eps, depth) {
      const m = (a + b) / 2, lm = (a + m) / 2, rm = (m + b) / 2, fl = F(lm), fr = F(rm);
      const L = simp(m - a, fa, fl, fm), R = simp(b - m, fm, fr, fb);
      let err = 0; for (let i = 0; i < n; i++) err = Math.max(err, Math.abs(L[i] + R[i] - wh[i]));
      const out_ = evals > cap || (bud && bud.left <= 0);
      if (depth <= 0 || err <= 15 * eps || out_) {
        if (out_) { capped = true; if (bud) bud.capped = true; }
        const out = new Array(n); for (let i = 0; i < n; i++) out[i] = L[i] + R[i] + (L[i] + R[i] - wh[i]) / 15; return out;
      }
      const A = rec(a, m, fa, fl, fm, L, eps / 2, depth - 1), B = rec(m, b, fm, fr, fb, R, eps / 2, depth - 1);
      for (let i = 0; i < n; i++) A[i] += B[i]; return A;
    }
    // split the interval at the quarter points first so a narrow peak is not missed
    const pts = [a, a + (b - a) / 4, m, a + 3 * (b - a) / 4, b], fv = [fa, q1, fm, q3, fb], out = new Array(n).fill(0);
    for (let s = 0; s < 4; s++) {
      const x0 = pts[s], x1 = pts[s + 1], xm = (x0 + x1) / 2, f0 = fv[s], f1 = F(xm), f2 = fv[s + 1];
      const r = rec(x0, x1, f0, f1, f2, simp(x1 - x0, f0, f1, f2), eps0 / 4, maxDepth);
      for (let i = 0; i < n; i++) out[i] += r[i];
    }
    return { v: out, evals, capped };
  }
  // Gauss–Legendre nodes/weights on [-1,1] by Newton iteration on P_n
  function glRule(n) {
    const x = new Float64Array(n), w = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      let z = Math.cos(PI * (i + 0.75) / (n + 0.5)), pp = 1;
      for (let it = 0; it < 100; it++) {
        let p1 = 1, p2 = 0; for (let j = 1; j <= n; j++) { const p3 = p2; p2 = p1; p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j; }
        pp = n * (z * p1 - p2) / (z * z - 1); const z1 = z; z = z1 - p1 / pp; if (Math.abs(z - z1) < 1e-15) break;
      }
      x[i] = z; w[i] = 2 / ((1 - z * z) * pp * pp);
    }
    return { x, w };
  }
  const GL16 = glRule(16);
  // composite Gauss–Legendre (16-point per panel), vector integrand
  function glComp(f, a, b, panels) {
    let out = null; const h = (b - a) / panels;
    for (let p = 0; p < panels; p++) {
      const c = a + (p + 0.5) * h, hh = h / 2;
      for (let i = 0; i < 16; i++) { const v = f(c + hh * GL16.x[i]), wgt = GL16.w[i] * hh; if (!out) out = new Array(v.length).fill(0); for (let k = 0; k < v.length; k++) out[k] += wgt * v[k]; }
    }
    return out;
  }

  /* ======================= 1. Coulomb: charge distributions ======================= */
  // P = [x,y,z] in m; charges in C (C/m, C/m²). Returns {E:[Ex,Ey,Ez] V/m, V (volts), ...}
  const vlen = v => Math.hypot(v[0], v[1], v[2]);
  function pointCharges(qs, P) {
    const E = [0, 0, 0]; let V = 0; const parts = [];
    for (const c of qs) {
      const R = [P[0] - c.x, P[1] - c.y, P[2] - c.z], r = vlen(R);
      if (!(r > 1e-12)) return { error: 'field point coincides with a point charge (|R| = 0): E and V are infinite there' };
      const k = KE * c.q / (r * r * r), Ei = [k * R[0], k * R[1], k * R[2]], Vi = KE * c.q / r;
      E[0] += Ei[0]; E[1] += Ei[1]; E[2] += Ei[2]; V += Vi; parts.push({ R, r, Ei, Vi });
    }
    return { E, V, parts };
  }
  // in-plane field for maps (no error object; returns NaN at a charge)
  function pointEV(qs, x, y, z) {
    let ex = 0, ey = 0, ez = 0, V = 0;
    for (let i = 0; i < qs.length; i++) { const c = qs[i], dx = x - c.x, dy = y - c.y, dz = z - c.z, r2 = dx * dx + dy * dy + dz * dz, r = Math.sqrt(r2), k = KE * c.q / (r2 * r); ex += k * dx; ey += k * dy; ez += k * dz; V += KE * c.q / r; }
    return [ex, ey, ez, V];
  }
  // finite line charge on the z-axis from z1 to z2 (z2 > z1), closed form
  function lineFinite(rhoL, z1, z2, P) {
    const rho = Math.hypot(P[0], P[1]), z = P[2];
    if (!(z2 > z1)) return { error: 'line charge needs z2 > z1' };
    if (rho < 1e-12) {
      if (z >= z1 - 1e-12 && z <= z2 + 1e-12) return { error: 'field point lies on the line charge (ρ = 0): E and V are infinite there' };
      // on the axis outside the segment: |Ez| = kρL (1/near − 1/far), pointing away from the line for ρL > 0
      const u1 = Math.abs(z - z1), u2 = Math.abs(z - z2), near = Math.min(u1, u2), far = Math.max(u1, u2), mag = KE * rhoL * (1 / near - 1 / far), s = z > z2 ? 1 : -1;
      return { E: [0, 0, s * mag], V: KE * rhoL * Math.log(far / near), rho, z, Erho: 0, Ez: s * mag, sinA1: -s, sinA2: -s, R1: u1, R2: u2 };
    }
    const R1 = Math.hypot(rho, z - z1), R2 = Math.hypot(rho, z - z2);
    const sinA1 = (z1 - z) / R1, sinA2 = (z2 - z) / R2;           // angles measured from the perpendicular at P, +ve toward +z
    const Erho = KE * rhoL / rho * (sinA2 - sinA1), Ez = KE * rhoL * (1 / R2 - 1 / R1);
    const V = KE * rhoL * (Math.asinh((z2 - z) / rho) - Math.asinh((z1 - z) / rho));
    const ux = P[0] / rho, uy = P[1] / rho;
    return { E: [Erho * ux, Erho * uy, Ez], V, rho, z, Erho, Ez, sinA1, sinA2, R1, R2 };
  }
  // generic numeric Coulomb integrand: element dQ at r' -> [dEx,dEy,dEz,dV] per unit dQ
  function coulombKernel(P, sx, sy, sz, dq) {
    const Rx = P[0] - sx, Ry = P[1] - sy, Rz = P[2] - sz, r2 = Rx * Rx + Ry * Ry + Rz * Rz, r = Math.sqrt(r2), k = KE * dq / (r2 * r);
    return [k * Rx, k * Ry, k * Rz, KE * dq / r];
  }
  function lineFiniteNum(rhoL, z1, z2, P, tol, budget) {
    const r = asimpV(zp => coulombKernel(P, 0, 0, zp, rhoL), z1, z2, { tol: tol || 1e-10, budget });
    return { E: r.v.slice(0, 3), V: r.v[3], evals: r.evals, capped: r.capped };
  }
  // infinite line on z-axis: Eρ = ρL/(2πε0ρ); V referenced to ρref (V(ρref) = 0)
  function lineInfinite(rhoL, P, rhoRef) {
    const rho = Math.hypot(P[0], P[1]); if (!(rho > 1e-12)) return { error: 'field point lies on the line charge (ρ = 0)' };
    const Erho = rhoL / (2 * PI * EPS0 * rho); rhoRef = rhoRef || 1;
    return { E: [Erho * P[0] / rho, Erho * P[1] / rho, 0], V: rhoL / (2 * PI * EPS0) * Math.log(rhoRef / rho), Erho, rho };
  }
  // infinite line, numeric: z' = z + ρ tanθ maps (−∞,∞) onto (−π/2, π/2); V is not integrable (needs a reference) so only E
  function lineInfiniteNum(rhoL, P, budget) {
    const rho = Math.hypot(P[0], P[1]), lim = PI / 2 * (1 - 1e-12);
    const r = asimpV(th => { const zp = P[2] + rho * Math.tan(th), jac = rho / (Math.cos(th) * Math.cos(th)); const k = coulombKernel(P, 0, 0, zp, rhoL); return [k[0] * jac, k[1] * jac, k[2] * jac]; }, -lim, lim, { tol: 1e-11, budget });
    return { E: r.v, evals: r.evals };
  }
  // ring of radius a in z = 0, centred on the z-axis
  function ringAxis(rhoL, a, z) {
    const d = Math.pow(a * a + z * z, 1.5);
    return { Ez: rhoL * a * z / (2 * EPS0 * d), V: rhoL * a / (2 * EPS0 * Math.sqrt(a * a + z * z)) };
  }
  function ringNum(rhoL, a, P, tol, budget) {
    const r = asimpV(ph => coulombKernel(P, a * Math.cos(ph), a * Math.sin(ph), 0, rhoL * a), 0, 2 * PI, { tol: tol || 1e-10, budget });
    return { E: r.v.slice(0, 3), V: r.v[3], evals: r.evals, capped: r.capped };
  }
  // disk of radius a in z = 0
  function diskAxis(rhoS, a, z) {
    const s = z > 0 ? 1 : z < 0 ? -1 : 0, R = Math.sqrt(a * a + z * z);
    return { Ez: rhoS / (2 * EPS0) * (s - z / R), V: rhoS / (2 * EPS0) * (R - Math.abs(z)) };
  }
  function diskNum(rhoS, a, P, tol, budget) {
    tol = tol || 1e-9;
    let ev = 0, cap = false;
    const r = asimpV(rp => { const q = asimpV(ph => coulombKernel(P, rp * Math.cos(ph), rp * Math.sin(ph), 0, rhoS * rp), 0, 2 * PI, { tol: tol * 0.03, cap: 60000, budget }); ev += q.evals; cap = cap || q.capped; return q.v; }, 0, a, { tol: tol, cap: 4000, budget });
    return { E: r.v.slice(0, 3), V: r.v[3], evals: ev, capped: cap || r.capped };
  }
  // infinite sheet z = 0: E = ρs/(2ε0) sgn(z) az; V referenced to the sheet (V = 0 at z = 0)
  function sheetInf(rhoS, z) { const s = z > 0 ? 1 : z < 0 ? -1 : 0; return { Ez: rhoS / (2 * EPS0) * s, V: -rhoS * Math.abs(z) / (2 * EPS0) }; }
  // numeric: by translational symmetry put the origin under P; ρ' = |z| tanθ, θ ∈ [0, π/2), φ' ∈ [0, 2π)
  function sheetInfNum(rhoS, z, tol, budget) {
    tol = tol || 1e-9;
    const az = Math.abs(z), P = [0, 0, z], lim = PI / 2 * (1 - 1e-12);
    const r = asimpV(th => { const rp = az * Math.tan(th), jac = az / (Math.cos(th) * Math.cos(th)); const q = asimpV(ph => coulombKernel(P, rp * Math.cos(ph), rp * Math.sin(ph), 0, rhoS * rp), 0, 2 * PI, { tol: tol * 0.03, cap: 20000, budget }); return [q.v[0] * jac, q.v[1] * jac, q.v[2] * jac]; }, 0, lim, { tol: tol, cap: 3000, budget });
    return { E: r.v };
  }
  // rectangular sheet x1<x'<x2, y1<y'<y2 in z = 0, uniform ρs: exact closed form (corner double-difference)
  function lnPlus(t, R, a2) { return t >= 0 ? Math.log(t + R) : Math.log(a2) - Math.log(R - t); } // ln(t + R), stable for t < 0
  function rectSheet(rhoS, x1, x2, y1, y2, P) {
    if (!(x2 > x1 && y2 > y1)) return { error: 'sheet needs x2 > x1 and y2 > y1' };
    const z = P[2]; if (Math.abs(z) < 1e-12) return { error: 'field point lies in the plane of the sheet (z = 0): Ez is discontinuous there; choose z ≠ 0' };
    let Ex = 0, Ey = 0, Ez = 0, V = 0; const corners = [];
    [[x2, 1], [x1, -1]].forEach(([xp, sx]) => [[y2, 1], [y1, -1]].forEach(([yp, sy]) => {
      const X = xp - P[0], Y = yp - P[1], R = Math.sqrt(X * X + Y * Y + z * z), s = sx * sy;
      const lY = lnPlus(Y, R, X * X + z * z), lX = lnPlus(X, R, Y * Y + z * z), at = Math.atan(X * Y / (z * R));
      Ex += s * lY; Ey += s * lX; Ez += s * at; V += s * (X * lY + Y * lX - z * at);
      corners.push({ X, Y, R, s, at });
    }));
    const k = KE * rhoS;
    return { E: [k * Ex, k * Ey, k * Ez], V: k * V, corners };
  }
  function rectSheetNum(rhoS, x1, x2, y1, y2, P, tol, budget) {
    tol = tol || 1e-9;
    let ev = 0, cap = false;
    const r = asimpV(xp => { const q = asimpV(yp => coulombKernel(P, xp, yp, 0, rhoS), y1, y2, { tol: tol * 0.03, cap: 60000, budget }); ev += q.evals; cap = cap || q.capped; return q.v; }, x1, x2, { tol: tol, cap: 4000, budget });
    return { E: r.v.slice(0, 3), V: r.v[3], evals: ev, capped: cap || r.capped };
  }
  // tutorial special case: square of side L centred under the z-axis, field on the axis
  function squareSheetAxis(rhoS, L, z) {
    const b = L / 2, R = Math.sqrt(2 * b * b + z * z), quad = Math.atan(b * b / (z * R));
    return { Ez: rhoS / (PI * EPS0) * quad, quad, R };     // = 4 quadrants × (ρs/4πε0)·atan(b²/(zR))
  }

  /* ======================= 2. Magnetostatics ======================= */
  const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  // finite filament on z-axis z1..z2, current I in +az. Hayt: H = I/(4πρ)(sin α2 − sin α1) aφ,
  // α measured at P from the perpendicular to the wire, positive toward +z (sin αi = (zi − z)/Ri).
  function filamentFinite(I, z1, z2, P) {
    const rho = Math.hypot(P[0], P[1]); if (!(rho > 1e-12)) return { error: 'field point lies on the filament axis (ρ = 0)' };
    if (!(z2 > z1)) return { error: 'filament needs z2 > z1' };
    const z = P[2], R1 = Math.hypot(rho, z1 - z), R2 = Math.hypot(rho, z2 - z), s1 = (z1 - z) / R1, s2 = (z2 - z) / R2;
    const Hphi = I / (4 * PI * rho) * (s2 - s1), ux = -P[1] / rho, uy = P[0] / rho;
    return { H: [Hphi * ux, Hphi * uy, 0], Hphi, rho, sinA1: s1, sinA2: s2, a1: Math.asin(s1) * R2D, a2: Math.asin(s2) * R2D, R1, R2 };
  }
  // numeric infinite filament: z' = z + ρ tanθ, Gauss–Legendre over θ ∈ (−π/2, π/2)
  function filamentInfNum(I, P) { const rho = Math.hypot(P[0], P[1]); return glComp(th => { const c = Math.cos(th), jac = rho / (c * c); return bsKernel(I, P, [0, 0, P[2] + rho * Math.tan(th)], [0, 0, jac]); }, -PI / 2, PI / 2, 8); }
  function filamentInf(I, P) { const rho = Math.hypot(P[0], P[1]); if (!(rho > 1e-12)) return { error: 'field point on the filament' }; const Hphi = I / (2 * PI * rho); return { H: [-Hphi * P[1] / rho, Hphi * P[0] / rho, 0], Hphi, rho }; }
  // numeric Biot–Savart, straight segment A→B: H = (I/4π) ∫ dL × R / |R|³
  function bsKernel(I, P, s, dl) { const R = [P[0] - s[0], P[1] - s[1], P[2] - s[2]], r = vlen(R), c = cross(dl, R), k = I / (4 * PI * r * r * r); return [k * c[0], k * c[1], k * c[2]]; }
  function bsSegment(I, A, B, P, panels) {
    const d = [B[0] - A[0], B[1] - A[1], B[2] - A[2]];
    return glComp(t => bsKernel(I, P, [A[0] + d[0] * t, A[1] + d[1] * t, A[2] + d[2] * t], d), 0, 1, panels || 64);
  }
  function bsPolygon(I, pts, P, panels) { const H = [0, 0, 0]; for (let i = 0; i < pts.length; i++) { const h = bsSegment(I, pts[i], pts[(i + 1) % pts.length], P, panels); H[0] += h[0]; H[1] += h[1]; H[2] += h[2]; } return H; }
  // circular loop radius a in z = 0, current counter-clockwise seen from +z (so H along +az at the centre)
  function bsCircle(I, a, P, panels, z0) {
    z0 = z0 || 0;
    return glComp(ph => bsKernel(I, P, [a * Math.cos(ph), a * Math.sin(ph), z0], [-a * Math.sin(ph), a * Math.cos(ph), 0]), 0, 2 * PI, panels || 32);
  }
  function loopAxis(I, a, z) { return { Hz: I * a * a / (2 * Math.pow(a * a + z * z, 1.5)) }; }
  // square loop of side w in z = 0 centred on the axis, CCW from +z
  function squareAxis(I, w, z) { const b = w / 2; return { Hz: 2 * I * b * b / (PI * (b * b + z * z) * Math.sqrt(2 * b * b + z * z)), b }; }
  const squarePts = w => { const b = w / 2; return [[b, -b, 0], [b, b, 0], [-b, b, 0], [-b, -b, 0]]; }; // CCW from +z
  // finite solenoid: N turns, length L centred at z = 0, radius a, on axis. Hz = (nI/2)(sin α2 − sin α1), same α convention as the filament
  function solenoidAxis(I, N, L, a, z) {
    const n = N / L, z1 = -L / 2, z2 = L / 2, s1 = (z1 - z) / Math.hypot(a, z1 - z), s2 = (z2 - z) / Math.hypot(a, z2 - z);
    return { Hz: n * I / 2 * (s2 - s1), n, sinA1: s1, sinA2: s2 };
  }
  // numeric: current sheet K = nI aφ on the cylinder, 2-D Gauss–Legendre over (z', φ')
  function solenoidNum(I, N, L, a, P, pz, pp) {
    const K = N * I / L; pz = pz || 48; pp = pp || 12;
    return glComp(zp => glComp(ph => bsKernel(K, P, [a * Math.cos(ph), a * Math.sin(ph), zp], [-a * Math.sin(ph), a * Math.cos(ph), 0]), 0, 2 * PI, pp), -L / 2, L / 2, pz);
  }
  // coaxial cable (Hayt §7.2): inner radius a carries +I (uniform), outer conductor b..c carries −I
  function coaxH(I, a, b, c, rho) {
    if (!(a > 0 && b > a && c > b)) return { error: 'coax needs 0 < a < b < c' };
    if (rho < 0) return { error: 'ρ must be ≥ 0' };
    let H, region, Ienc;
    if (rho < a) { Ienc = I * rho * rho / (a * a); region = 1; }
    else if (rho < b) { Ienc = I; region = 2; }
    else if (rho < c) { Ienc = I * (c * c - rho * rho) / (c * c - b * b); region = 3; }
    else { Ienc = 0; region = 4; }
    H = rho > 0 ? Ienc / (2 * PI * rho) : 0;
    return { Hphi: H, Ienc, region };
  }
  // numeric enclosed current: ∫ J·dS over 0..ρ using the (uniform) current densities
  function coaxIencNum(I, a, b, c, rho) {
    const J1 = I / (PI * a * a), J3 = -I / (PI * (c * c - b * b));
    const J = r => (r < a ? J1 : r < b ? 0 : r < c ? J3 : 0);
    const segs = [0, a, b, c].filter(x => x < rho).concat([rho]); let s = 0;
    for (let i = 0; i + 1 < segs.length; i++) s += glComp(r => [J(r) * 2 * PI * r], segs[i], segs[i + 1], 4)[0];
    return s;
  }
  // infinite current sheet K ax in z = 0: H = ½ K × aN  ->  Hy = −K/2 sgn(z)
  function sheetK(K, z) { const s = z > 0 ? 1 : z < 0 ? -1 : 0; return { Hy: -K / 2 * s }; }
  // numeric: sum of infinite filaments K dy' at y', each H = K dy'/(2π d) aφ; y' = |z| tanθ
  function sheetKNum(K, z) {
    const az = Math.abs(z), lim = PI / 2 * (1 - 1e-12);
    const r = asimpV(th => { const yp = az * Math.tan(th), jac = az / (Math.cos(th) * Math.cos(th)), d2 = yp * yp + z * z; // filament along +x at (y',0): H = K/(2π d²) (ax × R)
      const Ry = -yp, Rz = z; return [K * jac / (2 * PI * d2) * (-Rz), K * jac / (2 * PI * d2) * Ry]; }, -lim, lim, { tol: 1e-11 });
    return { Hy: r.v[0], Hz: r.v[1] };
  }

  /* ======================= 3. Faraday & displacement current ======================= */
  // (a) fixed N-turn loop of area A, B = B0 sin(ωt) normal to it: Φ = B0 A sin ωt, emf = −N dΦ/dt
  function faradayFixed(N, A, B0, f, t) { const w = 2 * PI * f; return { flux: B0 * A * Math.sin(w * t), emf: -N * B0 * A * w * Math.cos(w * t), peak: N * B0 * A * w }; }
  // (b) sliding bar on rails: separation l, bar at x(t) = x0 + v t, B uniform along +az, load R
  function slidingBar(B, l, v, R, x0, t) {
    const x = x0 + v * t, flux = B * l * x, emf = -B * l * v;   // emf around the circuit taken CCW (+az by the right-hand rule)
    const I = R > 0 ? Math.abs(emf) / R : NaN, F = B * I * l;
    return { x, flux, emf, mag: Math.abs(emf), I, F, Pmech: F * Math.abs(v), Pelec: R > 0 ? emf * emf / R : NaN };
  }
  // (c) N-turn loop of area A rotating at ω in uniform B; normal at angle ωt to B: Φ = B A cos ωt, emf = N B A ω sin ωt
  function rotatingLoop(N, A, B, f, t) { const w = 2 * PI * f; return { flux: B * A * Math.cos(w * t), emf: N * B * A * w * Math.sin(w * t), peak: N * B * A * w }; }
  // conduction vs displacement current, loss tangent σ/(ωε)
  const LT_GOOD_COND = 100, LT_GOOD_DIEL = 0.01;
  function classifyLT(lt) { return lt === 0 ? 'lossless' : lt <= LT_GOOD_DIEL ? 'good dielectric (low-loss)' : lt >= LT_GOOD_COND ? 'good conductor' : 'quasi-conductor (lossy)'; }
  function currentRatio(sigma, er, f, E0) {
    if (!(sigma >= 0 && er > 0 && f > 0)) return { error: 'need σ ≥ 0, εr > 0, f > 0' };
    const w = 2 * PI * f, eps = er * EPS0, lt = sigma / (w * eps);
    return { Jc: sigma * E0, Jd: w * eps * E0, lt, cls: classifyLT(lt), fCross: sigma > 0 ? sigma / (2 * PI * eps) : 0, omega: w };
  }

  /* ======================= 4. Plane waves & polarization ======================= */
  function medium(er, mr, sigma, f) {
    if (!(er > 0 && mr > 0 && sigma >= 0 && f > 0)) return { error: 'need εr > 0, μr > 0, σ ≥ 0, f > 0' };
    const w = 2 * PI * f, eps = er * EPS0, mu = mr * MU0;
    const jwmu = cx(0, w * mu), ysh = cx(sigma, w * eps);
    const gamma = csqrt(cmul(jwmu, ysh)), eta = csqrt(cdiv(jwmu, ysh));
    const alpha = gamma.re, beta = gamma.im, lt = sigma / (w * eps), cls = classifyLT(lt);
    const m = { er, mr, sigma, f, omega: w, eps, mu, gamma, alpha, beta, eta, lt, cls, lambda: 2 * PI / beta, vp: w / beta,
      delta: alpha > 0 ? 1 / alpha : Infinity, dBpm: 20 * Math.LOG10E * alpha, etaMag: cabs(eta), etaDeg: carg(eta) * R2D, pec: false };
    // alternative exact closed forms (Hayt eqs. for α and β in terms of the loss tangent)
    const k0 = w * Math.sqrt(mu * eps / 2), sq = Math.sqrt(1 + lt * lt);
    m.alphaAlt = k0 * Math.sqrt(sq - 1); m.betaAlt = k0 * Math.sqrt(sq + 1);
    // approximations by regime
    const etaL = Math.sqrt(mu / eps), betaL = w * Math.sqrt(mu * eps);
    if (lt === 0) m.approx = { name: 'lossless', alpha: 0, beta: betaL, eta: cx(etaL, 0) };
    else if (lt <= LT_GOOD_DIEL) m.approx = { name: 'low-loss', alpha: sigma / 2 * etaL, beta: betaL * (1 + lt * lt / 8), eta: cx(etaL, etaL * lt / 2) };
    else if (lt >= LT_GOOD_COND) { const ab = Math.sqrt(PI * f * mu * sigma), en = Math.sqrt(PI * f * mu / sigma); m.approx = { name: 'good conductor', alpha: ab, beta: ab, eta: cx(en, en) }; }
    else m.approx = null;
    return m;
  }
  const PEC = { pec: true, eta: cx(0, 0), er: Infinity, mr: 1, sigma: Infinity };
  // Polarization of E = ax Ex0 cos(ωt − βz + φx) + ay Ey0 cos(ωt − βz + φy), wave travelling in +az.
  // Handedness: IEEE (= Hayt): thumb of the hand along the direction of propagation, fingers curl the way E rotates in time.
  function polarization(Ex0, Ey0, phxDeg, phyDeg) {
    if (![Ex0, Ey0, phxDeg, phyDeg].every(fin)) return { error: 'non-numeric input' };
    if (Ex0 < 0) { Ex0 = -Ex0; phxDeg += 180; } if (Ey0 < 0) { Ey0 = -Ey0; phyDeg += 180; }
    if (Ex0 === 0 && Ey0 === 0) return { error: 'both amplitudes are zero: no wave' };
    let d = ((phyDeg - phxDeg) % 360 + 360) % 360; if (d > 180) d -= 360;   // δ = φy − φx in (−180, 180]
    const dr = d * D2R, sd = Math.sin(dr), cd = Math.cos(dr), mx = Math.max(Ex0, Ey0), tiny = 1e-9;
    const s = Ex0 * Ex0 + Ey0 * Ey0, q = Math.sqrt(Math.max(0, Ex0 ** 4 + Ey0 ** 4 + 2 * Ex0 * Ex0 * Ey0 * Ey0 * Math.cos(2 * dr)));
    const major = Math.sqrt((s + q) / 2), minor = Math.sqrt(Math.max(0, (s - q) / 2));
    let type, hand = null, tilt = 0.5 * Math.atan2(2 * Ex0 * Ey0 * cd, Ex0 * Ex0 - Ey0 * Ey0) * R2D;
    if (Ex0 <= tiny * mx || Ey0 <= tiny * mx || Math.abs(sd) < tiny) { type = 'linear'; tilt = Math.atan2(Ey0 * (cd < 0 ? -1 : 1), Ex0) * R2D; if (tilt > 90) tilt -= 180; if (tilt <= -90) tilt += 180; }
    else if (Math.abs(Ex0 - Ey0) <= tiny * mx && Math.abs(cd) < tiny) { type = 'circular'; tilt = NaN; }
    else type = 'elliptical';
    if (type !== 'linear') hand = sd < 0 ? 'right' : 'left';    // sin δ < 0: E turns from +x toward +y, i.e. CCW about +az
    const AR = type === 'linear' ? Infinity : type === 'circular' ? 1 : major / minor;
    return { type, hand, deltaDeg: d, AR, ARdB: 20 * Math.log10(AR), tiltDeg: tilt, major, minor, Ex0, Ey0, phx: phxDeg, phy: phyDeg };
  }
  // instantaneous E at plane z (in the medium m), time t
  function waveE(pol, m, z, t) {
    const att = Math.exp(-m.alpha * z), ph = m.omega * t - m.beta * z;
    return [pol.Ex0 * att * Math.cos(ph + pol.phx * D2R), pol.Ey0 * att * Math.cos(ph + pol.phy * D2R)];
  }
  // S_avg = ½ Re{E × H*} = (Ex0² + Ey0²)/2 · Re(1/η*) · e^{−2αz} az ; power through area A
  function poyntingWave(Ex0, Ey0, m, z, area) {
    const re = cdiv(cx(1, 0), m.eta).re, S0 = (Ex0 * Ex0 + Ey0 * Ey0) / 2 * re, Sz = S0 * Math.exp(-2 * m.alpha * z);
    return { S0, Sz, P: Sz * area, cosEta: Math.cos(carg(m.eta)), lossDb: 20 * Math.LOG10E * m.alpha * z };
  }

  /* ======================= 5. Normal incidence & standing waves ======================= */
  // medium 1 at z < 0 (incident travelling +z), medium 2 at z > 0 (or PEC). E0 = incident amplitude at z = 0.
  function normalIncidence(m1, m2) {
    const eta2 = m2.pec ? cx(0, 0) : m2.eta, G = cdiv(csub(eta2, m1.eta), cadd(eta2, m1.eta)), tau = cadd(cx(1, 0), G);
    const Gm = cabs(G); let thG = carg(G); if (thG < 0) thG += 2 * PI; if (Gm < 1e-14) thG = 0;
    const SWR = Gm < 1 - 1e-12 ? (1 + Gm) / (1 - Gm) : Infinity;
    const lam1 = m1.lambda; let lmax = NaN, lmin = NaN;
    if (Gm > 1e-12 && m1.alpha === 0) { lmax = thG / (4 * PI) * lam1; if (lmax >= lam1 / 2 - 1e-15 * lam1) lmax -= lam1 / 2; lmin = lmax >= lam1 / 4 ? lmax - lam1 / 4 : lmax + lam1 / 4; }
    const y1 = cdiv(cx(1, 0), m1.eta).re, R = Gm * Gm, T = m2.pec ? 0 : cabs(tau) ** 2 * cdiv(cx(1, 0), m2.eta).re / y1;
    return { G, Gm, thG, tau, SWR, lmax, lmin, R, T, lam1, eta1: m1.eta, eta2 };
  }
  // total phasor field in medium 1 (z ≤ 0) and medium 2 (z ≥ 0), per unit incident amplitude
  function fieldNormal(m1, m2, ni, z) {
    if (z <= 0) return cadd(cexp(cscale(m1.gamma, -z)), cmul(ni.G, cexp(cscale(m1.gamma, z))));
    if (m2.pec) return cx(0, 0);
    return cmul(ni.tau, cexp(cscale(m2.gamma, -z)));
  }
  // three-media problem: slab (m2, thickness d) between m1 (z<0) and m3 (z>d, may be PEC)
  function slab(m1, m2, d, m3) {
    const e1 = m1.eta, e2 = m2.eta, e3 = m3.pec ? cx(0, 0) : m3.eta, gd = cscale(m2.gamma, d);
    // ηin = η2 (η3 cosh γd + η2 sinh γd)/(η2 cosh γd + η3 sinh γd)  (= tanh form, but finite at βd = π/2)
    let ch, sh; if (gd.re > 20) { ch = cx(1, 0); sh = ctanh(gd); } else { ch = ccosh(gd); sh = csinh(gd); }
    const etaIn = cmul(e2, cdiv(cadd(cmul(e3, ch), cmul(e2, sh)), cadd(cmul(e2, ch), cmul(e3, sh))));
    const Gin = cdiv(csub(etaIn, e1), cadd(etaIn, e1)), G23 = cdiv(csub(e3, e2), cadd(e3, e2));
    const E0 = cadd(cx(1, 0), Gin), E2p = cdiv(E0, cadd(cx(1, 0), cmul(G23, cexp(cscale(m2.gamma, -2 * d)))));
    const Ed = cmul(cmul(E2p, cexp(cscale(m2.gamma, -d))), cadd(cx(1, 0), G23));
    const Gm = cabs(Gin), SWR = Gm < 1 - 1e-12 ? (1 + Gm) / (1 - Gm) : Infinity;
    const y1 = cdiv(cx(1, 0), e1).re, T3 = m3.pec ? 0 : cabs(Ed) ** 2 * cdiv(cx(1, 0), e3).re / y1;
    return { etaIn, Gin, Gm, SWR, G23, E2p, Ed, R: Gm * Gm, T3, absorbed: Math.max(0, 1 - Gm * Gm - T3) };
  }
  function fieldSlab(m1, m2, d, m3, s, z) {
    if (z <= 0) return cadd(cexp(cscale(m1.gamma, -z)), cmul(s.Gin, cexp(cscale(m1.gamma, z))));
    if (z <= d) return cmul(s.E2p, cadd(cexp(cscale(m2.gamma, -z)), cmul(s.G23, cexp(cscale(m2.gamma, z - 2 * d)))));
    if (m3.pec) return cx(0, 0);
    return cmul(s.Ed, cexp(cscale(m3.gamma, -(z - d))));
  }

  FSP.math.fields = {
    EPS0, MU0, C0, KE, ETA0, cx, cabs, carg, csqrt, ctanh, asimpV, glRule, glComp,
    pointCharges, pointEV, lineFinite, lineFiniteNum, lineInfinite, lineInfiniteNum, ringAxis, ringNum, diskAxis, diskNum,
    sheetInf, sheetInfNum, rectSheet, rectSheetNum, squareSheetAxis,
    filamentFinite, filamentInf, filamentInfNum, bsSegment, bsPolygon, bsCircle, loopAxis, squareAxis, squarePts, solenoidAxis, solenoidNum,
    coaxH, coaxIencNum, sheetK, sheetKNum,
    faradayFixed, slidingBar, rotatingLoop, currentRatio, classifyLT, LT_GOOD_COND, LT_GOOD_DIEL,
    medium, PEC, polarization, waveE, poyntingWave, normalIncidence, fieldNormal, slab, fieldSlab,
  };

  /* ======================= tests ======================= */

  FSP.math.fields._shared = { D2R, KE, PI, R2D, cabs, carg, cdiv, cmul, cx, fin, vlen };
})();
