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
  FSP.registerTests('fields', t => {
    const M = FSP.math.fields, nC = 1e-9, relv = (a, b) => vlen([a[0] - b[0], a[1] - b[1], a[2] - b[2]]) / Math.max(vlen(b), 1e-300);
    // --- Coulomb ---
    // 1 nC at origin, P = (1,0,0): Ex = Q/(4πε0·1²) = 8.987551787 V/m (k = 8.987551787e9), V = 8.987551787 V
    const p1 = M.pointCharges([{ q: nC, x: 0, y: 0, z: 0 }], [1, 0, 0]);
    t.check('point charge 1 nC at 1 m: E = 8.98755 V/m, V = 8.98755 V', t.rel(p1.E[0], 8.987551787, 1e-8) && t.rel(p1.V, 8.987551787, 1e-8) && p1.E[1] === 0);
    // dipole ±1 nC at (±1,0,0), P = (0,1,0): E = kq(−2,0,0)/2^1.5 = −6.355150 x̂ V/m (hand), V = 0
    const p2 = M.pointCharges([{ q: nC, x: 1, y: 0, z: 0 }, { q: -nC, x: -1, y: 0, z: 0 }], [0, 1, 0]);
    t.check('dipole on bisector: E = −6.35515 ax V/m, V = 0', t.rel(p2.E[0], -8.987551787 * 2 / Math.pow(2, 1.5), 1e-9) && t.near(p2.E[1], 0, 1e-12) && t.near(p2.V, 0, 1e-12), p2.E[0].toFixed(6));
    t.check('field point on a point charge -> error, not Infinity', !!M.pointCharges([{ q: nC, x: 0, y: 0, z: 0 }], [0, 0, 0]).error);
    // infinite line ρL = 1 nC/m at ρ = 0.5 m: ρL/(2πε0ρ) = 2k·ρL/ρ = 35.95020715 V/m
    const inf = M.lineInfinite(nC, [0.5, 0, 0], 1), lng = M.lineFiniteNum(nC, -1000, 1000, [0.5, 0, 0]);
    t.check('infinite line E = ρL/(2πε0ρ) = 35.9502 V/m', t.rel(inf.Erho, 2 * 8.987551787 / 0.5, 1e-8), inf.Erho.toFixed(5));
    t.check('numeric 2 km finite line (adaptive Simpson) → infinite-line value (rel < 1e-5)', t.rel(lng.E[0], inf.Erho, 1e-5) && Math.abs(lng.E[2]) < 1e-6 * inf.Erho, lng.E[0].toFixed(6));
    t.check('infinite line, mapped numeric integral z\' = z + ρ tanθ agrees (1e-7)', t.rel(M.lineInfiniteNum(nC, [0.3, 0.4, 2]).E[0], M.lineInfinite(nC, [0.3, 0.4, 2]).E[0], 1e-7));
    const lf = M.lineFinite(2 * nC, -0.5, 1.5, [0.3, -0.4, 0.2]), lfn = M.lineFiniteNum(2 * nC, -0.5, 1.5, [0.3, -0.4, 0.2]);
    t.check('finite line: closed form (Eρ, Ez, V) = numeric quadrature at off-centre point (1e-8)', relv(lfn.E, lf.E) < 1e-8 && t.rel(lfn.V, lf.V, 1e-8), 'rel ' + relv(lfn.E, lf.E).toExponential(2));
    const lfa = M.lineFinite(nC, -1, 1, [0, 0, 3]), lfan = M.lineFiniteNum(nC, -1, 1, [0, 0, 3]);
    t.check('finite line, field point on the axis beyond the end: kρL(1/2 − 1/4) and V = kρL ln 2', t.rel(lfa.E[2], 8.987551787 * 0.25, 1e-9) && t.rel(lfa.V, 8.987551787 * Math.LN2, 1e-9) && t.rel(lfan.E[2], lfa.E[2], 1e-8));
    // ring a = 1, z = 1, ρL = 1 nC/m: Ez = ρL a z / (2ε0 (a²+z²)^1.5) = 1e-9/(2·8.8541878128e-12·2^1.5) = 19.96498 V/m
    const ra = M.ringAxis(nC, 1, 1), rn = M.ringNum(nC, 1, [0, 0, 1]);
    t.check('ring on axis: closed form 19.965 V/m = numeric', t.rel(ra.Ez, 1e-9 / (2 * 8.8541878128e-12 * Math.pow(2, 1.5)), 1e-9) && t.rel(rn.E[2], ra.Ez, 1e-9) && Math.abs(rn.E[0]) < 1e-9 * ra.Ez && t.rel(rn.V, ra.V, 1e-9), ra.Ez.toFixed(5));
    const ro = M.ringNum(nC, 1, [0.6, 0.2, 0.4]), rog = M.glComp(ph => { const R = [0.6 - Math.cos(ph), 0.2 - Math.sin(ph), 0.4], r = vlen(R), k = KE * nC / r ** 3; return [k * R[0], k * R[1], k * R[2]]; }, 0, 2 * PI, 32);
    t.check('ring off axis: adaptive Simpson = composite Gauss–Legendre (1e-8)', relv(ro.E, rog) < 1e-8);
    // disk → infinite sheet: a = 1e4 z gives ρs/(2ε0)(1 − 1/√(1+1e8)) ≈ ρs/(2ε0)(1 − 1e-4)
    const big = M.diskAxis(nC, 1e4, 1), sh = M.sheetInf(nC, 1);
    t.check('disk a ≫ z → infinite sheet ρs/(2ε0) = 56.4705 V/m', t.rel(sh.Ez, 1e-9 / (2 * 8.8541878128e-12), 1e-9) && t.rel(big.Ez, sh.Ez, 1.1e-4) && big.Ez < sh.Ez, big.Ez.toFixed(4));
    const da = M.diskAxis(nC, 1, 0.5), dn = M.diskNum(nC, 1, [0, 0, 0.5], 1e-8);
    t.check('disk on axis: closed form = 2-D numeric (E and V)', t.rel(dn.E[2], da.Ez, 1e-7) && t.rel(dn.V, da.V, 1e-7), dn.E[2].toFixed(6) + ' vs ' + da.Ez.toFixed(6));
    t.check('infinite sheet numeric (mapped ρ\' = |z|tanθ) = ρs/(2ε0), below sheet negative', t.rel(M.sheetInfNum(nC, 0.7, 1e-8).E[2], sh.Ez, 1e-7) && t.rel(M.sheetInfNum(nC, -0.7, 1e-8).E[2], -sh.Ez, 1e-7));
    // square sheet side 2 m (b = 1), z = 1: Ez = (ρs/πε0) atan(1/√3) = (ρs/πε0)(π/6) = ρs/(6ε0) = 18.82350 V/m
    const sq = M.squareSheetAxis(nC, 2, 1), sqr = M.rectSheet(nC, -1, 1, -1, 1, [0, 0, 1]), sqn = M.rectSheetNum(nC, -1, 1, -1, 1, [0, 0, 1], 1e-8);
    t.check('square sheet tutorial: Ez = ρs/(6ε0) = 18.8235 V/m (4 quadrants, arctan form)', t.rel(sq.Ez, 1e-9 / (6 * 8.8541878128e-12), 1e-12) && t.rel(sqr.E[2], sq.Ez, 1e-12), sq.Ez.toFixed(5));
    t.check('square sheet: numeric quadrature = closed form, Ex = Ey = 0 by symmetry', t.rel(sqn.E[2], sq.Ez, 1e-7) && Math.abs(sqn.E[0]) < 1e-8 && Math.abs(sqr.E[0]) < 1e-12);
    const rg = M.rectSheet(nC, -0.5, 1.5, -1, 0.3, [0.9, 0.8, 0.25]), rgn = M.rectSheetNum(nC, -0.5, 1.5, -1, 0.3, [0.9, 0.8, 0.25], 1e-8);
    t.check('rectangular sheet, general point: closed-form E and V = numeric (1e-6)', relv(rgn.E, rg.E) < 1e-6 && t.rel(rgn.V, rg.V, 1e-6), 'E rel ' + relv(rgn.E, rg.E).toExponential(2) + ', V ' + rg.V.toFixed(5));
    t.check('square sheet z → 0+ tends to ρs/(2ε0); field point in the sheet plane -> error', t.rel(M.squareSheetAxis(nC, 2, 1e-6).Ez, sh.Ez, 1e-5) && !!M.rectSheet(nC, -1, 1, -1, 1, [0, 0, 0]).error);
    // --- magnetostatics ---
    // finite filament −1..1, I = 1 A, ρ = 1, z = 0: sinα2 = 1/√2, sinα1 = −1/√2 -> H = √2/(4π) = 0.1125395 A/m
    const ff = M.filamentFinite(1, -1, 1, [1, 0, 0]), ffn = M.bsSegment(1, [0, 0, -1], [0, 0, 1], [1, 0, 0]);
    t.check('finite filament H = I/(4πρ)(sinα2 − sinα1) = √2/(4π) and points along +aφ', t.rel(ff.Hphi, Math.SQRT2 / (4 * PI), 1e-12) && t.rel(ff.H[1], ff.Hphi, 1e-12) && t.rel(ffn[1], ff.Hphi, 1e-10), ff.Hphi.toFixed(7));
    const fa = M.filamentFinite(3, 0.2, 2.5, [-0.3, 0.5, 1.9]), fan = M.bsSegment(3, [0, 0, 0.2], [0, 0, 2.5], [-0.3, 0.5, 1.9]);
    t.check('finite filament, asymmetric: angle formula = numeric Biot–Savart (1e-9)', relv(fan, fa.H) < 1e-9);
    t.check('infinite filament: mapped numeric Biot–Savart = I/(2πρ)', t.rel(M.filamentInfNum(3, [0.3, 0.4, 1])[1], 3 / (2 * PI * 0.5) * 0.3 / 0.5, 1e-10));
    t.check('long filament → I/(2πρ)', t.rel(M.filamentFinite(5, -1e4, 1e4, [0.2, 0, 0]).Hphi, 5 / (2 * PI * 0.2), 1e-8));
    // loop I = 2 A, a = 0.5 m: centre H = I/(2a) = 2 A/m
    const lc = M.bsCircle(2, 0.5, [0, 0, 0]);
    t.check('loop centre H = I/(2a) = 2 A/m (formula and numeric)', t.rel(M.loopAxis(2, 0.5, 0).Hz, 2, 1e-12) && t.rel(lc[2], 2, 1e-10));
    t.check('loop on axis: Ia²/(2(a²+z²)^1.5) = numeric Biot–Savart', t.rel(M.bsCircle(2, 0.5, [0, 0, 0.3])[2], M.loopAxis(2, 0.5, 0.3).Hz, 1e-10));
    // square loop side w = 1, I = 1: centre H = 2√2 I/(π w) = 0.9003163 A/m
    t.check('square loop centre 2√2I/(πw) = 0.90032 A/m; on-axis formula = numeric', t.rel(M.squareAxis(1, 1, 0).Hz, 2 * Math.SQRT2 / PI, 1e-12) && t.rel(M.bsPolygon(1, M.squarePts(1), [0, 0, 0.4])[2], M.squareAxis(1, 1, 0.4).Hz, 1e-9));
    const so = M.solenoidAxis(1, 1000, 1, 0.01, 0);
    t.check('long solenoid centre → nI = 1000 A/m', t.rel(so.Hz, 1000, 1e-3), so.Hz.toFixed(3));
    const so2 = M.solenoidAxis(2, 200, 0.4, 0.05, 0.12), so2n = M.solenoidNum(2, 200, 0.4, 0.05, [0, 0, 0.12]);
    t.check('finite solenoid on axis: angle formula = numeric current-sheet Biot–Savart', t.rel(so2n[2], so2.Hz, 1e-8), so2.Hz.toFixed(4));
    // coax a=1mm, b=4mm, c=5mm, I=10 A
    const cA = 1e-3, cB = 4e-3, cC = 5e-3, h = r => M.coaxH(10, cA, cB, cC, r).Hphi;
    t.check('coax continuity at ρ = a and ρ = b; H(c) = 0; H = 0 outside', t.rel(h(cA * (1 - 1e-12)), 10 / (2 * PI * cA), 1e-9) && t.rel(h(cA), 10 / (2 * PI * cA), 1e-12) && t.rel(h(cB * (1 - 1e-12)), h(cB), 1e-9) && Math.abs(h(cC)) < 1e-9 && h(6e-3) === 0);
    t.check('coax H(ρ) = I_enc/(2πρ) with I_enc from numeric ∫J·dS (all 4 regions)', [0.5e-3, 2e-3, 4.5e-3, 7e-3].every(r => t.near(M.coaxIencNum(10, cA, cB, cC, r) / (2 * PI * r), h(r), 1e-9 * 10 / (2 * PI * cA))));
    t.check('coax with b ≤ a -> error', !!M.coaxH(1, 2, 1, 3, 0.5).error);
    t.check('current sheet K ax: H = −K/2 ay above, numeric filament sum agrees', t.near(M.sheetK(4, 1).Hy, -2, 1e-15) && t.rel(M.sheetKNum(4, 0.3).Hy, -2, 1e-8) && t.rel(M.sheetKNum(4, -0.3).Hy, 2, 1e-8));
    // --- Faraday ---
    const sb = M.slidingBar(0.5, 0.2, 3, 2, 0.1, 0);
    t.check('sliding bar: |emf| = Bvl = 0.3 V, I = 0.15 A, F = BIl = 0.015 N, Fv = I²R = 0.045 W', t.near(sb.mag, 0.3, 1e-15) && t.near(sb.I, 0.15, 1e-15) && t.near(sb.F, 0.015, 1e-15) && t.near(sb.Pmech, 0.045, 1e-15) && t.near(sb.Pelec, 0.045, 1e-15));
    const fx = M.faradayFixed(10, 0.01, 0.1, 50, 0.0013), dt = 1e-7, dphi = (M.faradayFixed(10, 0.01, 0.1, 50, 0.0013 + dt).flux - M.faradayFixed(10, 0.01, 0.1, 50, 0.0013 - dt).flux) / (2 * dt);
    t.check('fixed loop: emf = −N dΦ/dt (vs finite difference), peak NB0Aω = 10·0.1·0.01·100π = π V', t.rel(fx.emf, -10 * dphi, 1e-7) && t.rel(fx.peak, PI, 1e-12));
    t.check('rotating loop: peak emf NBAω, emf = 0 when flux is maximum', t.rel(M.rotatingLoop(50, 0.02, 0.4, 60, 0).peak, 50 * 0.4 * 0.02 * 120 * PI, 1e-12) && M.rotatingLoop(50, 0.02, 0.4, 60, 0).emf === 0);
    // seawater σ = 4 S/m, εr = 81, 1 MHz: σ/(ωε) = 4/(2π·1e6·81·8.8541878e-12) = 887.65
    const sw = M.currentRatio(4, 81, 1e6, 1);
    t.check('seawater 1 MHz: σ/(ωε) = 887.6 → good conductor; Jc/Jd = loss tangent', t.rel(sw.lt, 887.65, 1e-4) && sw.cls === 'good conductor' && t.rel(sw.Jc / sw.Jd, sw.lt, 1e-12), sw.lt.toFixed(2));
    t.check('seawater crossover f = σ/(2πε) = 887.6 MHz; at 10 GHz it is quasi-conductor', t.rel(sw.fCross, 887.65e6, 1e-4) && M.currentRatio(4, 81, 1e10, 1).cls.indexOf('quasi') === 0);
    // --- plane waves ---
    const sea = M.medium(81, 1, 4, 1e6);
    t.check('seawater 1 MHz: δ ≈ 1/√(πfμσ) = 0.2516 m, classification good conductor', t.rel(sea.delta, 0.2516, 2e-3) && sea.cls === 'good conductor', sea.delta.toFixed(5) + ' m');
    t.check('γ from √(jωμ(σ+jωε)) = Hayt α/β loss-tangent formulas (1e-10)', t.rel(sea.alpha, sea.alphaAlt, 1e-10) && t.rel(sea.beta, sea.betaAlt, 1e-10));
    const d4 = M.medium(4, 1, 0, 1e9);
    t.check('lossless εr = 4: η = η0/2 = 188.37 Ω (real), α = 0, β = 4πf/c = 41.917 rad/m', t.rel(d4.eta.re, 376.7303 / 2, 1e-6) && d4.eta.im === 0 && d4.alpha === 0 && t.rel(d4.beta, 4 * PI * 1e9 / 299792458, 1e-9) && d4.cls === 'lossless', d4.eta.re.toFixed(3));
    const ll = M.medium(2.25, 1, 1e-4, 1e9);
    t.check('low-loss approximations (α ≈ ση/2, β, η) agree with exact (lt = 8e-4)', ll.cls.indexOf('good dielectric') === 0 && t.rel(ll.approx.alpha, ll.alpha, 1e-6) && t.rel(ll.approx.beta, ll.beta, 1e-9) && t.rel(ll.approx.eta.im, ll.eta.im, 1e-5));
    t.check('copper 1 MHz: good-conductor α = β = 1/δ, η ≈ (1+j)/(σδ), δ = 66.09 µm', (() => { const cu = M.medium(1, 1, 5.8e7, 1e6); return t.rel(cu.delta, 66.09e-6, 1e-3) && t.rel(cu.approx.alpha, cu.alpha, 1e-9) && t.rel(cu.eta.re, 1 / (5.8e7 * cu.delta), 1e-9); })());
    t.check('medium with εr ≤ 0 or f ≤ 0 -> error', !!M.medium(0, 1, 0, 1e9).error && !!M.medium(1, 1, 0, 0).error && !!M.medium(1, 1, -1, 1).error);
    // polarization (IEEE/Hayt handedness, propagation +az)
    const P_ = (a, b, x, y) => M.polarization(a, b, x, y);
    t.check('pol: Ey0 = 0 -> linear along x (tilt 0)', P_(1, 0, 0, 37).type === 'linear' && t.near(P_(1, 0, 0, 37).tiltDeg, 0, 1e-12));
    t.check('pol: Ex0 = Ey0, δ = 0 -> linear at 45°; δ = 180° -> linear at −45°', P_(1, 1, 10, 10).type === 'linear' && t.near(P_(1, 1, 10, 10).tiltDeg, 45, 1e-9) && t.near(P_(1, 1, 0, 180).tiltDeg, -45, 1e-9));
    t.check('pol: Ex0 = Ey0, δ = −90° (Ey = sin(ωt−βz)) -> right-hand circular, AR = 1', P_(1, 1, 0, -90).type === 'circular' && P_(1, 1, 0, -90).hand === 'right' && P_(1, 1, 0, -90).AR === 1);
    t.check('pol: Ex0 = Ey0, δ = +90° -> left-hand circular', P_(2, 2, 30, 120).type === 'circular' && P_(2, 2, 30, 120).hand === 'left');
    t.check('pol: (2, 1, δ = −90°) -> RH elliptical, AR = 2, major axis along x', P_(2, 1, 0, -90).type === 'elliptical' && P_(2, 1, 0, -90).hand === 'right' && t.near(P_(2, 1, 0, -90).AR, 2, 1e-12) && t.near(P_(2, 1, 0, -90).tiltDeg, 0, 1e-12));
    // (1,1,δ=−45°): axes² = (2 ± √2)/2 -> AR = √((2+√2)/(2−√2)) = 1+√2, tilt 45°
    t.check('pol: (1, 1, δ = −45°) -> elliptical AR = 1+√2, tilt 45°', t.rel(P_(1, 1, 0, -45).AR, 1 + Math.SQRT2, 1e-12) && t.near(P_(1, 1, 0, -45).tiltDeg, 45, 1e-9));
    // independent handedness check: sample E(t) at z = 0 and take sign of (E × dE/dt)·az (CCW about +az = right-hand for +z travel)
    const sense = (a, b, x, y) => { const e = tt => [a * Math.cos(tt + x * D2R), b * Math.cos(tt + y * D2R)]; let s = 0; for (let i = 0; i < 16; i++) { const t0 = i * PI / 8, e0 = e(t0), e1 = e(t0 + 1e-4); s += e0[0] * e1[1] - e0[1] * e1[0]; } return s > 0 ? 'right' : 'left'; };
    t.check('pol: handedness agrees with sampled rotation sense for 6 random cases', [[1, 0.4, 0, 70], [0.3, 1, 20, -50], [1, 1, 0, -120], [2, 0.5, 200, 10], [1, 3, -30, 100], [0.7, 0.7, 45, -10]].every(c => P_.apply(null, c).hand === sense.apply(null, c)));
    t.check('pol: zero amplitudes -> error', !!P_(0, 0, 0, 0).error);
    const air = M.medium(1, 1, 0, 1e9);
    t.check('Poynting: 1 V/m in air -> S_avg = 1/(2η0) = 1.32721 mW/m²; 1 m² -> 1.32721 mW', t.rel(M.poyntingWave(1, 0, air, 5, 1).P, 1 / (2 * 376.7303), 1e-6));
    t.check('Poynting in lossy medium decays as e^{−2αz} (8.686α dB/m)', (() => { const p = M.poyntingWave(1, 1, sea, 0.5, 1); return t.rel(p.Sz / p.S0, Math.exp(-2 * sea.alpha * 0.5), 1e-12) && t.rel(p.lossDb, 8.685889638 * sea.alpha * 0.5, 1e-9); })());
    // --- normal incidence ---
    const ni = M.normalIncidence(air, d4);
    t.check('air → εr = 4: Γ = −1/3, τ = 2/3, SWR = 2', t.near(ni.G.re, -1 / 3, 1e-12) && t.near(ni.G.im, 0, 1e-15) && t.near(ni.tau.re, 2 / 3, 1e-12) && t.near(ni.SWR, 2, 1e-12));
    t.check('air → εr = 4: R = 1/9, T = 8/9; |E1| minimum at the boundary, maximum λ/4 away', t.near(ni.R, 1 / 9, 1e-12) && t.near(ni.T, 8 / 9, 1e-12) && t.near(ni.lmin, 0, 1e-15) && t.rel(ni.lmax, air.lambda / 4, 1e-12));
    t.check('|E1(z)| reaches 1+|Γ| and 1−|Γ| at those points', t.rel(M.cabs(M.fieldNormal(air, d4, ni, -ni.lmax)), 4 / 3, 1e-12) && t.rel(M.cabs(M.fieldNormal(air, d4, ni, 0)), 2 / 3, 1e-12));
    const np = M.normalIncidence(air, M.PEC), lam = air.lambda;
    t.check('PEC: Γ = −1, SWR = ∞, minima (zeros) at 0, λ/2, λ; |E1| = 2 at λ/4', t.near(np.G.re, -1, 1e-15) && np.SWR === Infinity && np.lmin === 0 && [0, 0.5, 1].every(k => M.cabs(M.fieldNormal(air, M.PEC, np, -k * lam)) < 1e-9) && t.rel(M.cabs(M.fieldNormal(air, M.PEC, np, -lam / 4)), 2, 1e-12));
    t.check('εr = 4 → air: Γ = +1/3, maximum at boundary', (() => { const r = M.normalIncidence(d4, air); return t.near(r.G.re, 1 / 3, 1e-12) && t.near(r.lmax, 0, 1e-15) && t.rel(r.lmin, d4.lambda / 4, 1e-12); })());
    // quarter-wave matching: air → εr2 = 4 (λ2/4) → εr3 = 16: η2 = √(η1η3) -> Γin = 0
    const e16 = M.medium(16, 1, 0, 1e9), qw = M.slab(air, d4, d4.lambda / 4, e16);
    t.check('quarter-wave slab √(η1η3), d = λ2/4 -> Γin = 0, all power to medium 3', qw.Gm < 1e-12 && t.rel(qw.T3, 1, 1e-12));
    t.check('half-wave slab is transparent: ηin = η3', (() => { const s = M.slab(air, M.medium(9, 1, 0, 1e9), M.medium(9, 1, 0, 1e9).lambda / 2, d4); return t.rel(s.etaIn.re, d4.eta.re, 1e-9) && Math.abs(s.etaIn.im) < 1e-9; })());
    t.check('slab of zero thickness = direct interface', (() => { const s = M.slab(air, e16, 0, d4); return t.near(s.Gin.re, -1 / 3, 1e-12); })());
  });
  /* ============================================================ UI */
  const css = n => { try { return getComputedStyle(document.documentElement).getPropertyValue(n).trim() || null; } catch (e) { return null; } };
  function theme() {
    return { bg: css('--panel2') || '#0f141f', text: css('--text') || '#e6edf3', muted: css('--muted') || '#8b98ab', border: css('--border') || '#243047',
      a: css('--accent2') || '#f9a03f', b: css('--ok') || '#2ecc71', c: css('--accent') || '#4cc9f0', pink: css('--pink') || '#f72585', warn: css('--warn') || '#ffb347', bad: css('--bad') || '#ff5c5c', mono: css('--mono') || 'monospace' };
  }
  const fe = (x, u, d) => FSP.fmtEng(x, u, d === undefined ? 4 : d);
  // plain number for worked solutions: 5 significant figures, ×10^n for very large/small, true minus sign
  function g(x, d) {
    if (x === Infinity) return '∞'; if (x === -Infinity) return '−∞'; if (!fin(x)) return '—';
    d = d || 5; if (x === 0) return '0'; const a = Math.abs(x);
    let s = (a >= 1e5 || a < 1e-3) ? x.toExponential(d - 1).replace(/\.?0+e/, 'e') : String(+x.toPrecision(d));
    return s.replace(/e\+?(-?\d+)/, '×10^$1').replace(/-/g, '−');
  }
  const clean = z => { const m = cabs(z); if (m < 1e-12) return cx(0, 0); return cx(Math.abs(z.re) < 1e-12 * m ? 0 : z.re, Math.abs(z.im) < 1e-12 * m ? 0 : z.im); };   // drop round-off dust such as 1e-17
  const z0 = x => (Math.abs(x) < 1e-12 ? 0 : x);
  const gc = (z, d) => { z = clean(z); return g(z.re, d) + (z.im < 0 ? ' − j' : ' + j') + g(Math.abs(z.im), d); };
  const gp = (z, d) => { z = clean(z); return g(cabs(z), d) + ' ∠ ' + g(carg(z) * R2D, 4) + '°'; };
  const vec3 = (v, d) => { const m = Math.max(Math.abs(v[0]), Math.abs(v[1]), Math.abs(v[2])); return '(' + v.map(x => g(Math.abs(x) < 1e-12 * m ? 0 : x, d)).join(', ') + ')'; };
  const tickStr = v => (Math.abs(v) < 1e-12 ? '0' : String(+v.toPrecision(4)));
  function niceTicks(lo, hi, n) {
    const span = hi - lo; if (!(span > 0)) return [lo];
    const raw = span / n, mag = Math.pow(10, Math.floor(Math.log10(raw))), f = raw / mag, step = (f < 1.5 ? 1 : f < 3 ? 2 : f < 7 ? 5 : 10) * mag, out = [];
    for (let v = Math.ceil(lo / step - 1e-9) * step; v <= hi + step * 1e-9; v += step) out.push(Math.abs(v) < step * 1e-9 ? 0 : +v.toPrecision(12));
    return out;
  }
  const logTicks = (lo, hi) => { const out = []; for (let e = Math.ceil(Math.log10(lo) - 1e-9); e <= Math.floor(Math.log10(hi) + 1e-9); e++) out.push(Math.pow(10, e)); return out; };
  const logStr = v => { const e = Math.round(Math.log10(v)); return e >= -2 && e <= 3 ? String(v) : '1e' + e; };
  const logspace = (a, b, n) => { const out = []; for (let i = 0; i < n; i++) out.push(a * Math.pow(b / a, i / (n - 1))); return out; };
  const linspace = (a, b, n) => { const out = []; for (let i = 0; i < n; i++) out.push(a + (b - a) * i / (n - 1)); return out; };
  // decade-snapped y range for a log axis from positive data
  function decadeRange(vals, span) {
    let lo = Infinity, hi = 0; vals.forEach(v => { if (fin(v) && v > 0) { lo = Math.min(lo, v); hi = Math.max(hi, v); } });
    if (!(hi > 0)) return [1, 10]; lo = Math.max(lo, hi * (span || 1e-8));
    let a = Math.pow(10, Math.floor(Math.log10(lo) + 1e-9)), b = Math.pow(10, Math.ceil(Math.log10(hi) - 1e-9));
    if (b <= a) { a /= 10; b *= 10; } if (b / a < 10) b = a * 10; return [a, b];
  }
  function linRange(vals, padFrac) {
    let lo = Infinity, hi = -Infinity; vals.forEach(v => { if (fin(v)) { lo = Math.min(lo, v); hi = Math.max(hi, v); } });
    if (!fin(lo)) return [-1, 1]; if (hi - lo < 1e-300 * Math.max(1, Math.abs(hi))) { const m = Math.abs(hi) || 1; return [lo - m, hi + m]; }
    const p = (hi - lo) * (padFrac === undefined ? 0.08 : padFrac); return [lo - p, hi + p];
  }

  // Generic axes + series (same engine as the other modules' plots). Returns mapping helpers.
  function plot(c, o) {
    const gg = c.prep(), ctx = gg.ctx, w = gg.w, h = gg.h, T = theme();
    const L = 54, R = 12, Tp = o.topPad || 12, B = 32, pw = w - L - R, ph = h - Tp - B;
    ctx.clearRect(0, 0, w, h); ctx.fillStyle = T.bg; ctx.fillRect(0, 0, w, h);
    if (pw < 20 || ph < 20) return null;
    const fx = o.xlog ? Math.log10 : v => v, fy = o.ylog ? Math.log10 : v => v;
    const x0 = fx(o.xmin), x1 = fx(o.xmax), y0 = fy(o.ymin), y1 = fy(o.ymax);
    const X = x => L + (fx(x) - x0) / (x1 - x0) * pw, Y = y => Tp + (1 - (fy(Math.max(y, o.ylog ? o.ymin * 1e-3 : -Infinity)) - y0) / (y1 - y0)) * ph;
    ctx.font = '11px ' + T.mono; ctx.lineWidth = 1;
    ctx.textBaseline = 'middle'; ctx.textAlign = 'right';
    (o.ylog ? logTicks(o.ymin, o.ymax) : niceTicks(o.ymin, o.ymax, 5)).forEach(v => { const yy = Y(v); ctx.strokeStyle = T.border; ctx.beginPath(); ctx.moveTo(L, yy); ctx.lineTo(L + pw, yy); ctx.stroke(); ctx.fillStyle = T.muted; ctx.fillText(o.ylog ? logStr(v) : tickStr(v), L - 4, yy); });
    ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    (o.xlog ? logTicks(o.xmin, o.xmax) : niceTicks(o.xmin, o.xmax, Math.max(3, Math.floor(pw / 70)))).forEach(v => { const xx = X(v); ctx.strokeStyle = T.border; ctx.beginPath(); ctx.moveTo(xx, Tp); ctx.lineTo(xx, Tp + ph); ctx.stroke(); ctx.fillStyle = T.muted; ctx.fillText(o.xlog ? logStr(v) : tickStr(v), xx, Tp + ph + 3); });
    ctx.fillStyle = T.muted; ctx.textAlign = 'center'; ctx.textBaseline = 'bottom'; ctx.fillText(o.xlabel || '', L + pw / 2, h - 1);
    ctx.save(); ctx.translate(10, Tp + ph / 2); ctx.rotate(-PI / 2); ctx.textBaseline = 'top'; ctx.fillText(o.ylabel || '', 0, 0); ctx.restore();
    ctx.save(); ctx.beginPath(); ctx.rect(L, Tp, pw, ph); ctx.clip();
    (o.vlines || []).forEach(v => { if (!fin(v.x)) return; ctx.strokeStyle = v.color || T.muted; ctx.setLineDash([4, 4]); ctx.beginPath(); ctx.moveTo(X(v.x), Tp); ctx.lineTo(X(v.x), Tp + ph); ctx.stroke(); ctx.setLineDash([]); if (v.label) { ctx.fillStyle = v.color || T.muted; ctx.textAlign = 'left'; ctx.textBaseline = 'top'; ctx.fillText(v.label, X(v.x) + 3, Tp + 2 + (v.dy || 0)); } });
    (o.hlines || []).forEach(v => { if (!fin(v.y)) return; ctx.strokeStyle = v.color || T.muted; ctx.setLineDash([4, 4]); ctx.beginPath(); ctx.moveTo(L, Y(v.y)); ctx.lineTo(L + pw, Y(v.y)); ctx.stroke(); ctx.setLineDash([]); if (v.label) { ctx.fillStyle = v.color || T.muted; ctx.textAlign = 'right'; ctx.textBaseline = 'bottom'; ctx.fillText(v.label, L + pw - 3, Y(v.y) - 2); } });
    (o.series || []).forEach(s => {
      ctx.strokeStyle = s.color; ctx.lineWidth = s.width || 2; ctx.setLineDash(s.dash || []); ctx.beginPath(); let pen = false;
      for (let i = 0; i < s.x.length; i++) { const xv = s.x[i], yv = s.y[i]; if (!fin(xv) || !fin(yv)) { pen = false; continue; } const px = X(xv), py = Y(yv); if (pen) ctx.lineTo(px, py); else { ctx.moveTo(px, py); pen = true; } }
      ctx.stroke(); ctx.setLineDash([]);
      if (s.dots) { ctx.fillStyle = s.color; for (let i = 0; i < s.x.length; i++) if (fin(s.y[i]) && fin(s.x[i])) { ctx.beginPath(); ctx.arc(X(s.x[i]), Y(s.y[i]), 3, 0, 2 * PI); ctx.fill(); } }
    });
    if (o.extra) o.extra(ctx, X, Y, T, { L, Tp, pw, ph });
    ctx.restore();
    (o.points || []).forEach(p => {
      if (!fin(p.x) || !fin(p.y)) return; const px = X(p.x), py = Y(p.y);
      ctx.fillStyle = p.color; ctx.strokeStyle = T.bg; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.arc(px, py, 5, 0, 2 * PI); ctx.fill(); ctx.stroke();
      if (p.label) { const right = px < L + pw * 0.6; ctx.textAlign = right ? 'left' : 'right'; ctx.textBaseline = p.below ? 'top' : 'bottom'; ctx.fillStyle = T.text; ctx.fillText(p.label, px + (right ? 8 : -8), py + (p.below ? 6 : -6)); }
    });
    ctx.font = '11px ' + T.mono; ctx.textAlign = 'right'; ctx.textBaseline = 'top';
    (o.legend || []).forEach((lg, i) => { ctx.fillStyle = lg.color; ctx.fillText(lg.text, L + pw - 4, Tp + 4 + i * 14); });
    return { ctx, X, Y, T, L, Tp, pw, ph, w, h };
  }
  function blank(c, text) { const gg = c.prep(), T = theme(); gg.ctx.clearRect(0, 0, gg.w, gg.h); gg.ctx.fillStyle = T.bg; gg.ctx.fillRect(0, 0, gg.w, gg.h); gg.ctx.fillStyle = T.muted; gg.ctx.font = '12px ' + T.mono; gg.ctx.textAlign = 'center'; gg.ctx.textBaseline = 'middle'; gg.ctx.fillText(text || 'Fix the inputs', gg.w / 2, gg.h / 2); }
  function arrow(ctx, x0, y0, x1, y1, hs) {
    ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1, y1); ctx.stroke();
    const a = Math.atan2(y1 - y0, x1 - x0), s = hs || 6; ctx.beginPath(); ctx.moveTo(x1, y1); ctx.lineTo(x1 - s * Math.cos(a - 0.4), y1 - s * Math.sin(a - 0.4)); ctx.lineTo(x1 - s * Math.cos(a + 0.4), y1 - s * Math.sin(a + 0.4)); ctx.closePath(); ctx.fill();
  }

  /* ---------- small DOM helpers ---------- */
  const el = function () { return FSP.ui.el.apply(null, arguments); };
  function layout(root) { const lay = el('div', { class: 'layout' }), ctl = el('div', { class: 'controls' }), stage = el('div', { class: 'stage' }); lay.appendChild(ctl); lay.appendChild(stage); root.appendChild(lay); return { ctl, stage }; }
  function msgBox(parent) { const m = el('div', { class: 'msg bad', hidden: '', role: 'alert' }); parent.appendChild(m); return { el: m, show(t) { m.hidden = false; m.textContent = t; }, hide() { m.hidden = true; } }; }
  function hud(parent) { const h = el('div', { class: 'hud' }); parent.appendChild(h); return h; }
  function hudSet(h, rows) { while (h.firstChild) h.removeChild(h.firstChild); rows.forEach(r => h.appendChild(el('div', null, el('span', { text: r[0] }), el('span', { text: r[1] })))); }
  function canvasIn(parent, label, height) { const wrap = el('div', { class: 'canvas-wrap' }); parent.appendChild(wrap); const c = FSP.ui.canvas(wrap, { height }); c.cv.setAttribute('role', 'img'); c.cv.setAttribute('aria-label', label); c.wrap = wrap; return c; }
  // sliders from definitions {k,l,u,min,max,step,v,log,m:[modes]}; the same store is shared by every mode of a sub-topic
  function makeSliders(parent, defs, st, onChange) {
    defs.forEach(d => { const s = FSP.ui.slider(parent, { label: d.l, min: d.min, max: d.max, step: d.step, value: d.v, unit: d.u, log: d.log, digits: d.digits, onInput: () => { onChange(); FSP.state.touch(); } }); s.def = d; st[d.k] = s; });
  }
  function modeVis(st, mode) { Object.keys(st).forEach(k => { const m = st[k].def && st[k].def.m; st[k].el.hidden = !!m && m.indexOf(mode) < 0; }); }
  const outVals = st => { const o = {}; Object.keys(st).forEach(k => { const v = st[k].get(); o[k] = fin(v) ? +v.toPrecision(7) : 0; }); return o; };
  function loadSt(st, o) { Object.keys(st).forEach(k => { if (o[k] !== undefined) { const v = parseFloat(o[k]); if (fin(v)) st[k].set(v, true); } }); }
  const oneOf = (v, list) => (list.indexOf(v) >= 0 ? v : null);
  // "σ = 0" checkbox for lossless media (a log slider cannot reach zero)
  function zeroBox(parent, label, checked, onChange) {
    const chk = el('input', { type: 'checkbox', 'aria-label': label }); chk.checked = !!checked;
    chk.addEventListener('change', () => { onChange(); FSP.state.touch(); });
    const row = el('label', { class: 'row' }, chk, el('span', { class: 'note', text: label })); parent.appendChild(row); return { chk, row };
  }
  const mkNote = (parent, text) => parent.appendChild(el('div', { class: 'note', text }));
  const relErr = (a, b) => Math.abs(a - b) / Math.max(Math.abs(b), 1e-300);
  const hyp3 = v => Math.hypot(v[0], v[1], v[2]);

  /* ---------- 1. charge distributions ---------- */
  function buildCoul(root) {
    const M = FSP.math.fields, { ctl, stage } = layout(root), st = {}, S = { mode: 'pts', view: 'lines' }, ui = {}, nC = 1e-9;
    const MODES = [['pts', 'Point charges (up to 4)'], ['linef', 'Finite line charge'], ['linei', 'Infinite line charge'], ['ring', 'Ring (on its axis)'], ['disk', 'Disk (on its axis)'], ['sheeti', 'Infinite sheet'], ['rect', 'Rectangular sheet']];
    const f0 = FSP.ui.fieldset(ctl, 'Source');
    ui.mode = FSP.ui.select(f0, 'Source', MODES, S.mode, v => { S.mode = v; update(); FSP.state.touch(); });
    ui.view = FSP.ui.select(f0, 'Map shows', [['lines', 'Field lines + equipotentials'], ['arrows', 'E arrows + equipotentials']], S.view, v => { S.view = v; drawMap(); FSP.state.touch(); });
    ui.viewRow = ui.view.parentNode; ui.msg = msgBox(f0);
    const f1 = FSP.ui.fieldset(ctl, 'Parameters');
    const D = [{ k: 'n', l: 'Charges', min: 1, max: 4, step: 1, v: 2, m: ['pts'] }];
    const dq = [[2, -1, 0], [-2, 1, 0], [1, 0, 1.5], [-1, 0, -1.5]];
    for (let i = 1; i <= 4; i++) { D.push({ k: 'q' + i, l: 'Q' + i, u: 'nC', min: -10, max: 10, step: 0.1, v: dq[i - 1][0], m: ['pts'] }, { k: 'x' + i, l: 'x' + i, u: 'm', min: -5, max: 5, step: 0.05, v: dq[i - 1][1], m: ['pts'] }, { k: 'y' + i, l: 'y' + i, u: 'm', min: -5, max: 5, step: 0.05, v: dq[i - 1][2], m: ['pts'] }); }
    D.push({ k: 'px', l: 'Probe x', u: 'm', min: -5, max: 5, step: 0.05, v: 0.5, m: ['pts'] }, { k: 'py', l: 'Probe y', u: 'm', min: -5, max: 5, step: 0.05, v: 1, m: ['pts'] }, { k: 'vw', l: 'View ±', u: 'm', min: 0.5, max: 10, step: 0.1, v: 3, m: ['pts'] },
      { k: 'rl', l: 'ρL', u: 'nC/m', min: -10, max: 10, step: 0.1, v: 2, m: ['linef', 'linei', 'ring'] },
      { k: 'rs', l: 'ρs', u: 'nC/m²', min: -10, max: 10, step: 0.1, v: 2, m: ['disk', 'sheeti', 'rect'] },
      { k: 'a', l: 'Radius a', u: 'm', min: 0.01, max: 10, log: true, v: 1, m: ['ring', 'disk'] },
      { k: 'z1', l: 'z1 (end)', u: 'm', min: -5, max: 5, step: 0.05, v: -1, m: ['linef'] }, { k: 'z2', l: 'z2 (end)', u: 'm', min: -5, max: 5, step: 0.05, v: 2, m: ['linef'] }, { k: 'zp', l: 'P at z', u: 'm', min: -5, max: 5, step: 0.05, v: 0.5, m: ['linef'] },
      { k: 'rx1', l: 'x1', u: 'm', min: -3, max: 3, step: 0.05, v: -1, m: ['rect'] }, { k: 'rx2', l: 'x2', u: 'm', min: -3, max: 3, step: 0.05, v: 1, m: ['rect'] }, { k: 'ry1', l: 'y1', u: 'm', min: -3, max: 3, step: 0.05, v: -1, m: ['rect'] }, { k: 'ry2', l: 'y2', u: 'm', min: -3, max: 3, step: 0.05, v: 1, m: ['rect'] },
      { k: 'rpx', l: 'P at x', u: 'm', min: -3, max: 3, step: 0.05, v: 0, m: ['rect'] }, { k: 'rpy', l: 'P at y', u: 'm', min: -3, max: 3, step: 0.05, v: 0, m: ['rect'] },
      { k: 'd', l: 'Distance d', u: 'm', min: 0.01, max: 20, log: true, v: 1, m: ['linef', 'linei', 'ring', 'disk', 'sheeti', 'rect'] },
      { k: 'dmax', l: 'Plot to', u: 'm', min: 0.1, max: 100, log: true, v: 10, m: ['linef', 'linei', 'ring', 'disk', 'sheeti', 'rect'] });
    makeSliders(f1, D, st, () => update());
    const note = el('div', { class: 'note' }); f1.appendChild(note);
    ui.note = note;
    // stage
    ui.mapWrap = el('div'); stage.appendChild(ui.mapWrap); ui.contWrap = el('div', { class: 'stage' }); stage.appendChild(ui.contWrap);
    ui.map = canvasIn(ui.mapWrap, 'Electric field map of the point charges in the z = 0 plane', 380); ui.map.onResize(() => drawMap());
    mkNote(ui.mapWrap, 'Plane z = 0 (x right, y up). Orange contours: V > 0, cyan: V < 0, grey: V = 0; levels are 0, ±V₀·2ⁿ. Click or drag on the map to move the probe P.');
    ui.cv = canvasIn(ui.contWrap, '|E| against distance, closed form and numerical integration', 300); ui.cv.onResize(() => drawCont());
    ui.hud = hud(stage); ui.work = FSP.ui.working(stage);
    const mapCv = ui.map.cv; let dragging = false;
    function pick(ev) {
      if (!ui.mapT) return; const r = mapCv.getBoundingClientRect(), x = ui.mapT.wx(ev.clientX - r.left), y = ui.mapT.wy(ev.clientY - r.top);
      st.px.set(x); st.py.set(y);
    }
    mapCv.addEventListener('pointerdown', ev => { dragging = true; try { mapCv.setPointerCapture(ev.pointerId); } catch (e) { /* ignore */ } pick(ev); });
    mapCv.addEventListener('pointermove', ev => { if (dragging) pick(ev); });
    mapCv.addEventListener('pointerup', () => { dragging = false; }); mapCv.addEventListener('pointercancel', () => { dragging = false; });

    let mapData = null, contData = null, dotTimer = 0; const dotCache = {};
    function chargesNow() { const n = Math.round(st.n.get()), qs = []; for (let i = 1; i <= n; i++) qs.push({ q: st['q' + i].get() * nC, x: st['x' + i].get(), y: st['y' + i].get(), z: 0, i }); return qs; }
    function params() {
      const v = {}; Object.keys(st).forEach(k => { v[k] = st[k].get(); });
      return { rl: v.rl * nC, rs: v.rs * nC, a: v.a, z1: v.z1, z2: v.z2, zp: v.zp, x1: v.rx1, x2: v.rx2, y1: v.ry1, y2: v.ry2, px: v.rpx, py: v.rpy, d: v.d, dmax: v.dmax };
    }
    // closed form at distance d -> {E:[3], V}; numeric -> {E:[3], V, evals}
    function closed(mode, p, d) {
      switch (mode) {
        case 'linef': { const r = M.lineFinite(p.rl, p.z1, p.z2, [d, 0, p.zp]); return r.error ? r : { E: r.E, V: r.V, r }; }
        case 'linei': { const r = M.lineInfinite(p.rl, [d, 0, 0], 1); return { E: r.E, V: r.V, r }; }
        case 'ring': { const r = M.ringAxis(p.rl, p.a, d); return { E: [0, 0, r.Ez], V: r.V }; }
        case 'disk': { const r = M.diskAxis(p.rs, p.a, d); return { E: [0, 0, r.Ez], V: r.V }; }
        case 'sheeti': { const r = M.sheetInf(p.rs, d); return { E: [0, 0, r.Ez], V: r.V }; }
        default: { const r = M.rectSheet(p.rs, p.x1, p.x2, p.y1, p.y2, [p.px, p.py, d]); return r.error ? r : { E: r.E, V: r.V, r }; }
      }
    }
    function numeric(mode, p, d, tol, bud) {
      switch (mode) {
        case 'linef': return M.lineFiniteNum(p.rl, p.z1, p.z2, [d, 0, p.zp], tol, bud);
        case 'linei': return M.lineInfiniteNum(p.rl, [d, 0, 0], bud);
        case 'ring': return M.ringNum(p.rl, p.a, [0, 0, d], tol, bud);
        case 'disk': return M.diskNum(p.rs, p.a, [0, 0, d], tol, bud);
        case 'sheeti': return M.sheetInfNum(p.rs, d, tol, bud);
        default: return M.rectSheetNum(p.rs, p.x1, p.x2, p.y1, p.y2, [p.px, p.py, d], tol, bud);
      }
    }
    function idealised(mode, p, d) { // the infinite-line / infinite-sheet limit drawn dashed for comparison
      if (mode === 'linef') return Math.abs(p.rl) / (2 * PI * M.EPS0 * d);
      if (mode === 'disk' || mode === 'rect') return Math.abs(p.rs) / (2 * M.EPS0);
      return NaN;
    }

    function update() {
      modeVis(st, S.mode); const pts = S.mode === 'pts';
      const n = Math.round(st.n.get()); for (let i = 1; i <= 4; i++) ['q', 'x', 'y'].forEach(c => { st[c + i].el.hidden = !pts || i > n; });
      ui.viewRow.hidden = !pts; ui.mapWrap.hidden = !pts; ui.contWrap.hidden = pts;
      ui.note.textContent = pts ? 'Charges lie in the z = 0 plane. Q in nC. Potential reference: V = 0 at infinity.' :
        { linef: 'Line charge on the z-axis from z1 to z2; P = (ρ = d, φ = 0, z = "P at z"). The plot sweeps ρ at that z.', linei: 'Infinite line charge on the z-axis; P at ρ = d. V is referenced to ρ = 1 m.', ring: 'Ring of radius a in the z = 0 plane, centred on the z-axis; P on the axis at z = d.', disk: 'Disk of radius a in z = 0; P on the axis at z = d.', sheeti: 'Infinite sheet in z = 0; P at z = d. V is referenced to the sheet (V = 0 at z = 0).', rect: 'Rectangle x1<x<x2, y1<y<y2 in z = 0; P = (P at x, P at y, z = d). The plot sweeps z.' }[S.mode];
      if (pts) updatePts(); else updateCont();
    }

    /* ----- point charges ----- */
    function updatePts() {
      const qs = chargesNow(), P = [st.px.get(), st.py.get(), 0], r = M.pointCharges(qs, P);
      mapData = { qs, P, res: r };
      ui.msg.hide();
      if (qs.every(c => c.q === 0)) ui.msg.show('All charges are zero: there is no field.');
      if (r.error) { ui.msg.show('Probe P: ' + r.error + '. Move the probe.'); hudSet(ui.hud, [['Probe', 'on a charge']]); ui.work.set('Move the probe away from the charges.'); drawMap(); return; }
      const E = r.E, Em = hyp3(E);
      hudSet(ui.hud, [['Ex', fe(E[0], 'V/m')], ['Ey', fe(E[1], 'V/m')], ['|E|', fe(Em, 'V/m')], ['Angle of E', g(Math.atan2(E[1], E[0]) * R2D, 4) + '°'], ['V at P', fe(r.V, 'V')], ['Charges', String(qs.length)]]);
      const L = ['Coulomb superposition (Hayt): E = Σ Qi (r − ri) / (4πε0 |r − ri|³),  V = Σ Qi / (4πε0 |r − ri|),  1/(4πε0) = ' + g(M.KE) + ' m/F', 'Probe P = (' + g(P[0]) + ', ' + g(P[1]) + ', 0) m', ''];
      r.parts.forEach((pt, i) => {
        const c = qs[i];
        L.push('Charge ' + (i + 1) + ':  Q = ' + g(c.q / nC) + ' nC at (' + g(c.x) + ', ' + g(c.y) + ') m');
        L.push('   R = r − r' + (i + 1) + ' = ' + vec3(pt.R) + ' m,   |R| = ' + g(pt.r) + ' m');
        L.push('   Ei = kQ R/|R|³ = ' + g(M.KE * c.q / (pt.r * pt.r * pt.r)) + ' × ' + vec3(pt.R) + ' = ' + vec3(pt.Ei) + ' V/m');
        L.push('   Vi = kQ/|R| = ' + g(pt.Vi) + ' V');
      });
      L.push('', 'Sum:  E = ' + vec3(E) + ' V/m,   |E| = ' + g(Em) + ' V/m at ' + g(Math.atan2(E[1], E[0]) * R2D, 4) + '° from +x', '      V = ' + g(r.V) + ' V');
      if (qs.length === 2 && Math.abs(qs[0].q + qs[1].q) < 1e-18 * 1e3 && qs[0].q !== 0) L.push('(Equal and opposite charges: this is a dipole; V = 0 on the perpendicular bisector plane.)');
      ui.work.set(L); drawMap();
    }
    function contours(ctx, V, nx, ny, cell, levels, colorFor) {
      const tbl = { 1: [[3, 2]], 2: [[2, 1]], 3: [[3, 1]], 4: [[0, 1]], 5: [[0, 3], [1, 2]], 6: [[0, 2]], 7: [[0, 3]], 8: [[0, 3]], 9: [[0, 2]], 10: [[0, 1], [2, 3]], 11: [[0, 1]], 12: [[3, 1]], 13: [[2, 1]], 14: [[3, 2]] };
      levels.forEach(Lv => {
        ctx.strokeStyle = colorFor(Lv); ctx.lineWidth = Lv === 0 ? 1.4 : 1; ctx.setLineDash(Lv === 0 ? [5, 3] : []); ctx.beginPath();
        for (let j = 0; j < ny - 1; j++) for (let i = 0; i < nx - 1; i++) {
          const a = V[j * nx + i], b = V[j * nx + i + 1], c = V[(j + 1) * nx + i + 1], d = V[(j + 1) * nx + i];
          if (!(fin(a) && fin(b) && fin(c) && fin(d))) continue;
          const idx = (a > Lv ? 8 : 0) | (b > Lv ? 4 : 0) | (c > Lv ? 2 : 0) | (d > Lv ? 1 : 0); if (idx === 0 || idx === 15) continue;
          const x0 = i * cell, y0 = j * cell;
          const pt = e => { switch (e) { case 0: return [x0 + cell * (Lv - a) / (b - a), y0]; case 1: return [x0 + cell, y0 + cell * (Lv - b) / (c - b)]; case 2: return [x0 + cell * (Lv - d) / (c - d), y0 + cell]; default: return [x0, y0 + cell * (Lv - a) / (d - a)]; } };
          tbl[idx].forEach(sg => { const p = pt(sg[0]), q = pt(sg[1]); ctx.moveTo(p[0], p[1]); ctx.lineTo(q[0], q[1]); });
        }
        ctx.stroke();
      }); ctx.setLineDash([]); ctx.lineWidth = 1;
    }
    function drawMap() {
      const gg = ui.map.prep(), ctx = gg.ctx, w = gg.w, h = gg.h, T = theme();
      ctx.clearRect(0, 0, w, h); ctx.fillStyle = T.bg; ctx.fillRect(0, 0, w, h);
      if (!mapData || S.mode !== 'pts') return;
      const qs = mapData.qs, vw = st.vw.get(), sc = w / (2 * vw), X = x => w / 2 + x * sc, Y = y => h / 2 - y * sc, wx = px => (px - w / 2) / sc, wy = py => (h / 2 - py) / sc;
      ui.mapT = { wx, wy };
      const hy = h / (2 * sc);
      // grid
      ctx.font = '10px ' + T.mono; ctx.lineWidth = 1; ctx.fillStyle = T.muted;
      niceTicks(-vw, vw, Math.max(3, Math.floor(w / 80))).forEach(v => { ctx.strokeStyle = T.border; ctx.globalAlpha = v === 0 ? 0.9 : 0.4; ctx.beginPath(); ctx.moveTo(X(v), 0); ctx.lineTo(X(v), h); ctx.stroke(); ctx.globalAlpha = 1; ctx.textAlign = 'center'; ctx.textBaseline = 'bottom'; ctx.fillText(tickStr(v), X(v), h - 2); });
      niceTicks(-hy, hy, Math.max(3, Math.floor(h / 80))).forEach(v => { ctx.strokeStyle = T.border; ctx.globalAlpha = v === 0 ? 0.9 : 0.4; ctx.beginPath(); ctx.moveTo(0, Y(v)); ctx.lineTo(w, Y(v)); ctx.stroke(); ctx.globalAlpha = 1; ctx.textAlign = 'left'; ctx.textBaseline = 'middle'; ctx.fillText(tickStr(v), 2, Y(v)); });
      const qmax = qs.reduce((m, c) => Math.max(m, Math.abs(c.q)), 0);
      if (qmax > 0) {
        // equipotentials
        const cell = 6, nx = Math.ceil(w / cell) + 1, ny = Math.ceil(h / cell) + 1, V = new Float64Array(nx * ny);
        for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) { const v = M.pointEV(qs, wx(i * cell), wy(j * cell), 0)[3]; V[j * nx + i] = fin(v) ? v : NaN; }
        const Vref = M.KE * qmax / (vw * 0.5), levels = [0]; for (let k = -3; k <= 3; k++) levels.push(Vref * Math.pow(2, k), -Vref * Math.pow(2, k));
        ctx.save(); ctx.globalAlpha = 0.85; contours(ctx, V, nx, ny, cell, levels, Lv => (Lv === 0 ? T.muted : Lv > 0 ? T.a : T.c)); ctx.restore();
        // field representation
        const Eat = (x, y) => { const e = M.pointEV(qs, x, y, 0); return [e[0], e[1]]; };
        ctx.strokeStyle = T.text; ctx.fillStyle = T.text;
        if (S.view === 'arrows') {
          const nxA = Math.max(8, Math.floor(w / 40)), sp = w / nxA, items = []; let lo = Infinity, hi = -Infinity;
          for (let py = sp / 2; py < h; py += sp) for (let px = sp / 2; px < w; px += sp) {
            const e = Eat(wx(px), wy(py)), m = Math.hypot(e[0], e[1]); if (!(m > 0) || !fin(m)) continue;
            if (qs.some(c => Math.hypot(wx(px) - c.x, wy(py) - c.y) < 0.12 * vw)) continue;
            const lm = Math.log10(m); lo = Math.min(lo, lm); hi = Math.max(hi, lm); items.push([px, py, e[0] / m, -e[1] / m, lm]);
          }
          items.forEach(it => { ctx.globalAlpha = 0.25 + 0.75 * (hi > lo ? (it[4] - lo) / (hi - lo) : 1); const L2 = sp * 0.36; arrow(ctx, it[0] - it[2] * L2, it[1] - it[3] * L2, it[0] + it[2] * L2, it[1] + it[3] * L2, 5); });
          ctx.globalAlpha = 1;
        } else {
          const pos = qs.filter(c => c.q > 0), src = pos.length ? pos : qs.filter(c => c.q < 0), dir = pos.length ? 1 : -1, ds = vw * 0.012, lim = [vw * 1.15, hy * 1.15];
          ctx.lineWidth = 1.2; ctx.globalAlpha = 0.8;
          src.forEach(c => {
            const cnt = Math.max(6, Math.round(16 * Math.abs(c.q) / qmax));
            for (let k = 0; k < cnt; k++) {
              const ang = (k + 0.5) / cnt * 2 * PI; let x = c.x + 2.5 * ds * Math.cos(ang), y = c.y + 2.5 * ds * Math.sin(ang); const path = [[x, y]];
              for (let s = 0; s < 2500; s++) {
                const e = Eat(x, y), m = Math.hypot(e[0], e[1]); if (!(m > 0) || !fin(m)) break;
                const xm = x + dir * e[0] / m * ds / 2, ym = y + dir * e[1] / m * ds / 2, e2 = Eat(xm, ym), m2 = Math.hypot(e2[0], e2[1]); if (!(m2 > 0) || !fin(m2)) break;
                x += dir * e2[0] / m2 * ds; y += dir * e2[1] / m2 * ds; path.push([x, y]);
                if (Math.abs(x) > lim[0] || Math.abs(y) > lim[1]) break;
                if (qs.some(o => o !== c && Math.hypot(x - o.x, y - o.y) < 1.6 * ds)) break;
              }
              ctx.beginPath(); path.forEach((p, i) => { if (i) ctx.lineTo(X(p[0]), Y(p[1])); else ctx.moveTo(X(p[0]), Y(p[1])); }); ctx.stroke();
              let acc = 0, next = 60; for (let i = 1; i < path.length; i++) {
                acc += Math.hypot(path[i][0] - path[i - 1][0], path[i][1] - path[i - 1][1]) * sc;
                if (acc >= next) { next += 140; const x1 = X(path[i][0]), y1 = Y(path[i][1]), x0 = X(path[i - 1][0]), y0 = Y(path[i - 1][1]); const a = Math.atan2(y1 - y0, x1 - x0) + (dir < 0 ? PI : 0); ctx.beginPath(); ctx.moveTo(x1 + 5 * Math.cos(a), y1 + 5 * Math.sin(a)); ctx.lineTo(x1 - 4 * Math.cos(a - 0.5), y1 - 4 * Math.sin(a - 0.5)); ctx.lineTo(x1 - 4 * Math.cos(a + 0.5), y1 - 4 * Math.sin(a + 0.5)); ctx.closePath(); ctx.fill(); }
              }
            }
          });
          ctx.globalAlpha = 1;
        }
      }
      // charges
      qs.forEach(c => { const rr = 8 + Math.min(8, 6 * Math.abs(c.q) / (Math.max(1e-30, qs.reduce((m, o) => Math.max(m, Math.abs(o.q)), 0)))); ctx.fillStyle = c.q > 0 ? T.a : c.q < 0 ? T.c : T.muted; ctx.beginPath(); ctx.arc(X(c.x), Y(c.y), rr, 0, 2 * PI); ctx.fill(); ctx.strokeStyle = T.bg; ctx.lineWidth = 1.5; ctx.stroke(); ctx.fillStyle = T.bg; ctx.font = 'bold 12px ' + T.mono; ctx.textAlign = 'center'; ctx.textBaseline = 'middle'; ctx.fillText(c.q > 0 ? '+' : c.q < 0 ? '−' : '0', X(c.x), Y(c.y) + 0.5); ctx.fillStyle = T.text; ctx.font = '10px ' + T.mono; ctx.fillText('Q' + c.i, X(c.x), Y(c.y) - rr - 7); });
      // probe
      const P = mapData.P, px = X(P[0]), py = Y(P[1]);
      ctx.strokeStyle = T.pink; ctx.fillStyle = T.pink; ctx.lineWidth = 2; ctx.beginPath(); ctx.moveTo(px - 6, py); ctx.lineTo(px + 6, py); ctx.moveTo(px, py - 6); ctx.lineTo(px, py + 6); ctx.stroke();
      ctx.font = 'bold 11px ' + T.mono; ctx.textAlign = 'left'; ctx.textBaseline = 'bottom'; ctx.fillText('P', px + 7, py - 5);
      if (mapData.res && mapData.res.E) { const E = mapData.res.E, m = Math.hypot(E[0], E[1]); if (m > 0) { ctx.lineWidth = 2.5; arrow(ctx, px, py, px + E[0] / m * 46, py - E[1] / m * 46, 9); } }
    }

    /* ----- continuous distributions ----- */
    function updateCont() {
      const mode = S.mode, p = params(); ui.msg.hide(); contData = null;
      let bad = null;
      if (mode === 'linef' && !(p.z2 > p.z1)) bad = 'Finite line needs z2 > z1.';
      if (mode === 'rect' && !(p.x2 > p.x1 && p.y2 > p.y1)) bad = 'Rectangle needs x2 > x1 and y2 > y1.';
      if (!bad && !(p.d > 0 && p.dmax > 0)) bad = 'Distance d and the plot range must be positive.';
      if (bad) { ui.msg.show(bad); blank(ui.cv, 'Fix the inputs'); hudSet(ui.hud, []); ui.work.set('Fix the input above.'); return; }
      const c = closed(mode, p, p.d); if (c.error) { ui.msg.show(c.error); blank(ui.cv, 'Fix the inputs'); hudSet(ui.hud, []); ui.work.set('Fix the input above.'); return; }
      const nm = numeric(mode, p, p.d, 1e-9, { left: 400000 }), Ec = hyp3(c.E), En = hyp3(nm.E), diff = relErr(En, Ec);
      const rows = [['E closed form', '(' + c.E.map(x => fe(x, '', 4).trim()).join(', ') + ') V/m'], ['|E| closed form', fe(Ec, 'V/m')], ['|E| numerical ∫', fe(En, 'V/m')], ['Relative difference', Ec === 0 ? (En < 1e-9 ? '0' : '—') : diff.toExponential(2)]];
      if (c.V !== undefined && fin(c.V)) rows.push(['V closed form', fe(c.V, 'V')]);
      if (nm.V !== undefined && fin(nm.V) && c.V !== undefined && mode !== 'sheeti' && mode !== 'linei') rows.push(['V numerical ∫', fe(nm.V, 'V')]);
      if (nm.capped) rows.push(['Numeric note', 'evaluation cap reached: precision reduced']);
      hudSet(ui.hud, rows);
      ui.work.set(workCont(mode, p, c, nm, Ec, En));
      // curves
      const xmin = p.dmax / 1000, xs = logspace(xmin, p.dmax, 140), ys = xs.map(d => { const r = closed(mode, p, d); return r.error ? NaN : hyp3(r.E); }), ideal = xs.map(d => idealised(mode, p, d));
      contData = { mode, p, xs, ys, ideal, xmin, Ec, En, dots: null };
      const key = JSON.stringify([mode, p.rl, p.rs, p.a, p.z1, p.z2, p.zp, p.x1, p.x2, p.y1, p.y2, p.px, p.py, p.dmax]);
      if (dotCache.key === key) contData.dots = dotCache.dots;
      else { clearTimeout(dotTimer); dotTimer = setTimeout(() => { if (!contData || contData.mode !== mode) return; const dx = logspace(xmin * 3, p.dmax, 9), dy = dx.map(d => { try { return hyp3(numeric(mode, p, d, 1e-6, { left: 150000 }).E); } catch (e) { return NaN; } }); dotCache.key = key; dotCache.dots = { x: dx, y: dy }; if (contData) { contData.dots = dotCache.dots; drawCont(); } }, 90); }
      drawCont();
    }
    function drawCont() {
      if (!contData || S.mode === 'pts') return;
      const D2 = contData, p = D2.p, series = [], all = D2.ys.concat(D2.ideal, [D2.Ec]);
      const yr = decadeRange(all);
      if (D2.ideal.some(fin)) series.push({ x: D2.xs, y: D2.ideal, color: theme().muted, dash: [5, 4], width: 1.5 });
      series.push({ x: D2.xs, y: D2.ys, color: theme().c, width: 2 });
      if (D2.dots) series.push({ x: D2.dots.x, y: D2.dots.y, color: theme().a, dots: true, width: 0.001 });
      const lg = [{ text: '— closed form', color: theme().c }, { text: '● numerical integral', color: theme().a }];
      if (D2.ideal.some(fin)) lg.push({ text: '- - ' + (S.mode === 'linef' ? 'infinite line' : 'infinite sheet'), color: theme().muted });
      plot(ui.cv, { xmin: D2.xmin, xmax: p.dmax, ymin: yr[0], ymax: yr[1], xlog: true, ylog: true, xlabel: (S.mode === 'linef' || S.mode === 'linei' ? 'ρ' : 'z') + ' (m)', ylabel: '|E| (V/m)', series, legend: lg, vlines: [{ x: p.d, color: theme().pink, label: 'd' }], points: fin(D2.Ec) && D2.Ec > 0 ? [{ x: p.d, y: D2.Ec, color: theme().pink }] : [] });
    }
    function workCont(mode, p, c, nm, Ec, En) {
      const L = [], k = M.KE, cmp = ['Numerical check (adaptive Simpson on the same integrand, no symmetry used): E = ' + vec3(nm.E) + ' V/m,  |E| = ' + g(En) + ' V/m', '   relative difference from the closed form = ' + (Ec > 0 ? relErr(En, Ec).toExponential(2) : 'n/a') + (nm.evals !== undefined ? '   (' + nm.evals + ' integrand evaluations)' : '')];
            if (mode === 'linef' || mode === 'linei') {
        const rho = p.d, ZZ = mode === 'linef' ? p.zp : 0;
        if (mode === 'linef') {
          const r = c.r;
          L.push('Finite line charge ρL = ' + g(p.rl / nC) + ' nC/m on the z-axis, z\' from ' + g(p.z1) + ' to ' + g(p.z2) + ' m.  Field point P(ρ = ' + g(rho) + ', 0, z = ' + g(ZZ) + ') m.', '',
            'Step 1 - source element:  dQ = ρL dz\'   at   r\' = z\' az',
            'Step 2 - R vector:       R = r − r\' = ρ aρ + (z − z\') az,   |R| = √(ρ² + (z − z\')²)',
            'Step 3 - element field:  dE = dQ R /(4πε0 |R|³) = ρL dz\' [ρ aρ + (z − z\') az] / (4πε0 [ρ² + (z − z\')²]^(3/2))',
            'Step 4 - symmetry:       the line is axially symmetric, so E has no aφ component; only Eρ and Ez survive.' + (Math.abs(ZZ - (p.z1 + p.z2) / 2) < 1e-9 ? ' P is level with the midpoint, so the az parts of dE from the upper and lower halves cancel: Ez = 0.' : ' P is not level with the midpoint, so Ez does not cancel.'),
            'Step 5 - integrate over z\' (let u = z\' − z):',
            '   Eρ = (ρL ρ/4πε0) ∫ dz\'/[ρ² + (z\' − z)²]^(3/2) = (ρL/4πε0 ρ) [ (z\' − z)/√(ρ² + (z\' − z)²) ] from z1 to z2',
            '   Ez = (ρL/4πε0) ∫ (z − z\') dz\'/[…]^(3/2) = (ρL/4πε0) [ 1/R2 − 1/R1 ]',
            '   with  sin α = (z\' − z)/R  at each end (α measured at P from the perpendicular to the line, positive toward +z):',
            '   Eρ = ρL/(4πε0 ρ) (sin α2 − sin α1),   Ez = ρL/(4πε0) (1/R2 − 1/R1)', '',
            'Numbers:  k = 1/(4πε0) = ' + g(k) + ' m/F;  kρL = ' + g(k * p.rl) + ' V',
            '   R1 = √(' + g(rho) + '² + (' + g(p.z1 - ZZ) + ')²) = ' + g(r.R1) + ' m,   R2 = √(' + g(rho) + '² + (' + g(p.z2 - ZZ) + ')²) = ' + g(r.R2) + ' m',
            '   sin α1 = ' + g(p.z1 - ZZ) + '/' + g(r.R1) + ' = ' + g(r.sinA1) + ' (α1 = ' + g(Math.asin(r.sinA1) * R2D, 4) + '°),   sin α2 = ' + g(p.z2 - ZZ) + '/' + g(r.R2) + ' = ' + g(r.sinA2) + ' (α2 = ' + g(Math.asin(r.sinA2) * R2D, 4) + '°)',
            '   Eρ = ' + g(k * p.rl) + '/' + g(rho) + ' × (' + g(r.sinA2) + ' − (' + g(r.sinA1) + ')) = ' + g(r.Erho) + ' V/m',
            '   Ez = ' + g(k * p.rl) + ' × (1/' + g(r.R2) + ' − 1/' + g(r.R1) + ') = ' + g(r.Ez) + ' V/m',
            '   |E| = √(Eρ² + Ez²) = ' + g(Ec) + ' V/m', '   V = kρL [asinh((z2 − z)/ρ) − asinh((z1 − z)/ρ)] = ' + g(c.V) + ' V   (zero at infinity)',
            'Limit check:  z1 → −∞, z2 → +∞ gives sin α2 − sin α1 → 2, so E → ρL/(2πε0 ρ) = ' + g(Math.abs(p.rl) / (2 * PI * M.EPS0 * rho)) + ' V/m (dashed curve).', '');
        } else {
          const r = c.r;
          L.push('Infinite line charge ρL = ' + g(p.rl / nC) + ' nC/m along the z-axis;  P at ρ = ' + g(rho) + ' m.', '',
            'Step 1:  dQ = ρL dz\',   R = ρ aρ + (z − z\') az,   |R| = √(ρ² + (z − z\')²)',
            'Step 2:  dE = ρL dz\' [ρ aρ + (z − z\') az] / (4πε0 |R|³)',
            'Step 3 - symmetry: elements at z + s and z − s have equal and opposite az parts, so Ez = 0 (equivalently α1 = −α2).',
            'Step 4:  Eρ = (ρL ρ/4πε0) ∫ from −∞ to ∞ dz\'/[ρ² + (z\' − z)²]^(3/2) = (ρL ρ/4πε0)(2/ρ²) = ρL/(2πε0 ρ)', '',
            'Numbers:  Eρ = ' + g(p.rl) + ' / (2π × ' + g(M.EPS0) + ' × ' + g(rho) + ') = ' + g(r.Erho) + ' V/m',
            '   V(ρ) = (ρL/2πε0) ln(ρref/ρ) = ' + g(c.V) + ' V with ρref = 1 m (V cannot be referenced to infinity for an infinite line).', '');
        }
      } else if (mode === 'ring') {
        const a = p.a, z = p.d, d3 = Math.pow(a * a + z * z, 1.5);
        L.push('Ring of radius a = ' + g(a) + ' m, ρL = ' + g(p.rl / nC) + ' nC/m in z = 0, centred on the z-axis.  P on the axis at z = ' + g(z) + ' m.', '',
          'Step 1 - source element:  dQ = ρL a dφ\'   at   r\' = a aρ(φ\')',
          'Step 2 - R vector:       R = r − r\' = −a aρ(φ\') + z az,   |R| = √(a² + z²)  (same for every element)',
          'Step 3 - element field:  dE = ρL a dφ\' (−a aρ + z az) / (4πε0 (a² + z²)^(3/2))',
          'Step 4 - symmetry:       the element at φ\' + π has the opposite aρ term, so all radial parts cancel; only Ez survives.',
          'Step 5 - integrate:      Ez = ρL a z /(4πε0 (a² + z²)^(3/2)) ∫0..2π dφ\' = ρL a z /(2ε0 (a² + z²)^(3/2)) = Q z/(4πε0 (a² + z²)^(3/2)),  Q = 2π a ρL', '',
          'Numbers:  (a² + z²)^(3/2) = (' + g(a * a) + ' + ' + g(z * z) + ')^(3/2) = ' + g(d3) + ' m³',
          '   Ez = ' + g(p.rl) + ' × ' + g(a) + ' × ' + g(z) + ' / (2 × ' + g(M.EPS0) + ' × ' + g(d3) + ') = ' + g(c.E[2]) + ' V/m',
          '   V = ρL a /(2ε0 √(a² + z²)) = ' + g(c.V) + ' V', '   Q = ' + g(2 * PI * a * p.rl / nC) + ' nC;  for z ≫ a, Ez → Q/(4πε0 z²) (point charge).', '');
      } else if (mode === 'disk') {
        const a = p.a, z = p.d, R = Math.hypot(a, z);
        L.push('Disk of radius a = ' + g(a) + ' m, ρs = ' + g(p.rs / nC) + ' nC/m² in z = 0.  P on the axis at z = ' + g(z) + ' m.', '',
          'Step 1 - source element:  dQ = ρs ρ\' dρ\' dφ\'   at   r\' = ρ\' aρ(φ\')',
          'Step 2 - R vector:       R = −ρ\' aρ(φ\') + z az,   |R| = √(ρ\'² + z²)',
          'Step 3 - symmetry:       integrating φ\' first (a ring of radius ρ\'), the aρ parts cancel in pairs, leaving only Ez.',
          'Step 4 - ring result:    dEz = ρs ρ\' dρ\' z /(2ε0 (ρ\'² + z²)^(3/2))',
          'Step 5 - integrate ρ\' from 0 to a:',
          '   Ez = (ρs z/2ε0) ∫0..a ρ\' dρ\'/(ρ\'² + z²)^(3/2) = (ρs z/2ε0) [ 1/|z| − 1/√(a² + z²) ] = (ρs/2ε0) [ sgn(z) − z/√(a² + z²) ]', '',
          'Numbers:  √(a² + z²) = √(' + g(a * a) + ' + ' + g(z * z) + ') = ' + g(R) + ' m',
          '   Ez = ' + g(p.rs) + '/(2 × ' + g(M.EPS0) + ') × (1 − ' + g(z) + '/' + g(R) + ') = ' + g(p.rs / (2 * M.EPS0)) + ' × ' + g(1 - z / R) + ' = ' + g(c.E[2]) + ' V/m',
          '   V = (ρs/2ε0)(√(a² + z²) − |z|) = ' + g(c.V) + ' V',
          'Limits:  a → ∞ gives ρs/(2ε0) = ' + g(p.rs / (2 * M.EPS0)) + ' V/m (dashed, infinite sheet);  z ≫ a gives Q/(4πε0 z²).', '');
      } else if (mode === 'sheeti') {
        const z = p.d;
        L.push('Infinite sheet ρs = ' + g(p.rs / nC) + ' nC/m² in z = 0;  P at z = ' + g(z) + ' m.', '',
          'Step 1:  dQ = ρs ρ\' dρ\' dφ\',   R = −ρ\' aρ + z az,   |R| = √(ρ\'² + z²)',
          'Step 2 - symmetry: for every element there is one diametrically opposite, so the aρ parts cancel and E = Ez az.',
          'Step 3:  Ez = (ρs z/2ε0) ∫0..∞ ρ\' dρ\'/(ρ\'² + z²)^(3/2) = (ρs z/2ε0)(1/|z|) = (ρs/2ε0) sgn(z)   (independent of distance!)', '',
          'Numbers:  Ez = ' + g(p.rs) + '/(2 × ' + g(M.EPS0) + ') = ' + g(c.E[2]) + ' V/m,   directed away from the sheet for ρs > 0.',
          '   V = −ρs |z|/(2ε0) = ' + g(c.V) + ' V relative to the sheet (a potential at infinity does not exist).', '');
      } else {
        const r = c.r;
        L.push('Rectangular sheet ρs = ' + g(p.rs / nC) + ' nC/m², x from ' + g(p.x1) + ' to ' + g(p.x2) + ' m, y from ' + g(p.y1) + ' to ' + g(p.y2) + ' m, in z = 0.  P = (' + g(p.px) + ', ' + g(p.py) + ', ' + g(p.d) + ') m.', '',
          'Step 1 - source element:  dQ = ρs dx\' dy\'   at   r\' = x\' ax + y\' ay',
          'Step 2 - R vector:       R = (x − x\') ax + (y − y\') ay + z az,   |R| = √((x − x\')² + (y − y\')² + z²)',
          'Step 3 - element field:  dE = ρs dx\' dy\' R /(4πε0 |R|³)',
          'Step 4 - symmetry:       ' + (Math.abs(p.px - (p.x1 + p.x2) / 2) < 1e-9 && Math.abs(p.py - (p.y1 + p.y2) / 2) < 1e-9 ? 'P lies over the centre, so Ex = Ey = 0 (opposite elements cancel).' : 'P is off the centre, so Ex and Ey do not cancel in general.'),
          'Step 5 - integrate:      Ez = (ρs z/4πε0) ∬ dx\' dy\' /(X² + Y² + z²)^(3/2) with X = x\' − x, Y = y\' − y.  The inner integral gives',
          '   Ez = kρs Σ s·atan( X Y /(z R) ) over the four corners (s = +1 for (x2,y2),(x1,y1); −1 for the other two), R = √(X² + Y² + z²).',
          '   Ex = kρs Σ s·ln(Y + R),   Ey = kρs Σ s·ln(X + R)   (with the corner differences taken in the same order).', '',
          'Numbers:  kρs = ' + g(k * p.rs) + ' V');
        (r.corners || []).forEach((cn, i) => L.push('   corner ' + (i + 1) + ': X = ' + g(cn.X) + ', Y = ' + g(cn.Y) + ', R = ' + g(cn.R) + ', s = ' + (cn.s > 0 ? '+1' : '−1') + ', atan(XY/zR) = ' + g(cn.at) + ' rad'));
        L.push('   Ez = ' + g(k * p.rs) + ' × (Σ s·atan) = ' + g(c.E[2]) + ' V/m,   Ex = ' + vec3(c.E).split(', ')[0].slice(1) + ',   Ey = ' + vec3(c.E).split(', ')[1] + ' V/m,   |E| = ' + g(Ec) + ' V/m',
          '   V = ' + g(c.V) + ' V.   Limits: z → 0⁺ gives ρs/(2ε0) = ' + g(p.rs / (2 * M.EPS0)) + ' V/m (dashed); z ≫ size gives Q/(4πε0 z²) with Q = ρs × area = ' + g(p.rs * (p.x2 - p.x1) * (p.y2 - p.y1) / nC) + ' nC.', '');
      }
      return L.concat(cmp);
    }
    return {
      update, draw() { drawMap(); drawCont(); },
      get: () => Object.assign({ mode: S.mode, view: S.view }, outVals(st)),
      set(o) {
        const m = oneOf(o.mode, MODES.map(x => x[0])); if (m) S.mode = m; const v = oneOf(o.view, ['lines', 'arrows']); if (v) S.view = v;
        loadSt(st, o); ui.mode.value = S.mode; ui.view.value = S.view;
      },
    };
  }

  /* ---------- 2. magnetostatics ---------- */
  function buildMag(root) {
    const M = FSP.math.fields, { ctl, stage } = layout(root), st = {}, S = { mode: 'fil' }, ui = {};
    const MODES = [['fil', 'Finite straight filament'], ['filinf', 'Infinite filament'], ['loop', 'Circular loop (on axis)'], ['sq', 'Square loop (on axis)'], ['sol', 'Finite solenoid (on axis)'], ['coax', 'Coaxial cable (4 regions)'], ['sheet', 'Infinite current sheet']];
    const f0 = FSP.ui.fieldset(ctl, 'Source');
    ui.mode = FSP.ui.select(f0, 'Source', MODES, S.mode, v => { S.mode = v; update(); FSP.state.touch(); });
    ui.msg = msgBox(f0);
    const f1 = FSP.ui.fieldset(ctl, 'Parameters');
    makeSliders(f1, [
      { k: 'I', l: 'Current I', u: 'A', min: -50, max: 50, step: 0.1, v: 5, m: ['fil', 'filinf', 'loop', 'sq', 'sol', 'coax'] },
      { k: 'K', l: 'Sheet K', u: 'A/m', min: -100, max: 100, step: 0.1, v: 10, m: ['sheet'] },
      { k: 'z1', l: 'z1 (end)', u: 'm', min: -5, max: 5, step: 0.05, v: -1, m: ['fil'] }, { k: 'z2', l: 'z2 (end)', u: 'm', min: -5, max: 5, step: 0.05, v: 2, m: ['fil'] }, { k: 'zp', l: 'P at z', u: 'm', min: -5, max: 5, step: 0.05, v: 0.5, m: ['fil'] },
      { k: 'rho', l: 'ρ (P)', u: 'm', min: 0.01, max: 20, log: true, v: 0.5, m: ['fil', 'filinf'] },
      { k: 'dmax', l: 'Plot to ρ', u: 'm', min: 0.1, max: 100, log: true, v: 10, m: ['fil', 'filinf'] },
      { k: 'a', l: 'Radius a', u: 'm', min: 0.01, max: 5, log: true, v: 0.5, m: ['loop'] },
      { k: 'w', l: 'Side w', u: 'm', min: 0.02, max: 10, log: true, v: 1, m: ['sq'] },
      { k: 'N', l: 'Turns N', min: 1, max: 5000, step: 1, v: 200, m: ['sol'] },
      { k: 'L', l: 'Length L', u: 'm', min: 0.05, max: 5, log: true, v: 0.4, m: ['sol'] },
      { k: 'sa', l: 'Radius a', u: 'm', min: 0.005, max: 1, log: true, v: 0.05, m: ['sol'] },
      { k: 'zz', l: 'P at z', u: 'm', min: -5, max: 5, step: 0.01, v: 0.3, m: ['loop', 'sq', 'sol', 'sheet'] },
      { k: 'zsp', l: 'Plot ± z', u: 'm', min: 0.1, max: 10, log: true, v: 2, m: ['loop', 'sq', 'sol', 'sheet'] },
      { k: 'ca', l: 'a (inner)', u: 'mm', min: 0.1, max: 20, log: true, v: 1, m: ['coax'] },
      { k: 'cb', l: 'b (shield in)', u: 'mm', min: 0.1, max: 30, log: true, v: 4, m: ['coax'] },
      { k: 'cc', l: 'c (shield out)', u: 'mm', min: 0.1, max: 40, log: true, v: 5, m: ['coax'] },
      { k: 'cr', l: 'ρ (P)', u: 'mm', min: 0, max: 60, step: 0.05, v: 2, m: ['coax'] },
    ], st, () => update());
    ui.note = el('div', { class: 'note' }); f1.appendChild(ui.note);
    ui.cv = canvasIn(stage, 'H against position, closed form and numerical Biot-Savart', 300); ui.cv.onResize(() => draw());
    ui.hud = hud(stage); ui.work = FSP.ui.working(stage);
    let data = null;

    function pr() { const v = {}; Object.keys(st).forEach(k => { v[k] = st[k].get(); }); return v; }
    // returns {H: signed component, name, num: signed numeric}
    function evalAt(mode, v, x) { // x = the swept coordinate
      switch (mode) {
        case 'fil': { const P = [x, 0, v.zp], c = M.filamentFinite(v.I, v.z1, v.z2, P); return { H: c.error ? NaN : c.Hphi, num: () => M.bsSegment(v.I, [0, 0, v.z1], [0, 0, v.z2], P)[1], r: c }; }
        case 'filinf': { const P = [x, 0, 0], c = M.filamentInf(v.I, P); return { H: c.error ? NaN : c.Hphi, num: () => M.filamentInfNum(v.I, P)[1], r: c }; }
        case 'loop': return { H: M.loopAxis(v.I, v.a, x).Hz, num: () => M.bsCircle(v.I, v.a, [0, 0, x])[2] };
        case 'sq': return { H: M.squareAxis(v.I, v.w, x).Hz, num: () => M.bsPolygon(v.I, M.squarePts(v.w), [0, 0, x])[2] };
        case 'sol': return { H: M.solenoidAxis(v.I, v.N, v.L, v.sa, x).Hz, num: (fast) => (fast ? M.solenoidNum(v.I, v.N, v.L, v.sa, [0, 0, x], 16, 8) : M.solenoidNum(v.I, v.N, v.L, v.sa, [0, 0, x]))[2] };
        case 'coax': { const c = M.coaxH(v.I, v.ca * 1e-3, v.cb * 1e-3, v.cc * 1e-3, x * 1e-3); return { H: c.error ? NaN : c.Hphi, num: () => (x > 0 ? M.coaxIencNum(v.I, v.ca * 1e-3, v.cb * 1e-3, v.cc * 1e-3, x * 1e-3) / (2 * PI * x * 1e-3) : 0), r: c }; }
        default: return { H: M.sheetK(v.K, x).Hy, num: () => M.sheetKNum(v.K, x).Hy };
      }
    }
    const NOTES = {
      fil: 'I flows along +az from z1 to z2; P at φ = 0. Angles α1, α2 are measured at P from the perpendicular (aρ) to each end, positive toward +z.',
      filinf: 'Infinite filament on the z-axis, I along +az; H = I/(2πρ) aφ.',
      loop: 'Loop in the z = 0 plane, current counter-clockwise seen from +z, so Hz > 0 on the axis for I > 0.',
      sq: 'Square loop in z = 0 centred on the z-axis, current counter-clockwise seen from +z.',
      sol: 'Solenoid of N turns, length L, centred at z = 0, current counter-clockwise seen from +z; Hz = (nI/2)(sin α2 − sin α1).',
      coax: 'Inner conductor (radius a) carries +I along +az, uniform density; outer conductor (b to c) carries −I uniformly. Hφ(ρ) = Ienc/(2πρ).',
      sheet: 'Sheet in z = 0 with K = K ax (A/m). H = ½ K × aN, aN pointing from the sheet to P. Hy is shown.',
    };
    function update() {
      const mode = S.mode, v = pr(); modeVis(st, mode); ui.note.textContent = NOTES[mode]; ui.msg.hide(); data = null;
      const x = mode === 'fil' || mode === 'filinf' ? v.rho : mode === 'coax' ? v.cr : v.zz;
      let bad = null;
      if (mode === 'fil' && !(v.z2 > v.z1)) bad = 'Filament needs z2 > z1.';
      if (mode === 'coax' && !(v.ca < v.cb && v.cb < v.cc)) bad = 'Coax needs a < b < c.';
      if (mode === 'sheet' && Math.abs(x) < 1e-9) bad = 'P is on the sheet (z = 0): H jumps by K there. Choose z ≠ 0.';
      if (mode === 'sol' && !(v.N >= 1 && v.L > 0 && v.sa > 0)) bad = 'Solenoid needs N ≥ 1, L > 0, a > 0.';
      if (bad) { ui.msg.show(bad); blank(ui.cv, 'Fix the inputs'); hudSet(ui.hud, []); ui.work.set('Fix the input above.'); return; }
      const ev = evalAt(mode, v, x); if (!fin(ev.H)) { ui.msg.show('Field point is not valid for this geometry.'); blank(ui.cv, 'Fix the inputs'); hudSet(ui.hud, []); ui.work.set('Fix the input above.'); return; }
      const num = ev.num(), diff = Math.abs(ev.H) > 1e-12 ? relErr(num, ev.H) : Math.abs(num - ev.H);
      const comp = { fil: 'Hφ', filinf: 'Hφ', loop: 'Hz', sq: 'Hz', sol: 'Hz', coax: 'Hφ', sheet: 'Hy' }[mode];
      const rows = [[comp + ' closed form', fe(ev.H, 'A/m')], [comp + ' numerical (Biot–Savart / Ampère integral)', fe(num, 'A/m')], ['Relative difference', Math.abs(ev.H) > 1e-12 ? diff.toExponential(2) : 'abs ' + diff.toExponential(2)], ['B = μ0 H', fe(M.MU0 * ev.H, 'T')]];
      if (mode === 'fil') rows.push(['α1, α2', g(ev.r.a1, 4) + '°, ' + g(ev.r.a2, 4) + '°']);
      if (mode === 'coax') rows.push(['Region', String(ev.r.region)], ['Ienc', fe(ev.r.Ienc, 'A')]);
      hudSet(ui.hud, rows);
      ui.work.set(working(mode, v, x, ev, num));
      // curves
      let xs, ys, dots = null, vl = [], hl = [], xlog = false, ylog = false, xl, yl, dashed = null;
      if (mode === 'fil' || mode === 'filinf') {
        xlog = ylog = true; xs = logspace(v.dmax / 1000, v.dmax, 140); ys = xs.map(r => Math.abs(evalAt(mode, v, r).H)); xl = 'ρ (m)'; yl = '|Hφ| (A/m)';
        if (mode === 'fil') dashed = { x: xs, y: xs.map(r => Math.abs(v.I) / (2 * PI * r)) };
        const dx = logspace(v.dmax / 300, v.dmax, 10); dots = { x: dx, y: dx.map(r => Math.abs(evalAt(mode, v, r).num())) };
      } else if (mode === 'coax') {
        const xm = 1.5 * v.cc; xs = linspace(0, xm, 400); ys = xs.map(r => evalAt(mode, v, r).H); xl = 'ρ (mm)'; yl = 'Hφ (A/m)';
        const dx = linspace(xm / 24, xm, 16); dots = { x: dx, y: dx.map(r => evalAt(mode, v, r).num()) };
        vl = [{ x: v.ca, label: 'a' }, { x: v.cb, label: 'b' }, { x: v.cc, label: 'c' }];
      } else {
        const zs = v.zsp; xs = linspace(-zs, zs, 321); if (mode === 'sheet') xs = xs.map(z => (Math.abs(z) < 1e-9 ? NaN : z));
        ys = xs.map(z => (fin(z) ? evalAt(mode, v, z).H : NaN)); xl = 'z (m)'; yl = (mode === 'sheet' ? 'Hy' : 'Hz') + ' (A/m)';
        const dx = linspace(-zs, zs, 13).map(z => (Math.abs(z) < 1e-9 ? 1e-3 * zs : z)); dots = { x: dx, y: dx.map(z => evalAt(mode, v, z).num(true)) };
        if (mode === 'sol') { vl = [{ x: -v.L / 2, label: '−L/2' }, { x: v.L / 2, label: 'L/2' }]; hl = [{ y: v.N / v.L * v.I, label: 'nI (infinite solenoid)' }]; }
      }
      data = { mode, xs, ys, dots, dashed, vl, hl, xlog, ylog, xl, yl, x, H: ev.H };
      draw();
    }
    function draw() {
      if (!data) return; const D2 = data, T = theme(), series = [];
      if (D2.dashed) series.push({ x: D2.dashed.x, y: D2.dashed.y, color: T.muted, dash: [5, 4], width: 1.5 });
      series.push({ x: D2.xs, y: D2.ys, color: T.c, width: 2 });
      if (D2.dots) series.push({ x: D2.dots.x, y: D2.dots.y, color: T.a, dots: true, width: 0.001 });
      const all = D2.ys.concat(D2.dots ? D2.dots.y : [], D2.dashed ? D2.dashed.y : [], D2.hl.map(h => h.y));
      let yr, xr = [D2.xs.filter(fin)[0], D2.xs[D2.xs.length - 1]];
      if (D2.ylog) yr = decadeRange(all); else { yr = linRange(all.concat([0])); }
      const lg = [{ text: '— closed form', color: T.c }, { text: '● numerical', color: T.a }]; if (D2.dashed) lg.push({ text: '- - infinite filament', color: T.muted });
      plot(ui.cv, { xmin: xr[0], xmax: xr[1], ymin: yr[0], ymax: yr[1], xlog: D2.xlog, ylog: D2.ylog, xlabel: D2.xl, ylabel: D2.yl, series, legend: lg, vlines: D2.vl.map(l => ({ x: l.x, label: l.label })).concat([{ x: D2.x, color: T.pink }]), hlines: D2.hl.map(h => ({ y: h.y, label: h.label, color: T.muted })), points: D2.ylog ? (Math.abs(D2.H) > 0 ? [{ x: D2.x, y: Math.abs(D2.H), color: T.pink }] : []) : [{ x: D2.x, y: D2.H, color: T.pink }] });
    }
    function working(mode, v, x, ev, num) {
      const L = [], cmp = ['', 'Numerical cross-check: ' + (mode === 'coax' ? 'Ienc = ∫J·dS (J uniform in each conductor), Hφ = Ienc/(2πρ)' : mode === 'sheet' ? 'the sheet as a sum of infinite filaments K dy\'' : mode === 'filinf' ? 'Biot–Savart integral with z\' = z + ρ tanθ' : 'Biot–Savart integral H = (I/4π)∮ dL × R/|R|³ evaluated by Gauss–Legendre quadrature') + ' → ' + g(num) + ' A/m;  relative difference ' + (Math.abs(ev.H) > 1e-12 ? relErr(num, ev.H).toExponential(2) : 'n/a'), 'Biot–Savart (Hayt): H = (I/4π) ∮ dL × aR / R² .   B = μ0 H = ' + g(M.MU0 * ev.H) + ' T  (μ0 = 4π×10^−7 H/m).'];
      if (mode === 'fil' || mode === 'filinf') {
        const rho = x, I = v.I;
        if (mode === 'fil') {
          const r = ev.r;
          L.push('Finite filament on the z-axis from z1 = ' + g(v.z1) + ' m to z2 = ' + g(v.z2) + ' m, I = ' + g(I) + ' A in +az.  P at ρ = ' + g(rho) + ' m, z = ' + g(v.zp) + ' m.', '',
            'Step 1 - source element:  dL = dz\' az  at  r\' = z\' az',
            'Step 2 - R vector:       R = r − r\' = ρ aρ + (z − z\') az,   |R| = √(ρ² + (z − z\')²)',
            'Step 3 - cross product:  dL × R = dz\' ρ (az × aρ) = ρ dz\' aφ,  so every dH points along +aφ (right-hand rule: thumb along I, fingers along aφ). No cancellation, no other components.',
            'Step 4 - integrate:      H = (I ρ/4π) aφ ∫ dz\'/[ρ² + (z − z\')²]^(3/2) = (I/4πρ) [ (z\' − z)/√(ρ² + (z\' − z)²) ] from z1 to z2  aφ',
            'Angle convention:        α_i is measured at P between the perpendicular to the wire (the aρ direction) and the line from P to each end, positive toward +z, so sin α_i = (z_i − z)/R_i.',
            '                         H = I/(4πρ) (sin α2 − sin α1) aφ.   A point level with the middle has α1 = −α2.', '',
            'Numbers:  R1 = √(' + g(rho) + '² + (' + g(v.z1 - v.zp) + ')²) = ' + g(r.R1) + ' m;   R2 = √(' + g(rho) + '² + (' + g(v.z2 - v.zp) + ')²) = ' + g(r.R2) + ' m',
            '   sin α1 = ' + g(v.z1 - v.zp) + '/' + g(r.R1) + ' = ' + g(r.sinA1) + '  (α1 = ' + g(r.a1, 4) + '°);   sin α2 = ' + g(v.z2 - v.zp) + '/' + g(r.R2) + ' = ' + g(r.sinA2) + '  (α2 = ' + g(r.a2, 4) + '°)',
            '   Hφ = ' + g(I) + '/(4π × ' + g(rho) + ') × (' + g(r.sinA2) + ' − (' + g(r.sinA1) + ')) = ' + g(r.Hphi) + ' A/m   (H = Hφ aφ; at φ = 0 that is along +ay)',
            'Limit: z1 → −∞, z2 → +∞ gives α2 = 90°, α1 = −90°, H → I/(2πρ) = ' + g(Math.abs(I) / (2 * PI * rho)) + ' A/m (dashed curve).');
        } else {
          L.push('Infinite filament along the z-axis, I = ' + g(I) + ' A in +az.  P at ρ = ' + g(rho) + ' m.', '',
            'Step 1:  dL = dz\' az,   R = ρ aρ + (z − z\') az,   dL × R = ρ dz\' aφ  (all elements add along +aφ).',
            'Step 2:  Hφ = (I ρ/4π) ∫ from −∞ to ∞ dz\'/[ρ² + (z − z\')²]^(3/2) = (I ρ/4π)(2/ρ²) = I/(2πρ)    (α2 = 90°, α1 = −90°).',
            '         Ampère check: ∮H·dL = Hφ (2πρ) = I.', '',
            'Numbers:  Hφ = ' + g(I) + '/(2π × ' + g(rho) + ') = ' + g(ev.H) + ' A/m');
        }
      } else if (mode === 'loop') {
        const a = v.a, z = x, d3 = Math.pow(a * a + z * z, 1.5);
        L.push('Circular loop radius a = ' + g(a) + ' m in z = 0, I = ' + g(v.I) + ' A (counter-clockwise from +z).  P on the axis at z = ' + g(z) + ' m.', '',
          'Step 1 - source element:  dL = a dφ\' aφ\'   at  r\' = a aρ(φ\')',
          'Step 2 - R vector:       R = −a aρ(φ\') + z az,   |R| = √(a² + z²)',
          'Step 3 - cross product:  dL × R = a dφ\' (aφ\' × (−a aρ\' + z az)) = a dφ\' (a az + z aρ\')   (since aφ × aρ = −az, aφ × az = aρ)',
          'Step 4 - symmetry:       the aρ\' parts of diametrically opposite elements cancel; only Hz survives.',
          'Step 5 - integrate:      Hz = (I/4π) a² ∫0..2π dφ\' /(a² + z²)^(3/2) = I a²/(2 (a² + z²)^(3/2))', '',
          'Numbers:  (a² + z²)^(3/2) = (' + g(a * a) + ' + ' + g(z * z) + ')^(3/2) = ' + g(d3) + ' m³',
          '   Hz = ' + g(v.I) + ' × ' + g(a * a) + ' / (2 × ' + g(d3) + ') = ' + g(ev.H) + ' A/m.   At the centre: H = I/(2a) = ' + g(v.I / (2 * a)) + ' A/m.');
      } else if (mode === 'sq') {
        const b = v.w / 2, z = x, dp = Math.hypot(b, z), s = b / Math.hypot(Math.SQRT2 * b, z);
        L.push('Square loop of side w = ' + g(v.w) + ' m (b = w/2 = ' + g(b) + ' m) in z = 0, I = ' + g(v.I) + ' A counter-clockwise from +z.  P on the axis at z = ' + g(z) + ' m.', '',
          'Step 1 - each side is a finite filament. The perpendicular distance from P to a side is d = √(b² + z²) = ' + g(dp) + ' m.',
          'Step 2 - the side subtends α = ±α0 about its midpoint foot with sin α0 = b/√(b² + d²) = b/√(2b² + z²) = ' + g(s) + ',  so |H_side| = I/(4π d) × 2 sin α0 = ' + g(v.I / (4 * PI * dp) * 2 * s) + ' A/m.',
          'Step 3 - direction: H_side is perpendicular to the plane through P and the side. Its component along the axis is |H_side| × (b/d); the in-plane components of opposite sides cancel.',
          'Step 4 - four sides:   Hz = 4 × I/(4π d) × 2 sin α0 × (b/d) = 2 I b² / (π (b² + z²) √(2b² + z²))', '',
          'Numbers:  Hz = 2 × ' + g(v.I) + ' × ' + g(b * b) + ' / (π × ' + g(b * b + z * z) + ' × ' + g(Math.sqrt(2 * b * b + z * z)) + ') = ' + g(ev.H) + ' A/m.   At the centre: H = 2√2 I/(π w) = ' + g(2 * Math.SQRT2 * v.I / (PI * v.w)) + ' A/m.');
      } else if (mode === 'sol') {
        const r = M.solenoidAxis(v.I, v.N, v.L, v.sa, x);
        L.push('Solenoid N = ' + g(v.N) + ' turns, length L = ' + g(v.L) + ' m, radius a = ' + g(v.sa) + ' m, I = ' + g(v.I) + ' A.  P on the axis at z = ' + g(x) + ' m (solenoid spans −L/2 to L/2).', '',
          'Step 1 - model:         n = N/L = ' + g(r.n) + ' turns/m; a slice dz\' carries n dz\' turns, i.e. a loop current n I dz\'.',
          'Step 2 - loop result:   dHz = (n I dz\') a² /(2 (a² + (z − z\')²)^(3/2))',
          'Step 3 - integrate:     Hz = (n I/2) ∫ a² dz\'/(a² + (z\' − z)²)^(3/2) = (n I/2) [ (z\' − z)/√(a² + (z\' − z)²) ] from −L/2 to L/2 = (n I/2)(sin α2 − sin α1)',
          '                        with sin α_i = (z_i − z)/√(a² + (z_i − z)²) (same angle convention as the filament, measured from the radial direction).', '',
          'Numbers:  sin α1 = ' + g(-v.L / 2 - x) + '/√(' + g(v.sa * v.sa) + ' + ' + g(Math.pow(-v.L / 2 - x, 2)) + ') = ' + g(r.sinA1) + ';   sin α2 = ' + g(v.L / 2 - x) + '/√(' + g(v.sa * v.sa) + ' + ' + g(Math.pow(v.L / 2 - x, 2)) + ') = ' + g(r.sinA2),
          '   Hz = ' + g(r.n * v.I / 2) + ' × (' + g(r.sinA2) + ' − (' + g(r.sinA1) + ')) = ' + g(ev.H) + ' A/m.   Long-solenoid limit: H = nI = ' + g(r.n * v.I) + ' A/m at the middle (ends: nI/2).');
      } else if (mode === 'coax') {
        const a = v.ca * 1e-3, b = v.cb * 1e-3, c = v.cc * 1e-3, rho = x * 1e-3, r = ev.r, I = v.I;
        const ex = ['0 < ρ < a:  Ienc = I (πρ²)/(πa²) = I ρ²/a²', 'a < ρ < b:  Ienc = I', 'b < ρ < c:  Ienc = I − I (ρ² − b²)/(c² − b²) = I (c² − ρ²)/(c² − b²)', 'ρ > c:  Ienc = I − I = 0'];
        L.push('Coax: inner conductor radius a = ' + g(v.ca) + ' mm carries +I = ' + g(I) + ' A (uniform, along +az); outer conductor b = ' + g(v.cb) + ' mm to c = ' + g(v.cc) + ' mm carries −I (uniform).  P at ρ = ' + g(x) + ' mm.', '',
          'Step 1 - symmetry:       the field is Hφ(ρ) aφ only and constant on a circle of radius ρ.',
          'Step 2 - Ampère\'s law:   ∮H·dL = Hφ 2πρ = Ienc,   so Hφ = Ienc/(2πρ).',
          'Step 3 - enclosed current in each region (J1 = I/πa², J3 = −I/π(c² − b²)):', '   ' + ex.join('\n   '), '',
          'This P is in region ' + r.region + ':  Ienc = ' + g(r.Ienc) + ' A',
          '   Hφ = ' + g(r.Ienc) + '/(2π × ' + g(rho) + ') = ' + g(ev.H) + ' A/m.   (Continuous at ρ = a, b and c; H = 0 outside, the cable is shielded.)',
          '   Peak H at ρ = a: I/(2πa) = ' + g(I / (2 * PI * a)) + ' A/m.');
      } else {
        const z = x;
        L.push('Infinite current sheet in z = 0 with K = ' + g(v.K) + ' A/m along +ax.  P at z = ' + g(z) + ' m.', '',
          'Step 1 - model:         split the sheet into filaments of current K dy\' along x (a strip of width dy\' carries K dy\').',
          'Step 2 - symmetry:      the filaments at +y\' and −y\' give equal and opposite az components of H, so only Hy survives, and it is the same for all positive z.',
          'Step 3 - integrate:     each filament at offset y\' is a distance √(y\'² + z²) away: |dH| = K dy\'/(2π√(y\'² + z²)), and its Hy part is |dH| × |z|/√(y\'² + z²).',
          '                        Hy = −sgn(z) (K/2π) ∫ from −∞ to ∞ |z| dy\'/(y\'² + z²) = −sgn(z) (K/2π)(π) = −sgn(z) K/2',
          'Step 4 - compact form:  H = ½ K × aN, aN = +az above the sheet (ax × az = −ay) and −az below (ax × (−az) = +ay).', '',
          'Numbers:  Hy = −sgn(z) K/2 = ' + (z > 0 ? '−' : '+') + g(v.K) + '/2 = ' + g(ev.H) + ' A/m, independent of |z|.');
      }
      return L.concat(cmp);
    }
    return {
      update, draw,
      get: () => Object.assign({ mode: S.mode }, outVals(st)),
      set(o) { const m = oneOf(o.mode, MODES.map(x => x[0])); if (m) S.mode = m; loadSt(st, o); ui.mode.value = S.mode; },
    };
  }

  /* ---------- 3. Faraday & displacement current ---------- */
  function buildFar(root) {
    const M = FSP.math.fields, { ctl, stage } = layout(root), st = {}, S = { mode: 'fixed' }, ui = {};
    const MODES = [['fixed', 'Fixed loop, B(t) = B0 sin ωt'], ['bar', 'Sliding bar on rails'], ['rot', 'Rotating loop in uniform B'], ['disp', 'Displacement vs conduction current']];
    const f0 = FSP.ui.fieldset(ctl, 'Case');
    ui.mode = FSP.ui.select(f0, 'Case', MODES, S.mode, v => { S.mode = v; update(); FSP.state.touch(); });
    ui.msg = msgBox(f0);
    const f1 = FSP.ui.fieldset(ctl, 'Parameters');
    makeSliders(f1, [
      { k: 'N', l: 'Turns N', min: 1, max: 1000, step: 1, v: 10, m: ['fixed', 'rot'] },
      { k: 'A', l: 'Loop area A', u: 'm²', min: 1e-4, max: 10, log: true, v: 0.01, m: ['fixed', 'rot'] },
      { k: 'B0', l: 'B', u: 'T', min: 1e-3, max: 10, log: true, v: 0.1, m: ['fixed', 'rot'] },
      { k: 'f', l: 'Frequency', u: 'Hz', min: 0.1, max: 1e4, log: true, v: 50, m: ['fixed', 'rot'] },
      { k: 'tt', l: 't / T', min: 0, max: 2, step: 0.01, v: 0.25, m: ['fixed', 'rot'] },
      { k: 'Bb', l: 'B (+az)', u: 'T', min: 0.01, max: 5, log: true, v: 0.5, m: ['bar'] },
      { k: 'l', l: 'Rail gap l', u: 'm', min: 0.01, max: 2, log: true, v: 0.2, m: ['bar'] },
      { k: 'v', l: 'Speed v (+ax)', u: 'm/s', min: -20, max: 20, step: 0.1, v: 3, m: ['bar'] },
      { k: 'R', l: 'Load R', u: 'Ω', min: 0.01, max: 1000, log: true, v: 2, m: ['bar'] },
      { k: 'x0', l: 'x(0)', u: 'm', min: 0.01, max: 5, log: true, v: 0.1, m: ['bar'] },
      { k: 'tsp', l: 'Plot to t', u: 's', min: 0.1, max: 20, log: true, v: 2, m: ['bar'] },
      { k: 'tb', l: 'Time t', u: 's', min: 0, max: 20, step: 0.01, v: 0.5, m: ['bar'] },
      { k: 'sig', l: 'σ', u: 'S/m', min: 1e-15, max: 1e8, log: true, v: 4, m: ['disp'] },
      { k: 'er', l: 'εr', min: 1, max: 1000, log: true, v: 81, m: ['disp'] },
      { k: 'fd', l: 'Frequency', u: 'Hz', min: 1, max: 1e15, log: true, v: 1e6, m: ['disp'] },
      { k: 'E0', l: 'E amplitude', u: 'V/m', min: 1e-3, max: 1e3, log: true, v: 1, m: ['disp'] },
    ], st, () => update());
    ui.note = el('div', { class: 'note' }); f1.appendChild(ui.note);
    ui.cvA = canvasIn(stage, 'Flux or current density plot', 230); ui.cvB = canvasIn(stage, 'Induced emf or current ratio plot', 230);
    ui.cvA.onResize(() => draw()); ui.cvB.onResize(() => draw());
    ui.hud = hud(stage); ui.work = FSP.ui.working(stage);
    let data = null;
    const NOTES = {
      fixed: 'B = B0 sin ωt along +az, perpendicular to a stationary N-turn loop. Positive emf/current = counter-clockwise seen from +z (the reference normal is +az).',
      bar: 'Rails along x, gap l along y, B along +az; the bar at x(t) = x(0) + v t closes the circuit through R. Positive emf = counter-clockwise seen from +z.',
      rot: 'Loop of area A spinning at ω = 2πf in uniform B; Φ = B A cos ωt where ωt is the angle between the loop normal and B. Positive emf = counter-clockwise seen from the tip of the loop normal.',
      disp: 'In a medium with J = σE + jωεE: conduction current σE against displacement current jωεE. Their ratio is the loss tangent σ/(ωε).',
    };
    const tUnit = span => (span < 1e-3 ? [1e-6, 'µs'] : span < 1 ? [1e-3, 'ms'] : [1, 's']);
    function lenz(dphi, what, normal, scale) {
      if (Math.abs(dphi) < 1e-9 * (scale || 1)) return 'dΦ/dt = 0 at this instant (flux is at an extremum or steady), so emf = 0 and there is no induced current.';
      const inc = dphi > 0;
      return 'Lenz: the flux through the loop (along ' + normal + ') is ' + (inc ? 'INCREASING' : 'DECREASING') + (what ? ' (' + what + ')' : '') + ', so the induced current opposes the change: its own field points ' + (inc ? 'against ' : 'along ') + normal + ' inside the loop, which means the current flows ' + (inc ? 'CLOCKWISE' : 'COUNTER-CLOCKWISE') + ' seen from the tip of ' + normal + '.';
    }
    function update() {
      const mode = S.mode; modeVis(st, mode); ui.note.textContent = NOTES[mode]; ui.msg.hide(); data = null;
      const v = {}; Object.keys(st).forEach(k => { v[k] = st[k].get(); });
      if (mode === 'fixed' || mode === 'rot') {
        const T = 1 / v.f, t = v.tt * T, fn = mode === 'fixed' ? M.faradayFixed : M.rotatingLoop, r0 = fn(v.N, v.A, v.B0, v.f, t), w = 2 * PI * v.f;
        const r = { flux: Math.abs(r0.flux) < 1e-12 * v.B0 * v.A ? 0 : r0.flux, emf: Math.abs(r0.emf) < 1e-12 * r0.peak ? 0 : r0.emf, peak: r0.peak };
        const ts = linspace(0, 2 * T, 361), ser = ts.map(x => fn(v.N, v.A, v.B0, v.f, x)), tu = tUnit(2 * T);
        const dphi = mode === 'fixed' ? v.B0 * v.A * w * Math.cos(w * t) : -v.B0 * v.A * w * Math.sin(w * t);
        const peakPhi = v.B0 * v.A;
        data = { kind: 'time', tx: ts.map(x => x / tu[0]), tu: tu[1], A: ser.map(x => x.flux), B: ser.map(x => x.emf), t: t / tu[0], a: r.flux, b: r.emf, ylA: 'Φ per turn (Wb)', ylB: 'emf (V)', xl: 't (' + tu[1] + ')' };
        hudSet(ui.hud, [['ω = 2πf', fe(w, 'rad/s')], ['Period T', fe(T, 's')], ['Φ(t) per turn', fe(r.flux, 'Wb')], ['emf(t)', fe(r.emf, 'V')], ['Peak emf', fe(r.peak, 'V')], ['Peak flux linkage NΦ', fe(v.N * peakPhi, 'Wb')]]);
        const L = mode === 'fixed' ? [
          'Fixed loop: N = ' + g(v.N) + ', A = ' + g(v.A) + ' m², B(t) = B0 sin ωt with B0 = ' + g(v.B0) + ' T, f = ' + g(v.f) + ' Hz.  Take t = ' + g(v.tt) + ' T = ' + g(t) + ' s.', '',
          'Step 1 - angular frequency:  ω = 2πf = ' + g(w) + ' rad/s',
          'Step 2 - flux through one turn (B ⟂ loop, uniform):  Φ(t) = B0 A sin ωt = ' + g(peakPhi) + ' sin ωt Wb;  at t: Φ = ' + g(r.flux) + ' Wb',
          'Step 3 - Faraday:  emf = −N dΦ/dt = −N B0 A ω cos ωt   (reference direction: counter-clockwise seen from +z, i.e. positive normal +az by the right-hand rule)',
          'Step 4 - peak:  emf_max = N B0 A ω = ' + g(v.N) + ' × ' + g(v.B0) + ' × ' + g(v.A) + ' × ' + g(w) + ' = ' + g(r.peak) + ' V',
          'Step 5 - at t:  emf = −' + g(r.peak) + ' × cos(' + g(w * t) + ' rad) = ' + g(r.emf) + ' V',
          'The emf is a cosine, 90° out of phase with the flux: it is largest when the flux is passing through zero (B changing fastest).', '', lenz(dphi, 'B = B0 sin ωt', '+az', v.B0 * v.A * w)] : [
          'Rotating loop: N = ' + g(v.N) + ', A = ' + g(v.A) + ' m², B = ' + g(v.B0) + ' T uniform, f = ' + g(v.f) + ' Hz.  Take t = ' + g(v.tt) + ' T = ' + g(t) + ' s.', '',
          'Step 1 - angular speed:  ω = 2πf = ' + g(w) + ' rad/s;  angle between loop normal and B:  θ = ωt = ' + g(w * t) + ' rad = ' + g((w * t * R2D) % 360, 4) + '°',
          'Step 2 - flux through one turn:  Φ(t) = B A cos ωt;  at t: Φ = ' + g(r.flux) + ' Wb',
          'Step 3 - Faraday:  emf = −N dΦ/dt = N B A ω sin ωt   (reference direction: counter-clockwise seen from the tip of the loop normal)',
          'Step 4 - peak:  emf_max = N B A ω = ' + g(v.N) + ' × ' + g(v.B0) + ' × ' + g(v.A) + ' × ' + g(w) + ' = ' + g(r.peak) + ' V',
          'Step 5 - at t:  emf = ' + g(r.peak) + ' × sin(' + g(w * t) + ' rad) = ' + g(r.emf) + ' V.   The emf is zero when the loop plane is perpendicular to B (flux maximum) and maximum when the plane is parallel to B.', '', lenz(dphi, 'loop turning in B', 'the loop normal a_n', v.B0 * v.A * w)];
        ui.work.set(L);
      } else if (mode === 'bar') {
        let tmax = v.tsp, note = null;
        if (v.v < 0) { const th = v.x0 / -v.v; if (th < tmax) { tmax = th; note = 'The bar reaches x = 0 at t = ' + g(th) + ' s (end of the rails); plot stops there.'; } }
        const tcur = Math.min(v.tb, tmax), r = M.slidingBar(v.Bb, v.l, v.v, v.R, v.x0, tcur);
        if (v.tb > tmax + 1e-12) note = (note || '') + ' Time t is beyond the plotted interval: values shown at t = ' + g(tmax) + ' s.';
        if (note) ui.msg.show(note.trim());
        const ts = linspace(0, tmax, 200), ser = ts.map(x => M.slidingBar(v.Bb, v.l, v.v, v.R, v.x0, x));
        data = { kind: 'time', tx: ts, tu: 's', A: ser.map(x => x.flux), B: ser.map(x => x.emf), t: tcur, a: r.flux, b: r.emf, ylA: 'Φ (Wb)', ylB: 'emf (V)', xl: 't (s)' };
        hudSet(ui.hud, [['Bar position x(t)', fe(r.x, 'm')], ['Flux Φ = B l x', fe(r.flux, 'Wb')], ['emf (ccw ref.)', fe(r.emf, 'V')], ['|emf| = B l v', fe(r.mag, 'V')], ['Current I', fe(r.I, 'A')], ['Magnetic force on bar', fe(r.F, 'N')], ['Mechanical power F·v', fe(r.Pmech, 'W')], ['Electrical power I²R', fe(r.Pelec, 'W')]]);
        const dphi = v.Bb * v.l * v.v;
        ui.work.set(['Sliding bar: B = ' + g(v.Bb) + ' T along +az, rail gap l = ' + g(v.l) + ' m, v = ' + g(v.v) + ' m/s along +ax, R = ' + g(v.R) + ' Ω, x(0) = ' + g(v.x0) + ' m.  Evaluate at t = ' + g(tcur) + ' s.', '',
          'Step 1 - position:   x(t) = x(0) + v t = ' + g(v.x0) + ' + ' + g(v.v) + ' × ' + g(tcur) + ' = ' + g(r.x) + ' m',
          'Step 2 - flux:       Φ = B × (l x) = ' + g(v.Bb) + ' × ' + g(v.l) + ' × ' + g(r.x) + ' = ' + g(r.flux) + ' Wb   (normal +az, counter-clockwise reference)',
          'Step 3 - Faraday:    emf = −dΦ/dt = −B l v = −' + g(v.Bb) + ' × ' + g(v.l) + ' × ' + g(v.v) + ' = ' + g(r.emf) + ' V',
          '   Check with the motional emf: the force per unit charge on the bar is u × B = (v ax) × (B az) = −vB ay, along the bar, so ∮(u × B)·dL = B l v in magnitude = ' + g(r.mag) + ' V (for v > 0 positive charge is pushed toward −y).',
          'Step 4 - current:     I = |emf|/R = ' + g(r.mag) + '/' + g(v.R) + ' = ' + g(r.I) + ' A',
          'Step 5 - force:       F = I l B = ' + g(r.I) + ' × ' + g(v.l) + ' × ' + g(v.Bb) + ' = ' + g(r.F) + ' N, directed OPPOSITE to the motion (an external agent must supply it to keep v constant).',
          'Step 6 - power check: P_mech = F v = ' + g(r.Pmech) + ' W;   P_elec = emf²/R = I²R = ' + g(r.Pelec) + ' W.   Equal: all the mechanical work becomes heat in R.', '',
          lenz(dphi, 'the loop area is ' + (v.v > 0 ? 'growing' : 'shrinking'), '+az', v.Bb * v.l)]);
      } else {
        const c = M.currentRatio(v.sig, v.er, v.fd, v.E0);
        if (c.error) { ui.msg.show(c.error); blank(ui.cvA, 'Fix the inputs'); blank(ui.cvB, ''); hudSet(ui.hud, []); ui.work.set('Fix the input above.'); return; }
        const fs = logspace(1, 1e15, 160), rr = fs.map(f => M.currentRatio(v.sig, v.er, f, v.E0));
        data = { kind: 'disp', fs, Jc: rr.map(x => x.Jc), Jd: rr.map(x => x.Jd), lt: rr.map(x => x.lt), f: v.fd, c };
        hudSet(ui.hud, [['|Jc| = σE', fe(c.Jc, 'A/m²')], ['|Jd| = ωεE', fe(c.Jd, 'A/m²')], ['Jc / Jd = σ/ωε', g(c.lt, 4)], ['Class', c.cls], ['Crossover f = σ/2πε', fe(c.fCross, 'Hz')], ['ω', fe(c.omega, 'rad/s')]]);
        ui.work.set(['Medium: σ = ' + g(v.sig) + ' S/m, εr = ' + g(v.er) + ' (ε = εr ε0 = ' + g(v.er * M.EPS0) + ' F/m), μr = 1;  f = ' + g(v.fd) + ' Hz, |E| = ' + g(v.E0) + ' V/m.', '',
          'Step 1 - Ampère-Maxwell:   ∇×H = J + ∂D/∂t = σE + ε ∂E/∂t;  in phasors  ∇×H = (σ + jωε) E.',
          'Step 2 - conduction current density:    |Jc| = σ E = ' + g(v.sig) + ' × ' + g(v.E0) + ' = ' + g(c.Jc) + ' A/m²',
          'Step 3 - displacement current density:  |Jd| = ω ε E = 2π × ' + g(v.fd) + ' × ' + g(v.er * M.EPS0) + ' × ' + g(v.E0) + ' = ' + g(c.Jd) + ' A/m²   (Jd = jωεE leads Jc by 90°)',
          'Step 4 - ratio (loss tangent):  |Jc|/|Jd| = σ/(ωε) = ' + g(c.lt, 5),
          'Step 5 - classification (rule of thumb used here): ratio ≥ 100 → good conductor;  ratio ≤ 0.01 → good dielectric (low-loss);  in between → quasi-conductor (lossy medium);  σ = 0 → lossless.',
          '   Result: ' + c.cls + '.',
          'Step 6 - crossover where |Jc| = |Jd|:  f_c = σ/(2π ε) = ' + g(v.sig) + '/(2π × ' + g(v.er * M.EPS0) + ') = ' + g(c.fCross) + ' Hz.  Below f_c the medium behaves as a conductor, above it as a dielectric.']);
      }
      draw();
    }
    function draw() {
      if (!data) return; const D2 = data, T = theme();
      if (D2.kind === 'time') {
        const mk = (c, y, yl, col, pt) => { const yr = linRange(y.concat([0])); plot(c, { xmin: D2.tx[0], xmax: D2.tx[D2.tx.length - 1], ymin: yr[0], ymax: yr[1], xlabel: D2.xl, ylabel: yl, series: [{ x: D2.tx, y, color: col }], vlines: [{ x: D2.t, color: T.pink }], points: [{ x: D2.t, y: pt, color: T.pink }], hlines: [{ y: 0, color: T.muted }] }); };
        mk(ui.cvA, D2.A, D2.ylA, T.c, D2.a); mk(ui.cvB, D2.B, D2.ylB, T.a, D2.b);
      } else {
        const yr = decadeRange(D2.Jc.concat(D2.Jd, [D2.c.Jc, D2.c.Jd]), 1e-14), xr = [1, 1e15], fc = D2.c.fCross;
        plot(ui.cvA, { xmin: 1, xmax: 1e15, ymin: yr[0], ymax: yr[1], xlog: true, ylog: true, xlabel: 'f (Hz)', ylabel: '|J| (A/m²)', series: [{ x: D2.fs, y: D2.Jc, color: T.c }, { x: D2.fs, y: D2.Jd, color: T.a }], legend: [{ text: '— |Jc| = σE', color: T.c }, { text: '— |Jd| = ωεE', color: T.a }], vlines: [{ x: D2.f, color: T.pink, label: 'f' }].concat(fc > 1 && fc < 1e15 ? [{ x: fc, color: T.muted, label: 'f_c', dy: 14 }] : []), points: [{ x: D2.f, y: D2.c.Jc, color: T.c }, { x: D2.f, y: D2.c.Jd, color: T.a }] });
        const yr2 = decadeRange(D2.lt.concat([1e-3, 1e3]));
        plot(ui.cvB, { xmin: xr[0], xmax: xr[1], ymin: yr2[0], ymax: yr2[1], xlog: true, ylog: true, xlabel: 'f (Hz)', ylabel: 'Jc / Jd = σ/ωε', series: [{ x: D2.fs, y: D2.lt, color: T.pink }], hlines: [{ y: 100, color: T.b, label: '≥ 100 good conductor' }, { y: 1, color: T.muted, label: '= 1' }, { y: 0.01, color: T.warn, label: '≤ 0.01 good dielectric' }], vlines: [{ x: D2.f, color: T.pink }], points: [{ x: D2.f, y: Math.max(D2.c.lt, 1e-300), color: T.pink, label: D2.c.cls }] });
      }
    }
    return {
      update, draw,
      get: () => Object.assign({ mode: S.mode }, outVals(st)),
      set(o) { const m = oneOf(o.mode, MODES.map(x => x[0])); if (m) S.mode = m; loadSt(st, o); ui.mode.value = S.mode; },
    };
  }

  /* ---------- 4. plane waves & polarization ---------- */
  function buildPW(root, host) {
    const M = FSP.math.fields, { ctl, stage } = layout(root), st = {}, ui = {}, S = { preset: 'sea' };
    const PRE = [['air', 'Air / vacuum', 1, 1, 0], ['glass', 'Lossless dielectric (εr = 4)', 4, 1, 0], ['fresh', 'Fresh water', 81, 1, 1e-3], ['sea', 'Sea water', 81, 1, 4], ['soil', 'Wet soil', 15, 1, 0.01], ['cu', 'Copper', 1, 1, 5.8e7]];
    const f1 = FSP.ui.fieldset(ctl, 'Medium');
    ui.pre = FSP.ui.select(f1, 'Preset', PRE.map(p => [p[0], p[1]]).concat([['custom', 'Custom']]), S.preset, v => {
      S.preset = v; const p = PRE.filter(x => x[0] === v)[0];
      if (p) { st.er.set(p[2], true); st.mr.set(p[3], true); ui.zero.chk.checked = p[4] === 0; if (p[4] > 0) st.sg.set(p[4], true); }
      update(); FSP.state.touch();
    });
    makeSliders(f1, [
      { k: 'er', l: 'εr', min: 1, max: 1000, log: true, v: 81 }, { k: 'mr', l: 'μr', min: 1, max: 1000, log: true, v: 1 },
      { k: 'sg', l: 'σ', u: 'S/m', min: 1e-9, max: 1e8, log: true, v: 4 }, { k: 'f', l: 'Frequency', u: 'Hz', min: 1e3, max: 1e12, log: true, v: 1e6 },
    ], st, () => { S.preset = 'custom'; ui.pre.value = 'custom'; update(); });
    ui.zero = zeroBox(f1, 'σ = 0 (lossless medium)', false, () => { S.preset = 'custom'; ui.pre.value = 'custom'; update(); });
    const f2 = FSP.ui.fieldset(ctl, 'Polarization (E at z = 0, wave travels in +z)');
    makeSliders(f2, [
      { k: 'ex', l: 'Ex0', u: 'V/m', min: 0, max: 10, step: 0.1, v: 1 }, { k: 'ey', l: 'Ey0', u: 'V/m', min: 0, max: 10, step: 0.1, v: 1 },
      { k: 'px', l: 'φx', u: '°', min: -180, max: 180, step: 1, v: 0 }, { k: 'py', l: 'φy', u: '°', min: -180, max: 180, step: 1, v: -90 },
    ], st, () => update());
    mkNote(f2, 'Ex = Ex0 cos(ωt − βz + φx),  Ey = Ey0 cos(ωt − βz + φy).  Handedness is the IEEE (= Hayt) convention: right-hand thumb along the propagation direction (+z), fingers curl the way E rotates in time. Viewed head-on (wave coming toward you, as on the canvas) a right-handed wave turns counter-clockwise. Optics texts use the opposite naming.');
    const f3 = FSP.ui.fieldset(ctl, 'Power (Poynting)');
    makeSliders(f3, [{ k: 'zd', l: 'Depth z', u: 'm', min: 1e-4, max: 1e3, log: true, v: 0.25 }, { k: 'ar', l: 'Area A', u: 'm²', min: 1e-4, max: 100, log: true, v: 1 }], st, () => update());
    ui.msg = msgBox(f1);
    ui.hudM = hud(stage);
    ui.cv = canvasIn(stage, 'Tip of the electric field vector in the xy plane', 320); ui.cv.onResize(() => drawPol(anim.th));
    const bar = el('div', { class: 'row' }); stage.appendChild(bar);
    ui.play = FSP.ui.button(bar, 'Pause', () => { anim.playing = !anim.playing; ui.play.textContent = anim.playing ? 'Pause' : 'Play'; if (anim.playing) start(); else stop(); });
    ui.play.className = 'btn';
    ui.hudP = hud(stage); ui.work = FSP.ui.working(stage);
    const anim = { raf: 0, th: 0, last: 0, playing: !(typeof window !== 'undefined' && window.matchMedia && window.matchMedia('(prefers-reduced-motion: reduce)').matches) };
    if (!anim.playing) ui.play.textContent = 'Play';
    let cur = null;
    const sigmaNow = () => (ui.zero.chk.checked ? 0 : st.sg.get());
    function running() { return host.active && host.sub === 'pw' && !(typeof document !== 'undefined' && document.hidden); }
    function start() { if (anim.raf || !anim.playing || !running() || !cur || !cur.pol || cur.pol.error) return; anim.last = 0; anim.raf = requestAnimationFrame(loop); }
    function stop() { if (anim.raf) { cancelAnimationFrame(anim.raf); anim.raf = 0; } }
    function loop(ts) {
      anim.raf = 0; if (!running() || !anim.playing) return;
      if (anim.last) anim.th = (anim.th + (ts - anim.last) / 1000 * 2 * PI * 0.25) % (2 * PI); anim.last = ts;
      drawPol(anim.th); anim.raf = requestAnimationFrame(loop);
    }
    function update() {
      ui.msg.hide();
      const sg = sigmaNow(), m = M.medium(st.er.get(), st.mr.get(), sg, st.f.get());
      const pol = M.polarization(st.ex.get(), st.ey.get(), st.px.get(), st.py.get());
      cur = { m, pol };
      if (m.error) { ui.msg.show(m.error); hudSet(ui.hudM, []); } else hudSet(ui.hudM, [['α (attenuation)', fe(m.alpha, 'Np/m')], ['β (phase)', fe(m.beta, 'rad/m')], ['η', gp(m.eta, 4) + ' Ω'], ['λ = 2π/β', fe(m.lambda, 'm')], ['vp = ω/β', fe(m.vp, 'm/s') + ' (' + g(m.vp / M.C0, 4) + ' c)'], ['Skin depth δ = 1/α', m.alpha > 0 ? fe(m.delta, 'm') : '∞ (lossless)'], ['Loss 8.686 α', g(m.dBpm, 4) + ' dB/m'], ['σ/(ωε)', g(m.lt, 4)], ['Class', m.cls]]);
      if (pol.error) { ui.msg.show(pol.error); hudSet(ui.hudP, []); blank(ui.cv, 'Set a non-zero amplitude'); stop(); ui.work.set(m.error ? 'Fix the input above.' : mediumWork(m)); return; }
      const rows = [['Polarization', pol.type + (pol.hand ? ', ' + pol.hand + '-handed' : '')], ['δ = φy − φx', g(pol.deltaDeg, 4) + '°'], ['Axial ratio', pol.AR === Infinity ? '∞' : g(pol.AR, 4) + ' (' + g(pol.ARdB, 4) + ' dB)'], ['Tilt angle τ', fin(pol.tiltDeg) ? g(pol.tiltDeg, 4) + '° from +x' : 'not defined (circular)'], ['Semi-axes a, b', g(pol.major, 4) + ', ' + g(pol.minor, 4) + ' V/m']];
      let pw = null;
      if (!m.error) {
        pw = M.poyntingWave(pol.Ex0, pol.Ey0, m, st.zd.get(), st.ar.get());
        rows.push(['⟨S⟩ at z = 0', fe(pw.S0, 'W/m²')], ['⟨S⟩ at depth z', fe(pw.Sz, 'W/m²')], ['Power through A', fe(pw.P, 'W')], ['Attenuation to z', g(pw.lossDb, 4) + ' dB']);
      }
      hudSet(ui.hudP, rows);
      ui.work.set((m.error ? [] : mediumWork(m)).concat(['', '─────────────────────────────────────────────', ''], polWork(pol), pw ? ['', '─────────────────────────────────────────────', ''].concat(powerWork(pol, m, pw)) : []));
      drawPol(anim.th); start();
    }
    function mediumWork(m) {
      const w = m.omega, jwmu = cx(0, w * m.mu), ys = cx(m.sigma, w * m.eps), prod = cmul(jwmu, ys), div = cdiv(jwmu, ys), L = [];
      L.push('PART A - medium (Hayt notation: fields ∝ e^{jωt − γz}, γ = α + jβ)', '',
        'Step 1 - constants:  ω = 2πf = ' + g(w) + ' rad/s;  ε = εr ε0 = ' + g(m.eps) + ' F/m;  μ = μr μ0 = ' + g(m.mu) + ' H/m',
        'Step 2 - loss tangent:  σ/(ωε) = ' + g(m.sigma) + ' / (' + g(w) + ' × ' + g(m.eps) + ') = ' + g(m.lt, 5) + '  →  ' + m.cls + '  (≤ 0.01 low-loss dielectric, ≥ 100 good conductor)',
        'Step 3 - propagation constant:  γ = √( jωμ (σ + jωε) )',
        '   jωμ = j' + g(w * m.mu) + ';   σ + jωε = ' + g(m.sigma) + ' + j' + g(w * m.eps) + ';   product = ' + gc(prod),
        '   γ = ' + gc(m.gamma) + ' m⁻¹   →   α = ' + g(m.alpha) + ' Np/m,   β = ' + g(m.beta) + ' rad/m',
        '   check with α = ω√(με/2) [√(1 + (σ/ωε)²) − 1]^½ = ' + g(m.alphaAlt) + ',   β = ω√(με/2) [√(1 + (σ/ωε)²) + 1]^½ = ' + g(m.betaAlt),
        'Step 4 - intrinsic impedance:  η = √( jωμ/(σ + jωε) ),  jωμ/(σ + jωε) = ' + gc(div) + ' Ω²,   η = ' + gc(m.eta) + ' Ω = ' + gp(m.eta, 5) + ' Ω',
        'Step 5 - derived quantities:  λ = 2π/β = ' + g(m.lambda) + ' m;   vp = ω/β = ' + g(m.vp) + ' m/s;   skin depth δ = 1/α = ' + (m.alpha > 0 ? g(m.delta) + ' m' : '∞') + ';   attenuation = 20 log10(e) α = 8.686 × ' + g(m.alpha) + ' = ' + g(m.dBpm) + ' dB/m');
      if (m.approx) {
        const ap = m.approx;
        L.push('Step 6 - ' + ap.name + ' approximation:');
        if (ap.name === 'lossless') L.push('   α = 0, β = ω√(με) = ' + g(ap.beta) + ' rad/m, η = √(μ/ε) = ' + g(ap.eta.re) + ' Ω  (exact here).');
        else if (ap.name === 'low-loss') L.push('   α ≈ (σ/2)√(μ/ε) = ' + g(ap.alpha) + ' Np/m;   β ≈ ω√(με)[1 + (σ/ωε)²/8] = ' + g(ap.beta) + ' rad/m;   η ≈ √(μ/ε)[1 + j σ/(2ωε)] = ' + gc(ap.eta) + ' Ω');
        else L.push('   α ≈ β ≈ √(π f μ σ) = ' + g(ap.alpha) + ' (= 1/δ);   η ≈ (1 + j)√(π f μ/σ) = ' + gc(ap.eta) + ' Ω  (45° impedance angle).');
      } else L.push('Step 6 - quasi-conductor: neither the low-loss nor the good-conductor approximation is valid; use the exact expressions above.');
      return L;
    }
    function polWork(p) {
      const L = ['PART B - polarization', '',
        'E(z,t) = ax Ex0 cos(ωt − βz + φx) + ay Ey0 cos(ωt − βz + φy),   Ex0 = ' + g(p.Ex0) + ', Ey0 = ' + g(p.Ey0) + ' V/m, φx = ' + g(p.phx) + '°, φy = ' + g(p.phy) + '°',
        'Step 1 - phase difference:  δ = φy − φx = ' + g(p.phy - p.phx) + '° ≡ ' + g(p.deltaDeg, 4) + '° (wrapped into −180°…180°)',
        'Step 2 - classify:  linear if Ex0 = 0, or Ey0 = 0, or δ = 0, ±180°;   circular if Ex0 = Ey0 and δ = ±90°;   otherwise elliptical.   → ' + p.type.toUpperCase()];
      if (p.type === 'linear') L.push('Step 3 - linear: the tip oscillates along a line at angle τ = ' + g(p.tiltDeg, 4) + '° to +x with amplitude ' + g(p.major, 4) + ' V/m. No handedness.');
      else {
        L.push('Step 3 - handedness:  sin δ = ' + g(Math.sin(p.deltaDeg * D2R), 4) + ' → ' + (p.hand === 'right' ? 'sin δ < 0: E turns from +x toward +y, counter-clockwise about +z' : 'sin δ > 0: E turns from +y toward +x, clockwise about +z') + ' → ' + p.hand.toUpperCase() + '-HANDED (IEEE/Hayt: thumb along +z, fingers along the rotation).',
          '         Check: Ey ' + (p.hand === 'right' ? 'lags' : 'leads') + ' Ex by |δ| = ' + g(Math.abs(p.deltaDeg), 4) + '°, so the tip swings from ' + (p.hand === 'right' ? '+x toward +y' : '+y toward +x') + '.');
        const s = p.Ex0 * p.Ex0 + p.Ey0 * p.Ey0, q = Math.sqrt(Math.max(0, Math.pow(p.Ex0, 4) + Math.pow(p.Ey0, 4) + 2 * p.Ex0 * p.Ex0 * p.Ey0 * p.Ey0 * Math.cos(2 * p.deltaDeg * D2R)));
        L.push('Step 4 - ellipse:  s = Ex0² + Ey0² = ' + g(s) + ',   q = √(Ex0⁴ + Ey0⁴ + 2Ex0²Ey0² cos 2δ) = ' + g(q),
          '         semi-major a = √((s + q)/2) = ' + g(p.major) + ' V/m,   semi-minor b = √((s − q)/2) = ' + g(p.minor) + ' V/m',
          '         axial ratio AR = a/b = ' + (p.AR === Infinity ? '∞' : g(p.AR, 5)) + ' (' + g(p.ARdB, 4) + ' dB);  AR = 1 for circular.',
          '         tilt τ = ½ atan2(2 Ex0 Ey0 cos δ, Ex0² − Ey0²) = ' + (fin(p.tiltDeg) ? g(p.tiltDeg, 4) + '°' : 'not defined (circular: any axis)') + ' from +x (major axis).');
      }
      return L;
    }
    function powerWork(p, m, pw) {
      const z = st.zd.get();
      return ['PART C - average Poynting vector', '',
        'For a wave in +z, H = (az × E)/η.  ⟨S⟩ = ½ Re{ E × H* } = ½ (Ex0² + Ey0²) Re(1/η*) e^{−2αz} az.   (Polarization shape does not change the average, only the total amplitude.)',
        'Step 1 - Re(1/η*) = cos θη/|η| = cos(' + g(m.etaDeg, 4) + '°)/' + g(m.etaMag, 5) + ' = ' + g(pw.cosEta / m.etaMag, 5) + ' S',
        'Step 2 - at z = 0:   ⟨S⟩ = ½ (' + g(p.Ex0) + '² + ' + g(p.Ey0) + '²) × ' + g(pw.cosEta / m.etaMag, 5) + ' = ' + g(pw.S0) + ' W/m²',
        'Step 3 - at depth z = ' + g(z) + ' m:   ⟨S⟩ = ⟨S⟩0 e^{−2αz} = ' + g(pw.S0) + ' × e^{−2 × ' + g(m.alpha) + ' × ' + g(z) + '} = ' + g(pw.Sz) + ' W/m²',
        'Step 4 - power through A = ' + g(st.ar.get()) + ' m²:   P = ⟨S⟩ A = ' + g(pw.P) + ' W;   loss to depth z = 8.686 α z = ' + g(pw.lossDb, 4) + ' dB (power), i.e. field amplitude falls to e^{−αz} = ' + g(Math.exp(-m.alpha * z), 4) + ' of its surface value.'];
    }
    function drawPol(th) {
      const gg = ui.cv.prep(), ctx = gg.ctx, w = gg.w, h = gg.h, T = theme(); ctx.clearRect(0, 0, w, h); ctx.fillStyle = T.bg; ctx.fillRect(0, 0, w, h);
      if (!cur || !cur.pol || cur.pol.error) return;
      const p = cur.pol, cxp = w / 2, cyp = h / 2, amp = Math.max(p.Ex0, p.Ey0), sc = Math.min(w, h) * 0.4 / amp, X = x => cxp + x * sc, Y = y => cyp - y * sc;
      const E = a => [p.Ex0 * Math.cos(a + p.phx * D2R), p.Ey0 * Math.cos(a + p.phy * D2R)];
      ctx.font = '11px ' + T.mono; ctx.lineWidth = 1; ctx.strokeStyle = T.border; ctx.fillStyle = T.muted;
      ctx.beginPath(); ctx.moveTo(cxp - amp * sc * 1.15, cyp); ctx.lineTo(cxp + amp * sc * 1.15, cyp); ctx.moveTo(cxp, cyp - amp * sc * 1.15); ctx.lineTo(cxp, cyp + amp * sc * 1.15); ctx.stroke();
      ctx.textAlign = 'left'; ctx.textBaseline = 'middle'; ctx.fillText('x', cxp + amp * sc * 1.15 + 3, cyp); ctx.textAlign = 'center'; ctx.textBaseline = 'bottom'; ctx.fillText('y', cxp, cyp - amp * sc * 1.15 - 2);
      // bounding box of the ellipse
      ctx.setLineDash([3, 3]); ctx.strokeRect(X(-p.Ex0), Y(p.Ey0), 2 * p.Ex0 * sc, 2 * p.Ey0 * sc); ctx.setLineDash([]);
      // path
      ctx.strokeStyle = T.muted; ctx.lineWidth = 1.5; ctx.beginPath(); for (let i = 0; i <= 120; i++) { const e = E(i / 120 * 2 * PI); if (i) ctx.lineTo(X(e[0]), Y(e[1])); else ctx.moveTo(X(e[0]), Y(e[1])); } ctx.stroke();
      // trail: the last ~100° of phase, fading
      const Et = a => E(a);
      for (let i = 0; i < 24; i++) { const a0 = th - (i + 1) * 0.075, a1 = th - i * 0.075, e0 = Et(a0), e1 = Et(a1); ctx.strokeStyle = T.c; ctx.globalAlpha = 0.9 * (1 - i / 24); ctx.beginPath(); ctx.moveTo(X(e0[0]), Y(e0[1])); ctx.lineTo(X(e1[0]), Y(e1[1])); ctx.stroke(); }
      ctx.globalAlpha = 1;
      const e = Et(th); ctx.strokeStyle = T.a; ctx.fillStyle = T.a; ctx.lineWidth = 2; arrow(ctx, cxp, cyp, X(e[0]), Y(e[1]), 8);
      ctx.fillStyle = T.pink; ctx.beginPath(); ctx.arc(X(e[0]), Y(e[1]), 5, 0, 2 * PI); ctx.fill();
      ctx.fillStyle = T.text; ctx.font = '12px ' + T.mono; ctx.textAlign = 'left'; ctx.textBaseline = 'top';
      ctx.fillText(p.type + (p.hand ? ', ' + p.hand + '-handed' : ''), 8, 6);
      ctx.fillStyle = T.muted; ctx.font = '10px ' + T.mono; ctx.fillText('wave travels out of the screen (+z); viewed head-on', 8, 22);
    }
    return {
      update, draw: () => drawPol(anim.th), stop, start,
      get: () => Object.assign(outVals(st), { sz: ui.zero.chk.checked ? '1' : '0' }),
      set(o) { loadSt(st, o); if (o.sz === '1' || o.sz === '0') ui.zero.chk.checked = o.sz === '1'; S.preset = 'custom'; ui.pre.value = 'custom'; },
    };
  }

  /* ---------- 5. normal incidence & standing waves ---------- */
  function buildNI(root) {
    const M = FSP.math.fields, { ctl, stage } = layout(root), st = {}, ui = {}, S = { lay: 'int', t2: 'mat', t3: 'mat' };
    const PRE = [
      ['p1', 'Air → dielectric (εr = 4)', { f: 1e9, m1: [1, 1, 0], m2: [4, 1, 0], lay: 'int' }],
      ['p2', 'Air → perfect conductor', { f: 1e9, m1: [1, 1, 0], m2: 'pec', lay: 'int' }],
      ['p3', 'Air → sea water (1 MHz)', { f: 1e6, m1: [1, 1, 0], m2: [81, 1, 4], lay: 'int' }],
      ['p4', 'Dielectric (εr = 4) → air', { f: 1e9, m1: [4, 1, 0], m2: [1, 1, 0], lay: 'int' }],
      ['p5', 'Quarter-wave layer (air | εr = 4 | εr = 16)', { f: 1e9, m1: [1, 1, 0], m2: [4, 1, 0], lay: 'slab', dq: 0.25, m3: [16, 1, 0] }],
      ['p6', 'Half-wave window (air | εr = 9 | air)', { f: 1e9, m1: [1, 1, 0], m2: [9, 1, 0], lay: 'slab', dq: 0.5, m3: [1, 1, 0] }],
    ];
    const f0 = FSP.ui.fieldset(ctl, 'Setup');
    ui.pre = FSP.ui.select(f0, 'Preset', [['', '- choose -']].concat(PRE.map(p => [p[0], p[1]])), '', v => { const p = PRE.filter(x => x[0] === v)[0]; if (p) applyPreset(p[2]); });
    ui.lay = FSP.ui.select(f0, 'Layout', [['int', 'Single interface'], ['slab', 'Slab (layer) between 1 and 3']], S.lay, v => { S.lay = v; update(); FSP.state.touch(); });
    makeSliders(f0, [{ k: 'f', l: 'Frequency', u: 'Hz', min: 1e3, max: 1e12, log: true, v: 1e9 }], st, () => update());
    ui.msg = msgBox(f0);
    ui.fs1 = FSP.ui.fieldset(ctl, 'Medium 1 (z < 0), incident from the left');
    makeSliders(ui.fs1, [{ k: 'e1', l: 'εr1', min: 1, max: 1000, log: true, v: 1 }, { k: 'u1', l: 'μr1', min: 1, max: 1000, log: true, v: 1 }, { k: 's1', l: 'σ1', u: 'S/m', min: 1e-9, max: 1e8, log: true, v: 1 }], st, () => update());
    ui.zero1 = zeroBox(ui.fs1, 'σ1 = 0 (lossless)', true, () => update());
    ui.fs2 = FSP.ui.fieldset(ctl, 'Medium 2 (z > 0)');
    ui.t2 = FSP.ui.select(ui.fs2, 'Type', [['mat', 'Material'], ['pec', 'Perfect conductor (PEC)']], S.t2, v => { S.t2 = v; update(); FSP.state.touch(); });
    makeSliders(ui.fs2, [{ k: 'e2', l: 'εr2', min: 1, max: 1000, log: true, v: 4 }, { k: 'u2', l: 'μr2', min: 1, max: 1000, log: true, v: 1 }, { k: 's2', l: 'σ2', u: 'S/m', min: 1e-9, max: 1e8, log: true, v: 1 }], st, () => update());
    ui.zero2 = zeroBox(ui.fs2, 'σ2 = 0 (lossless)', true, () => update());
    ui.fsd = FSP.ui.fieldset(ctl, 'Slab thickness');
    makeSliders(ui.fsd, [{ k: 'd', l: 'Thickness d', u: 'm', min: 1e-4, max: 10, log: true, v: 0.0375 }], st, () => update());
    const brow = el('div', { class: 'row' }); ui.fsd.appendChild(brow);
    FSP.ui.button(brow, 'd = λ2/4', () => setD(0.25)).className = 'btn'; FSP.ui.button(brow, 'd = λ2/2', () => setD(0.5)).className = 'btn';
    ui.fs3 = FSP.ui.fieldset(ctl, 'Medium 3 (z > d)');
    ui.t3 = FSP.ui.select(ui.fs3, 'Type', [['mat', 'Material'], ['pec', 'Perfect conductor (PEC)']], S.t3, v => { S.t3 = v; update(); FSP.state.touch(); });
    makeSliders(ui.fs3, [{ k: 'e3', l: 'εr3', min: 1, max: 1000, log: true, v: 16 }, { k: 'u3', l: 'μr3', min: 1, max: 1000, log: true, v: 1 }, { k: 's3', l: 'σ3', u: 'S/m', min: 1e-9, max: 1e8, log: true, v: 1 }], st, () => update());
    ui.zero3 = zeroBox(ui.fs3, 'σ3 = 0 (lossless)', true, () => update());
    ui.hud = hud(stage);
    ui.cv = canvasIn(stage, '|E| against position z for normal incidence', 320); ui.cv.onResize(() => draw());
    ui.pw = canvasIn(stage, 'Power fractions reflected, transmitted and absorbed', 70); ui.pw.onResize(() => draw());
    ui.work = FSP.ui.working(stage);
    let data = null;
    const sgm = (k, z) => (z.chk.checked ? 0 : st[k].get());
    function lam2() { const m = M.medium(st.e2.get(), st.u2.get(), sgm('s2', ui.zero2), st.f.get()); return m.error ? NaN : m.lambda; }
    function setD(frac) { const l = lam2(); if (fin(l)) { st.d.set(l * frac); } }
    function applyPreset(p) {
      st.f.set(p.f, true); [['e1', 'u1', 's1', ui.zero1, p.m1], ['e2', 'u2', 's2', ui.zero2, p.m2], ['e3', 'u3', 's3', ui.zero3, p.m3]].forEach(x => {
        const v = x[4]; if (!v) return; if (v === 'pec') return; st[x[0]].set(v[0], true); st[x[1]].set(v[1], true); x[3].chk.checked = v[2] === 0; if (v[2] > 0) st[x[2]].set(v[2], true);
      });
      S.t2 = p.m2 === 'pec' ? 'pec' : 'mat'; ui.t2.value = S.t2; S.lay = p.lay; ui.lay.value = p.lay;
      if (p.dq) { const l = lam2(); if (fin(l)) st.d.set(l * p.dq, true); }
      update(); FSP.state.touch();
    }
    const medOf = (k, z) => M.medium(st['e' + k].get(), st['u' + k].get(), sgm('s' + k, z), st.f.get());
    function maxPos(G, m1) {
      const Gm = cabs(G); if (!(Gm > 1e-12) || m1.alpha !== 0) return null;
      let th = carg(G); if (th < 0) th += 2 * PI; let lmax = th / (4 * PI) * m1.lambda; if (lmax >= m1.lambda / 2 - 1e-15 * m1.lambda) lmax -= m1.lambda / 2;
      return { lmax, lmin: lmax >= m1.lambda / 4 ? lmax - m1.lambda / 4 : lmax + m1.lambda / 4, th };
    }
    const ext = (m, k) => (m.pec ? NaN : m.alpha > 0 ? Math.min(k * m.lambda, 3 * m.delta) : k * m.lambda);
    const medDesc = (name, m) => name + ': εr = ' + g(m.er) + ', μr = ' + g(m.mr) + ', σ = ' + g(m.sigma) + ' S/m, σ/ωε = ' + g(m.lt, 4) + ' (' + m.cls + ')  →  η = ' + gp(m.eta, 5) + ' Ω,  β = ' + g(m.beta) + ' rad/m,  α = ' + g(m.alpha) + ' Np/m,  λ = ' + g(m.lambda) + ' m';
    function update() {
      ui.msg.hide(); data = null; const slab = S.lay === 'slab';
      ui.fs2.firstChild.textContent = slab ? 'Slab material (0 < z < d)' : 'Medium 2 (z > 0)'; ui.t2.parentNode.hidden = slab; ui.fsd.hidden = !slab; ui.fs3.hidden = !slab;
      const pec2 = !slab && S.t2 === 'pec', pec3 = slab && S.t3 === 'pec';
      ['e2', 'u2', 's2'].forEach(k => { st[k].el.hidden = pec2; }); ui.zero2.row.hidden = pec2; ['e3', 'u3', 's3'].forEach(k => { st[k].el.hidden = pec3; }); ui.zero3.row.hidden = pec3;
      const m1 = medOf(1, ui.zero1), m2 = pec2 ? M.PEC : medOf(2, ui.zero2), m3 = slab ? (pec3 ? M.PEC : medOf(3, ui.zero3)) : null;
      const bad = [m1, slab ? medOf(2, ui.zero2) : m2, m3].filter(m => m && m.error)[0];
      if (bad) { ui.msg.show(bad.error); blank(ui.cv, 'Fix the inputs'); hudSet(ui.hud, []); ui.work.set('Fix the input above.'); return; }
      const d = st.d.get(); let pts = [], hud_ = [], L = [], R, T, Gm, G, SWR, mp, segs = [], hl = [];
      const left = ext(m1, 2.5), segStyle = [];
      if (!slab) {
        const ni = M.normalIncidence(m1, m2); G = ni.G; Gm = ni.Gm; SWR = ni.SWR; R = ni.R; T = ni.T;
        mp = fin(ni.lmax) ? { lmax: ni.lmax, lmin: ni.lmin, th: ni.thG } : null;
        const right = m2.pec ? 0 : ext(m2, 1.5);
        segs.push({ z: linspace(-left, 0, 500), c: 'a', f: z => cabs(M.fieldNormal(m1, m2, ni, z)) });
        if (!m2.pec) segs.push({ z: linspace(0, right, 300), c: 'b', f: z => cabs(M.fieldNormal(m1, m2, ni, z)) }); else segs.push({ z: [0, 0.4 * m1.lambda || 1], c: 'b', f: () => 0 });
        hud_ = [['η1', gp(m1.eta, 4) + ' Ω'], ['η2', m2.pec ? '0 (PEC)' : gp(m2.eta, 4) + ' Ω'], ['Γ = (η2−η1)/(η2+η1)', gc(G, 4) + '  =  ' + gp(G, 4)], ['τ = 1 + Γ', gp(ni.tau, 4)], ['SWR', SWR === Infinity ? '∞' : g(SWR, 4)]];
        if (mp) hud_.push(['First |E| max (from z = 0)', fe(mp.lmax, 'm') + ' = ' + g(mp.lmax / m1.lambda, 4) + ' λ1'], ['First |E| min', fe(mp.lmin, 'm') + ' = ' + g(mp.lmin / m1.lambda, 4) + ' λ1']);
        else hud_.push(['|E| max / min positions', Gm <= 1e-12 ? 'none (Γ = 0: travelling wave only)' : 'only defined for lossless medium 1']);
        L = workInt(m1, m2, ni, mp);
      } else {
        const s = M.slab(m1, m2, d, m3); G = s.Gin; Gm = s.Gm; SWR = s.SWR; R = s.R; T = s.T3; mp = maxPos(G, m1);
        const right = m3.pec ? 0 : ext(m3, 1.5), nS = Math.min(1200, Math.max(120, Math.round(d / (m2.lambda || d) * 160)));
        segs.push({ z: linspace(-left, 0, 500), c: 'a', f: z => cabs(M.fieldSlab(m1, m2, d, m3, s, z)) });
        segs.push({ z: linspace(0, d, nS), c: 'c', f: z => cabs(M.fieldSlab(m1, m2, d, m3, s, z)) });
        if (!m3.pec) segs.push({ z: linspace(d, d + right, 300), c: 'b', f: z => cabs(M.fieldSlab(m1, m2, d, m3, s, z)) });
        hud_ = [['η1, η2 (slab), η3', gp(m1.eta, 3) + ', ' + gp(m2.eta, 3) + ', ' + (m3.pec ? '0' : gp(m3.eta, 3)) + ' Ω'], ['Slab thickness', fe(d, 'm') + ' = ' + g(d / m2.lambda, 4) + ' λ2 (β2 d = ' + g(m2.beta * d, 4) + ' rad)'], ['Input impedance ηin', gp(s.etaIn, 4) + ' Ω'], ['Γin = (ηin−η1)/(ηin+η1)', gc(G, 4) + '  =  ' + gp(G, 4)], ['SWR in medium 1', SWR === Infinity ? '∞' : g(SWR, 4)]];
        if (mp) hud_.push(['First |E1| max (from z = 0)', fe(mp.lmax, 'm') + ' = ' + g(mp.lmax / m1.lambda, 4) + ' λ1'], ['First |E1| min', fe(mp.lmin, 'm') + ' = ' + g(mp.lmin / m1.lambda, 4) + ' λ1']);
        L = workSlab(m1, m2, m3, d, s, mp);
      }
      const absorbed = z0(Math.max(0, 1 - R - T));
      hud_.push(['Reflected R = |Γ|²', g(100 * z0(R), 5) + ' %'], ['Transmitted T' + (slab ? ' (into medium 3)' : ''), g(100 * z0(T), 5) + ' %'], ['Absorbed (1 − R − T)', g(100 * z0(absorbed), 4) + ' %']);
      hudSet(ui.hud, hud_); ui.work.set(L);
      data = { segs, slab, d, m1, G, Gm, mp, left, R, T, absorbed, pec: slab ? pec3 : pec2 };
      draw();
    }
    function workInt(m1, m2, ni, mp) {
      const L = ['Normal incidence of a uniform plane wave from medium 1 (z < 0) onto medium 2 (z > 0); incident amplitude E0⁺ = 1 (all fields are normalised to it). Phasors e^{jωt}, +z waves e^{−γz}.', '',
        'Step 1 - media at f = ' + g(st.f.get()) + ' Hz  (ω = ' + g(2 * PI * st.f.get()) + ' rad/s):', '   ' + medDesc('Medium 1', m1), m2.pec ? '   Medium 2: perfect conductor, η2 = 0 and E2 = 0 (all power is reflected).' : '   ' + medDesc('Medium 2', m2), '',
        'Step 2 - boundary conditions (tangential E and H continuous at z = 0):  1 + Γ = τ,   (1 − Γ)/η1 = τ/η2.',
        'Step 3 - reflection coefficient:   Γ = (η2 − η1)/(η2 + η1) = (' + (m2.pec ? '0' : gc(m2.eta, 4)) + ' − (' + gc(m1.eta, 4) + ')) / (' + (m2.pec ? '0' : gc(m2.eta, 4)) + ' + ' + gc(m1.eta, 4) + ') = ' + gc(ni.G, 5) + ' = ' + gp(ni.G, 5),
        'Step 4 - transmission coefficient: τ = 1 + Γ = 2η2/(η2 + η1) = ' + gc(ni.tau, 5) + ' = ' + gp(ni.tau, 5),
        'Step 5 - standing-wave ratio:       SWR = (1 + |Γ|)/(1 − |Γ|) = (1 + ' + g(ni.Gm) + ')/(1 − ' + g(ni.Gm) + ') = ' + (ni.SWR === Infinity ? '∞' : g(ni.SWR, 5)) + '    (|E1|max = 1 + |Γ| = ' + g(1 + ni.Gm) + ',  |E1|min = 1 − |Γ| = ' + g(1 - ni.Gm) + ')'];
      if (mp) {
        L.push('Step 6 - positions (lossless medium 1): E1(z) = E0⁺ e^{−jβ1 z}[1 + Γ e^{2jβ1 z}],  so |E1| = √(1 + |Γ|² + 2|Γ| cos(2β1 z + θΓ)).',
          '   Maxima occur where 2β1 z + θΓ = 2nπ, i.e. at distance  d = −z = θΓ/(2β1) − nπ/β1  from the boundary (θΓ taken in 0…2π).  The first one: d_max = θΓ/(2β1) = θΓ λ1/(4π), reduced by λ1/2 if that is ≥ λ1/2:  d_max = ' + g(mp.th * R2D, 5) + '° × λ1/720° → ' + g(mp.lmax) + ' m = ' + g(mp.lmax / m1.lambda, 4) + ' λ1 in front of the boundary (z = −d_max).',
          '   Minima are λ1/4 from the maxima: first at d_min = ' + g(mp.lmin) + ' m = ' + g(mp.lmin / m1.lambda, 4) + ' λ1. Successive maxima (or minima) are λ1/2 = ' + g(m1.lambda / 2) + ' m apart.',
          '   (Rule of thumb: Γ real and negative → minimum at the boundary; Γ real and positive → maximum at the boundary.)');
      } else L.push('Step 6 - positions: ' + (ni.Gm <= 1e-12 ? 'Γ = 0, no reflection, |E1| = 1 everywhere (matched).' : 'maxima and minima are given here only for a lossless medium 1 (in a lossy medium the standing-wave pattern is not periodic).'));
      L.push('', 'Step 7 - power (normal incidence):  R = |Γ|² = ' + g(z0(ni.R), 6) + ';   T = |τ|² Re(1/η2)/Re(1/η1) = ' + g(z0(ni.T), 6) + ';   R + T = ' + g(ni.R + ni.T, 6) + (m2.pec ? '   (PEC: nothing is transmitted.)' : (Math.abs(1 - ni.R - ni.T) < 1e-9 ? '   (= 1: power is conserved.)' : '   (< 1: the difference is dissipated at the interface region / in medium 1 losses.)')));
      if (!m2.pec) L.push('   Transmitted field in medium 2: E2(z) = τ e^{−γ2 z}; at z = 0: |E2| = ' + g(cabs(ni.tau), 5) + ' (can exceed 1, but T ≤ 1 because power also depends on 1/η).');
      return L;
    }
    function workSlab(m1, m2, m3, d, s, mp) {
      const gd = cx(m2.alpha * d, m2.beta * d);
      const L = ['Three media: 1 (z < 0) | slab 2 (0 < z < d) | 3 (z > d). Incident amplitude 1 from medium 1. Phasors e^{jωt}.', '',
        'Step 1 - media at f = ' + g(st.f.get()) + ' Hz:', '   ' + medDesc('Medium 1', m1), '   ' + medDesc('Slab 2', m2), m3.pec ? '   Medium 3: perfect conductor (η3 = 0).' : '   ' + medDesc('Medium 3', m3), '',
        'Step 2 - electrical thickness:  γ2 d = (' + g(m2.alpha) + ' + j' + g(m2.beta) + ') × ' + g(d) + ' = ' + gc(gd, 5) + ';   β2 d = ' + g(m2.beta * d, 5) + ' rad = ' + g(m2.beta * d * R2D, 5) + '° = ' + g(d / m2.lambda, 5) + ' λ2.',
        'Step 3 - input impedance at z = 0 looking into the slab:',
        '   ηin = η2 (η3 + η2 tanh γ2 d)/(η2 + η3 tanh γ2 d)  [computed as η2(η3 cosh γd + η2 sinh γd)/(η2 cosh γd + η3 sinh γd)]',
        '       = ' + gc(s.etaIn, 5) + ' Ω = ' + gp(s.etaIn, 5) + ' Ω.   Special cases: d = λ2/4 → ηin = η2²/η3;   d = λ2/2 → ηin = η3.',
        'Step 4 - reflection seen in medium 1:  Γin = (ηin − η1)/(ηin + η1) = ' + gp(s.Gin, 5) + ';   SWR = (1 + |Γin|)/(1 − |Γin|) = ' + (s.SWR === Infinity ? '∞' : g(s.SWR, 5)),
        'Step 5 - power: R = |Γin|² = ' + g(z0(s.R), 6) + ';   power into medium 3: T3 = |E(d)|² Re(1/η3)/Re(1/η1) = ' + g(z0(s.T3), 6) + ';   dissipated in the slab/medium 3 losses = 1 − R − T3 = ' + g(z0(Math.max(0, 1 - s.R - s.T3)), 5) + '.'];
      if (mp) L.push('Step 6 - standing wave in medium 1 (lossless): first maximum at ' + g(mp.lmax) + ' m = ' + g(mp.lmax / m1.lambda, 4) + ' λ1 from z = 0, first minimum at ' + g(mp.lmin) + ' m; |E1|max = 1 + |Γin| = ' + g(1 + s.Gm, 5) + ', |E1|min = 1 − |Γin| = ' + g(1 - s.Gm, 5) + '.');
      const lossless = m1.alpha === 0 && m2.alpha === 0 && !m3.pec && m3.alpha === 0;
      if (lossless) {
        const eq = Math.sqrt(m1.eta.re * m3.eta.re);
        L.push('Quarter-wave matching: with d = λ2/4 the layer is matched (Γin = 0) when η2 = √(η1 η3) = √(' + g(m1.eta.re, 5) + ' × ' + g(m3.eta.re, 5) + ') = ' + g(eq, 5) + ' Ω.   This slab has η2 = ' + g(m2.eta.re, 5) + ' Ω (' + (Math.abs(eq - m2.eta.re) < 1e-3 * eq ? 'matched value' : 'not the matched value') + ').');
      }
      return L;
    }
    function draw() {
      if (!data) { return; }
      const D2 = data, T = theme(), col = { a: T.a, b: T.b, c: T.c };
      const series = D2.segs.map(sg => ({ x: sg.z, y: sg.z.map(sg.f), color: col[sg.c], width: 2 })), all = [];
      series.forEach(s => s.y.forEach(v => all.push(v)));
      const ymax = Math.max(1.2, 1.1 * Math.max.apply(null, all.filter(fin))), xmin = -D2.left, xmax = D2.segs[D2.segs.length - 1].z[D2.segs[D2.segs.length - 1].z.length - 1];
      const hl = [{ y: 1, color: T.muted, label: '|E0⁺| = 1' }], vl = [{ x: 0, color: T.pink, label: 'z = 0' }];
      if (D2.slab) vl.push({ x: D2.d, color: T.pink, label: 'z = d', dy: 14 });
      const dots = [];
      if (D2.mp && D2.m1.alpha === 0) {
        hl.push({ y: 1 + D2.Gm, color: T.warn, label: '|E|max = 1 + |Γ|' }, { y: 1 - D2.Gm, color: T.b, label: '|E|min = 1 − |Γ|' });
        const half = D2.m1.lambda / 2;
        for (let k = 0; k < 8; k++) { const a = -(D2.mp.lmax + k * half), b = -(D2.mp.lmin + k * half); if (a >= xmin) dots.push({ x: a, y: 1 + D2.Gm, color: T.warn }); if (b >= xmin) dots.push({ x: b, y: 1 - D2.Gm, color: T.b }); }
      }
      plot(ui.cv, { xmin, xmax, ymin: 0, ymax, xlabel: 'z (m)', ylabel: '|E| / |E0⁺|', series, hlines: hl, vlines: vl, points: dots });
      // power bar
      const gg = ui.pw.prep(), ctx = gg.ctx, w = gg.w, h = gg.h; ctx.clearRect(0, 0, w, h); ctx.fillStyle = T.bg; ctx.fillRect(0, 0, w, h);
      const parts = [['R', D2.R, T.pink], ['T', D2.T, T.b], ['abs', D2.absorbed, T.warn]], x0 = 10, bw = w - 20, y0 = 22, bh = 26; let x = x0;
      ctx.font = '11px ' + T.mono; ctx.textBaseline = 'middle'; ctx.textAlign = 'left'; ctx.fillStyle = T.muted; ctx.fillText('Power fractions of the incident power', x0, 10);
      parts.forEach(p => { const wd = bw * Math.max(0, Math.min(1, p[1])); ctx.fillStyle = p[2]; ctx.fillRect(x, y0, wd, bh); if (wd > 54) { ctx.fillStyle = T.bg; ctx.textAlign = 'center'; ctx.fillText(p[0] + ' ' + g(100 * p[1], 3) + '%', x + wd / 2, y0 + bh / 2); } x += wd; });
      ctx.strokeStyle = T.border; ctx.strokeRect(x0, y0, bw, bh);
    }
    return {
      update, draw,
      get: () => Object.assign(outVals(st), { lay: S.lay, t2: S.t2, t3: S.t3, zero1: ui.zero1.chk.checked ? '1' : '0', zero2: ui.zero2.chk.checked ? '1' : '0', zero3: ui.zero3.chk.checked ? '1' : '0' }),
      set(o) {
        loadSt(st, o); const l = oneOf(o.lay, ['int', 'slab']); if (l) S.lay = l; ['t2', 't3'].forEach(k => { const v = oneOf(o[k], ['mat', 'pec']); if (v) S[k] = v; });
        ['zero1', 'zero2', 'zero3'].forEach(k => { if (o[k] === '0' || o[k] === '1') ui[k].chk.checked = o[k] === '1'; });
        ui.lay.value = S.lay; ui.t2.value = S.t2; ui.t3.value = S.t3;
      },
    };
  }

  /* ---------- tab ---------- */
  const host = { active: false, sub: 'coul' };
  let subs = null; const panels = {}, btns = {};
  const DEFS = [['coul', 'Charge distributions', 'cb_'], ['mag', 'Magnetostatics', 'mg_'], ['far', 'Faraday & displacement current', 'fd_'], ['pw', 'Plane waves & polarization', 'pw_'], ['ni', 'Normal incidence & standing waves', 'ni_']];
  function showSub(id, silent) {
    if (!panels[id]) id = 'coul';
    if (host.sub === 'pw' && id !== 'pw' && subs) subs.pw.stop();
    host.sub = id;
    Object.keys(panels).forEach(k => { panels[k].hidden = k !== id; btns[k].setAttribute('aria-pressed', String(k === id)); btns[k].classList.toggle('active', k === id); });
    if (host.active && subs) subs[id].update();
    if (!silent) FSP.state.touch();
  }
  FSP.registerTab({
    id: 'fields', title: 'Fields & Polarization',
    init(panel) {
      if (typeof document !== 'undefined' && !document.getElementById('fl-style')) {
        const s = document.createElement('style'); s.id = 'fl-style';
        s.textContent = '#panel-fields .ctl select{width:100%;min-width:0;max-width:100%}#panel-fields .ctl{grid-template-columns:76px minmax(0,1fr) auto}#panel-fields fieldset{overflow:hidden}#panel-fields .ctl input[type=text]{max-width:100%}';
        document.head.appendChild(s);
      }
      const bar = FSP.ui.el('div', { class: 'seg row', role: 'group', 'aria-label': 'Fields topic' }); panel.appendChild(bar);
      DEFS.forEach(d => { btns[d[0]] = FSP.ui.el('button', { type: 'button', class: 'seg-btn', 'aria-pressed': 'false', text: d[1], onclick: () => showSub(d[0]) }); bar.appendChild(btns[d[0]]); });
      DEFS.forEach(d => { panels[d[0]] = FSP.ui.el('div', { hidden: '' }); panel.appendChild(panels[d[0]]); });
      subs = { coul: buildCoul(panels.coul), mag: buildMag(panels.mag), far: buildFar(panels.far), pw: buildPW(panels.pw, host), ni: buildNI(panels.ni) };
      FSP.state.bind('fields', {
        get() { const o = { sub: host.sub }; DEFS.forEach(d => { const gt = subs[d[0]].get(); Object.keys(gt).forEach(k => { o[d[2] + k] = gt[k]; }); }); return o; },
        set(o) {
          o = o || {};
          DEFS.forEach(d => { const r = {}; Object.keys(o).forEach(k => { if (k.indexOf(d[2]) === 0) r[k.slice(d[2].length)] = o[k]; }); try { subs[d[0]].set(r); } catch (e) { /* ignore bad state */ } });
          showSub(panels[o.sub] ? o.sub : host.sub, true);
        },
      });
      if (typeof document !== 'undefined') document.addEventListener('visibilitychange', () => { if (!subs) return; if (document.hidden) subs.pw.stop(); else if (host.active && host.sub === 'pw') subs.pw.start(); });
      showSub(host.sub, true);
    },
    activate() { host.active = true; if (subs) subs[host.sub].update(); },
    deactivate() { host.active = false; if (subs) subs.pw.stop(); },
  });

})();
