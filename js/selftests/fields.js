/* Numerical regression checks for fields; also loaded by ?selftest. */
(function () {
  'use strict';
  const { D2R, KE, PI, vlen } = FSP.math.fields._shared;
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

})();
