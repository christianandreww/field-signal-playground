/* Numerical regression checks for edm2; also loaded by ?selftest. */
(function () {
  'use strict';
  const E = FSP.math.edm2;
  const { MU0, PI, S3 } = E._shared;
  FSP.registerTests('edm2', t => {
    const tbl = E.parseTable(E.tableToText(E.DEFAULT_TABLE)), mat = E.makeMaterial({ mode: 'table', H: tbl.H, B: tbl.B });
    t.check('past-paper B–H table parses (7 points)', tbl.ok && tbl.H.length === 7);
    t.check('parseTable rejects non-monotone B, duplicate H, junk, too few points',
      !E.parseTable('100 0.5\n200 0.4').ok && !E.parseTable('100 0.5;100 0.6;300 0.9').ok && !E.parseTable('100 abc').ok && !E.parseTable('100 0.5').ok);
    // hand: B=0.9 lies between (350,0.85) and (500,1.0): H = 350 + 0.05/0.15*150 = 400
    t.check('linear interpolation: H(0.75)=286 (node), H(0.9)=400, B(600)=1.08−… = 1.0457 (hand)',
      t.near(mat.Hof(0.75), 286, 1e-9) && t.near(mat.Hof(0.9), 400, 1e-9) && t.near(mat.Bof(600), 1.0 + 100 / 175 * 0.08, 1e-12));
    // hand: below the first point, line from origin: H(0.2) = 125*0.2/0.4 = 62.5
    t.check('below first point: straight line from origin, H(0.2)=62.5; odd symmetry', t.near(mat.Hof(0.2), 62.5, 1e-12) && t.near(mat.Hof(-0.9), -400, 1e-9));
    t.check('beyond table flagged as extrapolation', mat.outside(1.2) && !mat.outside(1.175));

    // ---- three-limb worked example (default state). Coil on centre limb N=500 (↑); gap 0.2 mm in the right limb;
    // depth 2 cm; widths L 2 cm, C 4 cm, R 2 cm → A_L=A_R=4 cm², A_C=8 cm²; l_L=l_R=0.30 m, l_C=0.10 m.
    // Given B_R = 0.75 T (↓) → Φ_R = 3.0e-4 Wb.
    // Right limb: H = 286 A/m (table); iron length 0.30−0.0002 = 0.2998 m → 85.7428 At;
    //   gap H_g = 0.75/μ0 = 596 831 A/m; F_g = 596831×2e-4 = 119.366 At → D_R = 205.109 At.
    // Left limb (no coil) in parallel: H_L = 205.109/0.30 = 683.697 A/m → B_L = 1.08 + 8.697/325×0.07 = 1.081873 T
    //   → Φ_L = 1.081873×4e-4 = 4.32749e-4 Wb.
    // Flux law: Φ_C = 4.32749e-4 + 3e-4 = 7.32749e-4 Wb → B_C = 0.915937 T → H_C = 350 + 0.065937/0.15×150 = 415.937 A/m
    //   → F_C = 41.5937 At.  NI = 41.5937 + 205.109 = 246.703 At → I = 0.493405 A ≈ 0.493 A.
    const lim = (on, l, w, lg, N, I, s) => E.makeLimb({ on, l, w, d: 0.02, lg, fringe: false, N, I, sense: s });
    const c0 = { mat, limbs: [lim(true, 0.3, 0.02, 0, 0, 0, 1), lim(true, 0.1, 0.04, 0, 500, 0, 1), lim(true, 0.3, 0.02, 2e-4, 0, 0, 1)] };
    const r0 = E.mcSolveForCurrent(c0, 2, -0.75 * 4e-4, 1);
    t.check('three-limb past-paper chain: I = 0.493 A (hand 0.493405)', r0.ok && t.rel(r0.I, 0.493405328, 1e-6) && r0.I.toPrecision(3) === '0.493', r0.ok ? 'I=' + r0.I : r0.msg);
    t.check('three-limb chain intermediates: B_L=1.08187 T, B_C=0.915937 T, H_C=415.937 A/m, U=205.109 At',
      t.rel(Math.abs(r0.drops[0].B), 1.0818731, 1e-6) && t.rel(r0.drops[1].B, 0.9159366, 1e-6) && t.rel(r0.drops[1].H, 415.93657, 1e-6) && t.rel(r0.U, 205.10901, 1e-6));
    // round trip: put I back and solve for fluxes (Newton) → B_R = 0.75 T ↓
    const c0b = { mat, limbs: [c0.limbs[0], Object.assign({}, c0.limbs[1], { I: r0.I }), c0.limbs[2]] };
    const rb = E.mcSolveCurrents(c0b);
    t.check('given currents (Newton) recovers B_R = −0.75 T, converged', rb.ok && rb.converged && rb.method.indexOf('Newton') === 0 && t.rel(rb.drops[2].B, -0.75, 1e-8), rb.ok ? rb.method + ' iter=' + rb.iter : rb.msg);
    t.check('solution satisfies flux law (ΣΦ≈0) and MMF law (same U on every limb)', Math.abs(rb.kcl) < 1e-15 && rb.mmfResid < 1e-7);
    // symmetric: no gap, given B_C = 1.0 T ↑ (Φ_C = 8e-4) → outer limbs 4e-4 each → B = 1.0 T, H = 500 everywhere.
    // NI = 500×0.1 + 500×0.3 = 200 At → I = 0.4 A
    const cs = { mat, limbs: [lim(true, 0.3, 0.02, 0, 0, 0, 1), lim(true, 0.1, 0.04, 0, 500, 0, 1), lim(true, 0.3, 0.02, 0, 0, 0, 1)] };
    const rs = E.mcSolveForCurrent(cs, 1, 8e-4, 1);
    t.check('symmetric core, B_C = 1.0 T: I = 0.400 A (hand: 200 At / 500)', rs.ok && t.rel(rs.I, 0.4, 1e-12) && rs.method.indexOf('symmetry') > 0);
    // asymmetric (gap in right limb) with given B in centre: x = g → U by bisection; verify by Newton round trip
    const ra = E.mcSolveForCurrent(c0, 1, 8e-4, 1), rra = ra.ok && E.mcSolveCurrents({ mat, limbs: [c0.limbs[0], Object.assign({}, c0.limbs[1], { I: ra.I }), c0.limbs[2]] });
    t.check('given B_C with asymmetric outer limbs: current reproduces Φ_C under the forward solve', ra.ok && rra.ok && t.rel(rra.phi[1], 8e-4, 1e-8), ra.ok ? 'I=' + ra.I : ra.msg);
    // coil on left limb, flux given in right limb: other branch of the chain
    const cl = { mat, limbs: [lim(true, 0.3, 0.02, 0, 400, 0, 1), lim(true, 0.1, 0.04, 0, 0, 0, 1), lim(true, 0.3, 0.02, 2e-4, 0, 0, 1)] };
    const rl = E.mcSolveForCurrent(cl, 2, -2e-4, 0), rlb = rl.ok && E.mcSolveCurrents({ mat, limbs: [Object.assign({}, cl.limbs[0], { I: rl.I }), cl.limbs[1], cl.limbs[2]] });
    t.check('coil on left, given Φ_R: forward solve reproduces Φ_R', rl.ok && rlb.ok && t.rel(rlb.phi[2], -2e-4, 1e-8));
    // two coils, Newton vs bisection-on-U (independent algorithm)
    const c2c = { mat, limbs: [lim(true, 0.3, 0.02, 0, 300, 0.6, 1), lim(true, 0.1, 0.04, 0, 0, 0, 1), lim(true, 0.25, 0.03, 3e-4, 200, 1.5, -1)] };
    const n2 = E.mcSolveCurrents(c2c), b2 = E.mcSolveBisect(c2c);
    t.check('two coils: Newton fluxes == bisection-on-U fluxes (1e-9 rel)', n2.ok && b2.ok && [0, 1, 2].every(k => t.near(n2.phi[k], b2.phi[k], 1e-9 * 1e-3)), n2.method);
    // linear (μr=1000) three-limb vs closed-form reluctance network: Φ_C = NI/(R_C + R_L‖R_R)
    // R_L = 0.3/(μ0·1000·4e-4), R_C = 0.1/(μ0·1000·8e-4), R_R = 0.2998/(μ0·1000·4e-4) + 2e-4/(μ0·4e-4); python: Φ_C = 1.058346e-3 Wb
    const lin = E.makeMaterial({ mode: 'const', mur: 1000 });
    const cL = { mat: lin, limbs: [c0.limbs[0], Object.assign({}, c0.limbs[1], { I: 1 }), c0.limbs[2]] }, rL = E.mcSolveCurrents(cL);
    t.check('constant μr: Φ_C = NI/(R_C + R_L‖R_R) = 1.058346 mWb, Φ_L/Φ_R = R_R/R_L', rL.ok && t.rel(rL.phi[1], 1.058346045e-3, 1e-8) && t.rel(rL.phi[0] / rL.phi[2], 0.6613670334 / 0.3969790117, 1e-8));
    // Chapman Ex. 1-2 style 2-limb series: μr=4000, iron 0.40 m, A=12 cm², gap 0.05 cm, N=400, B=0.5 T, no fringing.
    // R_c = 0.4/(4000·μ0·12e-4) = 66 315 A/Wb; R_g = 5e-4/(μ0·12e-4) = 331 573 A/Wb; I = 0.5·12e-4·397 887/400 = 0.59683 A
    const mu4k = E.makeMaterial({ mode: 'const', mur: 4000 });
    const cS = { mat: mu4k, limbs: [E.makeLimb({ on: true, l: 0.2, w: 0.03, d: 0.04, lg: 0, N: 400, I: 0, sense: 1 }), E.makeLimb({ on: false, l: 0.1, w: 0.03, d: 0.04, lg: 0 }), E.makeLimb({ on: true, l: 0.2005, w: 0.03, d: 0.04, lg: 5e-4, N: 0, I: 0, sense: 1 })] };
    const rS = E.mcSolveForCurrent(cS, 2, -0.5 * 12e-4, 0);
    t.check('2-limb series circuit with gap (textbook): I = 0.5968 A', rS.ok && t.rel(rS.I, 0.5968310366, 1e-8), rS.ok ? 'I=' + rS.I : rS.msg);
    // fringing: a=3 cm, b=4 cm, lg=1 mm → A_g = 3.1×4.1 cm²;  R_g = 1e-3/(μ0·0.031·0.041) = 626 101 A/Wb
    const Lf = E.makeLimb({ on: true, l: 0.2, w: 0.03, d: 0.04, lg: 1e-3, fringe: true, N: 0, I: 0, sense: 1 });
    t.check('fringing: A_g = (a+l_g)(b+l_g) = 12.71 cm², R_g = 626 101 A/Wb', t.rel(Lf.Ag, 0.031 * 0.041, 1e-12) && t.rel(E.drop(Lf, mu4k, 1e-4).Fg / 1e-4, 626101.27, 1e-7));
    t.check('invalid circuits rejected (one limb, gap ≥ length, unknown coil N=0)',
      !E.mcSolveCurrents({ mat, limbs: [c0.limbs[0], Object.assign({}, c0.limbs[1], { on: false }), Object.assign({}, c0.limbs[2], { on: false })] }).ok &&
      !E.mcSolveCurrents({ mat, limbs: [E.makeLimb({ on: true, l: 0.001, w: 0.02, d: 0.02, lg: 0.002 }), c0.limbs[1], c0.limbs[2]] }).ok &&
      !E.mcSolveForCurrent(c0, 2, -3e-4, 0).ok);

    // ---- inductance & energy. N=200, lc=0.2, A=4 cm², μr=2000, lg=1 mm:
    // R_c = 0.2/(2000·μ0·4e-4) = 198 944; R_g = 1e-3/(μ0·4e-4) = 1 989 437; L = 200²/2 188 381 = 0.0182784 H
    const mu2k = E.makeMaterial({ mode: 'const', mur: 2000 }), kc = { N: 200, lc: 0.2, A: 4e-4, lg: 1e-3, Ag: 4e-4, mat: mu2k }, ka = E.coreAnalyze(kc, 3);
    t.check('gapped core L = N²/(R_c+R_g) = 18.278 mH (hand)', t.rel(ka.Lsec, 0.0182783573, 1e-8));
    t.check('linear: W = ½LI², W_gap/W_core = R_g/R_c = 10, co-energy = energy', t.rel(ka.W, 0.5 * ka.Lsec * 9, 1e-10) && t.rel(ka.Wgap / ka.Wcore, 10, 1e-10) && t.rel(ka.Wco, ka.W, 1e-10));
    const kt = { N: 300, lc: 0.3, A: 4e-4, lg: 2e-4, Ag: 4e-4, mat }, kta = E.coreAnalyze(kt, 1.5);
    t.check('B–H table core: field energy (exact trapezoids) == ∫ i dλ (Simpson); W + W\' = λi', kta.B > 0.9 && t.rel(kta.W, E.energyAtLambda(kt, kta.lam, 2000), 1e-5) && t.rel(kta.W + kta.Wco, kta.lam * 1.5, 1e-12) && t.rel(kta.Wco, E.coEnergy(kt, 1.5, 2000), 1e-4), 'B=' + kta.B.toFixed(4));
    // ---- force. Ideal iron (μr→∞), N=1000, i=1 A, 2 gaps of 1 mm, A=4 cm²: B = μ0·1000/2e-3 = 0.6283 T;
    // F = 2·B²A/(2μ0) = 125.66 N (attractive) = ½ i² N² μ0 A/(2 x²)
    const rel1 = { N: 1000, lc: 0.2, A: 4e-4, ng: 2, mat: E.makeMaterial({ mode: 'const', mur: 1e12 }) };
    t.check('relay (ideal iron): F = ½ i² dL/dx = −125.66 N = −2·B²A/(2μ0)', t.rel(E.relayForce(rel1, 1, 1e-3), -0.5 * 1e6 * MU0 * 4e-4 / (2 * 1e-6), 1e-9) && t.rel(E.relayForce(rel1, 1, 1e-3), -2 * Math.pow(MU0 * 1000 / 2e-3, 2) * 4e-4 / (2 * MU0), 1e-9));
    const rel2 = { N: 1000, lc: 0.2, A: 4e-4, ng: 2, mat: mu2k };
    t.check('linear: ½ i² dL/dx == finite difference of co-energy ∂W\'/∂x (1e-5)', t.rel(E.relayForce(rel2, 0.8, 1.5e-3), E.relayForceFDco(rel2, 0.8, 1.5e-3), 1e-5));
    const rel3 = { N: 1000, lc: 0.2, A: 4e-4, ng: 2, mat };
    // Three routes to the same force. Exact: F = −n_g·B²A/(2μ0), B = Φ/A (flux linkage held fixed, so only the gap energy moves).
    // The two finite-difference routes are checked against it. The co-energy route integrates λ(i') (piecewise linear: kinks at the
    // table breakpoints) by Simpson, whose error at a kink is O(Δi²·Δslope); dividing by the small difference 2h amplifies it, so a
    // 0.1 % step and n=400 gave 5e-4 error. With h = 1 % of x and n=2000 the observed error is 3e-5, hence the 1e-4 tolerance.
    const fex = E.relayForce(rel3, 0.6, 0.3e-3), Bx = E.coreFlux(E.relayCore(rel3, 0.3e-3), 0.6) / 4e-4;
    const fco = E.relayForceFDco(rel3, 0.6, 0.3e-3, 3e-6, 2000), fen = E.relayForceFDenergy(rel3, 0.6, 0.3e-3);
    t.check('nonlinear (table): F = −n_g B²A/(2μ0) == ∂W\'/∂x|i == −∂W/∂x|λ (three routes, 1e-4)', fex < 0 && t.rel(fex, -2 * Bx * Bx * 4e-4 / (2 * MU0), 1e-12) && t.rel(fco, fex, 1e-4) && t.rel(fen, fex, 1e-4), fex.toFixed(4) + ' / ' + fco.toFixed(4) + ' / ' + fen.toFixed(4));
    const ip = E.relayPullIn(rel2, 2e-3, 30), ipn = E.relayPullIn(rel3, 2e-3, 30);
    t.check('pull-in current: |F(i_pull, x0)| = F_s (linear closed form and nonlinear bisection)', t.rel(Math.abs(E.relayForce(rel2, ip, 2e-3)), 30, 1e-9) && t.rel(Math.abs(E.relayForce(rel3, ipn, 2e-3)), 30, 1e-6));
    const th = 0.4, h = 1e-6, fdT = 0.5 * 4 * (E.rotL(0.2, 0.05, th + h) - E.rotL(0.2, 0.05, th - h)) / (2 * h);
    t.check('rotary: T = ½ i² dL/dθ = −½i²(Ld−Lq)sin2θ (finite-difference check); zero at aligned/unaligned', t.rel(E.rotTorque(0.2, 0.05, 2, th), fdT, 1e-7) && Math.abs(E.rotTorque(0.2, 0.05, 2, 0)) < 1e-15 && Math.abs(E.rotTorque(0.2, 0.05, 2, PI / 2)) < 1e-15);
    // 4.44 rule: 230 V, 50 Hz, 200 turns: Φmax = 230/(√2·π·50·200) = 5.17682 mWb; 4.44·f·N·Φ within 0.07 %
    const sf = E.sineFlux(230, 50, 200, 25e-4);
    t.check('E = 4.44 f N Φmax: Φmax = 5.1768 mWb, Bmax = 2.0707 T for 25 cm²', t.rel(sf.phim, 5.17681882e-3, 1e-8) && t.rel(4.44 * 50 * 200 * sf.phim, 230, 7e-4) && t.rel(sf.Bm, 5.17681882e-3 / 25e-4, 1e-8));
    // separation, hand: P1=100 W @50 Hz, P2=130 W @60 Hz: P/f = 2, 2.16667 → b = 0.0166667, a = 1.16667;
    // at 50 Hz: P_h = 58.333 W, P_e = 41.667 W
    const sp = E.separateLosses({ f1: 50, P1: 100, f2: 60, P2: 130, Bm: 1.2, n: 1.6, vol: 1e-3 });
    t.check('loss separation (hand): a = 1.16667 J, b = 0.0166667 J·s, P_h(50) = 58.333 W, P_e(50) = 41.667 W', sp.ok && t.rel(sp.a, 7 / 6, 1e-12) && t.rel(sp.b, 1 / 60, 1e-12) && t.rel(sp.Ph1, 175 / 3, 1e-12) && t.rel(sp.Pe1, 125 / 3, 1e-12));
    const kh0 = 120, ke0 = 0.05, gen = f => E.coreLoss({ kh: kh0, ke: ke0, n: 1.7, f, Bm: 1.3, vol: 2e-3 }).P;
    const sp2 = E.separateLosses({ f1: 25, P1: gen(25), f2: 70, P2: gen(70), Bm: 1.3, n: 1.7, vol: 2e-3 });
    t.check('loss separation recovers k_h, k_e from synthetic data', sp2.ok && t.rel(sp2.kh, kh0, 1e-10) && t.rel(sp2.ke, ke0, 1e-10));
    t.check('loss separation: equal frequencies / negative term rejected', !E.separateLosses({ f1: 50, P1: 100, f2: 50, P2: 120, Bm: 1, n: 1.6, vol: 1 }).ok && !E.separateLosses({ f1: 50, P1: 100, f2: 60, P2: 110, Bm: 1, n: 1.6, vol: 1 }).ok);

    // ---- autotransformer: 5 kVA 220/110 V → 330/220 V step-down (common = 220 V winding):
    // S_auto = 5·330/110 = 15 kVA; I_H = 15000/330 = 45.4545 A (= 5000/110, series winding rated);
    // I_L = 15000/220 = 68.1818 A; common winding I_L − I_H = 22.7273 A (= 5000/220); conducted 10 kVA, transformed 5 kVA
    const au = E.auto({ V1: 220, V2: 110, S: 5000, common: 'V1', dir: 'down', x: 1, pf: 0.8, Pcore: 50, Pcu: 100 });
    t.check('auto 5 kVA 220/110 → 330/220: 15 kVA, I_H = 45.45 A, I_L = 68.18 A, I_common = 22.73 A',
      au.ok && au.VH === 330 && au.VL === 220 && t.rel(au.Sa, 15000, 1e-12) && t.rel(au.IH, 5000 / 110, 1e-12) && t.rel(au.IL, 15000 / 220, 1e-12) && t.rel(au.Icom, 5000 / 220, 1e-12));
    t.check('auto power split: transformed 5 kVA + conducted 10 kVA', t.rel(au.Strans, 5000, 1e-12) && t.rel(au.Scond, 10000, 1e-12));
    // common = 110 V winding → 330/110 V: S = 5·330/220 = 7.5 kVA
    const au2 = E.auto({ V1: 220, V2: 110, S: 5000, common: 'V2', dir: 'up' });
    t.check('auto with the 110 V winding common: 330/110 V, 7.5 kVA', au2.ok && au2.VL === 110 && t.rel(au2.Sa, 7500, 1e-12));
    // efficiencies (hand): 2-w: 4000/(4000+150) = 0.963855; auto: 12000/12150 = 0.987654
    t.check('auto efficiency: 96.386 % (2-winding) → 98.765 % (auto), same losses', t.rel(au.eta2w, 4000 / 4150, 1e-12) && t.rel(au.etaAuto, 12000 / 12150, 1e-12));
    t.check('auto rejects zero/negative rating', !E.auto({ V1: 0, V2: 110, S: 5000 }).ok);

    // ---- three-phase: 150 kVA bank Y-Δ, 11 kV / 415 V.
    // HV winding = 11000/√3 = 6350.85 V, LV winding = 415 V → a = 15.3033; line ratio 11000/415 = √3·a.
    // I_HL = 150e3/(√3·11e3) = 7.8730 A (= HV winding); I_LL = 150e3/(√3·415) = 208.68 A; LV winding 120.482 A = 50e3/415
    const yd = E.threePhase({ conn: 'YD', S: 150e3, VHL: 11000, VLL: 415, conv: 'ANSI', Req: 12, Xeq: 30 });
    t.check('Y-Δ: line-voltage ratio = √3·a (a = winding ratio)', yd.ok && t.rel(yd.lineRatio, S3 * yd.a, 1e-12) && t.rel(yd.a, 11000 / S3 / 415, 1e-12));
    t.check('Y-Δ currents: I_HL = 7.873 A, I_LL = 208.68 A, LV winding 120.48 A (= 50 kVA/415 V)', t.rel(yd.IHL, 150e3 / (S3 * 11000), 1e-12) && t.rel(yd.ILw, 50e3 / 415, 1e-12) && t.rel(yd.IHw, yd.IHL, 1e-12));
    t.check('Y-Δ phase shift: −30° (ANSI, LV lags) / +30° (IEC clock 11)', yd.shift === -30 && E.threePhase({ conn: 'YD', S: 150e3, VHL: 11000, VLL: 415, conv: 'IEC11' }).shift === 30);
    const dy = E.threePhase({ conn: 'DY', S: 150e3, VHL: 11000, VLL: 415 }), dd = E.threePhase({ conn: 'DD', S: 150e3, VHL: 11000, VLL: 415 });
    t.check('Δ-Y ratio = a/√3; Δ-Δ and Y-Y: ratio = a, no shift', t.rel(dy.lineRatio, dy.a / S3, 1e-12) && t.rel(dd.lineRatio, dd.a, 1e-12) && dd.shift === 0 && E.threePhase({ conn: 'YY', S: 1, VHL: 11000, VLL: 415 }).shift === 0);
    const vv = E.threePhase({ conn: 'VV', S: 150e3, VHL: 11000, VLL: 415 });
    t.check('open-Δ (V-V): capacity = √3·S_unit = 57.7 % of the Δ-Δ bank, 86.6 % unit utilisation', t.rel(vv.vvFraction, 1 / S3, 1e-12) && t.near(vv.vvFraction * 100, 57.735, 1e-3) && t.rel(vv.vvUtil, S3 / 2, 1e-12));
    // per-phase equivalent: Y on HV → Z_Y(HV) = Z_w = 12 + j30; referred to LV × (415/11000)² = 0.0170805 + j0.0427012
    t.check('per-phase equivalent: Y-HV → Z_Y = Z_w; referred to LV by (V_LL/V_HL)²; Δ-HV → Z_w/3',
      t.rel(yd.ZYH.re, 12, 1e-12) && t.rel(yd.ZYL.im, 30 * (415 / 11000) ** 2, 1e-12) && t.rel(E.threePhase({ conn: 'DY', S: 150e3, VHL: 11000, VLL: 415, Req: 12, Xeq: 30 }).ZYH.re, 4, 1e-12));

    // ---- all-day efficiency (hand): 100 kVA, P_core 1 kW, P_cu,FL 1.5 kW; 6 h 80 kW pf 0.8 (x=1), 10 h 40 kW pf 1 (x=0.4), 8 h idle.
    // E_out = 480+400 = 880 kWh; E_cu = 6·1.5 + 10·1.5·0.16 = 11.4 kWh; E_core = 24 kWh → η = 880/915.4 = 96.133 %
    const pr = E.parseProfile('6 80 0.8\n10 40 1\n8 0 1'), ad = pr.ok && E.allDay({ S: 100, Pcore: 1, Pcu: 1.5, rows: pr.rows });
    t.check('all-day efficiency = 880/915.4 = 96.133 % (hand)', ad && ad.ok && t.rel(ad.eta, 880 / 915.4, 1e-12) && t.rel(ad.Ecu, 11.4, 1e-12) && t.rel(ad.Ecore, 24, 1e-12));
    const me = E.maxEffPoint(100, 1, 1.5, 0.8);
    t.check('max-efficiency load fraction x = √(P_core/P_cu) = 0.8165 and η(x) is the maximum', t.rel(me.x, 0.8164965809, 1e-9) && t.rel(me.eta, E.etaAt(100, 1, 1.5, me.x, 0.8), 1e-12) && E.etaAt(100, 1, 1.5, me.x * 1.05, 0.8) < me.eta && E.etaAt(100, 1, 1.5, me.x * 0.95, 0.8) < me.eta);
    t.check('profile validation: > 24 h, pf > 1, junk rejected', !E.parseProfile('20 10 1\n5 10 1').ok && !E.parseProfile('2 10 1.2').ok && !E.parseProfile('two 10').ok);
  });


})();
