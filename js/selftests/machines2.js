/* Numerical regression checks for rotating; also loaded by ?selftest. */
(function () {
  'use strict';
  const { PI, S3, cabs, fin } = FSP.math.rotating._shared;
  FSP.registerTests('rotating', t => {
    const M = FSP.math.rotating;
    const P0 = { VL: 400, conn: 'Y', f: 50, poles: 4, R1: 0.7, X1: 1.2, R2: 0.6, X2: 1.2, Xm: 35, Prot: 150 };
    // --- induction ---
    t.check('ns = 1500 rpm for 4 poles, 50 Hz', M.syncRpm(50, 4) === 1500);
    t.check('ns = 1200 rpm for 6 poles, 60 Hz; 3000 rpm for 2 poles, 50 Hz', M.syncRpm(60, 6) === 1200 && M.syncRpm(50, 2) === 3000);
    t.check('ws = 2*pi*ns/60', t.rel(M.wsync(50, 4), 2 * PI * 1500 / 60, 1e-12));
    const th = M.imThevenin(P0);
    t.check('Thevenin values match independent calc (|Vth|=223.2429 V, Zth=0.654116+j1.172870)', t.rel(th.VthMag, 223.24290087793207, 1e-9) && t.near(th.Rth, 0.6541157804001739, 1e-9) && t.near(th.Xth, 1.172869642162434, 1e-9));
    t.check('T(s=0) = 0 (both methods)', M.imTorqueTh(P0, 0) === 0 && M.imTorqueFull(P0, 0) === 0);
    let worst = 0; [0.04, 0.3, 1, -0.04, 1.7].forEach(s => { const a = M.imTorqueTh(P0, s), b = M.imTorqueFull(P0, s); worst = Math.max(worst, Math.abs(a - b)); });
    t.check('Thevenin torque == direct full-circuit torque at 5 slips (<1e-9)', worst < 1e-9, 'max err ' + worst.toExponential(2));
    t.check('T(0.04) = 56.9541 N·m (independent numpy)', t.rel(M.imTorqueTh(P0, 0.04), 56.954108468870785, 1e-9));
    t.check('T(1) = 79.2822 N·m starting torque; T(2) = 43.6560', t.rel(M.imTorqueTh(P0, 1), 79.28215227042001, 1e-9) && t.rel(M.imTorqueTh(P0, 2), 43.65599982620026, 1e-9));
    t.check('T > 0 for 0<s<1 (50 points)', Array.from({ length: 50 }, (_, i) => (i + 1) / 51).every(s => M.imTorqueTh(P0, s) > 0));
    t.check('T < 0 for s<0 (50 points)', Array.from({ length: 50 }, (_, i) => -(i + 1) / 51).every(s => M.imTorqueTh(P0, s) < 0));
    t.check('T(-0.04) = -67.526 (generator)', t.rel(M.imTorqueTh(P0, -0.04), -67.52603779639693, 1e-9));
    const bd = M.imBreakdown(P0, th), smForm = P0.R2 / Math.sqrt(th.Rth * th.Rth + (th.Xth + P0.X2) ** 2);
    t.check('breakdown slip = R2/sqrt(Rth^2+(Xth+X2)^2) = 0.243766', t.rel(bd.sm, smForm, 1e-12) && t.near(bd.sm, 0.24376596826587157, 1e-12));
    t.check('breakdown torque = 152.7566 N·m and equals T(sm)', t.rel(bd.Tmax, 152.75657675282008, 1e-9) && t.rel(M.imTorqueTh(P0, bd.sm), bd.Tmax, 1e-12));
    t.check('generator breakdown at -sm: -263.333 N·m', t.rel(bd.Tgen, -263.33325981545613, 1e-9) && t.rel(M.imTorqueTh(P0, bd.sg), bd.Tgen, 1e-12));
    t.check('T(sm) is a maximum (neighbours smaller)', M.imTorqueTh(P0, bd.sm * 0.98) < bd.Tmax && M.imTorqueTh(P0, bd.sm * 1.02) < bd.Tmax);
    const op = M.imOp(P0, 0.04);
    t.check('op point: T and Pag/ws agree; Pmech = (1-s)Pag; Pcu2 = s*Pag', t.rel(op.T, 56.954108468870785, 1e-9) && t.rel(op.Pmech, (1 - 0.04) * op.Pag, 1e-12) && t.rel(op.Pcu2, 0.04 * op.Pag, 1e-12));
    t.check('op point: Pin = Pcu1 + Pag (energy balance)', t.rel(op.Pin, op.Pcu1 + op.Pag, 1e-12), 'Pin=' + op.Pin.toFixed(1));
    t.check('op point: |I1| = 15.7857 A, speed 1440 rpm, mode motor', t.rel(op.I1, 15.785727132098748, 1e-9) && t.near(op.nrpm, 1440, 1e-9) && op.mode === 'motor');
    t.check('efficiency in (0,1) for rated motor, NaN for braking', op.eff > 0.5 && op.eff < 1 && !fin(M.imOp(P0, 1.5).eff));
    t.check('mode labels', M.imMode(-0.1) === 'generator' && M.imMode(0.5) === 'motor' && M.imMode(1.5).indexOf('braking') === 0);
    t.check('generator op point: negative torque, negative Pin, efficiency finite', (() => { const g = M.imOp(P0, -0.04); return g.T < 0 && g.Pin < 0 && fin(g.eff) && g.eff > 0 && g.eff < 1; })());
    t.check('Delta connection uses VL as phase voltage', M.imVphase(Object.assign({}, P0, { conn: 'D' })) === 400 && t.rel(M.imVphase(P0), 400 / S3, 1e-15));
    t.check('torque scales with V^2 (halving V -> T/4)', t.rel(M.imTorqueTh(Object.assign({}, P0, { VL: 200 }), 0.1), M.imTorqueTh(P0, 0.1) / 4, 1e-9));
    t.check('sweep covers slip -1..2 and all finite', (() => { const sw = M.imSweep(P0, 300); return sw.s[0] === -1 && t.near(sw.s[300], 2, 1e-12) && sw.T.every(fin); })());
    // --- DC ---
    const D0 = { type: 'shunt', V: 240, Vf: 240, Ra: 0.5, Rf: 120, Rs: 0.1, Laf: 0.8, Vb: 2, Prot: 300 };
    const d1 = M.dcSolve(D0, 50);
    t.check('shunt at Ia=50: If=2, Ea=213, w=133.125, Tem=80 (hand calc)', t.near(d1.If, 2, 1e-12) && t.near(d1.Ea, 213, 1e-12) && t.near(d1.w, 133.125, 1e-12) && t.near(d1.Tem, 80, 1e-12));
    t.check('shunt: Pin = 12480 W, Pout = 10350 W, losses = 2130 W', t.near(d1.Pin, 12480, 1e-9) && t.near(d1.Pout, 10350, 1e-9) && t.near(d1.losses, 2130, 1e-9));
    const bal = [0, 5, 25, 50, 75].map(ia => ['sep', 'shunt', 'series'].map(ty => Math.abs(M.dcSolve(Object.assign({}, D0, { type: ty, Laf: ty === 'series' ? 0.032 : 0.8 }), ia).balance))).flat();
    t.check('power balance Pin = Pout + losses (<1e-9) for all three types', Math.max.apply(null, bal) < 1e-9, 'max ' + Math.max.apply(null, bal).toExponential(2));
    t.check('Pem = Ea*Ia = Tem*w', t.rel(d1.Pem, d1.Tem * d1.w, 1e-12));
    const sp = D0, pts = [10, 20, 30, 40].map(ia => M.dcSolve(Object.assign({}, sp, { type: 'sep' }), ia));
    const sl1 = (pts[1].w - pts[0].w) / (pts[1].Tem - pts[0].Tem), sl2 = (pts[3].w - pts[2].w) / (pts[3].Tem - pts[2].Tem);
    t.check('separately-excited speed falls linearly with torque', t.rel(sl1, sl2, 1e-9), 'slope ' + sl1.toFixed(5));
    t.check('separately-excited slope = -Ra/K^2', t.rel(sl1, -D0.Ra / ((D0.Laf * 240 / 120) ** 2), 1e-9));
    const sr = Object.assign({}, D0, { type: 'series', Laf: 0.032 }), s50 = M.dcSolve(sr, 50);
    t.check('series: Tem = Laf*Ia^2, w = (V-Vb-Ia(Ra+Rs))/(Laf*Ia) = 130 rad/s', t.near(s50.Tem, 0.032 * 2500, 1e-12) && t.near(s50.w, 130, 1e-12));
    const s0 = M.dcSolve(sr, 0);
    t.check('series no-load guarded: runaway flag, speed NaN (not Infinity)', s0.runaway === true && Number.isNaN(s0.w) && !fin(s0.rpm) && fin(s0.Pin));
    t.check('series speed rises as load falls (Ia 50 -> 5)', M.dcSolve(sr, 5).w > 5 * s50.w);
    t.check('shunt speed: no-load 240/1.6 = 150 rad/s, rated 133.125 (hand calc)', t.near(M.dcSolve(D0, 0).w, 150, 1e-12) && t.near(M.dcSolve(D0, 0).w - d1.w, 16.875, 1e-12));
    t.check('stalled flag when Ia beyond stall current', M.dcSolve(D0, M.dcStallCurrent(D0) * 1.1).stalled && !M.dcSolve(D0, 10).stalled);
    t.check('efficiency NaN at no-load (Pout<=0) and in (0,1) at rated', !fin(M.dcSolve(D0, 0).eff) && d1.eff > 0.8 && d1.eff < 1);
    // --- three-phase ---
    let maxSum = 0, maxMag = 0, maxErr = 0;
    for (let k = 0; k < 2000; k++) {
      const tt = k * 7.3e-5; [1, -1].forEach(seq => {
        const i = M.phases(tt, 3.7, 50, seq); maxSum = Math.max(maxSum, Math.abs(i[0] + i[1] + i[2]));
        const F = M.mmfVector(tt, 3.7, 50, seq, 1); maxMag = Math.max(maxMag, Math.abs(F.mag - 1.5 * 3.7));
      });
    }
    t.check('phase sum = 0 at every sample (<1e-12), both sequences', maxSum < 1e-12, maxSum.toExponential(2));
    t.check('resultant MMF magnitude = 1.5 x peak (<1e-12) with N=1', maxMag < 1e-12, maxMag.toExponential(2));
    t.check('resultant MMF scales with turns N', t.rel(M.mmfVector(0.003, 2, 50, 1, 40).mag, 1.5 * 2 * 40, 1e-12));
    const w1 = M.rotationRate(0.0013, 1e-5, 2, 50, 1), w2 = M.rotationRate(0.0013, 1e-5, 2, 50, -1);
    t.check('MMF rotates at supply frequency: +2*pi*50 rad/s (a-b-c)', t.rel(w1, 2 * PI * 50, 1e-9));
    t.check('swapping two phases reverses rotation: -2*pi*50 rad/s', t.rel(w2, -2 * PI * 50, 1e-9));
    t.check('rotation rate scales with frequency (60 Hz)', t.rel(M.rotationRate(0.0007, 1e-5, 1, 60, 1), 2 * PI * 60, 1e-9));
    t.check('MMF angle = omega*t exactly for a-b-c', t.near(M.mmfVector(0.002, 1, 50, 1, 1).ang, 2 * PI * 50 * 0.002 - 2 * PI, 1e-12) || t.near(M.mmfVector(0.002, 1, 50, 1, 1).ang, 2 * PI * 50 * 0.002, 1e-12));
    t.check('spatial MMF distribution = 1.5 Im cos(th - wt) and its peak is at the vector angle', (() => {
      let e = 0; for (let k = 0; k < 24; k++) { const th = k * 0.27, tt = 0.0031; e = Math.max(e, Math.abs(M.mmfAt(th, tt, 2, 50, 1, 1) - 3 * Math.cos(th - 2 * PI * 50 * tt))); } return e < 1e-12;
    })());
    t.check('single winding alone pulsates (no rotation): peak = Im, not 1.5 Im', (() => { let mx = 0; for (let k = 0; k < 400; k++) mx = Math.max(mx, Math.abs(Math.cos(2 * PI * 50 * k * 5e-5))); return t.near(mx, 1, 1e-3); })());
    const qY = M.threePhaseQuantities('Y', 400, 8, 6), qD = M.threePhaseQuantities('D', 400, 8, 6);
    t.check('Y: Vph = VL/sqrt3, IL = Iph', t.rel(qY.Vph, 400 / S3, 1e-12) && qY.IL === qY.Iph);
    t.check('Delta: Vph = VL, IL = sqrt3 Iph', qD.Vph === 400 && t.rel(qD.IL, S3 * qD.Iph, 1e-12));
    t.check('P = sqrt3 VL IL cos(phi) = 3 Vph Iph cos(phi) for Y and Delta', t.rel(qY.P, S3 * 400 * qY.IL * 0.8, 1e-12) && t.rel(qD.P, S3 * 400 * qD.IL * 0.8, 1e-12) && t.rel(qD.P, 3 * qD.Vph * qD.Iph * 0.8, 1e-12));
    t.check('Delta draws 3x the power of Y for the same load impedance', t.rel(qD.P, 3 * qY.P, 1e-12));
    t.check('S^2 = P^2 + Q^2', t.rel(qY.S * qY.S, qY.P * qY.P + qY.Q * qY.Q, 1e-12));
    t.check('zero-impedance load gives NaN (guarded), not Infinity', !fin(M.threePhaseQuantities('Y', 400, 0, 0).Iph));
    const sY = M.phasorSets('Y', 400, 8, 6, 1), sD = M.phasorSets('D', 400, 8, 6, 1);
    const ang = z => Math.atan2(z.im, z.re) * 180 / PI;
    t.check('Y: |Vab| = sqrt3 |Van| and Vab leads Van by 30 deg', t.rel(cabs(sY.line[0]), S3 * cabs(sY.phase[0]), 1e-12) && t.near(ang(sY.line[0]) - ang(sY.phase[0]), 30, 1e-9));
    t.check('Delta: |Ia| = sqrt3 |Iab| and Ia lags Iab by 30 deg', t.rel(cabs(sD.line[0]), S3 * cabs(sD.phase[0]), 1e-12) && t.near(ang(sD.line[0]) - ang(sD.phase[0]), -30, 1e-9));
    t.check('line phasors also sum to zero (balanced)', Math.abs(sY.line[0].re + sY.line[1].re + sY.line[2].re) < 1e-9 && Math.abs(sD.line[0].im + sD.line[1].im + sD.line[2].im) < 1e-9);
    t.check('swapped sequence: Vbn at +120 deg', t.near(ang(M.phasorSets('Y', 400, 8, 6, -1).phase[1]), 120, 1e-9));
    // --- induction speed control (tutorial-style hand example): delta 240 V, 4 poles, 50 Hz, R1 = 0, X1 = 1, Xm = 9, R2' = 0.2, X2' = 0.1 Ω
    //   Zth = jXm·jX1/(j(X1+Xm)) = j0.9 Ω, Vth = 240·9/10 = 216 V, Xth+X2' = 1.0 Ω, ωs = 50π = 157.08 rad/s
    //   Tmax = 3·216²/(2·157.08·1.0) = 139968/(100π) = 445.5 N·m, sm = R2'/(Xth+X2') = 0.2
    const H0 = { VL: 240, conn: 'D', f: 50, poles: 4, R1: 0, X1: 1, R2: 0.2, X2: 0.1, Xm: 9, Prot: 0 };
    const h50 = M.imBreakdown(H0);
    t.check('hand example: Vth = 216 V, Tmax = 139968/(100π) = 445.5 N·m, sm = 0.2', t.rel(M.imThevenin(H0).VthMag, 216, 1e-12) && t.rel(h50.Tmax, 139968 / (100 * PI), 1e-12) && t.near(h50.sm, 0.2, 1e-12));
    // V/f with R1 = 0:  Vth ∝ V' = aV, ωs' = aωs, Xth+X2' = a(Xth+X2')  =>  Tmax' = 3a²Vth²/(2aωs·a(X)) = Tmax  (sm' = sm/a)
    [0.5, 0.2, 0.8].forEach(a => {
      const c = M.imControl(H0, 'vf', { f: 50 * a, boost: 0 }), b = M.imBreakdown(c.p2);
      t.check('V/f, R1 = 0: breakdown torque identical at ' + (50 * a) + ' Hz and 50 Hz; sm = 0.2/' + a, c.ok && t.rel(b.Tmax, h50.Tmax, 1e-12) && t.rel(b.sm, 0.2 / a, 1e-12));
    });
    t.check('V/f, R1 = 0: speed drop at breakdown ns\'·sm is the same at 25 Hz (750·0.4) and 50 Hz (1500·0.2) = 300 rpm', (() => { const c = M.imControl(H0, 'vf', { f: 25, boost: 0 }); return t.rel(M.syncRpm(25, 4) * M.imBreakdown(c.p2).sm, 300, 1e-12) && t.rel(M.syncRpm(50, 4) * h50.sm, 300, 1e-12); })());
    t.check('V/f: T(a·s) at 50 Hz = T(s) at 25 Hz (R1 = 0), i.e. the curve is the same shape with slip scaled by 1/a', (() => { const c = M.imControl(H0, 'vf', { f: 25, boost: 0 }); return t.rel(M.imTorqueTh(c.p2, 0.1), M.imTorqueTh(H0, 0.05), 1e-12); })());
    const cv0 = M.imControl(P0, 'vf', { f: 25, boost: 0 }), cvb = M.imControl(P0, 'vf', { f: 25, boost: 0.15 });
    const r0 = M.imBreakdown(cv0.p2).Tmax / M.imBreakdown(P0).Tmax, rb = M.imBreakdown(cvb.p2).Tmax / M.imBreakdown(P0).Tmax;
    t.check('V/f with R1 > 0: Tmax at 25 Hz falls below the 50 Hz value (only approximately constant); boost recovers it', r0 < 0.99 && r0 > 0.5 && rb > r0, 'ratio ' + r0.toFixed(3) + ' → ' + rb.toFixed(3) + ' with 15 % boost');
    t.check('V/f voltage law: a = 0.5, b = 0.2 → 0.5 + 0.2·0.5 = 0.6; b = 0 → a; a ≥ 1 → rated; machine VL scaled and X ∝ f', (() => { const c = M.imControl(P0, 'vf', { f: 25, boost: 0.2 }); return t.near(M.imVfFrac(0.5, 0.2), 0.6, 1e-15) && t.near(M.imVfFrac(0.5, 0), 0.5, 1e-15) && M.imVfFrac(1.4, 0.3) === 1 && t.rel(c.p2.VL, 400 * 0.6, 1e-12) && t.rel(c.p2.X1, 0.6, 1e-12) && t.rel(c.p2.Xm, 17.5, 1e-12) && c.p2.R2 === P0.R2; })());
    t.check('V/f above base frequency: V constant (field weakening), Tmax falls ∝ 1/a² (R1 = 0): 75 Hz → 445.5/2.25', (() => { const c = M.imControl(H0, 'vf', { f: 75, boost: 0 }); return t.rel(M.imBreakdown(c.p2).Tmax, h50.Tmax / 2.25, 1e-12); })());
    // added rotor resistance: K = |Zth + jX2'| has no R2', so sm = R2'/K ∝ R2'; Tmax = 3Vth²/(2ωs(Rth+K)) unchanged
    const cr = M.imControl(P0, 'rot', { Rext: 0.9 }), bP = M.imBreakdown(P0), bR = M.imBreakdown(cr.p2);
    t.check('rotor resistance: sm scales by (R2\'+Rext)/R2\' = 1.5/0.6 = 2.5, Tmax unchanged (R1 > 0 case)', t.rel(bR.sm, 2.5 * bP.sm, 1e-12) && t.rel(bR.Tmax, bP.Tmax, 1e-12));
    // T depends on s only through R2'/s: R2'/s unchanged when s_new = s_old·(R2+Rext)/R2 = 2.5·s_old
    t.check('rotor resistance: T_new(2.5·s) = T_old(s) (same curve, slip axis stretched by (R2+Rext)/R2)', t.rel(M.imTorqueTh(cr.p2, 0.125), M.imTorqueTh(P0, 0.05), 1e-12) && t.rel(M.imTorqueTh(cr.p2, 2.5), M.imTorqueTh(P0, 1), 1e-12));
    t.check('rotor resistance of Rext = sm·K − R2 makes the starting torque the maximum: T(1) = Tmax', (() => { const K = M.imThevenin(P0).K, c = M.imControl(P0, 'rot', { Rext: K - P0.R2 }); return t.rel(M.imTorqueTh(c.p2, 1), bP.Tmax, 1e-9); })());
    // stator voltage control: Vth ∝ V so T ∝ V² at fixed slip, slip at Tmax unchanged
    const ck = M.imControl(P0, 'volt', { k: 0.8 });
    t.check('voltage control: T(s) scales by (V/V0)² = 0.64 at s = 0.03, 0.2, 1; sm unchanged; Tmax × 0.64', [0.03, 0.2, 1].every(s => t.rel(M.imTorqueTh(ck.p2, s), 0.64 * M.imTorqueTh(P0, s), 1e-12)) && t.rel(M.imBreakdown(ck.p2).sm, bP.sm, 1e-12) && t.rel(M.imBreakdown(ck.p2).Tmax, 0.64 * bP.Tmax, 1e-12));
    t.check('speed-control inputs validated (k ≤ 0, Rext < 0, f ≤ 0, boost > 1, unknown method) and p not mutated', !M.imControl(P0, 'volt', { k: 0 }).ok && !M.imControl(P0, 'rot', { Rext: -1 }).ok && !M.imControl(P0, 'vf', { f: -5, boost: 0 }).ok && !M.imControl(P0, 'vf', { f: 25, boost: 2 }).ok && !M.imControl(P0, 'xyz', {}).ok && P0.VL === 400 && P0.R2 === 0.6 && P0.f === 50);
    // operating point, constant torque TL = 200 N·m on the hand example.  With R1 = 0:  T = 3Vth²R2·s/(ωs(R2² + s²X²))  (X = Xth+X2' = 1)
    //   TL·ωs·X²·s² − 3Vth²R2·s + TL·ωs·R2² = 0   → a = 31 415.9, b = 27 993.6, c = 1256.64, s = (b − √(b²−4ac))/(2a) = 0.04741 (stable, smaller root)
    const wsH = 50 * PI, TLc = 200, Aq = TLc * wsH, Bq = 3 * 216 * 216 * 0.2, Cq = TLc * wsH * 0.04, sq = (Bq - Math.sqrt(Bq * Bq - 4 * Aq * Cq)) / (2 * Aq);
    const opc = M.imOperatingPoint(H0, { type: 'const', T: TLc, n: 1500 });
    t.check('constant load 200 N·m: slip from the quadratic (0.04741), n = 1500(1−s) = 1428.9 rpm', opc.ok && t.near(sq, 0.04741, 2e-5) && t.rel(opc.s, sq, 1e-9) && t.rel(opc.n, 1500 * (1 - sq), 1e-9) && t.rel(opc.T, 200, 1e-9));
    const op25 = M.imOperatingPoint(M.imControl(H0, 'vf', { f: 25, boost: 0 }).p2, { type: 'const', T: TLc, n: 1500 });
    t.check('V/f (R1 = 0), constant torque: the speed drop below ns\' is the same at 25 Hz and 50 Hz (Δn = 71.1 rpm); n = 750 − Δn', op25.ok && t.rel(750 - op25.n, 1500 - opc.n, 1e-9) && t.rel(op25.n, 750 - (1500 - opc.n), 1e-9));
    // fan load  TL = 200 (n/1450)²;  verify against a brute-force scan for the stable crossing (different algorithm)
    const fan = { type: 'fan', T: 200, n: 1450 }, fan2 = { type: 'fan', T: 60, n: 1450 }, opf = M.imOperatingPoint(P0, fan);
    let sScan = NaN; { const th2 = M.imThevenin(P0), smx = M.imBreakdown(P0, th2).sm; let prev = null; for (let i = 1; i <= 200000; i++) { const s = smx * i / 200000, gg = M.imTorqueTh(P0, s, th2) - 200 * Math.pow(1500 * (1 - s) / 1450, 2); if (prev !== null && prev < 0 && gg >= 0) { sScan = s; break; } prev = gg; } }
    t.check('fan load T_L = 200(n/1450)²: bisection root equals brute-force scan (to the 1e-5 grid) and T(s) = T_L(n)', opf.ok && Math.abs(opf.s - sScan) < 2e-6 && t.rel(opf.T, opf.TL, 1e-9) && t.rel(opf.TL, 200 * Math.pow(opf.n / 1450, 2), 1e-12), 's=' + opf.s);
    t.check('voltage control, fan load: lower voltage → lower speed (k = 1, 0.9, 0.8; fan 60 N·m at 1450 rpm), all stable', (() => { const n = [1, 0.9, 0.8].map(k => M.imOperatingPoint(M.imControl(P0, 'volt', { k }).p2, fan2)); return n.every(o => o.ok) && n[0].n > n[1].n && n[1].n > n[2].n; })());
    t.check('rotor resistance, constant torque: slip rises in proportion to total R2\' (s\' = s·(R2+Rext)/R2, same torque) and speed falls', (() => { const L = { type: 'const', T: 60, n: 1500 }, a = M.imOperatingPoint(P0, L), b = M.imOperatingPoint(cr.p2, L); return a.ok && b.ok && t.rel(b.s, a.s * 2.5, 1e-9) && b.n < a.n; })());
    t.check('load above breakdown torque: no stable point (stalled), no NaN speed; zero load → s = 0 at ns', (() => { const a = M.imOperatingPoint(P0, { type: 'const', T: 400, n: 1500 }), b = M.imOperatingPoint(P0, { type: 'const', T: 0, n: 1500 }); return !a.ok && a.stalled && !('n' in a) && b.ok && b.s === 0 && b.n === 1500; })());
    t.check('operating point rejects negative load torque and fan load with n_ref ≤ 0; start check T(1) vs load', !M.imOperatingPoint(P0, { type: 'const', T: -1, n: 1500 }).ok && !M.imOperatingPoint(P0, { type: 'fan', T: 10, n: 0 }).ok && M.imOperatingPoint(P0, { type: 'const', T: 100, n: 1500 }).canStart === false && M.imOperatingPoint(P0, { type: 'const', T: 50, n: 1500 }).canStart === true);

    // --- compound DC machine (motor).  V = 240 V, Ra = 0.5, Rs = 0.1, Rf = 120 Ω, Laf = 1 (per shunt-field ampere), Nf = 1000, Ns = 2 (r = 0.002), no brush drop, linear flux
    //   no load (Ia = 0), long shunt: If = 2 A, Ea = 240, K = 2.0 → ω0 = 120 rad/s (short shunt: series coil carries If, ω0 = Rf/(Laf(1+r)) = 119.76).
    //   long shunt, Ia = 50 A: Is = Ia; Ea = 240 − 50(0.5+0.1) = 210 V.  Cumulative: Ife = 2 + 0.002·50 = 2.1 → ω = 100 rad/s;  differential: Ife = 1.9 → ω = 110.526 rad/s
    //   plain shunt (no series winding): Ea = 240 − 25 = 215 V, K = 2 → 107.5 rad/s.   Speed falls: 20 (16.7 %), 12.5 (10.4 %), 9.47 (7.9 %) rad/s.
    const C0 = { V: 240, Ra: 0.5, Rf: 120, Rs: 0.1, Nf: 1000, Ns: 2, Laf: 1, Vb: 0, Isat: 0, conn: 'long', comp: 'cum' };
    const cum = M.dcCompoundMotor(C0, 50), dif = M.dcCompoundMotor(Object.assign({}, C0, { comp: 'diff' }), 50), shu = M.dcCompoundMotor(Object.assign({}, C0, { Ns: 0, Rs: 0 }), 50);
    t.check('compound motor, long shunt, Ia = 50 A: cumulative ω = 210/2.1 = 100 rad/s, Tem = 2.1·50 = 105 N·m; differential ω = 210/1.9 = 110.526', t.near(cum.w, 100, 1e-9) && t.near(cum.Tem, 105, 1e-9) && t.near(cum.Ea, 210, 1e-12) && t.near(dif.w, 210 / 1.9, 1e-9) && t.near(dif.Tem, 95, 1e-9));
    t.check('plain shunt at 50 A: ω = 215/2 = 107.5; no load: 120 rad/s for shunt, cumulative and differential', t.near(shu.w, 107.5, 1e-9) && [C0, Object.assign({}, C0, { comp: 'diff' }), Object.assign({}, C0, { Ns: 0, Rs: 0 })].every(c => t.near(M.dcCompoundMotor(c, 0).w, 120, 1e-9)));
    t.check('speed drop no-load → 50 A: cumulative (20) > shunt (12.5) > differential (9.47 rad/s)', t.near(120 - cum.w, 20, 1e-9) && t.near(120 - shu.w, 12.5, 1e-9) && t.near(120 - dif.w, 120 - 210 / 1.9, 1e-9) && (120 - cum.w) > (120 - shu.w) && (120 - shu.w) > (120 - dif.w));
    // short shunt, Ia = 50: If = (V − Ia·Rs)/(Rf + Rs) = 235/120.1 = 1.956703 A;  Va = If·Rf = 234.8043 V;  IL = 51.95670 A (= Is);  Ea = Va − Ia·Ra = 209.8043 V
    //   Ife = 1.956703 + 0.002·51.9567 = 2.060616 → ω = 209.8043/2.060616 = 101.8163 rad/s
    const shs = M.dcCompoundMotor(Object.assign({}, C0, { conn: 'short' }), 50), Ifs = 235 / 120.1, ILs = 50 + Ifs, Eas = Ifs * 120 - 25;
    t.check('short-shunt cumulative at 50 A: If = 235/120.1, Va = 234.804 V, Ea = 209.804 V, ω = Ea/(If + 0.002·IL) = 101.816 rad/s', t.rel(shs.If, Ifs, 1e-12) && t.rel(shs.Va, Ifs * 120, 1e-12) && t.rel(shs.Ea, Eas, 1e-12) && t.rel(shs.w, Eas / (Ifs + 0.002 * ILs), 1e-12) && t.near(shs.w, 101.8163, 1e-3) && t.rel(shs.Is, ILs, 1e-12));
    t.check('short shunt at no load: series coil carries If → ω0 = Rf/(Laf(1+r)) = 120/1.002 = 119.760; under load it differs from long shunt', t.near(M.dcCompoundMotor(Object.assign({}, C0, { conn: 'short' }), 0).w, 120 / 1.002, 1e-9) && Math.abs(shs.w - cum.w) > 0.5);
    t.check('compound with Ns = 0, Rs = 0 equals the module\'s shunt model dcSolve (w, Tem) incl. brush drop', (() => { const c = Object.assign({}, C0, { Ns: 0, Rs: 0, Vb: 2 }), a = M.dcCompoundMotor(c, 37), b = M.dcSolve({ type: 'shunt', V: 240, Vf: 240, Ra: 0.5, Rf: 120, Rs: 0.1, Laf: 1, Vb: 2, Prot: 0 }, 37); return t.rel(a.w, b.w, 1e-12) && t.rel(a.Tem, b.Tem, 1e-12); })());
    t.check('series-only (Rf = ∞, Laf·Ns/Nf per amp) equals dcSolve series with Laf\' = Laf·Ns/Nf = 0.002', (() => { const c = Object.assign({}, C0, { Rf: Infinity, Vb: 2 }), a = M.dcCompoundMotor(c, 20), b = M.dcSolve({ type: 'series', V: 240, Vf: 0, Ra: 0.5, Rf: 120, Rs: 0.1, Laf: 0.002, Vb: 2, Prot: 0 }, 20); return t.rel(a.w, b.w, 1e-12) && t.rel(a.Tem, b.Tem, 1e-12); })());
    t.check('power: Pin = V·IL, Pem = Ea·Ia = Tem·ω for all connections', ['long', 'short'].every(cn => ['cum', 'diff'].every(cm => { const r = M.dcCompoundMotor(Object.assign({}, C0, { conn: cn, comp: cm }), 30); return t.rel(r.Pin, 240 * r.IL, 1e-12) && t.rel(r.Pem, r.Ea * r.Ia, 1e-12) && t.rel(r.Pem, r.Tem * r.w, 1e-12) && t.rel(r.IL, r.Ia + r.If, 1e-12); })));
    t.check('differential compounding with a strong series winding: flux reverses → flagged, speed not a finite number', (() => { const r = M.dcCompoundMotor(Object.assign({}, C0, { comp: 'diff', Ns: 100 }), 50); return r.fluxRev && !r.runaway === (r.Ea <= 0) && !fin(r.w); })());

    // --- compound generator (Fröhlich saturation K = Laf·Ife/(1+Ife/Isat)), speed 150 rad/s, Laf = 1.2, Isat = 4 A, Rf = 120 Ω
    //   no load, short shunt (Is = IL = 0), Ra = Rs = Vb = 0:  ω·Laf/(1+If/Isat) = Rf  → If = Isat(ωLaf/Rf − 1) = 4(180/120 − 1) = 2 A, Vt = 240 V
    //   with Ra = 0.5: ωLaf/(1+If/Isat) = Rf + Ra → If = 4(180/120.5 − 1) = 1.975104 A, Vt = Rf·If = 237.0124 V
    const G0 = { V: 0, Ra: 0, Rf: 120, Rs: 0, Nf: 1000, Ns: 20, Laf: 1.2, Vb: 0, Isat: 4, w: 150, conn: 'short', comp: 'cum' };
    const g1 = M.dcCompoundGen(G0, 0), g2 = M.dcCompoundGen(Object.assign({}, G0, { Ra: 0.5 }), 0);
    t.check('generator no-load build-up: If = 2 A, Vt = 240 V (Ra = 0); with Ra = 0.5: Vt = 237.012 V', g1.ok && t.rel(g1.Vt, 240, 1e-9) && t.rel(g1.If, 2, 1e-9) && g2.ok && t.rel(g2.Vt, 120 * 4 * (180 / 120.5 - 1), 1e-9) && t.near(g2.Vt, 237.0124, 1e-3));
    // series-only generator (Rf = ∞): Ife = (Ns/Nf)·IL = 0.05·40 = 2 A → K = 1.2·2/(1+2/4) = 1.6, E = 150·1.6 = 240 V; Vt = E − IL(Ra + Rs) = 240 − 40·0.5 = 220 V
    const gs = M.dcCompoundGen(Object.assign({}, G0, { Rf: Infinity, Ns: 50, Ra: 0.3, Rs: 0.2, conn: 'long' }), 40);
    t.check('series generator: Vt = 150·1.6 − 40·(0.3+0.2) = 220 V', gs.ok && t.rel(gs.Vt, 220, 1e-9) && t.rel(gs.Ea, 240, 1e-9));
    // KVL/KCL check of a loaded compound solution, recomputed from first principles (long shunt, cumulative, IL = 40 A)
    const GL = Object.assign({}, G0, { Ra: 0.4, Rs: 0.1, Vb: 2, Ns: 20, conn: 'long', comp: 'cum' }), gl = M.dcCompoundGen(GL, 40);
    { const If = gl.Vt / 120, Ia = 40 + If, Ife = If + 0.02 * Ia, E = 150 * 1.2 * Ife / (1 + Ife / 4);
      t.check('loaded long-shunt cumulative generator: E = Vt + Ia(Ra+Rs) + Vb and E = ω·K(Ife) both satisfied (<1e-9 V)', gl.ok && !gl.collapsed && Math.abs(E - (gl.Vt + Ia * 0.5 + 2)) < 1e-9 && t.rel(gl.Ia, Ia, 1e-9), 'Vt=' + gl.Vt.toFixed(3)); }
    { const gS = Object.assign({}, GL, { conn: 'short' }), gk = M.dcCompoundGen(gS, 40), Va = gk.Vt + 40 * 0.1, If = Va / 120, Ia = 40 + If, Ife = If + 0.02 * 40, E = 150 * 1.2 * Ife / (1 + Ife / 4);
      t.check('loaded short-shunt generator: Va = Vt + IL·Rs, Ia = IL + If, E = Va + Ia·Ra + Vb = ω·K(If + r·IL)', gk.ok && Math.abs(E - (Va + Ia * 0.4 + 2)) < 1e-9); }
    const vOf = (cm, ns, rs) => M.dcCompoundGen(Object.assign({}, GL, { comp: cm, Ns: ns, Rs: rs }), 40).Vt;
    t.check('generator regulation at IL = 40 A: cumulative > plain shunt > differential terminal voltage', vOf('cum', 20, 0.1) > vOf('cum', 0, 0) && vOf('cum', 0, 0) > vOf('diff', 20, 0.1), [vOf('cum', 20, 0.1), vOf('cum', 0, 0), vOf('diff', 20, 0.1)].map(v => v.toFixed(1)).join(' / '));
    t.check('short-shunt compound generator at no load (Is = IL = 0) equals the plain shunt generator', t.rel(M.dcCompoundGen(Object.assign({}, GL, { conn: 'short' }), 0).Vt, M.dcCompoundGen(Object.assign({}, GL, { Ns: 0, conn: 'short' }), 0).Vt, 1e-9));
    t.check('shunt generator under extreme load: voltage collapses (flagged), no NaN; invalid inputs (Isat = 0, ω ≤ 0) rejected', (() => { const c = M.dcCompoundGen(Object.assign({}, G0, { Ns: 0, Ra: 0.5 }), 2000); return c.ok && c.collapsed && c.Vt === 0 && !M.dcCompoundGen(Object.assign({}, G0, { Isat: 0 }), 10).ok && !M.dcCompoundGen(Object.assign({}, G0, { w: 0 }), 10).ok; })());
    t.check('generator terminal voltage is monotone non-increasing with load for plain shunt (0, 20, 40, 60 A)', (() => { const v = [0, 20, 40, 60].map(i => M.dcCompoundGen(Object.assign({}, G0, { Ns: 0, Ra: 0.3, Rs: 0 }), i).Vt); return v[0] > v[1] && v[1] > v[2] && v[2] > v[3]; })());
  });


})();
