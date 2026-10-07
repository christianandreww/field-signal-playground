/* FSP module: Rotating Machines (Phase 8c induction motor, 8d DC machine, 8e three-phase & rotating field).
   Pure math lives in FSP.math.rotating (no DOM). All DOM work happens inside init/activate. */
(function () {
  'use strict';
  const PI = Math.PI, S3 = Math.sqrt(3);

  /* ======================= complex helpers (tiny, self-contained) ======================= */
  const cx = (re, im) => ({ re: re, im: im || 0 });
  const cadd = (a, b) => cx(a.re + b.re, a.im + b.im);
  const csub = (a, b) => cx(a.re - b.re, a.im - b.im);
  const cmul = (a, b) => cx(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
  const cdiv = (a, b) => { const d = b.re * b.re + b.im * b.im; return cx((a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d); };
  const cabs = a => Math.hypot(a.re, a.im);
  const cpol = (m, deg) => cx(m * Math.cos(deg * PI / 180), m * Math.sin(deg * PI / 180));
  const fin = x => Number.isFinite(x);

  /* ======================= 8c: induction motor ======================= */
  // p: {VL, conn:'Y'|'D', f, poles, R1, X1, R2, X2, Xm, Prot}. Per-phase, rotor referred to stator, core loss in Prot.
  const IM = {
    syncRpm: (f, poles) => 120 * f / poles,
    wsync: (f, poles) => 4 * PI * f / poles,
    vphase: p => (p.conn === 'D' ? p.VL : p.VL / S3),
    thevenin(p) {
      const Z1 = cx(p.R1, p.X1), Zm = cx(0, p.Xm), den = cadd(Z1, Zm);
      const Zth = cdiv(cmul(Zm, Z1), den), Vth = cdiv(cmul(cx(IM.vphase(p), 0), Zm), den);
      const Rth = Zth.re, Xth = Zth.im, K = Math.hypot(Rth, Xth + p.X2);
      return { Zth, Vth, Rth, Xth, VthMag: cabs(Vth), K };
    },
    // torque from the Thevenin form (N·m), sign convention: positive = motoring (in direction of rotating field)
    torqueTh(p, s, th) {
      if (s === 0) return 0; th = th || IM.thevenin(p);
      const r = p.R2 / s, ws = IM.wsync(p.f, p.poles);
      return 3 * th.VthMag * th.VthMag * r / (ws * ((th.Rth + r) * (th.Rth + r) + (th.Xth + p.X2) * (th.Xth + p.X2)));
    },
    // full per-phase circuit solved directly (no Thevenin reduction)
    solveFull(p, s) {
      const V = cx(IM.vphase(p), 0), Z1 = cx(p.R1, p.X1), Zm = cx(0, p.Xm);
      const Z2 = cx(p.R2 / s, p.X2), Zp = cdiv(cmul(Zm, Z2), cadd(Zm, Z2));
      const I1 = cdiv(V, cadd(Z1, Zp)), Vm = csub(V, cmul(I1, Z1)), I2 = cdiv(Vm, Z2);
      const Pag = 3 * (I2.re * I2.re + I2.im * I2.im) * p.R2 / s;
      return { V, I1, I2, Vm, Pag };
    },
    torqueFull(p, s) { if (s === 0) return 0; return IM.solveFull(p, s).Pag / IM.wsync(p.f, p.poles); },
    mode(s) { return s < 0 ? 'generator' : s === 0 ? 'synchronous (no torque)' : s <= 1 ? 'motor' : 'braking (plugging)'; },
    breakdown(p, th) {
      th = th || IM.thevenin(p); const ws = IM.wsync(p.f, p.poles), v2 = th.VthMag * th.VthMag;
      return {
        sm: p.R2 / th.K, Tmax: 3 * v2 / (2 * ws * (th.Rth + th.K)),
        sg: -p.R2 / th.K, Tgen: -3 * v2 / (2 * ws * (-th.Rth + th.K)),
      };
    },
    // operating point at slip s
    op(p, s) {
      const ns = IM.syncRpm(p.f, p.poles), ws = IM.wsync(p.f, p.poles);
      const wm = ws * (1 - s), full = IM.solveFull(p, s), T = full.Pag / ws;
      const Pag = full.Pag, Pmech = T * wm, Pcu2 = s * Pag, Pcu1 = 3 * (full.I1.re ** 2 + full.I1.im ** 2) * p.R1;
      const Pin = 3 * (full.V.re * full.I1.re + full.V.im * full.I1.im);
      const Pout = Pmech - p.Prot, mode = IM.mode(s);
      let eff = NaN;
      if (s > 0 && s < 1 && Pin > 0 && Pout > 0) eff = Pout / Pin;
      else if (s < 0 && Pin < 0 && (-Pmech + p.Prot) > 0) eff = (-Pin) / (-Pmech + p.Prot);
      const I1 = cabs(full.I1), pf = I1 > 0 ? Math.cos(Math.atan2(full.I1.im, full.I1.re)) : NaN;
      return { s, ns, ws, nrpm: ns * (1 - s), wm, T, Pag, Pmech, Pcu1, Pcu2, Pin, Pout, eff, I1, I2: cabs(full.I2), pf, mode };
    },
    sweep(p, n) { // slip -1..2
      const th = IM.thevenin(p), s = [], T = [];
      for (let i = 0; i <= n; i++) { const si = -1 + 3 * i / n; s.push(si); T.push(IM.torqueTh(p, Math.abs(si) < 1e-12 ? 0 : si, th)); }
      return { s, T };
    },
  };

  /* ======================= 8d: DC machine (motoring, armature reaction ignored, linear flux) ======================= */
  // p: {type:'sep'|'shunt'|'series', V, Vf, Ra, Rf, Rs, Laf, Vb, Prot}; flux linkage constant K = Laf*If, E = K*w, Tem = K*Ia
  const DC = {
    solve(p, Ia) {
      const series = p.type === 'series';
      const If = p.type === 'sep' ? p.Vf / p.Rf : p.type === 'shunt' ? p.V / p.Rf : Ia;
      const Raeff = p.Ra + (series ? p.Rs : 0);
      const Vdrop = Ia > 0 ? p.Vb : 0;
      const Ea = p.V - Vdrop - Ia * Raeff, K = p.Laf * If;
      const runaway = !(K > 1e-12) && Ea > 0;                 // no flux at all and still driven: speed -> infinity
      const w = runaway ? NaN : Ea / K;
      const Tem = K * Ia, Pem = Ea * Ia, Pshaft = Pem - p.Prot;
      const Pin = p.type === 'sep' ? p.V * Ia + p.Vf * If : p.type === 'shunt' ? p.V * (Ia + If) : p.V * Ia;
      const Pcua = Ia * Ia * Raeff, Pfield = series ? 0 : If * If * p.Rf, Pbrush = Vdrop * Ia;
      const Pout = Pshaft, losses = Pcua + Pfield + Pbrush + p.Prot;
      const eff = Pin > 0 && Pout > 0 ? Pout / Pin : NaN;
      const Tshaft = fin(w) && w > 0 ? Pshaft / w : NaN;
      return { Ia, If, Ea, K, w, rpm: fin(w) ? w * 30 / PI : NaN, Tem, Tshaft, Pem, Pshaft, Pin, Pout, Pcua, Pfield, Pbrush, Prot: p.Prot, losses, eff, runaway, stalled: Ea < 0, balance: Pin - Pout - losses };
    },
    stallCurrent: p => (p.V - p.Vb) / (p.Ra + (p.type === 'series' ? p.Rs : 0)),
  };

  /* ======================= 8e: three-phase & rotating field ======================= */
  const TP = {
    // instantaneous phase quantities; seq = +1 (a-b-c) or -1 (b and c swapped)
    phases(t, Am, f, seq) {
      const w = 2 * PI * f * t, a = 2 * PI / 3;
      return [Am * Math.cos(w), Am * Math.cos(w - seq * a), Am * Math.cos(w + seq * a)];
    },
    // resultant MMF vector (per-unit winding turns N): F = N * sum i_k e^{j k 120deg}
    mmfVector(t, Im, f, seq, N) {
      const i = TP.phases(t, Im, f, seq), a = 2 * PI / 3; N = N === undefined ? 1 : N;
      let x = 0, y = 0;
      for (let k = 0; k < 3; k++) { x += N * i[k] * Math.cos(k * a); y += N * i[k] * Math.sin(k * a); }
      return { x, y, mag: Math.hypot(x, y), ang: Math.atan2(y, x) };
    },
    // MMF at spatial angle th (rad) from the three pulsating windings
    mmfAt(th, t, Im, f, seq, N) {
      const i = TP.phases(t, Im, f, seq), a = 2 * PI / 3; N = N === undefined ? 1 : N;
      let s = 0; for (let k = 0; k < 3; k++) s += N * i[k] * Math.cos(th - k * a); return s;
    },
    // average angular speed (rad/s, + = ccw a->b->c) of the resultant, from unwrapped angle over dt
    rotationRate(t, dt, Im, f, seq) {
      const a0 = TP.mmfVector(t, Im, f, seq, 1).ang, a1 = TP.mmfVector(t + dt, Im, f, seq, 1).ang;
      let d = a1 - a0; while (d > PI) d -= 2 * PI; while (d < -PI) d += 2 * PI; return d / dt;
    },
    // balanced load: conn 'Y' or 'D', VL (rms line), load per phase Z = R + jX
    quantities(conn, VL, R, X) {
      const Zm = Math.hypot(R, X), Vph = conn === 'D' ? VL : VL / S3;
      const Iph = Zm > 1e-12 ? Vph / Zm : NaN, IL = conn === 'D' ? S3 * Iph : Iph, phi = Math.atan2(X, R);
      const S = 3 * Vph * Iph, P = S * Math.cos(phi), Q = S * Math.sin(phi);
      return { Vph, Iph, IL, VL, phi, S, P, Q, pf: Math.cos(phi), Zm };
    },
    // phasor sets (rms): Y -> phase voltages Van.. and line voltages; D -> phase currents Iab.. and line currents
    phasorSets(conn, VL, R, X, seq) {
      const q = TP.quantities(conn, VL, R, X), a = 120 * seq, ph = [], mag = conn === 'D' ? q.Iph : q.Vph;
      const base = conn === 'D' ? -q.phi * 180 / PI : 0;           // Vab = VL at 0 deg, phase current lags by phi
      [0, -a, a].forEach(o => ph.push(cpol(mag, base + o)));
      const line = ph.map((z, k) => (conn === 'D' ? csub(z, ph[(k + 2) % 3]) : csub(z, ph[(k + 1) % 3])));
      return { q, phase: ph, line, names: conn === 'D' ? { ph: ['Iab', 'Ibc', 'Ica'], ln: ['Ia', 'Ib', 'Ic'] } : { ph: ['Van', 'Vbn', 'Vcn'], ln: ['Vab', 'Vbc', 'Vca'] } };
    },
  };

  FSP.math.rotating = {
    syncRpm: IM.syncRpm, wsync: IM.wsync, imVphase: IM.vphase, imThevenin: IM.thevenin, imTorqueTh: IM.torqueTh,
    imSolveFull: IM.solveFull, imTorqueFull: IM.torqueFull, imBreakdown: IM.breakdown, imOp: IM.op, imMode: IM.mode, imSweep: IM.sweep,
    dcSolve: DC.solve, dcStallCurrent: DC.stallCurrent,
    phases: TP.phases, mmfVector: TP.mmfVector, mmfAt: TP.mmfAt, rotationRate: TP.rotationRate,
    threePhaseQuantities: TP.quantities, phasorSets: TP.phasorSets,
  };

  /* ======================= tests ======================= */
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
  });

  /* ======================= UI ======================= */
  const css = n => { try { return getComputedStyle(document.documentElement).getPropertyValue(n).trim() || null; } catch (e) { return null; } };
  function theme() {
    return { bg: css('--panel2') || '#0f141f', text: css('--text') || '#e6edf3', muted: css('--muted') || '#8b98ab', border: css('--border') || '#243047',
      a: css('--accent2') || '#f9a03f', b: css('--ok') || '#2ecc71', c: css('--accent') || '#4cc9f0', pink: css('--pink') || '#f72585', warn: css('--warn') || '#ffb347', mono: css('--mono') || 'monospace' };
  }
  const fmt = (x, d) => (fin(x) ? x.toFixed(d === undefined ? 2 : d) : '—');
  const tickStr = v => (Math.abs(v) < 1e-12 ? '0' : String(+v.toPrecision(4)));
  function niceTicks(lo, hi, n) {
    const span = hi - lo; if (!(span > 0)) return [lo];
    const raw = span / n, mag = Math.pow(10, Math.floor(Math.log10(raw))), f = raw / mag, step = (f < 1.5 ? 1 : f < 3 ? 2 : f < 7 ? 5 : 10) * mag, out = [];
    for (let v = Math.ceil(lo / step - 1e-9) * step; v <= hi + step * 1e-9; v += step) out.push(Math.abs(v) < step * 1e-9 ? 0 : +v.toPrecision(12));
    return out;
  }
  // generic line plot. o: {xmin,xmax,ymin,ymax,y2min,y2max,xlabel,ylabel,y2label,series:[{x,y,color,width,dash,axis}],points:[{x,y,label,color,axis}],vlines:[{x,label}],topNotes:[{x,text}]}
  function drawPlot(c, o) {
    const g = c.prep(), ctx = g.ctx, w = g.w, h = g.h, T = theme(), y2 = o.y2min !== undefined;
    const L = 54, R = y2 ? 40 : 10, Tp = o.topPad || 10, B = 32, pw = w - L - R, ph = h - Tp - B;
    ctx.clearRect(0, 0, w, h); ctx.fillStyle = T.bg; ctx.fillRect(0, 0, w, h);
    if (pw < 20 || ph < 20) return;
    const X = x => L + (x - o.xmin) / (o.xmax - o.xmin) * pw, Y = y => Tp + (1 - (y - o.ymin) / (o.ymax - o.ymin)) * ph;
    const Y2 = y => Tp + (1 - (y - o.y2min) / (o.y2max - o.y2min)) * ph;
    ctx.font = '11px ' + T.mono; ctx.lineWidth = 1;
    ctx.textBaseline = 'middle'; ctx.textAlign = 'right';
    niceTicks(o.ymin, o.ymax, 5).forEach(v => { const yy = Y(v); ctx.strokeStyle = T.border; ctx.beginPath(); ctx.moveTo(L, yy); ctx.lineTo(L + pw, yy); ctx.stroke(); ctx.fillStyle = T.muted; ctx.fillText(tickStr(v), L - 4, yy); });
    if (y2) { ctx.textAlign = 'left'; niceTicks(o.y2min, o.y2max, 5).forEach(v => { ctx.fillStyle = T.muted; ctx.fillText(tickStr(v), L + pw + 4, Y2(v)); }); }
    ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    niceTicks(o.xmin, o.xmax, Math.max(3, Math.floor(pw / 70))).forEach(v => { const xx = X(v); ctx.strokeStyle = T.border; ctx.beginPath(); ctx.moveTo(xx, Tp); ctx.lineTo(xx, Tp + ph); ctx.stroke(); ctx.fillStyle = T.muted; ctx.fillText(tickStr(v), xx, Tp + ph + 3); });
    // zero axes
    ctx.strokeStyle = T.muted; if (o.ymin < 0 && o.ymax > 0) { ctx.beginPath(); ctx.moveTo(L, Y(0)); ctx.lineTo(L + pw, Y(0)); ctx.stroke(); }
    if (o.xmin < 0 && o.xmax > 0) { ctx.beginPath(); ctx.moveTo(X(0), Tp); ctx.lineTo(X(0), Tp + ph); ctx.stroke(); }
    ctx.fillStyle = T.muted; ctx.textAlign = 'center'; ctx.textBaseline = 'bottom'; ctx.fillText(o.xlabel || '', L + pw / 2, h - 1);
    ctx.save(); ctx.translate(10, Tp + ph / 2); ctx.rotate(-PI / 2); ctx.textBaseline = 'top'; ctx.fillText(o.ylabel || '', 0, 0); ctx.restore();
    if (y2) { ctx.save(); ctx.translate(w - 3, Tp + ph / 2); ctx.rotate(-PI / 2); ctx.textBaseline = 'bottom'; ctx.fillText(o.y2label || '', 0, 0); ctx.restore(); }
    ctx.save(); ctx.beginPath(); ctx.rect(L, Tp, pw, ph); ctx.clip();
    (o.vlines || []).forEach(v => { ctx.strokeStyle = T.muted; ctx.setLineDash([4, 4]); ctx.beginPath(); ctx.moveTo(X(v.x), Tp); ctx.lineTo(X(v.x), Tp + ph); ctx.stroke(); ctx.setLineDash([]); });
    (o.series || []).forEach(s => {
      const fy = s.axis === 'y2' ? Y2 : Y; ctx.strokeStyle = s.color; ctx.lineWidth = s.width || 2; ctx.setLineDash(s.dash || []); ctx.beginPath(); let pen = false;
      for (let i = 0; i < s.x.length; i++) { const xv = s.x[i], yv = s.y[i]; if (!fin(xv) || !fin(yv)) { pen = false; continue; } const px = X(xv), py = fy(yv); if (pen) ctx.lineTo(px, py); else { ctx.moveTo(px, py); pen = true; } }
      ctx.stroke(); ctx.setLineDash([]);
    });
    ctx.restore();
    ctx.font = '11px ' + T.mono;
    (o.points || []).forEach(p => {
      if (!fin(p.x) || !fin(p.y)) return; const px = X(p.x), py = (p.axis === 'y2' ? Y2 : Y)(p.y);
      if (px < L - 1 || px > L + pw + 1) return;
      ctx.fillStyle = p.color; ctx.strokeStyle = T.bg; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.arc(px, py, 5, 0, 2 * PI); ctx.fill(); ctx.stroke();
      if (p.label) { const right = px < L + pw * 0.62; ctx.textAlign = right ? 'left' : 'right'; ctx.textBaseline = p.below ? 'top' : 'bottom'; ctx.fillStyle = T.text; ctx.fillText(p.label, px + (right ? 8 : -8), py + (p.below ? 6 : -6)); }
    });
    (o.topNotes || []).forEach(n => { ctx.fillStyle = n.color || T.muted; ctx.textAlign = n.align || 'center'; ctx.textBaseline = 'top'; ctx.font = (n.small ? '10px ' : '11px ') + T.mono; ctx.fillText(n.text, X(n.x), Tp + (n.dy || 3)); });
    (o.legend || []).forEach((lg, i) => { ctx.font = '11px ' + T.mono; ctx.fillStyle = lg.color; ctx.textAlign = 'right'; ctx.textBaseline = 'bottom'; ctx.fillText(lg.text, L + pw - 4, Tp + ph - 4 - i * 14); });
  }

  // slider factory: defs [{k,l,min,max,v,u,log,step}]
  function makeSliders(parent, defs, store, onChange) {
    defs.forEach(d => { store[d.k] = FSP.ui.slider(parent, { label: d.l, min: d.min, max: d.max, step: d.step, value: d.v, unit: d.u, log: d.log, onInput: () => { onChange(); FSP.state.touch(); } }); });
  }
  const val = (st, k) => st[k].get();
  function hudSet(el, rows) {
    while (el.firstChild) el.removeChild(el.firstChild);
    rows.forEach(r => el.appendChild(FSP.ui.el('div', null, FSP.ui.el('span', { text: r[0] }), FSP.ui.el('span', { text: r[1] }))));
  }
  const wrapCanvas = (parent, id, height) => { const wrap = FSP.ui.el('div', { class: 'canvas-wrap' }); parent.appendChild(wrap); const c = FSP.ui.canvas(wrap, { height }); c.cv.setAttribute('role', 'img'); if (id) c.cv.setAttribute('aria-label', id); return c; };

  /* ---------- 8c panel ---------- */
  function buildInduction(root) {
    const st = {}, ui = {}; const lay = FSP.ui.el('div', { class: 'layout' }); root.appendChild(lay);
    const ctl = FSP.ui.el('div', { class: 'controls' }), stage = FSP.ui.el('div', { class: 'stage' }); lay.appendChild(ctl); lay.appendChild(stage);
    const fs1 = FSP.ui.fieldset(ctl, 'Supply & construction');
    makeSliders(fs1, [
      { k: 'VL', l: 'V line', min: 50, max: 1000, v: 400, u: 'V' },
      { k: 'f', l: 'f', min: 1, max: 400, v: 50, u: 'Hz' }], st, update);
    ui.conn = FSP.ui.select(fs1, 'Connection', [['Y', 'Star (Y)'], ['D', 'Delta (Δ)']], 'Y', () => { update(); FSP.state.touch(); });
    ui.poles = FSP.ui.select(fs1, 'Poles', ['2', '4', '6', '8', '10', '12'], '4', () => { update(); FSP.state.touch(); });
    const fs2 = FSP.ui.fieldset(ctl, 'Per-phase circuit (referred to stator)');
    makeSliders(fs2, [
      { k: 'R1', l: 'R1', min: 0.01, max: 5, v: 0.7, u: 'Ω' }, { k: 'X1', l: 'X1', min: 0.01, max: 10, v: 1.2, u: 'Ω' },
      { k: 'R2', l: "R2'", min: 0.01, max: 5, v: 0.6, u: 'Ω' }, { k: 'X2', l: "X2'", min: 0.01, max: 10, v: 1.2, u: 'Ω' },
      { k: 'Xm', l: 'Xm', min: 2, max: 200, v: 35, u: 'Ω' }, { k: 'Prot', l: 'P rot.', min: 0, max: 2000, v: 150, u: 'W' }], st, update);
    const fs3 = FSP.ui.fieldset(ctl, 'Operating point');
    makeSliders(fs3, [{ k: 'sr', l: 'Slip s', min: -1, max: 2, v: 0.04, step: 0.001, u: '' }], st, update);
    fs3.appendChild(FSP.ui.el('div', { class: 'note', text: 'Core loss is lumped into rotational loss (no Rc branch). Slip < 0 generator, 0–1 motor, > 1 braking.' }));
    ui.msg = FSP.ui.el('div', { class: 'msg warn', hidden: '' }); ctl.appendChild(ui.msg);
    ui.c1 = wrapCanvas(stage, 'Induction motor torque versus slip', 300); ui.c1.onResize(update);
    ui.c2 = wrapCanvas(stage, 'Mechanical power and efficiency versus slip', 220); ui.c2.onResize(update);
    ui.hud = FSP.ui.el('div', { class: 'hud' }); stage.appendChild(ui.hud);
    ui.work = FSP.ui.working(stage);

    function params() { return { VL: val(st, 'VL'), conn: ui.conn.value, f: val(st, 'f'), poles: parseInt(ui.poles.value, 10) || 4, R1: val(st, 'R1'), X1: val(st, 'X1'), R2: val(st, 'R2'), X2: val(st, 'X2'), Xm: val(st, 'Xm'), Prot: val(st, 'Prot') }; }
    function update() {
      if (!ui.c1) return;
      const p = params(), th = IM.thevenin(p), sw = IM.sweep(p, 360), bd = IM.breakdown(p, th), sr = val(st, 'sr');
      const op = IM.op(p, Math.abs(sr) < 1e-9 ? 0 : sr), start = IM.op(p, 1), ns = IM.syncRpm(p.f, p.poles);
      const T = theme();
      const tmax = Math.max(bd.Tmax, Math.abs(bd.Tgen), 1e-9), ymin = -tmax * 1.2, ymax = tmax * 1.25;
      const nt = s => 'n=' + Math.round(ns * (1 - s));
      drawPlot(ui.c1, {
        xmin: -1, xmax: 2, ymin, ymax, xlabel: 'slip s  (speed n = ns(1−s))', ylabel: 'torque (N·m)', topPad: 34,
        series: [{ x: sw.s, y: sw.T, color: T.c, width: 2.2 }], vlines: [{ x: 0 }, { x: 1 }],
        topNotes: [{ x: -0.5, text: 'GENERATOR', color: T.b }, { x: 0.5, text: 'MOTOR', color: T.a }, { x: 1.5, text: 'BRAKING', color: T.pink },
          { x: -1, text: nt(-1), align: 'left', small: true, dy: 17 }, { x: 0, text: nt(0), small: true, dy: 17 }, { x: 1, text: nt(1), small: true, dy: 17 }, { x: 2, text: nt(2), align: 'right', small: true, dy: 17 }],
        points: [{ x: 1, y: start.T, label: 'start', color: T.pink, below: true }, { x: bd.sm, y: bd.Tmax, label: 'breakdown', color: T.warn }, { x: bd.sg, y: bd.Tgen, label: 'gen. breakdown', color: T.warn, below: true }, { x: op.s, y: op.T, label: 'op point', color: T.text }],
      });
      // second plot: motor region power & efficiency
      const sa = [], pm = [], pi = [], ef = []; for (let i = 1; i <= 200; i++) { const s = i / 200, o = IM.op(p, s); sa.push(s); pm.push(o.Pmech / 1000); pi.push(o.Pin / 1000); ef.push(fin(o.eff) ? o.eff * 100 : NaN); }
      const pmax = Math.max.apply(null, pi.concat(pm).filter(fin)) * 1.1;
      drawPlot(ui.c2, {
        xmin: 0, xmax: 1, ymin: 0, ymax: pmax > 0 ? pmax : 1, y2min: 0, y2max: 100, xlabel: 'slip s (motor region)', ylabel: 'power (kW)', y2label: 'efficiency (%)',
        series: [{ x: sa, y: pi, color: T.muted, width: 1.5, dash: [5, 3] }, { x: sa, y: pm, color: T.a, width: 2.2 }, { x: sa, y: ef, color: T.b, width: 2, axis: 'y2' }],
        points: op.s > 0 && op.s <= 1 ? [{ x: op.s, y: op.Pmech / 1000, color: T.text }] : [],
        legend: [{ text: 'Pmech', color: T.a }, { text: 'Pin', color: T.muted }, { text: 'efficiency', color: T.b }],
      });
      hudSet(ui.hud, [
        ['ns = 120f/P', fmt(ns, 1) + ' rpm  (ωs ' + fmt(op.ws, 2) + ' rad/s)'], ['phase voltage', fmt(IM.vphase(p), 1) + ' V'],
        ['Zth', fmt(th.Rth, 3) + ' + j' + fmt(th.Xth, 3) + ' Ω'], ['|Vth|', fmt(th.VthMag, 2) + ' V'],
        ['mode @ s=' + fmt(op.s, 3), op.mode], ['speed', fmt(op.nrpm, 1) + ' rpm'],
        ['torque', fmt(op.T, 2) + ' N·m'], ['air-gap power', fmt(op.Pag / 1000, 3) + ' kW'],
        ['Pmech (T·ωm)', fmt(op.Pmech / 1000, 3) + ' kW'], ['shaft out', fmt(op.Pout / 1000, 3) + ' kW'],
        ['Pin', fmt(op.Pin / 1000, 3) + ' kW'], ['efficiency', fin(op.eff) ? fmt(op.eff * 100, 1) + ' %' : '—'],
        ['stator current', fmt(op.I1, 2) + ' A'], ['power factor', fmt(op.pf, 3)],
        ['starting torque', fmt(start.T, 2) + ' N·m'], ['starting current', fmt(start.I1, 1) + ' A'],
        ['breakdown slip', fmt(bd.sm, 4)], ['breakdown torque', fmt(bd.Tmax, 2) + ' N·m'],
        ['Tmax / Tstart', fmt(start.T > 0 ? bd.Tmax / start.T : NaN, 2)], ['gen. breakdown', fmt(bd.Tgen, 2) + ' N·m @ s=' + fmt(bd.sg, 4)],
      ]);
      ui.work.set([
        'ns = 120 f / P = 120·' + p.f + ' / ' + p.poles + ' = ' + fmt(ns, 2) + ' rpm;  ωs = 2π ns/60 = ' + fmt(op.ws, 4) + ' rad/s',
        'Vph = ' + (p.conn === 'D' ? 'VL' : 'VL/√3') + ' = ' + fmt(IM.vphase(p), 3) + ' V (per-phase analysis)',
        'Thevenin seen from the rotor branch:',
        '  Zth = jXm(R1+jX1)/(R1+j(X1+Xm)) = ' + fmt(th.Rth, 5) + ' + j' + fmt(th.Xth, 5) + ' Ω',
        '  Vth = Vph·jXm/(R1+j(X1+Xm)),  |Vth| = ' + fmt(th.VthMag, 4) + ' V',
        'T(s) = 3|Vth|² (R2\'/s) / ( ωs [ (Rth + R2\'/s)² + (Xth + X2\')² ] )',
        '  at s = ' + fmt(op.s, 4) + ':  T = ' + fmt(op.T, 4) + ' N·m  (direct full-circuit solve gives the same)',
        'Pag = T·ωs = ' + fmt(op.Pag, 2) + ' W;  Pcu2 = s·Pag = ' + fmt(op.Pcu2, 2) + ' W;  Pmech = (1−s)Pag = ' + fmt(op.Pmech, 2) + ' W',
        'Pin = 3 Re(V I1*) = ' + fmt(op.Pin, 2) + ' W = Pcu1 (' + fmt(op.Pcu1, 2) + ') + Pag;   Pout = Pmech − Prot = ' + fmt(op.Pout, 2) + ' W',
        'Breakdown: sm = R2\'/√(Rth²+(Xth+X2\')²) = ' + fmt(bd.sm, 5) + ';  Tmax = 3|Vth|² / (2ωs (Rth + √(Rth²+(Xth+X2\')²))) = ' + fmt(bd.Tmax, 3) + ' N·m',
      ]);
      ui.msg.hidden = true;
      if (!(bd.Tmax > 0)) { ui.msg.hidden = false; ui.msg.textContent = 'Parameters give no valid torque curve.'; }
    }
    return {
      update,
      get: () => ({ VL: val(st, 'VL'), f: val(st, 'f'), conn: ui.conn.value, poles: ui.poles.value, R1: val(st, 'R1'), X1: val(st, 'X1'), R2: val(st, 'R2'), X2: val(st, 'X2'), Xm: val(st, 'Xm'), Prot: val(st, 'Prot'), sr: val(st, 'sr') }),
      set(o) {
        Object.keys(st).forEach(k => { if (o[k] !== undefined) { const v = parseFloat(o[k]); if (fin(v)) st[k].set(v, true); } });
        if (o.conn === 'Y' || o.conn === 'D') ui.conn.value = o.conn;
        if (['2', '4', '6', '8', '10', '12'].indexOf(o.poles) >= 0) ui.poles.value = o.poles;
      },
    };
  }

  /* ---------- 8d panel ---------- */
  function buildDC(root) {
    const st = {}, ui = {}; const lay = FSP.ui.el('div', { class: 'layout' }); root.appendChild(lay);
    const ctl = FSP.ui.el('div', { class: 'controls' }), stage = FSP.ui.el('div', { class: 'stage' }); lay.appendChild(ctl); lay.appendChild(stage);
    const LAF = { sep: 0.8, shunt: 0.8, series: 0.032 };
    const fs1 = FSP.ui.fieldset(ctl, 'Machine & excitation');
    ui.type = FSP.ui.select(fs1, 'Excitation', [['sep', 'Separately excited'], ['shunt', 'Shunt'], ['series', 'Series']], 'shunt', () => { st.Laf.set(LAF[ui.type.value], true); update(); FSP.state.touch(); });
    makeSliders(fs1, [
      { k: 'V', l: 'V term.', min: 12, max: 600, v: 240, u: 'V' }, { k: 'Vf', l: 'V field', min: 12, max: 600, v: 240, u: 'V' },
      { k: 'Laf', l: 'Laf', min: 0.005, max: 3, v: 0.8, u: 'H', log: true }], st, update);
    ui.vfRow = st.Vf.el; fs1.appendChild(FSP.ui.el('div', { class: 'note', text: 'Φ ∝ field current (no saturation): E = Laf·If·ω, Tem = Laf·If·Ia. For series, If = Ia.' }));
    const fs2 = FSP.ui.fieldset(ctl, 'Resistances & losses');
    makeSliders(fs2, [
      { k: 'Ra', l: 'Ra', min: 0.01, max: 5, v: 0.5, u: 'Ω' }, { k: 'Rf', l: 'Rf', min: 10, max: 500, v: 120, u: 'Ω' },
      { k: 'Rs', l: 'Rs', min: 0.01, max: 2, v: 0.1, u: 'Ω' }, { k: 'Vb', l: 'Vbrush', min: 0, max: 10, v: 2, u: 'V' },
      { k: 'Prot', l: 'P rot.', min: 0, max: 2000, v: 300, u: 'W' }], st, update);
    const fs3 = FSP.ui.fieldset(ctl, 'Load');
    makeSliders(fs3, [{ k: 'Irated', l: 'Ia rated', min: 5, max: 300, v: 50, u: 'A' }, { k: 'load', l: 'Load', min: 0, max: 150, v: 100, step: 1, u: '% Ia' }], st, update);
    ui.warn = FSP.ui.el('div', { class: 'msg warn', hidden: '' }); ctl.appendChild(ui.warn);
    ui.c1 = wrapCanvas(stage, 'DC machine speed versus torque', 280); ui.c1.onResize(update);
    ui.c2 = wrapCanvas(stage, 'DC machine efficiency and shaft power versus armature current', 220); ui.c2.onResize(update);
    stage.appendChild(FSP.ui.el('div', { class: 'note', text: 'Power flow at the operating point (width ∝ power):' }));
    ui.flow = FSP.ui.el('div', { class: 'row' }); ui.flow.style.cssText = 'display:flex;flex-wrap:nowrap;gap:0;border:1px solid var(--border);border-radius:6px;overflow:hidden;min-height:26px'; stage.appendChild(ui.flow);
    ui.flowLegend = FSP.ui.el('div', { class: 'note mono' }); stage.appendChild(ui.flowLegend);
    ui.hud = FSP.ui.el('div', { class: 'hud' }); stage.appendChild(ui.hud);
    ui.work = FSP.ui.working(stage);

    function params() { return { type: ui.type.value, V: val(st, 'V'), Vf: val(st, 'Vf'), Ra: val(st, 'Ra'), Rf: val(st, 'Rf'), Rs: val(st, 'Rs'), Laf: val(st, 'Laf'), Vb: val(st, 'Vb'), Prot: val(st, 'Prot') }; }
    function update() {
      if (!ui.c1) return;
      const p = params(), T = theme(), type = p.type, Ir = val(st, 'Irated'), Ia = Ir * val(st, 'load') / 100;
      st.Vf.el.hidden = type !== 'sep'; st.Rs.el.hidden = type !== 'series'; st.Rf.el.hidden = type === 'series';
      const r = DC.solve(p, Ia), ref = DC.solve(p, Ir), wref = fin(ref.w) && ref.w > 0 ? ref.w : NaN;
      const wcap = fin(wref) ? 4 * wref : Infinity;
      // curves vs armature current
      const n = 240, imax = Math.min(1.6 * Ir, DC.stallCurrent(p) * 0.999), ia = [], rpm = [], tem = [], eff = [], psh = [];
      for (let i = 0; i <= n; i++) {
        const I = imax * i / n, s = DC.solve(p, I);
        if (type === 'series' && !(fin(s.w) && s.w <= wcap)) continue;       // hide the runaway region
        if (!fin(s.w) || s.w < 0) continue;
        ia.push(I); rpm.push(s.rpm); tem.push(s.Tem); eff.push(fin(s.eff) ? s.eff * 100 : NaN); psh.push(s.Pshaft / 1000);
      }
      const rmax = Math.max.apply(null, rpm.concat([1])) * 1.1, tmax = Math.max.apply(null, tem.concat([1e-6])) * 1.05;
      drawPlot(ui.c1, {
        xmin: 0, xmax: tmax, ymin: 0, ymax: rmax, xlabel: 'electromagnetic torque (N·m)', ylabel: 'speed (rpm)',
        series: [{ x: tem, y: rpm, color: T.c, width: 2.2 }],
        points: [{ x: r.Tem, y: r.rpm, label: fin(r.rpm) ? 'op point' : '', color: T.text }, { x: ref.Tem, y: ref.rpm, label: 'rated', color: T.warn, below: true }],
      });
      const pmax = Math.max.apply(null, psh.concat([0.001])) * 1.1;
      drawPlot(ui.c2, {
        xmin: 0, xmax: Math.max(imax, 1e-6), ymin: Math.min(0, Math.min.apply(null, psh.filter(fin))), ymax: pmax, y2min: 0, y2max: 100, xlabel: 'armature current Ia (A)', ylabel: 'shaft power (kW)', y2label: 'efficiency (%)',
        series: [{ x: ia, y: psh, color: T.a, width: 2.2 }, { x: ia, y: eff, color: T.b, width: 2, axis: 'y2' }],
        points: [{ x: Ia, y: fin(r.eff) ? r.eff * 100 : NaN, color: T.text, axis: 'y2' }],
        legend: [{ text: 'shaft power', color: T.a }, { text: 'efficiency', color: T.b }],
      });
      // warnings
      const warns = [];
      if (r.runaway || (type === 'series' && fin(wref) && fin(r.w) && r.w > wcap)) warns.push('Series motor runaway: flux ∝ Ia, so at (near) zero load the speed → ∞. Speed is shown as — ; raise the load. A real series motor must never run unloaded.');
      else if (type === 'series' && fin(r.w) && fin(wref) && r.w > 2 * wref) warns.push('Light load on a series motor: speed is already more than 2× rated and rises rapidly as load falls.');
      if (r.stalled) warns.push('Armature current exceeds the stall current (V−Vb)/Ra = ' + fmt(DC.stallCurrent(p), 1) + ' A: back-EMF would be negative.');
      if (!(r.K > 1e-12) && !r.runaway) warns.push('No field flux (K = 0).');
      ui.warn.hidden = warns.length === 0; ui.warn.textContent = warns.join(' ');
      // power flow bar
      while (ui.flow.firstChild) ui.flow.removeChild(ui.flow.firstChild);
      const parts = [['Pout', Math.max(0, r.Pout), T.b], ['Cu(armature)', r.Pcua, T.a], ['field', r.Pfield, T.pink], ['brush', r.Pbrush, T.warn], ['rot.', r.Prot, T.muted]].filter(x => x[1] > 0 && fin(x[1]));
      const tot = parts.reduce((a, x) => a + x[1], 0);
      parts.forEach(x => { const d = FSP.ui.el('div', { title: x[0] + ' ' + fmt(x[1], 0) + ' W' }); d.style.cssText = 'height:26px;background:' + x[2] + ';flex:' + (x[1] / tot) + ' 1 0;min-width:2px'; ui.flow.appendChild(d); });
      ui.flowLegend.textContent = parts.map(x => x[0] + ' ' + fmt(x[1] / 1000, 3) + ' kW').join('  |  ') + (parts.length ? '' : '—');
      hudSet(ui.hud, [
        ['armature current', fmt(r.Ia, 2) + ' A'], ['field current', fmt(r.If, 3) + ' A'],
        ['back-EMF E', fmt(r.Ea, 2) + ' V'], ['speed', fin(r.rpm) ? fmt(r.rpm, 1) + ' rpm' : '— (runaway)'],
        ['Tem = K·Ia', fmt(r.Tem, 2) + ' N·m'], ['shaft torque', fmt(r.Tshaft, 2) + ' N·m'],
        ['Pin', fmt(r.Pin / 1000, 3) + ' kW'], ['Pem = E·Ia', fmt(r.Pem / 1000, 3) + ' kW'],
        ['Pout (shaft)', fmt(r.Pout / 1000, 3) + ' kW'], ['total losses', fmt(r.losses / 1000, 3) + ' kW'],
        ['efficiency', fin(r.eff) ? fmt(r.eff * 100, 1) + ' %' : '—'], ['balance Pin−Pout−losses', fmt(r.balance, 9) + ' W'],
      ]);
      ui.work.set([
        'Type: ' + type + '.  Flux linkage constant K = Laf·If = ' + fmt(p.Laf, 4) + ' × ' + fmt(r.If, 4) + ' = ' + fmt(r.K, 5) + ' V·s/rad (A·… per unit Ia gives N·m)',
        'If = ' + (type === 'sep' ? 'Vf/Rf' : type === 'shunt' ? 'V/Rf' : 'Ia') + ' = ' + fmt(r.If, 4) + ' A',
        'E = V − Vb − Ia·Ra' + (type === 'series' ? '(+Rs)' : '') + ' = ' + fmt(p.V, 2) + ' − ' + fmt(p.Vb, 2) + ' − ' + fmt(Ia, 2) + '·' + fmt(p.Ra + (type === 'series' ? p.Rs : 0), 3) + ' = ' + fmt(r.Ea, 3) + ' V',
        'ω = E / K = ' + (fin(r.w) ? fmt(r.w, 3) + ' rad/s = ' + fmt(r.rpm, 1) + ' rpm' : 'undefined (K→0, runaway)') + ';   Tem = K·Ia = ' + fmt(r.Tem, 3) + ' N·m',
        'Pem = E·Ia = ' + fmt(r.Pem, 2) + ' W = Tem·ω;   Pshaft = Pem − Prot = ' + fmt(r.Pshaft, 2) + ' W',
        'Pin = ' + (type === 'sep' ? 'V·Ia + Vf·If' : type === 'shunt' ? 'V·(Ia+If)' : 'V·Ia') + ' = ' + fmt(r.Pin, 2) + ' W',
        'Losses = Ia²Ra' + (type === 'series' ? '(+Rs)' : '') + ' ' + fmt(r.Pcua, 1) + ' + field ' + fmt(r.Pfield, 1) + ' + brush ' + fmt(r.Pbrush, 1) + ' + rot ' + fmt(r.Prot, 1) + ' = ' + fmt(r.losses, 2) + ' W',
        'Pin − Pout − losses = ' + fmt(r.balance, 12) + ' W',
        type === 'sep' ? 'Separately excited: ω = (V − Vb)/K − Ra·Tem/K²  (straight line, slope −Ra/K²)' : type === 'shunt' ? 'Shunt (constant V): same straight line as separately excited, nearly constant speed' : 'Series: ω = (V − Vb − Ia(Ra+Rs)) / (Laf·Ia) — hyperbolic, ∞ at Ia→0',
      ]);
    }
    return {
      update,
      get: () => ({ type: ui.type.value, V: val(st, 'V'), Vf: val(st, 'Vf'), Laf: val(st, 'Laf'), Ra: val(st, 'Ra'), Rf: val(st, 'Rf'), Rs: val(st, 'Rs'), Vb: val(st, 'Vb'), Prot: val(st, 'Prot'), Irated: val(st, 'Irated'), load: val(st, 'load') }),
      set(o) {
        if (['sep', 'shunt', 'series'].indexOf(o.type) >= 0) ui.type.value = o.type;
        Object.keys(st).forEach(k => { if (o[k] !== undefined) { const v = parseFloat(o[k]); if (fin(v)) st[k].set(v, true); } });
      },
    };
  }

  /* ---------- 8e panel ---------- */
  function buildThree(root, host) {
    const st = {}, ui = {}; const lay = FSP.ui.el('div', { class: 'layout' }); root.appendChild(lay);
    const ctl = FSP.ui.el('div', { class: 'controls' }), stage = FSP.ui.el('div', { class: 'stage' }); lay.appendChild(ctl); lay.appendChild(stage);
    let seq = 1, tsim = 0, paused = false, raf = 0, last = 0;
    const fs1 = FSP.ui.fieldset(ctl, 'Supply & load');
    makeSliders(fs1, [
      { k: 'VL', l: 'V line', min: 50, max: 1000, v: 400, u: 'V rms' }, { k: 'f', l: 'f', min: 1, max: 400, v: 50, u: 'Hz' },
      { k: 'R', l: 'R / phase', min: 0.1, max: 100, v: 8, u: 'Ω', log: true }, { k: 'X', l: 'X / phase', min: -100, max: 100, v: 6, u: 'Ω' }], st, redraw);
    ui.conn = FSP.ui.select(fs1, 'Connection', [['Y', 'Star (Y)'], ['D', 'Delta (Δ)']], 'Y', () => { redraw(); FSP.state.touch(); });
    fs1.appendChild(FSP.ui.el('div', { class: 'note', text: 'Y shows phase voltages and line voltages; Δ shows phase currents and line currents (balanced load, X > 0 inductive).' }));
    const fs2 = FSP.ui.fieldset(ctl, 'Rotating field');
    makeSliders(fs2, [{ k: 'slow', l: 'Slow-mo ×', min: 0.001, max: 0.5, v: 0.02, u: '', log: true }], st, () => { /* speed only */ });
    ui.row = FSP.ui.el('div', { class: 'row' }); fs2.appendChild(ui.row);
    ui.swap = FSP.ui.el('button', { type: 'button', class: 'btn', 'aria-pressed': 'false', text: 'Swap phases b ↔ c', onclick: () => { seq = -seq; redraw(); FSP.state.touch(); } }); ui.row.appendChild(ui.swap);
    ui.pause = FSP.ui.el('button', { type: 'button', class: 'btn', 'aria-pressed': 'false', text: 'Pause', onclick: () => { paused = !paused; ui.pause.textContent = paused ? 'Play' : 'Pause'; ui.pause.setAttribute('aria-pressed', String(paused)); refreshAnim(); } }); ui.row.appendChild(ui.pause);
    ui.dir = FSP.ui.el('div', { class: 'readout mono' }); fs2.appendChild(ui.dir);
    fs2.appendChild(FSP.ui.el('div', { class: 'note', text: 'Real rotation is ' + '50 Hz × 60 = 3000 rev/min for 2 poles; slowed so the eye can follow.' }));
    ui.cP = wrapCanvas(stage, 'Three-phase phasor diagram', 280); ui.cP.onResize(redraw);
    ui.cW = wrapCanvas(stage, 'Three-phase waveforms', 220); ui.cW.onResize(redraw);
    ui.cM = wrapCanvas(stage, 'Rotating MMF from three windings', 320); ui.cM.onResize(redraw);
    ui.hud = FSP.ui.el('div', { class: 'hud' }); stage.appendChild(ui.hud);
    ui.work = FSP.ui.working(stage);

    const Im = 1; // normalised winding current peak for the MMF picture
    function cfg() { return { conn: ui.conn.value, VL: val(st, 'VL'), f: val(st, 'f'), R: val(st, 'R'), X: val(st, 'X') }; }
    function draw() {
      if (!ui.cP.cv.offsetWidth || ui.cP.cv.offsetWidth < 60 || !ui.cM.cv.offsetWidth || !ui.cW.cv.offsetWidth) return null;
      const c = cfg(), T = theme(), cols = [T.a, T.b, T.c], ph = 2 * PI * c.f * tsim, sets = TP.phasorSets(c.conn, c.VL, c.R, c.X, seq), q = sets.q;
      // --- phasors (rotated by wt so projections on the horizontal axis are the instantaneous values)
      let g = ui.cP.prep(), ctx = g.ctx, w = g.w, h = g.h; ctx.clearRect(0, 0, w, h); ctx.fillStyle = T.bg; ctx.fillRect(0, 0, w, h);
      const cxp = w / 2, cyp = h / 2, rad = Math.min(w, h) / 2 - 22;
      const mx = Math.max.apply(null, sets.phase.concat(sets.line).map(cabs).concat([1e-9])), sc = rad / (mx * 1.05);
      ctx.strokeStyle = T.border; ctx.lineWidth = 1; ctx.beginPath(); ctx.arc(cxp, cyp, rad, 0, 2 * PI); ctx.stroke();
      ctx.beginPath(); ctx.moveTo(cxp - rad, cyp); ctx.lineTo(cxp + rad, cyp); ctx.moveTo(cxp, cyp - rad); ctx.lineTo(cxp, cyp + rad); ctx.stroke();
      const rot = (z) => cx(z.re * Math.cos(ph) - z.im * Math.sin(ph), z.re * Math.sin(ph) + z.im * Math.cos(ph));
      function arrow(z, color, label, dash, lw) {
        if (!fin(z.re) || !fin(z.im)) return; const r = rot(z), x1 = cxp + r.re * sc, y1 = cyp - r.im * sc;
        ctx.strokeStyle = color; ctx.fillStyle = color; ctx.lineWidth = lw; ctx.setLineDash(dash || []); ctx.beginPath(); ctx.moveTo(cxp, cyp); ctx.lineTo(x1, y1); ctx.stroke(); ctx.setLineDash([]);
        const a = Math.atan2(y1 - cyp, x1 - cxp); ctx.beginPath(); ctx.moveTo(x1, y1); ctx.lineTo(x1 - 9 * Math.cos(a - 0.35), y1 - 9 * Math.sin(a - 0.35)); ctx.lineTo(x1 - 9 * Math.cos(a + 0.35), y1 - 9 * Math.sin(a + 0.35)); ctx.closePath(); ctx.fill();
        ctx.font = '11px ' + T.mono; ctx.textAlign = 'center'; ctx.textBaseline = 'middle'; ctx.fillText(label, x1 + 14 * Math.cos(a), y1 + 14 * Math.sin(a));
      }
      sets.line.forEach((z, k) => arrow(z, cols[k], sets.names.ln[k], [5, 4], 1.4));
      sets.phase.forEach((z, k) => arrow(z, cols[k], sets.names.ph[k], null, 2.4));
      ctx.fillStyle = T.muted; ctx.font = '11px ' + T.mono; ctx.textAlign = 'left'; ctx.textBaseline = 'top'; ctx.fillText('solid: phase, dashed: line (rms)', 6, 5);
      // --- waveforms of the phase set
      g = ui.cW.prep(); ctx = g.ctx; w = g.w; h = g.h; ctx.clearRect(0, 0, w, h); ctx.fillStyle = T.bg; ctx.fillRect(0, 0, w, h);
      const L = 40, Rm = 10, Tp = 12, B = 22, pw = w - L - Rm, pht = h - Tp - B, per = 1 / c.f, X = tt => L + tt / (2 * per) * pw, Y = v => Tp + (1 - (v + 1.1) / 2.2) * pht;
      ctx.font = '11px ' + T.mono; ctx.strokeStyle = T.border; ctx.lineWidth = 1; ctx.fillStyle = T.muted; ctx.textBaseline = 'middle'; ctx.textAlign = 'right';
      [-1, 0, 1].forEach(v => { ctx.beginPath(); ctx.moveTo(L, Y(v)); ctx.lineTo(L + pw, Y(v)); ctx.stroke(); ctx.fillText(String(v), L - 4, Y(v)); });
      ctx.textAlign = 'center'; ctx.textBaseline = 'top'; [0, 1, 2].forEach(k => { ctx.beginPath(); ctx.moveTo(X(k * per), Tp); ctx.lineTo(X(k * per), Tp + pht); ctx.stroke(); ctx.fillText(k + 'T', X(k * per), Tp + pht + 3); });
      const npt = 240, sums = [];
      const wave = (k, color, lw, dash) => { ctx.strokeStyle = color; ctx.lineWidth = lw; ctx.setLineDash(dash || []); ctx.beginPath(); for (let i = 0; i <= npt; i++) { const tt = 2 * per * i / npt, v = k < 3 ? TP.phases(tt, 1, c.f, seq)[k] : (() => { const a = TP.phases(tt, 1, c.f, seq); return a[0] + a[1] + a[2]; })(); if (i) ctx.lineTo(X(tt), Y(v)); else ctx.moveTo(X(tt), Y(v)); } ctx.stroke(); ctx.setLineDash([]); };
      wave(0, cols[0], 2); wave(1, cols[1], 2); wave(2, cols[2], 2); wave(3, T.text, 1.5, [4, 3]);
      const tc = tsim % per; ctx.strokeStyle = T.pink; ctx.lineWidth = 1.2; ctx.beginPath(); ctx.moveTo(X(tc), Tp); ctx.lineTo(X(tc), Tp + pht); ctx.stroke();
      const inst = TP.phases(tc, 1, c.f, seq); inst.forEach((v, k) => { ctx.fillStyle = cols[k]; ctx.beginPath(); ctx.arc(X(tc), Y(v), 4, 0, 2 * PI); ctx.fill(); });
      ctx.fillStyle = T.muted; ctx.textAlign = 'left'; ctx.textBaseline = 'top'; ctx.fillText('p.u. of peak; dashed = a+b+c = 0', L + 4, Tp + 1);
      // --- rotating MMF
      g = ui.cM.prep(); ctx = g.ctx; w = g.w; h = g.h; ctx.clearRect(0, 0, w, h); ctx.fillStyle = T.bg; ctx.fillRect(0, 0, w, h);
      const mcx = w / 2, mcy = h / 2, R0 = Math.min(w, h) / 2 - 26, unit = R0 / 1.5 / Im * 0.98, i3 = TP.phases(tsim, Im, c.f, seq), F = TP.mmfVector(tsim, Im, c.f, seq, 1);
      ctx.strokeStyle = T.border; ctx.lineWidth = 1; ctx.beginPath(); ctx.arc(mcx, mcy, R0, 0, 2 * PI); ctx.stroke();
      ctx.beginPath(); ctx.arc(mcx, mcy, unit * 1.5, 0, 2 * PI); ctx.setLineDash([3, 4]); ctx.stroke(); ctx.setLineDash([]);
      const axn = ['A', 'B', 'C'];
      for (let k = 0; k < 3; k++) {
        const a = k * 2 * PI / 3, ux = Math.cos(a), uy = -Math.sin(a);
        ctx.strokeStyle = T.border; ctx.beginPath(); ctx.moveTo(mcx - ux * R0 * 0.2, mcy - uy * R0 * 0.2); ctx.lineTo(mcx + ux * R0, mcy + uy * R0); ctx.stroke();
        ctx.fillStyle = cols[k]; ctx.font = 'bold 12px ' + T.mono; ctx.textAlign = 'center'; ctx.textBaseline = 'middle'; ctx.fillText(axn[k], mcx + ux * (R0 + 13), mcy + uy * (R0 + 13));
        // pulsating component of this winding along its own axis (offset slightly perpendicular for visibility)
        const m = i3[k] * unit, nx = -uy * 5, ny = ux * 5; ctx.strokeStyle = cols[k]; ctx.lineWidth = 4; ctx.beginPath(); ctx.moveTo(mcx + nx, mcy + ny); ctx.lineTo(mcx + nx + ux * m, mcy + ny + uy * m); ctx.stroke();
      }
      const fx = mcx + F.x * unit, fy = mcy - F.y * unit; ctx.strokeStyle = T.pink; ctx.fillStyle = T.pink; ctx.lineWidth = 3.2; ctx.beginPath(); ctx.moveTo(mcx, mcy); ctx.lineTo(fx, fy); ctx.stroke();
      const aa = Math.atan2(fy - mcy, fx - mcx); ctx.beginPath(); ctx.moveTo(fx, fy); ctx.lineTo(fx - 12 * Math.cos(aa - 0.35), fy - 12 * Math.sin(aa - 0.35)); ctx.lineTo(fx - 12 * Math.cos(aa + 0.35), fy - 12 * Math.sin(aa + 0.35)); ctx.closePath(); ctx.fill();
      ctx.fillStyle = T.text; ctx.font = '11px ' + T.mono; ctx.textAlign = 'left'; ctx.textBaseline = 'top'; ctx.fillText('F = ' + fmt(F.mag, 3) + ' (1.5·Im)', 6, 5);
      ctx.fillStyle = T.muted; ctx.fillText('thick bars: windings A, B, C', 6, 19); ctx.fillText('arrow: resultant MMF', 6, 33);
      const deg = ((F.ang * 180 / PI) % 360 + 360) % 360;
      ui.dir.textContent = 'sequence ' + (seq > 0 ? 'a-b-c' : 'a-c-b (b,c swapped)') + ' → rotation ' + (seq > 0 ? 'counter-clockwise' : 'clockwise') + ', angle ' + fmt(deg, 0) + '°';
      ui.swap.setAttribute('aria-pressed', String(seq < 0));
      return { c, q, sets };
    }
    let lastCfgKey = '';
    function redraw() {
      if (!ui.cP) return; const r = draw(); if (!r) { lastCfgKey = ''; return; } const c = r.c, q = r.q, sets = r.sets;
      const key = JSON.stringify([c, seq]); if (key === lastCfgKey) return; lastCfgKey = key;
      const rpm2 = 60 * c.f;
      hudSet(ui.hud, [
        ['phase voltage Vph', fmt(q.Vph, 2) + ' V'], ['line voltage VL', fmt(q.VL, 2) + ' V'],
        ['phase current Iph', fmt(q.Iph, 3) + ' A'], ['line current IL', fmt(q.IL, 3) + ' A'],
        ['VL/Vph', fmt(q.VL / q.Vph, 4) + (c.conn === 'Y' ? '  (√3)' : '  (1)')], ['IL/Iph', fin(q.IL / q.Iph) ? fmt(q.IL / q.Iph, 4) + (c.conn === 'D' ? '  (√3)' : '  (1)') : '—'],
        ['power factor', fmt(q.pf, 3) + (q.phi >= 0 ? ' lagging' : ' leading')], ['P (total)', fmt(q.P / 1000, 3) + ' kW'],
        ['Q (total)', fmt(q.Q / 1000, 3) + ' kvar'], ['S (total)', fmt(q.S / 1000, 3) + ' kVA'],
        ['field speed (2-pole)', fmt(rpm2, 0) + ' rev/min'], ['direction', seq > 0 ? 'CCW (a→b→c)' : 'CW (reversed)'],
      ]);
      ui.work.set([
        'Phases: a = cos ωt, b = cos(ωt − ' + (seq > 0 ? '' : '(−)') + '120°), c = cos(ωt + ' + (seq > 0 ? '' : '(−)') + '120°); a + b + c = 0 at every instant.',
        c.conn === 'Y' ? 'Star: Vph = VL/√3 = ' + fmt(q.Vph, 3) + ' V;  IL = Iph. Vab = Van − Vbn = √3 Van ∠+30° (for a-b-c).' : 'Delta: Vph = VL;  IL = √3 Iph = ' + fmt(q.IL, 3) + ' A. Ia = Iab − Ica = √3 Iab ∠−30°.',
        '|Z| = √(R²+X²) = ' + fmt(q.Zm, 4) + ' Ω, φ = atan(X/R) = ' + fmt(q.phi * 180 / PI, 2) + '°; Iph = Vph/|Z| = ' + fmt(q.Iph, 4) + ' A',
        'S = 3 Vph Iph = √3 VL IL = ' + fmt(q.S, 2) + ' VA;  P = S cos φ = ' + fmt(q.P, 2) + ' W;  Q = S sin φ = ' + fmt(q.Q, 2) + ' var',
        'Windings are 120° apart in space: F(θ,t) = Σ N·ik(t)·cos(θ − k·120°) = (3/2) N Im cos(θ − ωt)  → constant magnitude 1.5·N·Im, angle ωt, i.e. one revolution per electrical cycle (ns = 120f/P).',
        'Swapping two phases makes the sequence a-c-b, so the angle becomes −ωt: the field (and motor) reverses.',
      ]);
    }
    function loop(ts) {
      raf = 0; if (!host.active || host.sub !== 'three' || document.hidden || paused) return;
      const dt = Math.min(0.1, Math.max(0, (ts - last) / 1000)); last = ts; tsim += dt * val(st, 'slow'); draw(); raf = requestAnimationFrame(loop);
    }
    function refreshAnim() {
      const want = host.active && host.sub === 'three' && !document.hidden && !paused;
      if (want && !raf) { last = (typeof performance !== 'undefined' ? performance.now() : 0); raf = requestAnimationFrame(loop); }
      else if (!want && raf) { cancelAnimationFrame(raf); raf = 0; }
      if (!want) redraw();
    }
    return {
      update: redraw, refreshAnim, stop() { if (raf) { cancelAnimationFrame(raf); raf = 0; } },
      isRunning: () => !!raf,
      get: () => ({ VL: val(st, 'VL'), f: val(st, 'f'), R: val(st, 'R'), X: val(st, 'X'), conn: ui.conn.value, seq: seq, slow: val(st, 'slow') }),
      set(o) {
        ['VL', 'f', 'R', 'X', 'slow'].forEach(k => { if (o[k] !== undefined) { const v = parseFloat(o[k]); if (fin(v)) st[k].set(v, true); } });
        if (o.conn === 'Y' || o.conn === 'D') ui.conn.value = o.conn;
        if (o.seq === '-1') seq = -1; else if (o.seq === '1') seq = 1;
      },
    };
  }

  /* ---------- tab ---------- */
  const host = { active: false, sub: 'induction' };
  let subs = null, panels = {}, btns = {};
  function showSub(id, silent) {
    host.sub = id;
    Object.keys(panels).forEach(k => { panels[k].hidden = k !== id; btns[k].setAttribute('aria-pressed', String(k === id)); btns[k].classList.toggle('active', k === id); });
    if (host.active) { subs[id].update(); subs.three.refreshAnim(); }
    if (!silent) FSP.state.touch();
  }
  FSP.registerTab({
    id: 'rotating', title: 'Rotating Machines',
    init(panel) {
      const bar = FSP.ui.el('div', { class: 'seg row', role: 'group', 'aria-label': 'Rotating machine topic' }); panel.appendChild(bar);
      const defs = [['induction', 'Induction motor'], ['dc', 'DC machine'], ['three', 'Three-phase & rotating field']];
      defs.forEach(d => { btns[d[0]] = FSP.ui.el('button', { type: 'button', class: 'seg-btn', 'aria-pressed': 'false', text: d[1], onclick: () => showSub(d[0]) }); bar.appendChild(btns[d[0]]); });
      defs.forEach(d => { panels[d[0]] = FSP.ui.el('div', { hidden: '' }); panel.appendChild(panels[d[0]]); });
      subs = { induction: buildInduction(panels.induction), dc: buildDC(panels.dc), three: buildThree(panels.three, host) };
      FSP.state.bind('rotating', {
        get() {
          const o = { sub: host.sub }, add = (pre, obj) => Object.keys(obj).forEach(k => { o[pre + k] = obj[k]; });
          add('im_', subs.induction.get()); add('dc_', subs.dc.get()); add('tp_', subs.three.get()); return o;
        },
        set(o) {
          const pick = pre => { const r = {}; Object.keys(o || {}).forEach(k => { if (k.indexOf(pre) === 0) r[k.slice(pre.length)] = o[k]; }); return r; };
          subs.induction.set(pick('im_')); subs.dc.set(pick('dc_')); subs.three.set(pick('tp_'));
          showSub(o && panels[o.sub] ? o.sub : 'induction', true);
        },
      });
      if (typeof document !== 'undefined') document.addEventListener('visibilitychange', () => { if (subs) subs.three.refreshAnim(); });
      showSub(host.sub, true);
    },
    activate() { host.active = true; subs[host.sub].update(); subs.three.refreshAnim(); },
    deactivate() { host.active = false; if (subs) subs.three.stop(); },
  });
})();
