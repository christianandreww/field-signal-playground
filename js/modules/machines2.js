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

  /* ======================= 8c+: induction-motor speed control (pure math) ======================= */
  // Every method returns a modified parameter set p2 for the SAME torque function IM.torqueTh / IM.thevenin.
  //  (a) stator voltage control : VL -> k·VL.   Vth ∝ V, so T(s) -> k²·T(s) at every slip; sm unchanged.
  //  (b) added rotor resistance : R2' -> R2' + Rext (Rext referred to the stator, = Rext_actual/(turns ratio)²).
  //      K = |Zth + jX2'| does not contain R2', so sm = R2'/K scales with R2' and Tmax is unchanged;
  //      the whole curve is the old one with s replaced by s·R2'/(R2'+Rext).
  //  (c) constant V/f : f -> a·f (a = f'/f), reactances X1, X2', Xm -> a·X, V -> V·(a + b(1−a)) for a<1 (b = low-frequency boost
  //      fraction, b = 0: V ∝ f), constant V above base frequency (a>1, field weakening).  Rth ≈ R1·(Xm/(X1+Xm))² does NOT scale with a.
  //      Exactly constant Tmax needs R1 = 0 (then Rth = 0, Zth = jXth with Xth ∝ a):
  //        Tmax = 3 Vth² / (2 ωs' (Xth + X2')),  Vth = V'·Xm/(X1+Xm) ∝ V' ∝ a,  ωs' ∝ a,  Xth+X2' ∝ a   =>   Tmax ∝ a²/(a·a) = const.
  //      and sm = R2'/(Xth+X2') ∝ 1/a, so the speed drop ns'·sm = ns·sm0 is the same for every frequency.
  const SC = {
    vfFrac: (a, b) => (a >= 1 ? 1 : a + (b || 0) * (1 - a)),
    control(p, mode, o) {
      if (mode === 'volt') {
        if (!(o.k > 0) || !(o.k <= 1.5)) return { ok: false, msg: 'Voltage ratio V/V0 must be in (0, 1.5].' };
        return { ok: true, p2: Object.assign({}, p, { VL: p.VL * o.k }), a: 1, vfrac: o.k };
      }
      if (mode === 'rot') {
        if (!(o.Rext >= 0) || !fin(o.Rext)) return { ok: false, msg: 'Added rotor resistance must be ≥ 0.' };
        return { ok: true, p2: Object.assign({}, p, { R2: p.R2 + o.Rext }), a: 1, vfrac: 1 };
      }
      if (mode === 'vf') {
        if (!(o.f > 0) || !fin(o.f)) return { ok: false, msg: 'Frequency must be > 0.' };
        if (!(o.boost >= 0) || !(o.boost <= 1)) return { ok: false, msg: 'Boost must be between 0 and 100 % of rated voltage.' };
        const a = o.f / p.f, vfrac = SC.vfFrac(a, o.boost);
        return { ok: true, p2: Object.assign({}, p, { VL: p.VL * vfrac, f: o.f, X1: p.X1 * a, X2: p.X2 * a, Xm: p.Xm * a }), a, vfrac };
      }
      return { ok: false, msg: 'Unknown control method.' };
    },
    // load torque at speed n (rpm): load = {type:'const'|'fan', T (N·m at n_ref), n (rpm, n_ref)}
    loadTorque(load, n) { return load.type === 'fan' ? load.T * (n / load.n) * (n / load.n) : load.T; },
    // Stable operating point: T(s) = T_L(ns(1−s)) on 0 < s <= sm.  T(s) increases and T_L(n(s)) does not increase with s there,
    // so g(s) = T − T_L is strictly increasing and the root is unique (bisection).
    operatingPoint(p, load) {
      if (!(load.T >= 0) || !fin(load.T)) return { ok: false, msg: 'Load torque must be ≥ 0.' };
      if (load.type === 'fan' && !(load.n > 0)) return { ok: false, msg: 'Reference speed must be > 0.' };
      const th = IM.thevenin(p), bd = IM.breakdown(p, th), ns = IM.syncRpm(p.f, p.poles);
      if (!(bd.Tmax > 0) || !fin(bd.Tmax) || !(bd.sm > 0)) return { ok: false, msg: 'Parameters give no valid torque curve.' };
      const g = s => IM.torqueTh(p, s, th) - SC.loadTorque(load, ns * (1 - s));
      const Tstart = IM.torqueTh(p, 1, th), TLstart = SC.loadTorque(load, 0);
      const res = { th, bd, ns, Tstart, TLstart, canStart: Tstart > TLstart };
      if (load.T === 0) return Object.assign(res, { ok: true, s: 0, n: ns, T: 0, TL: 0 });
      if (g(bd.sm) < 0) return Object.assign(res, { ok: false, stalled: true, msg: 'Load exceeds the breakdown torque at this setting: no stable operating point (motor stalls).' });
      let lo = 0, hi = bd.sm;
      for (let i = 0; i < 200; i++) { const m = 0.5 * (lo + hi); if (g(m) < 0) lo = m; else hi = m; if (hi - lo < 1e-15) break; }
      const s = 0.5 * (lo + hi), n = ns * (1 - s);
      return Object.assign(res, { ok: true, s, n, T: IM.torqueTh(p, s, th), TL: SC.loadTorque(load, n) });
    },
  };

  /* ======================= 8d+: compound DC machine ======================= */
  // Linear flux per ampere-turn (motor) or Fröhlich saturation (generator).  Field MMF in shunt-field-ampere units:
  //   I_fe = I_f + σ (Ns/Nf) I_s ,  σ = +1 cumulative, −1 differential,   K = Laf·I_fe  (/(1+|I_fe|/Isat) when Isat > 0)
  //   E = K ω,  Tem = K Ia.
  // long shunt : shunt field across the supply/load terminals, series field in series with the armature   (Is = Ia)
  // short shunt: shunt field across the armature only, series field in the line                              (Is = IL)
  // p: {V, Ra, Rf, Rs, Nf, Ns, Laf, Vb, conn:'long'|'short', comp:'cum'|'diff', Isat (0 = linear), w (generator speed, rad/s)}
  const CD = {
    sigma: p => (p.comp === 'diff' ? -1 : 1),
    kOf(p, Ife) { return p.Isat > 0 && fin(p.Isat) ? p.Laf * Ife / (1 + Math.abs(Ife) / p.Isat) : p.Laf * Ife; },
    motor(p, Ia) {
      const r = p.Ns / p.Nf, sg = CD.sigma(p), vb = Ia > 0 ? p.Vb : 0;
      let If, IL, Is, Va, Ea;
      if (p.conn === 'short') { If = (p.V - Ia * p.Rs) / (p.Rf + p.Rs); IL = Ia + If; Is = IL; Va = p.V - IL * p.Rs; Ea = Va - vb - Ia * p.Ra; }
      else { If = p.V / p.Rf; IL = Ia + If; Is = Ia; Va = p.V; Ea = p.V - vb - Ia * (p.Ra + p.Rs); }
      const Ife = If + sg * r * Is, K = CD.kOf(p, Ife), okK = K > 1e-9;
      const w = okK ? Ea / K : NaN;
      return { Ia, If, Is, IL, Va, Ea, Ife, K, w, rpm: fin(w) ? w * 30 / PI : NaN, Tem: K * Ia, Pem: Ea * Ia, Pin: p.V * IL, runaway: !okK && Ea > 0, fluxRev: Ife <= 0, stalled: Ea < 0, Vb: vb };
    },
    // residual h(Vt) = E produced − E required at terminal voltage Vt and load current IL (generator)
    genState(p, Vt, IL) {
      const r = p.Ns / p.Nf, sg = CD.sigma(p), vb = p.Vb;
      let If, Ia, Is, Va, Ereq;
      if (p.conn === 'short') { Va = Vt + IL * p.Rs; If = Va / p.Rf; Ia = IL + If; Is = IL; Ereq = Va + Ia * p.Ra + (Ia > 0 ? vb : 0); }
      else { If = Vt / p.Rf; Ia = IL + If; Is = Ia; Va = Vt; Ereq = Vt + Ia * (p.Ra + p.Rs) + (Ia > 0 ? vb : 0); }
      const Ife = If + sg * r * Is, K = CD.kOf(p, Ife), Eprod = p.w * K;
      return { Vt, IL, If, Ia, Is, Va, Ife, K, Ea: Eprod, Ereq, h: Eprod - Ereq };
    },
    // largest self-consistent terminal voltage for load current IL; collapsed = no solution (voltage collapses to ~0)
    generator(p, IL) {
      if (!(p.Isat > 0) || !fin(p.Isat)) return { ok: false, msg: 'A self-excited generator needs a magnetic saturation knee Isat > 0 (the linear model has no stable voltage).' };
      if (!(p.w > 0)) return { ok: false, msg: 'Speed must be > 0.' };
      const Vhi = p.w * p.Laf * p.Isat * (1 + 1e-9), N = 400;
      let prev = Vhi, hit = -1;
      for (let i = 0; i <= N; i++) {
        const v = Vhi * (1 - i / N), h = CD.genState(p, v, IL).h;
        if (h >= -1e-12) { hit = i; break; }
        prev = v;
      }
      if (hit < 0) return Object.assign(CD.genState(p, 0, IL), { ok: true, collapsed: true, Vt: 0 });
      if (hit === 0) return Object.assign(CD.genState(p, Vhi, IL), { ok: true });
      let lo = Vhi * (1 - hit / N), hi = prev;      // h(lo) >= 0 > h(hi)
      for (let i = 0; i < 200; i++) { const m = 0.5 * (lo + hi); if (CD.genState(p, m, IL).h >= 0) lo = m; else hi = m; if (hi - lo < 1e-13 * Vhi) break; }
      return Object.assign(CD.genState(p, lo, IL), { ok: true, collapsed: false });
    },
  };

  FSP.math.rotating = {
    syncRpm: IM.syncRpm, wsync: IM.wsync, imVphase: IM.vphase, imThevenin: IM.thevenin, imTorqueTh: IM.torqueTh,
    imSolveFull: IM.solveFull, imTorqueFull: IM.torqueFull, imBreakdown: IM.breakdown, imOp: IM.op, imMode: IM.mode, imSweep: IM.sweep,
    dcSolve: DC.solve, dcStallCurrent: DC.stallCurrent,
    imControl: SC.control, imLoadTorque: SC.loadTorque, imOperatingPoint: SC.operatingPoint, imVfFrac: SC.vfFrac,
    dcCompoundMotor: CD.motor, dcCompoundGen: CD.generator, dcCompoundGenState: CD.genState,
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

  /* ---------- view switch inside a sub-topic ---------- */
  function viewBar(root, defs, label, onPick) {
    const bar = FSP.ui.el('div', { class: 'seg row', role: 'group', 'aria-label': label }); bar.style.marginBottom = '8px'; const btn = {};
    defs.forEach(d => { btn[d[0]] = FSP.ui.el('button', { type: 'button', class: 'seg-btn', 'aria-pressed': 'false', text: d[1], onclick: () => onPick(d[0]) }); bar.appendChild(btn[d[0]]); });
    root.appendChild(bar);
    return { set(id) { Object.keys(btn).forEach(k => { btn[k].setAttribute('aria-pressed', String(k === id)); btn[k].classList.toggle('active', k === id); }); } };
  }
  const pick = (v, list, dflt) => (list.indexOf(v) >= 0 ? v : dflt);
  const setNum = (store, o) => Object.keys(store).forEach(k => { if (o[k] !== undefined) { const v = parseFloat(o[k]); if (fin(v)) store[k].set(v, true); } });

  /* ---------- 8c+: induction motor speed control ---------- */
  function buildSpeedControl(lay, getP) {
    const st = {}, ui = {};
    const ctl = FSP.ui.el('div', { class: 'controls' }), stage = FSP.ui.el('div', { class: 'stage' }); lay.appendChild(ctl); lay.appendChild(stage);
    const touch = () => { update(); FSP.state.touch(); };
    const fs0 = FSP.ui.fieldset(ctl, 'Speed-control method');
    ui.mode = FSP.ui.select(fs0, 'Method', [['volt', '(a) Stator voltage control'], ['rot', '(b) Added rotor resistance'], ['vf', '(c) Constant V/f (frequency)']], 'volt', touch);
    fs0.appendChild(FSP.ui.el('div', { class: 'note', text: 'Machine data (V, f, poles, R1, X1, R2′, X2′, Xm) come from the Torque–slip analysis view. The three methods keep their own settings; the plot shows the selected one against the base curve.' }));
    const fs1 = FSP.ui.fieldset(ctl, 'Method settings');
    makeSliders(fs1, [
      { k: 'kv', l: 'V / V0', min: 20, max: 120, v: 80, step: 1, u: '%' },
      { k: 'rx', l: 'Added R2′ (ext.)', min: 0, max: 10, v: 1.2, u: 'Ω' },
      { k: 'fn', l: 'Frequency f′', min: 1, max: 200, v: 25, u: 'Hz' },
      { k: 'bst', l: 'LF voltage boost', min: 0, max: 30, v: 0, step: 0.5, u: '% V0' }], st, update);
    ui.rnote = FSP.ui.el('div', { class: 'note', text: 'Added R2′ is the external rotor resistance already referred to the stator (R_ext,actual ÷ turns-ratio²). Boost: V = V0·[a + b(1−a)] for a = f′/f < 1, constant V above base frequency.' }); fs1.appendChild(ui.rnote);
    const fs2 = FSP.ui.fieldset(ctl, 'Load torque');
    ui.lt = FSP.ui.select(fs2, 'Load type', [['const', 'Constant torque'], ['fan', 'Fan / pump: T_L ∝ n²']], 'const', touch);
    makeSliders(fs2, [{ k: 'lT', l: 'T_L (at n_ref)', min: 0, max: 500, v: 50, u: 'N·m' }, { k: 'ln', l: 'n_ref', min: 50, max: 6000, v: 1440, u: 'rpm' }], st, update);
    ui.msg = FSP.ui.el('div', { class: 'msg bad', hidden: '' }); ctl.appendChild(ui.msg);
    ui.c1 = wrapCanvas(stage, 'Induction motor torque versus speed with load torque curve for the selected speed-control method', 340); ui.c1.onResize(update);
    ui.hud = FSP.ui.el('div', { class: 'hud' }); stage.appendChild(ui.hud);
    ui.work = FSP.ui.working(stage);
    const load = () => ({ type: ui.lt.value, T: st.lT.get(), n: st.ln.get() });
    const opts = () => ({ k: st.kv.get() / 100, Rext: st.rx.get(), f: st.fn.get(), boost: st.bst.get() / 100 });
    const NAMES_M = { volt: '(a) stator voltage', rot: '(b) rotor resistance', vf: '(c) V/f' };

    function update() {
      if (!ui.c1) return;
      const mode = ui.mode.value, p = getP(), L = load(), o = opts(), T = theme();
      st.kv.el.hidden = mode !== 'volt'; st.rx.el.hidden = mode !== 'rot'; st.fn.el.hidden = mode !== 'vf'; st.bst.el.hidden = mode !== 'vf'; ui.rnote.hidden = mode === 'volt';
      st.ln.el.hidden = L.type !== 'fan';
      const base = SC.operatingPoint(p, L), ctrl = {};
      ['volt', 'rot', 'vf'].forEach(m => { const c = SC.control(p, m, o); ctrl[m] = c.ok ? Object.assign({ c }, SC.operatingPoint(c.p2, L)) : { c, ok: false, msg: c.msg }; });
      const sel = ctrl[mode], msgs = [];
      if (!sel.c.ok) msgs.push(sel.c.msg); else if (!sel.ok) msgs.push(sel.msg);
      if (!base.ok) msgs.push('Base case: ' + base.msg);
      if (sel.c.ok && sel.th && !sel.canStart) msgs.push('Note: starting torque ' + fmt(sel.Tstart, 1) + ' N·m is below the load torque at standstill (' + fmt(sel.TLstart, 1) + ' N·m): the motor cannot start against this load with this setting.');
      ui.msg.hidden = msgs.length === 0; ui.msg.textContent = msgs.join(' ');
      if (!base.th || !sel.c.ok || !sel.th) { hudSet(ui.hud, [['status', 'fix the inputs']]); ui.work.set('Inputs are invalid for this setting: ' + msgs.join(' ')); drawPlot(ui.c1, { xmin: 0, xmax: 1, ymin: 0, ymax: 1, xlabel: '', ylabel: '' }); return; }
      const p2 = sel.c.p2, ns0 = base.ns, ns2 = sel.ns, nmax = Math.max(ns0, ns2) * 1.04;
      const curve = (pp, th, ns) => { const x = [], y = []; for (let i = 0; i <= 240; i++) { const s = 1 - i / 240; x.push(ns * (1 - s)); y.push(IM.torqueTh(pp, Math.abs(s) < 1e-12 ? 0 : s, th)); } return { x, y }; };
      const c0 = curve(p, base.th, ns0), c2 = curve(p2, sel.th, ns2);
      const lx = [], ly = []; for (let i = 0; i <= 160; i++) { const n = nmax * i / 160; lx.push(n); ly.push(SC.loadTorque(L, n)); }
      const tpk = Math.max(base.bd.Tmax, sel.bd.Tmax), ymax = Math.max(tpk, Math.min(Math.max.apply(null, ly), tpk * 1.5), 1e-6) * 1.12;
      const pts = [];
      if (base.ok) pts.push({ x: base.n, y: base.T, label: 'base ' + fmt(base.n, 0), color: T.muted, below: true });
      if (sel.ok) pts.push({ x: sel.n, y: sel.T, label: 'new ' + fmt(sel.n, 0) + ' rpm', color: T.text });
      drawPlot(ui.c1, {
        xmin: 0, xmax: nmax, ymin: 0, ymax, xlabel: 'speed n (rpm)', ylabel: 'torque (N·m)',
        series: [{ x: c0.x, y: c0.y, color: T.muted, width: 1.6, dash: [5, 3] }, { x: c2.x, y: c2.y, color: T.c, width: 2.4 }, { x: lx, y: ly, color: T.pink, width: 2 }],
        vlines: [{ x: ns0 }, { x: ns2 }], points: pts,
        legend: [{ text: 'load T_L(n)', color: T.pink }, { text: NAMES_M[mode], color: T.c }, { text: 'base curve', color: T.muted }],
      });
      const row = (m, lab) => { const r = ctrl[m]; return [lab, r.ok ? fmt(r.n, 1) + ' rpm, s = ' + fmt(r.s, 4) + ', T = ' + fmt(r.T, 1) + ' N·m' : (r.c && !r.c.ok ? 'invalid input' : 'stalls (no stable point)')]; };
      const opSel = sel.ok ? IM.op(p2, sel.s) : null;
      hudSet(ui.hud, [
        ['base (V0, f0, R2′)', base.ok ? fmt(base.n, 1) + ' rpm, s = ' + fmt(base.s, 4) + ', T = ' + fmt(base.T, 1) + ' N·m' : 'stalls'],
        row('volt', '(a) V = ' + fmt(st.kv.get(), 0) + ' % V0'), row('rot', '(b) R2′ + ' + fmt(st.rx.get(), 2) + ' Ω'), row('vf', '(c) f′ = ' + fmt(st.fn.get(), 1) + ' Hz'),
        ['selected: ns′', fmt(ns2, 1) + ' rpm'], ['breakdown slip sm', fmt(sel.bd.sm, 4)], ['breakdown torque', fmt(sel.bd.Tmax, 2) + ' N·m  (' + fmt(100 * sel.bd.Tmax / base.bd.Tmax, 1) + ' % of base)'],
        ['starting torque', fmt(sel.Tstart, 2) + ' N·m'], ['stator current', opSel ? fmt(opSel.I1, 2) + ' A' : '—'], ['power factor', opSel ? fmt(opSel.pf, 3) : '—'],
        ['efficiency', opSel && fin(opSel.eff) ? fmt(opSel.eff * 100, 1) + ' %' : '—'], ['speed drop ns′ − n', sel.ok ? fmt(ns2 - sel.n, 1) + ' rpm' : '—'],
      ]);
      // ---- working
      const W = [], a2 = sel.bd;
      W.push('Machine: ' + p.conn + ' connection, VL = ' + fmt(p.VL, 1) + ' V → Vph = ' + fmt(IM.vphase(p), 3) + ' V, ' + p.poles + ' poles, f = ' + fmt(p.f, 2) + ' Hz → ns = 120f/P = ' + fmt(ns0, 1) + ' rpm.');
      W.push('Torque (the same equivalent-circuit/Thevenin formula as the analysis view):');
      W.push('  Vth = Vph·jXm/(R1 + j(X1+Xm)),  Zth = jXm(R1+jX1)/(R1+j(X1+Xm)) = Rth + jXth,  ωs = 4πf/P');
      W.push('  T(s) = 3|Vth|²·(R2′/s) / ( ωs·[ (Rth + R2′/s)² + (Xth + X2′)² ] ),   sm = R2′/√(Rth² + (Xth+X2′)²),   Tmax = 3|Vth|² / (2ωs(Rth + √(Rth² + (Xth+X2′)²)))');
      W.push('Base: |Vth| = ' + fmt(base.th.VthMag, 3) + ' V, Zth = ' + fmt(base.th.Rth, 4) + ' + j' + fmt(base.th.Xth, 4) + ' Ω, sm = ' + fmt(base.bd.sm, 4) + ', Tmax = ' + fmt(base.bd.Tmax, 2) + ' N·m, T(start) = ' + fmt(base.Tstart, 2) + ' N·m.');
      W.push('Load: ' + (L.type === 'fan' ? 'T_L(n) = T_ref·(n/n_ref)² = ' + fmt(L.T, 2) + '·(n/' + fmt(L.n, 0) + ')² N·m' : 'T_L = ' + fmt(L.T, 2) + ' N·m (constant)') + ';  n = ns(1 − s).  Operating point: T(s) = T_L(n(s)) on 0 < s ≤ sm (bisection; T rises and T_L falls with s, so the root is unique).');
      if (base.ok) W.push('Base operating point: s = ' + fmt(base.s, 5) + ', n = ' + fmt(ns0, 1) + '·(1 − ' + fmt(base.s, 5) + ') = ' + fmt(base.n, 2) + ' rpm, T = ' + fmt(base.T, 3) + ' N·m.');
      W.push('');
      const c = sel.c;
      if (mode === 'volt') {
        W.push('(a) Stator voltage control: V = ' + fmt(o.k, 4) + '·V0 = ' + fmt(p2.VL, 2) + ' V line.');
        W.push('  Vth ∝ V, so T ∝ V² at the same slip:  T′(s) = (V/V0)²·T(s) = ' + fmt(o.k * o.k, 4) + '·T(s).  Slip at breakdown is unchanged (sm = ' + fmt(a2.sm, 4) + '); Tmax′ = ' + fmt(o.k * o.k, 4) + ' × ' + fmt(base.bd.Tmax, 2) + ' = ' + fmt(a2.Tmax, 2) + ' N·m.');
        W.push('  Speed range is small: the stable region is 0 < s < sm, so n can only fall to ns(1 − sm) = ' + fmt(ns0 * (1 - a2.sm), 1) + ' rpm (and a constant-torque load needs T_L < Tmax′). Rotor copper loss s·Pag grows with s.');
      } else if (mode === 'rot') {
        W.push('(b) Added rotor resistance: R2′_total = ' + fmt(base.th ? p.R2 : 0, 4) + ' + ' + fmt(o.Rext, 4) + ' = ' + fmt(p2.R2, 4) + ' Ω (referred to the stator).');
        W.push('  K = √(Rth² + (Xth+X2′)²) = ' + fmt(base.th.K, 4) + ' Ω does not contain R2′, so sm′ = R2′_total/K = ' + fmt(p2.R2, 4) + '/' + fmt(base.th.K, 4) + ' = ' + fmt(a2.sm, 4) + ' (× ' + fmt(p2.R2 / p.R2, 3) + ') and Tmax′ = ' + fmt(a2.Tmax, 2) + ' N·m is unchanged.');
        W.push('  T depends on s only through R2′/s, so the same torque occurs at s′ = s·R2′_total/R2′. Rotor loss is dissipated mostly in the external resistor; wound-rotor motors only.');
      } else {
        const a = c.a, vf = c.vfrac;
        W.push('(c) Constant V/f: a = f′/f = ' + fmt(o.f, 3) + '/' + fmt(p.f, 3) + ' = ' + fmt(a, 4) + ';  ns′ = 120f′/P = ' + fmt(ns2, 2) + ' rpm.');
        W.push('  Reactances scale with frequency: X1′ = ' + fmt(p2.X1, 4) + ', X2′′ = ' + fmt(p2.X2, 4) + ', Xm′ = ' + fmt(p2.Xm, 3) + ' Ω; R1, R2′ unchanged.');
        W.push('  Voltage: V′ = V0·[a + b(1 − a)]' + (a >= 1 ? ' → a ≥ 1: constant V0 (field weakening)' : ' = ' + fmt(p.VL, 1) + '·[' + fmt(a, 4) + ' + ' + fmt(o.boost, 3) + '·' + fmt(1 - a, 4) + '] = ' + fmt(p2.VL, 2) + ' V line') + ';  V′/f′ = ' + fmt(vf / a, 4) + ' × (V0/f).');
        W.push('  Zth′ = ' + fmt(sel.th.Rth, 4) + ' + j' + fmt(sel.th.Xth, 4) + ' Ω, |Vth′| = ' + fmt(sel.th.VthMag, 3) + ' V;  sm′ = ' + fmt(a2.sm, 4) + ' (speed at breakdown ns′(1 − sm′) = ' + fmt(ns2 * (1 - a2.sm), 1) + ' rpm);  Tmax′ = ' + fmt(a2.Tmax, 2) + ' N·m = ' + fmt(100 * a2.Tmax / base.bd.Tmax, 1) + ' % of base.');
        W.push('  Why Tmax is ~constant: with R1 = 0, Rth = 0 and Tmax = 3Vth²/(2ωs(Xth+X2′)); Vth ∝ V′ ∝ a, ωs ∝ a, (Xth+X2′) ∝ a, so Tmax ∝ a²/(a·a) = constant and sm ∝ 1/a. With R1 ≠ 0 (Rth does not scale) or with boost it is only approximately constant: it falls at low f unless boost compensates.');
      }
      W.push('');
      if (sel.ok) {
        W.push('New operating point: solve T′(s) = T_L(n):  s = ' + fmt(sel.s, 5) + ',  n = ns′(1 − s) = ' + fmt(ns2, 2) + '·(1 − ' + fmt(sel.s, 5) + ') = ' + fmt(sel.n, 2) + ' rpm.');
        W.push('  Check: T′(s) = ' + fmt(sel.T, 4) + ' N·m,  T_L(n) = ' + fmt(sel.TL, 4) + ' N·m.   Speed change vs base: ' + (base.ok ? fmt(sel.n - base.n, 2) + ' rpm.' : '—'));
        W.push('  Stator current ' + fmt(opSel.I1, 2) + ' A at pf ' + fmt(opSel.pf, 3) + (fin(opSel.eff) ? ', efficiency ' + fmt(opSel.eff * 100, 1) + ' %.' : '.'));
      } else W.push('No stable operating point: T_L exceeds the breakdown torque (or inputs invalid).');
      ui.work.set(W);
    }
    return {
      update,
      get: () => ({ mode: ui.mode.value, kv: st.kv.get(), rx: st.rx.get(), fn: st.fn.get(), bst: st.bst.get(), ltype: ui.lt.value, lT: st.lT.get(), ln: st.ln.get() }),
      set(o) { ui.mode.value = pick(o.mode, ['volt', 'rot', 'vf'], 'volt'); ui.lt.value = pick(o.ltype, ['const', 'fan'], 'const'); setNum(st, o); },
    };
  }

  /* ---------- 8d+: compound DC machine ---------- */
  function buildCompound(lay) {
    const st = {}, ui = {};
    const ctl = FSP.ui.el('div', { class: 'controls' }), stage = FSP.ui.el('div', { class: 'stage' }); lay.appendChild(ctl); lay.appendChild(stage);
    const touch = () => { update(); FSP.state.touch(); };
    const fs0 = FSP.ui.fieldset(ctl, 'Compound machine');
    ui.mode = FSP.ui.select(fs0, 'Operation', [['motor', 'Motor (speed–torque)'], ['gen', 'Generator (voltage–load)']], 'motor', touch);
    ui.conn = FSP.ui.select(fs0, 'Connection', [['long', 'Long shunt'], ['short', 'Short shunt']], 'long', touch);
    ui.comp = FSP.ui.select(fs0, 'Compounding', [['cum', 'Cumulative (series aids shunt)'], ['diff', 'Differential (series opposes)']], 'cum', touch);
    fs0.appendChild(FSP.ui.el('div', { class: 'note', text: 'Long shunt: shunt field across the terminals, series field in series with the armature. Short shunt: shunt field across the armature only, series field in the line.' }));
    const fs1 = FSP.ui.fieldset(ctl, 'Machine');
    makeSliders(fs1, [
      { k: 'V', l: 'Supply V', min: 12, max: 600, v: 240, u: 'V' }, { k: 'N', l: 'Speed', min: 100, max: 4000, v: 1450, u: 'rpm' },
      { k: 'Ra', l: 'Ra', min: 0.01, max: 5, v: 0.5, u: 'Ω' }, { k: 'Rf', l: 'Rf (shunt)', min: 10, max: 500, v: 120, u: 'Ω' }, { k: 'Rs', l: 'Rs (series)', min: 0.01, max: 2, v: 0.1, u: 'Ω' },
      { k: 'Nf', l: 'Shunt turns/pole', min: 100, max: 5000, v: 1000, u: '' }, { k: 'Ns', l: 'Series turns/pole', min: 1, max: 200, v: 6, u: '' },
      { k: 'Laf', l: 'Laf', min: 0.1, max: 5, v: 1.2, step: 0.01, u: 'H' }, { k: 'Isat', l: 'Saturation knee', min: 0.5, max: 50, v: 4, step: 0.1, u: 'A (If)' },
      { k: 'Vb', l: 'Vbrush', min: 0, max: 10, v: 2, u: 'V' }], st, update);
    fs1.appendChild(FSP.ui.el('div', { class: 'note', text: 'Field MMF in shunt-ampere units: I_fe = If ± (Ns/Nf)·I_series. Motor: flux linear in I_fe (K = Laf·I_fe). Generator: Fröhlich saturation K = Laf·I_fe/(1 + I_fe/Isat); it needs saturation to settle at a stable voltage. Armature reaction ignored.' }));
    const fs2 = FSP.ui.fieldset(ctl, 'Load');
    makeSliders(fs2, [{ k: 'Ir', l: 'Rated load current', min: 5, max: 300, v: 50, u: 'A' }, { k: 'ld', l: 'Load', min: 0, max: 150, v: 100, step: 1, u: '% rated' }], st, update);
    ui.warn = FSP.ui.el('div', { class: 'msg warn', hidden: '' }); ctl.appendChild(ui.warn);
    ui.c1 = wrapCanvas(stage, 'Compound DC machine characteristic compared with shunt and series machines', 340); ui.c1.onResize(update);
    ui.hud = FSP.ui.el('div', { class: 'hud' }); stage.appendChild(ui.hud);
    ui.work = FSP.ui.working(stage);

    function params() {
      return { V: st.V.get(), Ra: st.Ra.get(), Rf: st.Rf.get(), Rs: st.Rs.get(), Nf: st.Nf.get(), Ns: st.Ns.get(), Laf: st.Laf.get(), Vb: st.Vb.get(), conn: ui.conn.value, comp: ui.comp.value, Isat: 0, w: st.N.get() * PI / 30 };
    }
    function update() {
      if (!ui.c1) return;
      const T = theme(), gen = ui.mode.value === 'gen', Ir = st.Ir.get(), x = Ir * st.ld.get() / 100, p0 = params();
      st.V.el.hidden = gen; st.N.el.hidden = !gen; st.Isat.el.hidden = !gen;
      const p = Object.assign({}, p0, { Isat: gen ? st.Isat.get() : 0 });
      const warns = [];
      const variants = pp => ({
        shunt: Object.assign({}, pp, { Ns: 0, Rs: 0, comp: 'cum' }),
        cum: Object.assign({}, pp, { comp: 'cum' }), diff: Object.assign({}, pp, { comp: 'diff' }),
      });
      const V = variants(p);
      let serP;   // series machine for comparison: same armature, series winding sized so its full-load MMF equals the shunt no-load MMF
      const sel = ui.comp.value;
      const cols = { shunt: T.muted, series: T.pink, cum: T.c, diff: T.a };
      const series = [], points = [], lines = [];
      if (!gen) {
        const ifNL = p.V / p.Rf; serP = Object.assign({}, p, { Rf: Infinity, Ns: p.Nf * ifNL / Ir, conn: 'long', comp: 'cum' });
        const m0 = CD.motor(V.shunt, 0), wcap = fin(m0.w) ? m0.w * 2.4 : Infinity, imax = 1.5 * Ir;
        const sw = pp => { const X = [], Y = []; for (let i = 0; i <= 150; i++) { const r = CD.motor(pp, imax * i / 150); if (fin(r.w) && r.w >= 0 && r.w <= wcap && r.Ea >= 0) { X.push(r.Tem); Y.push(r.rpm); } } return { x: X, y: Y }; };
        const defs = [['series', serP, 'series (same armature)'], ['shunt', V.shunt, 'shunt'], [sel === 'cum' ? 'diff' : 'cum', sel === 'cum' ? V.diff : V.cum, sel === 'cum' ? 'differential compound' : 'cumulative compound'], [sel, V[sel], (sel === 'cum' ? 'cumulative' : 'differential') + ' compound, ' + p.conn + ' shunt']];
        defs.forEach(d => { const c = sw(d[1]); series.push({ x: c.x, y: c.y, color: cols[d[0]], width: d[0] === sel ? 2.8 : 1.7, dash: d[0] === sel ? [] : (d[0] === 'series' || d[0] === 'shunt' ? [5, 3] : []) }); lines.push({ text: d[2], color: cols[d[0]] }); });
        const r = CD.motor(V[sel], x);
        if (fin(r.rpm)) points.push({ x: r.Tem, y: r.rpm, label: fmt(r.rpm, 0) + ' rpm', color: T.text });
        const rmaxs = Math.max.apply(null, series.reduce((a, s) => a.concat(s.y), [1]));
        drawPlot(ui.c1, { xmin: 0, xmax: Math.max.apply(null, series.reduce((a, s) => a.concat(s.x), [1e-6])) * 1.05, ymin: 0, ymax: rmaxs * 1.1, xlabel: 'electromagnetic torque (N·m)', ylabel: 'speed (rpm)', series, points, legend: lines.reverse() });
        if (r.runaway) warns.push('K ≈ 0 at this load: the speed is unbounded (runaway).'); else if (r.fluxRev) warns.push('Differential compounding has driven the net field to zero or negative at this load: the speed rises without limit / the motor is unstable. Reduce the series turns or the load.');
        if (r.stalled) warns.push('Back-EMF would be negative: load current beyond the stall current.');
        const nl = CD.motor(V[sel], 0), nlS = V.shunt, rs = CD.motor(V.shunt, x);
        const reg = (a, b) => (fin(a) && fin(b) && b > 0 ? fmt(100 * (a - b) / a, 1) + ' % fall' : '—');
        hudSet(ui.hud, [
          ['armature current Ia', fmt(r.Ia, 2) + ' A'], ['shunt field If', fmt(r.If, 3) + ' A'], ['series current Is', fmt(r.Is, 2) + ' A'], ['line current IL', fmt(r.IL, 2) + ' A'],
          ['field MMF I_fe', fmt(r.Ife, 4) + ' A'], ['K = Laf·I_fe', fmt(r.K, 4) + ' V·s/rad'], ['back-EMF Ea', fmt(r.Ea, 2) + ' V'], ['speed', fin(r.rpm) ? fmt(r.rpm, 1) + ' rpm' : '—'],
          ['Tem = K·Ia', fmt(r.Tem, 2) + ' N·m'], ['Pem = Ea·Ia', fmt(r.Pem / 1000, 3) + ' kW'], ['Pin = V·IL', fmt(r.Pin / 1000, 3) + ' kW'],
          ['no-load speed (this)', fin(nl.rpm) ? fmt(nl.rpm, 1) + ' rpm' : '—'], ['speed change NL → load', reg(nl.rpm, r.rpm) + (fin(nl.rpm) && fin(r.rpm) && r.rpm > nl.rpm ? ' (speed RISES)' : '')],
          ['plain shunt NL → load', reg(CD.motor(nlS, 0).rpm, rs.rpm)],
        ]);
        const cs = sel === 'cum' ? '+' : '−', r0 = CD.motor(V[sel], 0);
        const W = [];
        W.push((p.conn === 'long' ? 'Long shunt' : 'Short shunt') + ', ' + (sel === 'cum' ? 'cumulative' : 'differential') + ' compound motor.  Ia = ' + fmt(x, 2) + ' A (' + fmt(st.ld.get(), 0) + ' % of ' + fmt(Ir, 1) + ' A).');
        if (p.conn === 'long') {
          W.push('If = V/Rf = ' + fmt(p.V, 2) + '/' + fmt(p.Rf, 2) + ' = ' + fmt(r.If, 4) + ' A;   Is = Ia = ' + fmt(r.Is, 3) + ' A;   IL = Ia + If = ' + fmt(r.IL, 3) + ' A');
          W.push('Ea = V − Vb − Ia(Ra + Rs) = ' + fmt(p.V, 2) + ' − ' + fmt(r.Vb, 2) + ' − ' + fmt(x, 2) + '·(' + fmt(p.Ra, 3) + ' + ' + fmt(p.Rs, 3) + ') = ' + fmt(r.Ea, 3) + ' V');
        } else {
          W.push('Shunt field across the armature: Va = V − IL·Rs, If = Va/Rf, IL = Ia + If  →  If = (V − Ia·Rs)/(Rf + Rs) = (' + fmt(p.V, 2) + ' − ' + fmt(x, 2) + '·' + fmt(p.Rs, 3) + ')/(' + fmt(p.Rf, 2) + ' + ' + fmt(p.Rs, 3) + ') = ' + fmt(r.If, 4) + ' A');
          W.push('IL = ' + fmt(r.IL, 3) + ' A = Is (series field carries the line current);  Va = V − IL·Rs = ' + fmt(r.Va, 3) + ' V');
          W.push('Ea = Va − Vb − Ia·Ra = ' + fmt(r.Va, 3) + ' − ' + fmt(r.Vb, 2) + ' − ' + fmt(x, 2) + '·' + fmt(p.Ra, 3) + ' = ' + fmt(r.Ea, 3) + ' V');
        }
        W.push('Net field MMF (shunt-ampere units): I_fe = If ' + cs + ' (Ns/Nf)·Is = ' + fmt(r.If, 4) + ' ' + cs + ' (' + fmt(p.Ns, 1) + '/' + fmt(p.Nf, 0) + ')·' + fmt(r.Is, 3) + ' = ' + fmt(r.Ife, 4) + ' A');
        W.push('K = Laf·I_fe = ' + fmt(p.Laf, 3) + '·' + fmt(r.Ife, 4) + ' = ' + fmt(r.K, 4) + ' V·s/rad;   ω = Ea/K = ' + (fin(r.w) ? fmt(r.w, 3) + ' rad/s = ' + fmt(r.rpm, 1) + ' rpm' : 'undefined (K ≤ 0)') + ';   Tem = K·Ia = ' + fmt(r.Tem, 3) + ' N·m');
        W.push('No load (Ia = 0): ω0 = ' + (fin(r0.w) ? fmt(r0.w, 3) + ' rad/s = ' + fmt(r0.rpm, 1) + ' rpm' : '—') + '.  Plain shunt (series winding removed): Ea = ' + fmt(rs.Ea, 3) + ' V, K = Laf·If = ' + fmt(rs.K, 4) + ' → ' + fmt(rs.rpm, 1) + ' rpm.');
        W.push('Cumulative: series MMF adds, flux rises with load, so speed falls MORE than a shunt motor (and torque is larger for the same Ia). Differential: series MMF opposes, flux falls with load, so speed falls LESS, can even rise, and is unstable if I_fe → 0.');
        W.push('Comparison curves: shunt = same machine with Ns = 0, Rs = 0.  Series = same armature (Ra, Rs), no shunt field, series turns chosen so I_fe at rated current equals the shunt no-load If (' + fmt(serP.Ns, 2) + ' turns/pole equivalent).  Curves are cut at 2.4× the shunt no-load speed (series runaway).');
        ui.work.set(W);
      } else {
        // ----- generator
        const ifNL = CD.generator(Object.assign({}, V.shunt, { conn: 'short' }), 0);
        const bad = !ifNL.ok ? ifNL.msg : null;
        if (bad) { ui.warn.hidden = false; ui.warn.textContent = bad; ui.work.set(bad); hudSet(ui.hud, [['status', 'fix inputs']]); drawPlot(ui.c1, { xmin: 0, xmax: 1, ymin: 0, ymax: 1 }); return; }
        const Ifnl = ifNL.collapsed ? 0.2 : ifNL.If;
        serP = Object.assign({}, p, { Rf: Infinity, Ns: p.Nf * Ifnl / Ir, conn: 'long', comp: 'cum' });
        const imax = 1.6 * Ir;
        const sw = pp => { const X = [], Y = []; for (let i = 0; i <= 100; i++) { const il = imax * i / 100, g = CD.generator(pp, il); if (!g.ok || g.collapsed) break; X.push(il); Y.push(g.Vt); } return { x: X, y: Y }; };
        const defs = [['series', serP, 'series generator'], ['shunt', V.shunt, 'shunt generator'], [sel === 'cum' ? 'diff' : 'cum', sel === 'cum' ? V.diff : V.cum, sel === 'cum' ? 'differential compound' : 'cumulative compound'], [sel, V[sel], (sel === 'cum' ? 'cumulative' : 'differential') + ' compound, ' + p.conn + ' shunt']];
        defs.forEach(d => { const c = sw(d[1]); series.push({ x: c.x, y: c.y, color: cols[d[0]], width: d[0] === sel ? 2.8 : 1.7, dash: d[0] === sel ? [] : (d[0] === 'series' || d[0] === 'shunt' ? [5, 3] : []) }); lines.push({ text: d[2], color: cols[d[0]] }); });
        const g = CD.generator(V[sel], x), g0 = CD.generator(V[sel], 0), gs = CD.generator(V.shunt, x), g0s = CD.generator(V.shunt, 0);
        if (g.ok && !g.collapsed) points.push({ x, y: g.Vt, label: fmt(g.Vt, 1) + ' V', color: T.text });
        const vmax = Math.max.apply(null, series.reduce((a, s) => a.concat(s.y), [1]));
        drawPlot(ui.c1, { xmin: 0, xmax: imax, ymin: 0, ymax: vmax * 1.1, xlabel: 'load current IL (A)', ylabel: 'terminal voltage Vt (V)', series, points, legend: lines.reverse() });
        if (g.collapsed) warns.push('Voltage collapse: at this load there is no self-consistent terminal voltage (the generator has stalled out). Reduce the load or raise the speed.');
        if (g0.collapsed) warns.push('No voltage build-up at no load: Laf·ω must exceed the shunt-circuit resistance (critical resistance) for self-excitation.');
        const reg = (a, b) => (a > 0 ? fmt(100 * (a - b) / a, 1) + ' %' : '—');
        hudSet(ui.hud, [
          ['load current IL', fmt(x, 2) + ' A'], ['terminal voltage Vt', g.collapsed ? '0 (collapsed)' : fmt(g.Vt, 2) + ' V'], ['shunt field If', fmt(g.If, 3) + ' A'], ['armature current Ia', fmt(g.Ia, 2) + ' A'],
          ['series current Is', fmt(g.Is, 2) + ' A'], ['field MMF I_fe', fmt(g.Ife, 4) + ' A'], ['K', fmt(g.K, 4) + ' V·s/rad'], ['generated EMF E = Kω', fmt(g.Ea, 2) + ' V'],
          ['no-load voltage (this)', fmt(g0.Vt, 2) + ' V'], ['voltage drop NL → load', g0.Vt > 0 ? reg(g0.Vt, g.Vt) + (g.Vt > g0.Vt ? ' (voltage RISES)' : '') : '—'],
          ['plain shunt: NL / load', fmt(g0s.Vt, 1) + ' V / ' + (gs.collapsed ? 'collapsed' : fmt(gs.Vt, 1) + ' V')], ['load power Vt·IL', fmt(g.Vt * x / 1000, 3) + ' kW'],
        ]);
        const cs = sel === 'cum' ? '+' : '−', W = [];
        W.push((p.conn === 'long' ? 'Long shunt' : 'Short shunt') + ', ' + (sel === 'cum' ? 'cumulative' : 'differential') + ' compound generator at ω = ' + fmt(p.w, 3) + ' rad/s (' + fmt(st.N.get(), 0) + ' rpm), IL = ' + fmt(x, 2) + ' A.');
        if (p.conn === 'long') {
          W.push('Shunt field across the terminals: If = Vt/Rf;  Ia = IL + If;  series field in the armature path: Is = Ia.');
          W.push('KVL: E = Vt + Ia(Ra + Rs) + Vb');
        } else {
          W.push('Shunt field across the armature: Va = Vt + IL·Rs;  If = Va/Rf;  Ia = IL + If;  Is = IL.');
          W.push('KVL: E = Va + Ia·Ra + Vb  = Vt + IL·Rs + Ia·Ra + Vb');
        }
        W.push('I_fe = If ' + cs + ' (Ns/Nf)·Is;   E = ω·K(I_fe) = ω·Laf·I_fe/(1 + I_fe/Isat)   (Laf = ' + fmt(p.Laf, 3) + ' H, Isat = ' + fmt(p.Isat, 2) + ' A, Ns/Nf = ' + fmt(p.Ns, 1) + '/' + fmt(p.Nf, 0) + ')');
        W.push('Both expressions for E must agree: solve for the largest Vt (bisection) — the stable upper intersection of the field-resistance line with the saturating magnetisation curve.');
        if (!g.collapsed) {
          W.push('Result: Vt = ' + fmt(g.Vt, 3) + ' V;  If = ' + fmt(g.If, 4) + ' A;  Ia = ' + fmt(g.Ia, 3) + ' A;  Is = ' + fmt(g.Is, 3) + ' A;  I_fe = ' + fmt(g.Ife, 4) + ' A;  K = ' + fmt(g.K, 4) + ';  E = ' + fmt(g.Ea, 3) + ' V');
          W.push('Check: E required = ' + fmt(g.Ereq, 3) + ' V, E produced = ' + fmt(g.Ea, 3) + ' V.');
        } else W.push('No solution for this load: the voltage collapses.');
        W.push('No-load voltage ' + fmt(g0.Vt, 2) + ' V (' + (p.conn === 'long' ? 'long shunt carries If in the series field at no load' : 'short shunt: series field carries IL = 0') + ').  Cumulative compounding keeps Vt up (the series MMF compensates the armature drop and can over-compound); differential makes it fall sharply (welding sets).');
        W.push('Comparison: shunt = Ns = 0, Rs = 0 (drooping curve, collapses at high load); series = no shunt field, series turns sized so I_fe at rated current equals the no-load If of the shunt machine; Vt rises with IL, then falls.');
        ui.work.set(W);
      }
      ui.warn.hidden = warns.length === 0; ui.warn.textContent = warns.join(' ');
    }
    return {
      update,
      get: () => ({ c_mode: ui.mode.value, c_conn: ui.conn.value, c_comp: ui.comp.value, c_V: st.V.get(), c_N: st.N.get(), c_Ra: st.Ra.get(), c_Rf: st.Rf.get(), c_Rs: st.Rs.get(), c_Nf: st.Nf.get(), c_Ns: st.Ns.get(), c_Laf: st.Laf.get(), c_Isat: st.Isat.get(), c_Vb: st.Vb.get(), c_Ir: st.Ir.get(), c_ld: st.ld.get() }),
      set(o) {
        ui.mode.value = pick(o.c_mode, ['motor', 'gen'], 'motor'); ui.conn.value = pick(o.c_conn, ['long', 'short'], 'long'); ui.comp.value = pick(o.c_comp, ['cum', 'diff'], 'cum');
        Object.keys(st).forEach(k => { if (o['c_' + k] !== undefined) { const v = parseFloat(o['c_' + k]); if (fin(v)) st[k].set(v, true); } });
      },
    };
  }

  /* ---------- 8c panel ---------- */
  function buildInduction(root) {
    const st = {}, ui = {}; let view = 'analysis';
    const vbar = viewBar(root, [['analysis', 'Torque–slip analysis'], ['speed', 'Speed control']], 'Induction motor view', v => setView(v));
    const lay = FSP.ui.el('div', { class: 'layout' }); root.appendChild(lay);
    const lay2 = FSP.ui.el('div', { class: 'layout', hidden: '' }); root.appendChild(lay2);
    const sc = buildSpeedControl(lay2, () => params());
    function setView(v, silent) { view = v; vbar.set(v); lay.hidden = v !== 'analysis'; lay2.hidden = v !== 'speed'; updateAll(); if (!silent) FSP.state.touch(); }
    function updateAll() { if (view === 'analysis') update(); else sc.update(); }
    vbar.set(view);
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
      update: updateAll,
      get: () => Object.assign({ view, VL: val(st, 'VL'), f: val(st, 'f'), conn: ui.conn.value, poles: ui.poles.value, R1: val(st, 'R1'), X1: val(st, 'X1'), R2: val(st, 'R2'), X2: val(st, 'X2'), Xm: val(st, 'Xm'), Prot: val(st, 'Prot'), sr: val(st, 'sr') }, sc.get()),
      set(o) {
        sc.set(o); setView(pick(o.view, ['analysis', 'speed'], 'analysis'), true);
        Object.keys(st).forEach(k => { if (o[k] !== undefined) { const v = parseFloat(o[k]); if (fin(v)) st[k].set(v, true); } });
        if (o.conn === 'Y' || o.conn === 'D') ui.conn.value = o.conn;
        if (['2', '4', '6', '8', '10', '12'].indexOf(o.poles) >= 0) ui.poles.value = o.poles;
      },
    };
  }

  /* ---------- 8d panel ---------- */
  function buildDC(root) {
    const st = {}, ui = {}; let view = 'single';
    const vbar = viewBar(root, [['single', 'Shunt / series / separate'], ['compound', 'Compound machine']], 'DC machine view', v => setView(v));
    const lay = FSP.ui.el('div', { class: 'layout' }); root.appendChild(lay);
    const lay2 = FSP.ui.el('div', { class: 'layout', hidden: '' }); root.appendChild(lay2);
    const cc = buildCompound(lay2);
    function setView(v, silent) { view = v; vbar.set(v); lay.hidden = v !== 'single'; lay2.hidden = v !== 'compound'; updateAll(); if (!silent) FSP.state.touch(); }
    function updateAll() { if (view === 'single') update(); else cc.update(); }
    vbar.set(view);
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
      update: updateAll,
      get: () => Object.assign({ view, type: ui.type.value, V: val(st, 'V'), Vf: val(st, 'Vf'), Laf: val(st, 'Laf'), Ra: val(st, 'Ra'), Rf: val(st, 'Rf'), Rs: val(st, 'Rs'), Vb: val(st, 'Vb'), Prot: val(st, 'Prot'), Irated: val(st, 'Irated'), load: val(st, 'load') }, cc.get()),
      set(o) {
        cc.set(o); setView(pick(o.view, ['single', 'compound'], 'single'), true);
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
