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

  FSP.math.rotating._shared = { CD, DC, IM, PI, S3, SC, TP, cabs, cx, fin };
})();
