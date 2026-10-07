/* Magnetics & Transformer module (Phase 8a + 8b).
   Pure math lives in FSP.math.magnetics (no DOM); UI is built only from init()/activate(). */
(function () {
'use strict';
const MU0 = 1.25663706212e-6;
const M = {};
FSP.math.magnetics = M;
M.MU0 = MU0;

/* =============================================================
   Part A: magnetic circuit math
   ============================================================= */

// Monotone cubic Hermite (Fritsch-Carlson / PCHIP). xs strictly increasing.
function pchip(xs, ys) {
  const n = xs.length, h = [], d = [], m = new Array(n);
  for (let i = 0; i < n - 1; i++) { h[i] = xs[i + 1] - xs[i]; d[i] = (ys[i + 1] - ys[i]) / h[i]; }
  function endSlope(h0, h1, d0, d1) {
    let s = ((2 * h0 + h1) * d0 - h0 * d1) / (h0 + h1);
    if (Math.sign(s) !== Math.sign(d0)) s = 0;
    else if (Math.sign(d0) !== Math.sign(d1) && Math.abs(s) > 3 * Math.abs(d0)) s = 3 * d0;
    return s;
  }
  if (n === 2) { m[0] = m[1] = d[0]; }
  else {
    for (let i = 1; i < n - 1; i++) {
      if (d[i - 1] * d[i] <= 0) m[i] = 0;
      else { const w1 = 2 * h[i] + h[i - 1], w2 = h[i] + 2 * h[i - 1]; m[i] = (w1 + w2) / (w1 / d[i - 1] + w2 / d[i]); }
    }
    m[0] = endSlope(h[0], h[1], d[0], d[1]);
    m[n - 1] = endSlope(h[n - 2], h[n - 3], d[n - 2], d[n - 3]);
  }
  function seg(x) { let lo = 0, hi = n - 2; while (lo < hi) { const mid = (lo + hi + 1) >> 1; if (xs[mid] <= x) lo = mid; else hi = mid - 1; } return lo; }
  return {
    eval(x) {
      const i = seg(x), t = (x - xs[i]) / h[i], t2 = t * t, t3 = t2 * t;
      return (2 * t3 - 3 * t2 + 1) * ys[i] + (t3 - 2 * t2 + t) * h[i] * m[i] + (-2 * t3 + 3 * t2) * ys[i + 1] + (t3 - t2) * h[i] * m[i + 1];
    },
    deriv(x) {
      const i = seg(x), t = (x - xs[i]) / h[i], t2 = t * t;
      return (6 * t2 - 6 * t) / h[i] * ys[i] + (3 * t2 - 4 * t + 1) * m[i] + (-6 * t2 + 6 * t) / h[i] * ys[i + 1] + (3 * t2 - 2 * t) * m[i + 1];
    },
  };
}
M.pchip = pchip;

// Default M-19-like silicon steel (H in A/m, B in T); leading (0,0) node is required.
M.DEFAULT_TABLE = {
  H: [0, 10, 20, 40, 60, 80, 100, 150, 200, 300, 500, 1000, 2000, 5000, 10000, 20000],
  B: [0, 0.10, 0.35, 0.80, 1.05, 1.20, 1.28, 1.38, 1.43, 1.50, 1.58, 1.70, 1.80, 1.95, 2.05, 2.15],
};

// Parse "H B" rows (whitespace/comma/semicolon separated). Returns {ok, H, B, msg}.
M.parseTable = function (text) {
  const rows = [];
  String(text == null ? '' : text).split(/\r?\n|;/).forEach(line => {
    line = line.replace(/#.*$/, '').trim(); if (!line) return;
    const p = line.split(/[\s,]+/).map(Number);
    if (p.length >= 2 && Number.isFinite(p[0]) && Number.isFinite(p[1])) rows.push([p[0], p[1]]);
  });
  if (rows.length && !(rows[0][0] === 0 && rows[0][1] === 0)) rows.unshift([0, 0]);
  rows.sort((a, b) => a[0] - b[0]);
  if (rows.length < 3) return { ok: false, msg: 'Need at least two (H, B) points with H > 0.' };
  for (let i = 0; i < rows.length; i++) {
    if (rows[i][0] < 0 || rows[i][1] < 0) return { ok: false, msg: 'H and B must be non-negative (first quadrant only).' };
    if (i > 0 && !(rows[i][0] > rows[i - 1][0])) return { ok: false, msg: 'H values must be strictly increasing (duplicate H at ' + rows[i][0] + ').' };
    if (i > 0 && !(rows[i][1] > rows[i - 1][1])) return { ok: false, msg: 'B must strictly increase with H (check row H = ' + rows[i][0] + ').' };
  }
  return { ok: true, H: rows.map(r => r[0]), B: rows.map(r => r[1]), msg: '' };
};

// Material: constant mu_r or B-H table. All functions are odd in their argument.
// Beyond the last table point the curve continues with slope mu0 (fully saturated).
M.makeMaterial = function (o) {
  if (!o || o.mode !== 'table') {
    const mur = o && o.mur > 0 ? o.mur : 1, mu = MU0 * mur;
    return { mode: 'const', mur, muRef: mu, Hof: B => B / mu, Bof: H => H * mu, dHdB: () => 1 / mu, Bsat: Infinity };
  }
  const Hs = o.H, Bs = o.B, n = Hs.length, Bl = Bs[n - 1], Hl = Hs[n - 1];
  const f = pchip(Bs, Hs);
  let muRef = 0; for (let i = 1; i < n; i++) muRef = Math.max(muRef, Bs[i] / Hs[i]);
  return {
    mode: 'table', muRef, Bsat: Bl, H: Hs, B: Bs,
    Hof(B) { const a = Math.abs(B), s = B < 0 ? -1 : 1; return s * (a <= Bl ? f.eval(a) : Hl + (a - Bl) / MU0); },
    dHdB(B) { const a = Math.abs(B); return a <= Bl ? Math.max(f.deriv(a), 1 / (1e4 * muRef)) : 1 / MU0; },
    Bof(H) {
      const a = Math.abs(H), s = H < 0 ? -1 : 1;
      if (a >= Hl) return s * (Bl + (a - Hl) * MU0);
      let lo = 0, hi = Bl;
      for (let i = 0; i < 100; i++) { const mid = 0.5 * (lo + hi); if (f.eval(mid) < a) lo = mid; else hi = mid; if (hi - lo <= 1e-16 * Bl) break; }
      return s * 0.5 * (lo + hi);
    },
  };
};

// c: {N, lc, Ac, lg, fringe, mat, par: null | {l2, A2, l3, A3}} (SI units)
M.makeCircuit = function (p) {
  return { N: p.N, lc: p.lc, Ac: p.Ac, lg: p.lg, Ag: p.Ac * (p.fringe > 0 ? p.fringe : 1), mat: p.mat, par: p.par || null };
};

// Two parallel core legs (after the winding leg+gap) sharing flux phi: find common mmf drop F.
function parSolve(c, phi) {
  const p = c.par, mat = c.mat;
  const leg = (F, A, l) => A * mat.Bof(F / l);
  const flux = F => leg(F, p.A2, p.l2) + leg(F, p.A3, p.l3);
  if (!(phi > 0)) return { F: 0, phi2: 0, phi3: 0, B2: 0, B3: 0, H2: 0, H3: 0 };
  let hi = 1, k = 0; while (flux(hi) < phi && k++ < 400) hi *= 2;
  let lo = 0;
  for (let i = 0; i < 200; i++) { const mid = 0.5 * (lo + hi); if (flux(mid) < phi) lo = mid; else hi = mid; if (hi - lo <= 1e-16 * hi) break; }
  const F = 0.5 * (lo + hi), H2 = F / p.l2, H3 = F / p.l3, B2 = mat.Bof(H2), B3 = mat.Bof(H3);
  return { F, H2, H3, B2, B3, phi2: p.A2 * B2, phi3: p.A3 * B3 };
}

// Forward evaluation: total NI required for flux phi (>= 0).
M.mmf = function (c, phi) {
  phi = Math.max(0, phi);
  const Bc = phi / c.Ac, Hc = c.mat.Hof(Bc), Fc = Hc * c.lc;
  const Bg = phi / c.Ag, Hg = Bg / MU0, Fg = Hg * c.lg;
  let par = null, Fp = 0;
  if (c.par) { par = parSolve(c, phi); Fp = par.F; }
  return { NI: Fc + Fg + Fp, Bc, Hc, Fc, Bg, Hg, Fg, Fp, par };
};

// Inverse: flux for a given NI by bracketed bisection. Reports convergence.
M.solvePhi = function (c, NI) {
  if (!(NI > 0)) return { phi: 0, iter: 0, bracket: 0, residual: Math.abs(NI) || 0, relResidual: 0, converged: NI === 0 || !Number.isFinite(NI) ? NI === 0 : false };
  let hi = 1e-9, k = 0;
  while (M.mmf(c, hi).NI < NI && k < 400) { hi *= 2; k++; }
  let lo = 0, it = 0;
  for (; it < 300; it++) {
    const mid = 0.5 * (lo + hi);
    if (M.mmf(c, mid).NI < NI) lo = mid; else hi = mid;
    if (hi - lo <= 4e-16 * hi) { it++; break; }
  }
  const phi = 0.5 * (lo + hi), residual = Math.abs(M.mmf(c, phi).NI - NI), rel = residual / NI;
  return { phi, iter: it, bracket: k, residual, relResidual: rel, converged: rel <= 1e-10 && k < 400 };
};

// Stored field energy W = integral_0^phi NI(phi') dphi' (Simpson, n even).
M.energy = function (c, phi, n) {
  n = n || 200; if (n % 2) n++;
  if (!(phi > 0)) return 0;
  const h = phi / n; let s = M.mmf(c, 0).NI + M.mmf(c, phi).NI;
  for (let i = 1; i < n; i++) s += (i % 2 ? 4 : 2) * M.mmf(c, i * h).NI;
  return s * h / 3;
};

// Full analysis. spec: {NI} or {phi}.
M.analyze = function (c, spec) {
  let sol = null, phi;
  if (spec.phi !== undefined) phi = Math.max(0, spec.phi);
  else { sol = M.solvePhi(c, spec.NI); phi = sol.phi; }
  const m = M.mmf(c, phi), NI = m.NI;
  const pr = Math.max(phi, 1e-10), mr = phi >= 1e-10 ? m : M.mmf(c, pr);
  const Rcore = mr.Fc / pr, Rgap = mr.Fg / pr, Rpar = mr.Fp / pr, Rtot = mr.NI / pr;
  const h = Math.max(1e-4 * phi, 1e-10), lo = Math.max(0, phi - h);
  const Rdiff = (M.mmf(c, phi + h).NI - M.mmf(c, lo).NI) / (phi + h - lo);
  const W = M.energy(c, phi);
  const Wgap = 0.5 * m.Bg * m.Bg / MU0 * c.Ag * c.lg;
  // saturation: differential permeability relative to the best secant permeability
  const legs = [{ name: 'core', B: m.Bc }];
  if (c.par) { legs.push({ name: 'leg 2', B: m.par.B2 }, { name: 'leg 3', B: m.par.B3 }); }
  let ratio = Infinity;
  legs.forEach(g => { g.muD = 1 / c.mat.dHdB(g.B); g.ratio = g.muD / c.mat.muRef; ratio = Math.min(ratio, g.ratio); });
  let sat = 'linear';
  if (c.mat.mode === 'table') sat = ratio < 0.03 ? 'saturated' : ratio < 0.3 ? 'knee' : 'linear';
  const N = c.N;
  return {
    sol, phi, NI, I: NI / N, lambda: N * phi, m,
    Rcore, Rgap, Rpar, Rtot, Rdiff, Lsec: N * N / Rtot, Linc: N * N / Rdiff,
    W, Wgap, ratio, sat, legs, murEff: m.Hc > 0 ? m.Bc / (MU0 * m.Hc) : (c.mat.mode === 'const' ? c.mat.mur : c.mat.muRef / MU0),
  };
};

// Closed-form linear total reluctance (constant mu_r only) for cross-checks.
M.linearReluctance = function (c) {
  const mu = c.mat.mur * MU0;
  let R = c.lc / (mu * c.Ac) + c.lg / (MU0 * c.Ag);
  if (c.par) { const R2 = c.par.l2 / (mu * c.par.A2), R3 = c.par.l3 / (mu * c.par.A3); R += R2 * R3 / (R2 + R3); }
  return R;
};

/* =============================================================
   Part B: single-phase transformer math (referred to primary)
   ============================================================= */

// s: {S (VA), V1, V2 rated, oc:{V,I,P,side:'pri'|'sec'}, sc:{V,I,P,side}}
M.xfmrFromTests = function (s) {
  const bad = msg => ({ ok: false, msg });
  const pos = v => Number.isFinite(v) && v > 0;
  if (![s.S, s.V1, s.V2, s.oc.V, s.oc.I, s.oc.P, s.sc.V, s.sc.I, s.sc.P].every(pos)) return bad('All ratings and test readings must be positive numbers.');
  const a = s.V1 / s.V2, refer = side => side === 'sec' ? a * a : 1;
  const Ic = s.oc.P / s.oc.V, Iphi2 = s.oc.I * s.oc.I - Ic * Ic;
  if (!(Iphi2 > 0)) return bad('Open-circuit data inconsistent: P_oc > V_oc·I_oc (power factor above 1).');
  const Zsc = s.sc.V / s.sc.I, Req0 = s.sc.P / (s.sc.I * s.sc.I), Xeq2 = Zsc * Zsc - Req0 * Req0;
  if (!(Xeq2 > 0)) return bad('Short-circuit data inconsistent: P_sc > V_sc·I_sc.');
  const Rc = s.oc.V * s.oc.V / s.oc.P * refer(s.oc.side), Xm = s.oc.V / Math.sqrt(Iphi2) * refer(s.oc.side);
  const Req = Req0 * refer(s.sc.side), Xeq = Math.sqrt(Xeq2) * refer(s.sc.side);
  const p = { ok: true, msg: '', S: s.S, V1r: s.V1, V2r: s.V2, a, Rc, Xm, Req, Xeq, R1: Req / 2, X1: Xeq / 2, R2p: Req / 2, X2p: Xeq / 2 };
  p.I1r = s.S / s.V1; p.I2r = s.S / s.V2;
  p.Pcore = s.V1 * s.V1 / Rc; p.PcuRated = Req * p.I1r * p.I1r;
  p.Zeq = Math.hypot(Req, Xeq);
  p.Rpu = Req * p.I1r / s.V1; p.Xpu = Xeq * p.I1r / s.V1;
  return p;
};

// Load impedance referred to primary for load fraction x, power factor pf, lead = true for leading current.
M.loadZ = function (p, x, pf, lead) {
  const mag = p.V1r * p.V1r / (x * p.S), th = Math.acos(Math.min(1, Math.max(0, pf))) * (lead ? -1 : 1);
  return C(mag * Math.cos(th), mag * Math.sin(th));
};

function finish(p, V1, I1, Vm, I2, V2) {
  const Ym = C(1 / p.Rc, -1 / p.Xm), Ie = Vm.mul(Ym);
  const Pin = reMulConj(V1, I1), Pout = reMulConj(V2, I2);
  const Pcu1 = I1.abs2() * p.R1, Pcu2 = I2.abs2() * p.R2p, Pcore = Vm.abs2() / p.Rc;
  return { V1, I1, Vm, I2, V2, Ie, Ic: Vm.scale(1 / p.Rc), Im: Vm.div(C(0, p.Xm)), Pin, Pout, Pcu1, Pcu2, Pcu: Pcu1 + Pcu2, Pcore, Qm: Vm.abs2() / p.Xm, eta: Pin > 0 ? Pout / Pin : NaN };
}

// Given primary voltage phasor V1 and load ZL (primary-referred, null = open circuit): exact T-model.
M.xfmrForward = function (p, V1, ZL) {
  const Z1 = C(p.R1, p.X1), Z2 = C(p.R2p, p.X2p), Ym = C(1 / p.Rc, -1 / p.Xm);
  let Y2 = C(0, 0), Zb2 = null;
  if (ZL) { Zb2 = Z2.add(ZL); Y2 = C(1, 0).div(Zb2); }
  const Zin = Z1.add(C(1, 0).div(Ym.add(Y2)));
  const I1 = V1.div(Zin), Vm = V1.sub(I1.mul(Z1));
  const I2 = ZL ? Vm.div(Zb2) : C(0, 0), V2 = ZL ? I2.mul(ZL) : Vm;
  return finish(p, V1, I1, Vm, I2, V2);
};

// Given secondary terminal voltage V2 (primary-referred phasor) and load ZL: back-calculate the source.
M.xfmrBack = function (p, V2, ZL) {
  const Z1 = C(p.R1, p.X1), Z2 = C(p.R2p, p.X2p), Ym = C(1 / p.Rc, -1 / p.Xm);
  const I2 = V2.div(ZL), Vm = V2.add(I2.mul(Z2)), Ie = Vm.mul(Ym), I1 = I2.add(Ie), V1 = Vm.add(I1.mul(Z1));
  return finish(p, V1, I1, Vm, I2, V2);
};

// Regulation at constant rated primary voltage: (|V2nl| - |V2|)/|V2nl|.
M.regulation = function (p, x, pf, lead) {
  const V1 = C(p.V1r, 0), nl = M.xfmrForward(p, V1, null), fl = M.xfmrForward(p, V1, M.loadZ(p, x, pf, lead));
  const a0 = nl.V2.abs(), a1 = fl.V2.abs();
  return { VR: (a0 - a1) / a0, V2nl: a0 / p.a, V2fl: a1 / p.a, nl, fl };
};

M.effExact = function (p, x, pf, lead) { return M.xfmrForward(p, C(p.V1r, 0), M.loadZ(p, x, pf, lead)).eta; };
// Conventional constant-core-loss, I^2R-copper efficiency.
M.effApprox = function (p, x, pf) { const po = x * p.S * pf; return po / (po + p.Pcore + x * x * p.PcuRated); };
M.maxEff = function (p, pf) { const x = Math.sqrt(p.Pcore / p.PcuRated), po = x * p.S * pf; return { x, eta: po / (po + 2 * p.Pcore), Pcore: p.Pcore, Pcu: x * x * p.PcuRated }; };

/* =============================================================
   Part C: UI (all DOM work inside init/activate)
   ============================================================= */
const TXT = '#8b98ab', GRID = '#243047', COL = { core: '#4cc9f0', gap: '#f9a03f', par: '#f72585', ok: '#2ecc71', fg: '#e6edf3' };
const fin = Number.isFinite;
const fE = (x, u, d) => fin(x) ? FSP.fmtEng(x, u, d || 4) : '—';
const fN = (x, d) => fin(x) ? FSP.fmtNum(x, d || 5) : '—';
const fP = (x, d) => fin(x) ? (100 * x).toFixed(d === undefined ? 3 : d) + ' %' : '—';
const monoFont = () => { try { return getComputedStyle(document.body).getPropertyValue('--mono') || 'monospace'; } catch (e) { return 'monospace'; } };

function niceTicks(lo, hi, n) {
  const span = hi - lo || 1, raw = span / n, p = Math.pow(10, Math.floor(Math.log10(raw))), f = raw / p;
  const step = (f < 1.5 ? 1 : f < 3 ? 2 : f < 7 ? 5 : 10) * p, out = [];
  for (let v = Math.ceil(lo / step - 1e-9) * step; v <= hi + step * 1e-9; v += step) out.push(+v.toPrecision(12));
  return out;
}
function tickLabel(v) { if (v === 0) return '0'; const a = Math.abs(v); if (a >= 1e5 || a < 1e-3) return v.toExponential(0).replace('e+', 'e'); return String(+v.toPrecision(3)); }
function arrow(ctx, x0, y0, x1, y1, hs) {
  const a = Math.atan2(y1 - y0, x1 - x0); ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1, y1); ctx.stroke();
  ctx.beginPath(); ctx.moveTo(x1, y1); ctx.lineTo(x1 - hs * Math.cos(a - 0.4), y1 - hs * Math.sin(a - 0.4)); ctx.lineTo(x1 - hs * Math.cos(a + 0.4), y1 - hs * Math.sin(a + 0.4)); ctx.closePath(); ctx.fill();
}
// Generic line plot. sp: {x0,x1,y0,y1,xlabel,ylabel,series:[{xs,ys,color,dash,width}],vlines:[{x,color,label}],marks:[{x,y,color,label}],clipLines:[...]}
function drawPlot(ctx, w, h, sp) {
  const L = 46, R = 10, T = 12, B = 34, pw = w - L - R, ph = h - T - B;
  const px = x => L + (x - sp.x0) / (sp.x1 - sp.x0) * pw, py = y => T + ph - (y - sp.y0) / (sp.y1 - sp.y0) * ph;
  ctx.clearRect(0, 0, w, h); ctx.font = '11px ' + monoFont(); ctx.lineWidth = 1;
  ctx.strokeStyle = GRID; ctx.fillStyle = TXT; ctx.textAlign = 'center'; ctx.textBaseline = 'top';
  niceTicks(sp.x0, sp.x1, Math.max(3, Math.floor(pw / 70))).forEach(v => { const x = px(v); ctx.beginPath(); ctx.moveTo(x, T); ctx.lineTo(x, T + ph); ctx.stroke(); ctx.fillText(tickLabel(v), x, T + ph + 3); });
  ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
  niceTicks(sp.y0, sp.y1, 5).forEach(v => { const y = py(v); ctx.beginPath(); ctx.moveTo(L, y); ctx.lineTo(L + pw, y); ctx.stroke(); ctx.fillText(tickLabel(v), L - 4, y); });
  ctx.strokeStyle = '#3a4a68'; ctx.strokeRect(L, T, pw, ph);
  ctx.fillStyle = TXT; ctx.textAlign = 'center'; ctx.textBaseline = 'bottom'; ctx.fillText(sp.xlabel || '', L + pw / 2, h - 1);
  if (sp.ylabel) { ctx.save(); ctx.translate(11, T + ph / 2); ctx.rotate(-Math.PI / 2); ctx.textBaseline = 'top'; ctx.fillText(sp.ylabel, 0, -4); ctx.restore(); }
  ctx.save(); ctx.beginPath(); ctx.rect(L, T, pw, ph); ctx.clip();
  (sp.series || []).forEach(s => {
    ctx.strokeStyle = s.color; ctx.lineWidth = s.width || 2; ctx.setLineDash(s.dash || []); ctx.beginPath(); let pen = false;
    for (let i = 0; i < s.xs.length; i++) {
      if (!fin(s.xs[i]) || !fin(s.ys[i])) { pen = false; continue; }
      const x = px(s.xs[i]), y = Math.max(-1e4, Math.min(1e4, py(s.ys[i]))); if (pen) ctx.lineTo(x, y); else { ctx.moveTo(x, y); pen = true; }
    }
    ctx.stroke(); ctx.setLineDash([]);
  });
  (sp.vlines || []).forEach(v => { if (!fin(v.x)) return; const x = px(v.x); ctx.strokeStyle = v.color; ctx.lineWidth = 1; ctx.setLineDash([4, 4]); ctx.beginPath(); ctx.moveTo(x, T); ctx.lineTo(x, T + ph); ctx.stroke(); ctx.setLineDash([]); });
  ctx.restore();
  ctx.lineWidth = 1;
  (sp.vlines || []).forEach(v => { if (!fin(v.x) || !v.label) return; const x = px(v.x); ctx.fillStyle = v.color; ctx.textBaseline = 'top'; ctx.textAlign = x > L + pw * 0.6 ? 'right' : 'left'; ctx.fillText(v.label, x + (x > L + pw * 0.6 ? -4 : 4), T + 3 + (v.row || 0) * 13); });
  (sp.marks || []).forEach(k => {
    if (!fin(k.x) || !fin(k.y)) return; const x = px(k.x), y = py(k.y); if (x < L || x > L + pw || y < T || y > T + ph) return;
    ctx.fillStyle = k.color; ctx.strokeStyle = '#0b0f17'; ctx.beginPath(); ctx.arc(x, y, 4.5, 0, 2 * Math.PI); ctx.fill(); ctx.stroke();
    if (k.label) { ctx.textBaseline = 'bottom'; ctx.textAlign = x > L + pw * 0.6 ? 'right' : 'left'; ctx.fillText(k.label, x + (x > L + pw * 0.6 ? -7 : 7), y - 5); }
  });
  (sp.legend || []).forEach((g, i) => { ctx.fillStyle = g.color; ctx.textAlign = 'right'; ctx.textBaseline = 'bottom'; ctx.fillText(g.text, L + pw - 4, T + ph - 3 - i * 13); });
}

function renderHud(hud, rows) {
  while (hud.firstChild) hud.removeChild(hud.firstChild);
  rows.forEach(r => { const v = FSP.ui.el('span', { text: r[1] }); if (r[2]) v.className = 'badge ' + r[2]; hud.appendChild(FSP.ui.el('div', null, FSP.ui.el('span', { text: r[0] }), v)); });
}

/* ---------- state ---------- */
const CTL = {
  N: { v: 500, min: 1, max: 20000, log: true, int: true, label: 'N', unit: 'turns', tip: 'Winding turns' },
  lc: { v: 0.3, min: 0.01, max: 5, log: true, label: 'l core', unit: 'm', tip: 'Mean magnetic path length of the winding leg' },
  Ac: { v: 10, min: 0.1, max: 1000, log: true, label: 'A core', unit: 'cm²', tip: 'Core cross-section (1e-3 m² = 10 cm²)' },
  lg: { v: 1, min: 0, max: 20, step: 0.01, label: 'gap', unit: 'mm', tip: 'Air gap length (0 = no gap)' },
  fr: { v: 1, min: 1, max: 2, step: 0.01, label: 'fringe', unit: '×A', tip: 'Fringing factor: gap area = factor × core area' },
  mur: { v: 2000, min: 1, max: 200000, log: true, label: 'μr', unit: '', tip: 'Constant relative permeability' },
  I: { v: 1, min: 0.001, max: 200, log: true, label: 'I', unit: 'A', tip: 'Winding current (solve mode: current → flux)' },
  Bc: { v: 0.5, min: 0.001, max: 2.6, log: true, label: 'B core', unit: 'T', tip: 'Core flux density (solve mode: flux → current)' },
  l2: { v: 0.2, min: 0.01, max: 5, log: true, label: 'l leg 2', unit: 'm', tip: 'Parallel leg 2 mean length' },
  A2: { v: 5, min: 0.1, max: 1000, log: true, label: 'A leg 2', unit: 'cm²', tip: 'Parallel leg 2 area' },
  l3: { v: 0.2, min: 0.01, max: 5, log: true, label: 'l leg 3', unit: 'm', tip: 'Parallel leg 3 mean length' },
  A3: { v: 5, min: 0.1, max: 1000, log: true, label: 'A leg 3', unit: 'cm²', tip: 'Parallel leg 3 area' },
  kva: { v: 10, min: 0.01, max: 100000, log: true, label: 'S rated', unit: 'kVA', tip: 'Rated apparent power' },
  V1: { v: 2400, min: 1, max: 1e6, log: true, label: 'V1 rated', unit: 'V', tip: 'Rated primary voltage' },
  V2: { v: 240, min: 1, max: 1e6, log: true, label: 'V2 rated', unit: 'V', tip: 'Rated secondary voltage' },
  f: { v: 50, min: 10, max: 400, log: true, label: 'f', unit: 'Hz', tip: 'Frequency (informational)' },
  ocV: { v: 240, min: 0.1, max: 1e6, log: true, label: 'V_oc', unit: 'V', tip: 'Open-circuit test voltage (usually rated, on the LV side)' },
  ocI: { v: 1.5, min: 1e-4, max: 1e4, log: true, label: 'I_oc', unit: 'A', tip: 'Open-circuit (no-load) current' },
  ocP: { v: 120, min: 0.01, max: 1e6, log: true, label: 'P_oc', unit: 'W', tip: 'Open-circuit power (core loss)' },
  scV: { v: 60, min: 0.01, max: 1e5, log: true, label: 'V_sc', unit: 'V', tip: 'Short-circuit test voltage (reduced, HV side)' },
  scI: { v: 4.1667, min: 1e-3, max: 1e4, log: true, label: 'I_sc', unit: 'A', tip: 'Short-circuit current (usually rated)' },
  scP: { v: 150, min: 0.01, max: 1e6, log: true, label: 'P_sc', unit: 'W', tip: 'Short-circuit power (copper loss)' },
  x: { v: 1, min: 0.01, max: 1.5, step: 0.01, label: 'load', unit: '× rated', tip: 'Load as fraction of rated kVA' },
  pf: { v: 0.8, min: 0.05, max: 1, step: 0.01, label: 'pf', unit: '', tip: 'Load power factor' },
};
const STR = { sub: 'a', mode: 'const', solve: 'I', par: '0', ocSide: 'sec', scSide: 'pri', lead: '0', table: '' };
const S = {}; Object.keys(CTL).forEach(k => { S[k] = CTL[k].v; }); Object.assign(S, STR);
S.table = tableToText(M.DEFAULT_TABLE);
function tableToText(t) { return t.H.map((h, i) => h + ' ' + t.B[i]).join('\n'); }

let inited = false, rootEl = null, active = false;
const widgets = []; // sync functions: push S into widgets
const draws = { a: [], b: [] };

/* ---------- circuit build from state ---------- */
let tableParsed = M.parseTable(S.table), tableMsg = '';
function currentCircuit() {
  const mat = S.mode === 'table' && tableParsed.ok ? M.makeMaterial({ mode: 'table', H: tableParsed.H, B: tableParsed.B }) : M.makeMaterial({ mode: 'const', mur: S.mur });
  const cm2 = 1e-4;
  return M.makeCircuit({
    N: S.N, lc: S.lc, Ac: S.Ac * cm2, lg: S.lg * 1e-3, fringe: S.fr, mat,
    par: S.par === '1' ? { l2: S.l2, A2: S.A2 * cm2, l3: S.l3, A3: S.A3 * cm2 } : null,
  });
}

let RA = null; // last magnetic analysis {c, r, working}
function computeA() {
  const c = currentCircuit();
  let r;
  if (S.solve === 'I') r = M.analyze(c, { NI: S.N * S.I });
  else r = M.analyze(c, { phi: S.Bc * c.Ac });
  RA = { c, r };
  return RA;
}

/* ---------- drawing: reluctance network ---------- */
function drawNetwork(cv) {
  if (!RA) return; const { c, r } = RA;
  const { ctx, w, h } = cv.prep(); if (w < 40) return;
  ctx.clearRect(0, 0, w, h); ctx.font = '11px ' + monoFont(); ctx.lineWidth = 2;
  const xL = 36, xR = w - 18, yT = 52, yB = h - 44, par = !!c.par, nar = w < 520;
  const fR = v => nar ? fE(v, '', 3) : fE(v, 'A/Wb', 3), fF = v => nar ? '' : 'F=' + fE(v, 'A', 3);
  const xP1 = par ? xR - Math.min(110, w * 0.3) : xR;
  const wire = (x0, y0, x1, y1, col) => { ctx.strokeStyle = col || '#6f819c'; ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1, y1); ctx.stroke(); };
  const boxH = (cx, y, name, val, col, mmf) => {
    ctx.fillStyle = '#0f141f'; ctx.strokeStyle = col; ctx.fillRect(cx - 28, y - 10, 56, 20); ctx.strokeRect(cx - 28, y - 10, 56, 20);
    ctx.fillStyle = col; ctx.textAlign = 'center'; ctx.textBaseline = 'bottom'; ctx.fillText(name, cx, y - 13);
    ctx.fillStyle = TXT; ctx.textBaseline = 'top'; ctx.fillText(val, cx, y + 13); ctx.fillText(mmf, cx, y + 25);
  };
  const boxV = (x, cy, name, col) => {
    ctx.fillStyle = '#0f141f'; ctx.strokeStyle = col; ctx.fillRect(x - 10, cy - 24, 20, 48); ctx.strokeRect(x - 10, cy - 24, 20, 48);
    ctx.fillStyle = col; ctx.textAlign = 'right'; ctx.textBaseline = 'middle'; ctx.fillText(name, x - 14, cy);
  };
  const ymid = (yT + yB) / 2;
  // loop wires
  const xs1 = xL + 30, xs2 = (par ? xP1 : xR) - 10, span = xs2 - xs1;
  const cx1 = xs1 + span * 0.27, cx2 = xs1 + span * 0.73;
  wire(xL, yT, cx1 - 28, yT); wire(cx1 + 28, yT, cx2 - 28, yT); wire(cx2 + 28, yT, par ? xR : xR, yT);
  wire(xL, yB, xR, yB); wire(xR, yT, xR, par ? ymid - 24 : yB); wire(xL, yT, xL, ymid - 16); wire(xL, ymid + 16, xL, yB);
  if (par) {
    wire(xP1, yT, xP1, ymid - 24); wire(xP1, ymid + 24, xP1, yB); wire(xR, ymid + 24, xR, yB);
    boxV(xP1, ymid, 'R₂', COL.par); boxV(xR, ymid, 'R₃', COL.par);
    ctx.fillStyle = TXT; ctx.textAlign = 'right'; ctx.textBaseline = 'middle'; ctx.fillText('R₂‖R₃ = ' + fR(r.Rpar), xP1 - 14, ymid + 16);
  }
  // MMF source
  ctx.strokeStyle = COL.fg; ctx.fillStyle = '#0f141f'; ctx.beginPath(); ctx.arc(xL, ymid, 16, 0, 2 * Math.PI); ctx.fill(); ctx.stroke();
  ctx.fillStyle = COL.fg; ctx.textAlign = 'center'; ctx.textBaseline = 'middle'; ctx.fillText('NI', xL, ymid);
  ctx.textBaseline = 'top'; ctx.fillStyle = TXT; ctx.textAlign = 'left'; ctx.fillText(fE(r.NI, 'A·t', 3), xL + 20, ymid - 6);
  boxH(cx1, yT, 'R_core', fR(r.Rcore), COL.core, fF(r.m.Fc));
  boxH(cx2, yT, 'R_gap', fR(r.Rgap), COL.gap, fF(r.m.Fg));
  // flux arrow on bottom rail
  ctx.strokeStyle = COL.ok; ctx.fillStyle = COL.ok; ctx.lineWidth = 2; arrow(ctx, xL + 30, yB, xL + 80, yB, 8);
  ctx.textAlign = 'left'; ctx.textBaseline = 'top'; ctx.fillText('φ = ' + fE(r.phi, 'Wb', 3) + (nar ? '' : '   R_tot = ' + fE(r.Rtot, 'A/Wb', 3)), xL + 88 > w - 150 ? xL : xL + 88, yB + 5);
  ctx.lineWidth = 1;
}

/* ---------- drawing: B-H ---------- */
function drawBH(cv) {
  if (!RA) return; const { c, r } = RA; const mat = c.mat;
  const { ctx, w, h } = cv.prep(); if (w < 40) return;
  const legsH = [r.m.Hc]; if (c.par) legsH.push(r.m.par.H2, r.m.par.H3);
  let x1 = Math.max(10, 1.6 * Math.max.apply(null, legsH)); if (mat.mode === 'table') x1 = Math.max(x1, 2 * mat.H[3] || 0);
  const n = 160, xs = [], ys = [];
  for (let i = 0; i <= n; i++) { const H = x1 * i / n; xs.push(H); ys.push(mat.Bof(H)); }
  const Bops = [r.m.Bc]; if (c.par) Bops.push(r.m.par.B2, r.m.par.B3);
  const y1 = Math.max(0.1, ys[n], 0) * 1.05;
  const series = [{ xs, ys, color: COL.core, width: 2.5 }];
  if (!c.par && c.lg > 0 && r.NI > 0) { // load line: Hc = NI/lc - k*Bc
    const k = c.lg * c.Ac / (c.lc * c.Ag * MU0), Hint = r.NI / c.lc;
    series.push({ xs: [0, Hint], ys: [Hint / k, 0], color: COL.gap, dash: [6, 4], width: 1.5 });
  }
  const marks = [{ x: r.m.Hc, y: r.m.Bc, color: COL.ok, label: 'core' }];
  if (c.par) { marks.push({ x: r.m.par.H2, y: r.m.par.B2, color: COL.par, label: 'leg2' }, { x: r.m.par.H3, y: r.m.par.B3, color: COL.par, label: 'leg3' }); }
  const legend = [{ color: COL.core, text: mat.mode === 'table' ? 'B–H table' : 'B = μ0 μr H' }];
  if (series.length > 1) legend.push({ color: COL.gap, text: 'gap load line' });
  drawPlot(ctx, w, h, { x0: 0, x1, y0: 0, y1, xlabel: 'H (A/m)', ylabel: 'B (T)', series, marks, legend });
}

/* ---------- working text (A) ---------- */
function workingA() {
  const { c, r } = RA, m = r.m;
  const L = [];
  L.push('Series path: winding-leg core (l=' + fN(c.lc) + ' m, A=' + fN(c.Ac) + ' m²) + air gap (l_g=' + fN(c.lg) + ' m, A_g=' + fN(c.Ag) + ' m²)' + (c.par ? ' + parallel legs R2‖R3' : ''));
  L.push('Ampère: NI = H_c·l_c + B_g·l_g/μ0' + (c.par ? ' + F_par' : '') + '   (Hopkinson: F = φ·R)');
  if (c.mat.mode === 'const') {
    L.push('R_core = l_c/(μ0 μr A_c) = ' + fN(c.lc) + '/(' + fN(MU0) + '·' + fN(c.mat.mur) + '·' + fN(c.Ac) + ') = ' + fE(r.Rcore, 'A/Wb', 6));
  } else {
    L.push('Core: H_c = H(B_c) by monotone cubic (PCHIP) interpolation of the table; B_c = φ/A_c = ' + fN(m.Bc) + ' T → H_c = ' + fN(m.Hc) + ' A/m; secant R_core = H_c·l_c/φ = ' + fE(r.Rcore, 'A/Wb', 6));
  }
  L.push('R_gap = l_g/(μ0 A_g) = ' + fE(r.Rgap, 'A/Wb', 6));
  if (c.par) L.push('Parallel legs share mmf drop F_par = ' + fN(m.Fp) + ' A; φ2 = ' + fE(m.par.phi2, 'Wb', 4) + ', φ3 = ' + fE(m.par.phi3, 'Wb', 4) + ' (sum = φ); R_par = F_par/φ = ' + fE(r.Rpar, 'A/Wb', 6));
  L.push('R_total = NI/φ = ' + fE(r.Rtot, 'A/Wb', 6) + '   →   L = N²/R_total = ' + fN(c.N) + '²/' + fN(r.Rtot) + ' = ' + fE(r.Lsec, 'H', 6) + ' (secant)');
  L.push('Incremental L_inc = N²/(dNI/dφ) = ' + fE(r.Linc, 'H', 6) + (c.mat.mode === 'table' ? '   (falls below the secant value in saturation)' : ''));
  L.push('Energy W = ∫ NI dφ = ' + fE(r.W, 'J', 6) + '  (gap: ' + fE(r.Wgap, 'J', 5) + ' = B_g²/(2μ0)·A_g·l_g);   linear check ½LI² = ' + fE(0.5 * r.Lsec * r.I * r.I, 'J', 6));
  if (r.sol) L.push('Solver: bracketed bisection on φ, ' + r.sol.iter + ' iterations (bracket doublings ' + r.sol.bracket + '), residual |NI(φ) − NI| = ' + r.sol.residual.toExponential(2) + ' A·t (relative ' + r.sol.relResidual.toExponential(2) + '), converged = ' + r.sol.converged);
  else L.push('Direct mode: φ = B·A_c, NI evaluated in closed form (no iteration).');
  L.push('Saturation: μ_diff/μ_secant,max = ' + fN(r.ratio, 3) + ' → ' + r.sat);
  return L;
}

function updateA(U) {
  if (!inited) return;
  computeA();
  const { c, r } = RA, m = r.m;
  const satCls = r.sat === 'saturated' ? 'bad' : r.sat === 'knee' ? 'warn' : 'ok';
  const rows = [
    ['NI (mmf)', fE(r.NI, 'A·t', 5)], ['Current I', fE(r.I, 'A', 5)], ['Flux φ', fE(r.phi, 'Wb', 5)], ['Flux linkage Nφ', fE(r.lambda, 'Wb·t', 5)],
    ['B core', fN(m.Bc, 5) + ' T'], ['H core', fE(m.Hc, 'A/m', 4)], ['B gap', fN(m.Bg, 5) + ' T'], ['H gap', fE(m.Hg, 'A/m', 4)],
    ['F core / F gap', fE(m.Fc, 'A', 4) + ' / ' + fE(m.Fg, 'A', 4)],
    ['R core', fE(r.Rcore, 'A/Wb', 5)], ['R gap', fE(r.Rgap, 'A/Wb', 5)], ['R total', fE(r.Rtot, 'A/Wb', 5)],
    ['L = N²/R (secant)', fE(r.Lsec, 'H', 5)], ['L incremental', fE(r.Linc, 'H', 5)],
    ['Stored energy W', fE(r.W, 'J', 5)], ['Gap share of W', fin(r.W) && r.W > 0 ? fP(r.Wgap / r.W, 2) : '—'],
    ['μr effective (core)', fN(r.murEff, 5)], ['Saturation', r.sat, satCls],
  ];
  if (c.par) { rows.splice(9, 0, ['φ2 / φ3', fE(m.par.phi2, 'Wb', 4) + ' / ' + fE(m.par.phi3, 'Wb', 4)], ['R par (R2‖R3)', fE(r.Rpar, 'A/Wb', 5)]); }
  if (r.sol) rows.push(['Solver', (r.sol.converged ? 'converged' : 'NOT converged') + ' · ' + r.sol.iter + ' it · res ' + r.sol.relResidual.toExponential(1), r.sol.converged ? 'ok' : 'bad']);
  else rows.push(['Solver', 'direct (φ given)', 'ok']);
  renderHud(U.hud, rows);
  U.working.set(workingA());
  U.warn.hidden = true;
  if (S.mode === 'table' && !tableParsed.ok) { U.warn.hidden = false; U.warn.textContent = 'Table invalid, using constant μr instead: ' + tableParsed.msg; }
  else if (r.sat === 'saturated') { U.warn.hidden = false; U.warn.textContent = 'Core is deep in saturation: the incremental permeability is under 3% of its peak, so L_inc << L_secant.'; }
  else if (r.sol && !r.sol.converged) { U.warn.hidden = false; U.warn.textContent = 'Solver did not converge to tolerance.'; }
  drawNetwork(U.net); drawBH(U.bh);
}

/* ---------- transformer compute/draw ---------- */
let RB = null;
function computeB() {
  const p = M.xfmrFromTests({ S: S.kva * 1000, V1: S.V1, V2: S.V2, oc: { V: S.ocV, I: S.ocI, P: S.ocP, side: S.ocSide }, sc: { V: S.scV, I: S.scI, P: S.scP, side: S.scSide } });
  RB = { p };
  if (!p.ok) return RB;
  const lead = S.lead === '1';
  RB.lead = lead;
  RB.ZL = M.loadZ(p, S.x, S.pf, lead);
  RB.reg = M.regulation(p, S.x, S.pf, lead);
  RB.fl = RB.reg.fl; RB.nl = RB.reg.nl;
  RB.me = M.maxEff(p, S.pf);
  RB.exactAtMax = M.effExact(p, RB.me.x, S.pf, lead);
  return RB;
}

function plotVR(cv) {
  if (!RB || !RB.p.ok) return; const p = RB.p;
  const { ctx, w, h } = cv.prep(); if (w < 40) return;
  const xs = [], lag = [], led = []; let lo = 0, hi = 0;
  for (let i = 0; i <= 98; i++) { const pf = 0.02 + i * 0.01; xs.push(pf); const a = M.regulation(p, S.x, pf, false).VR * 100, b = M.regulation(p, S.x, pf, true).VR * 100; lag.push(a); led.push(b); lo = Math.min(lo, a, b); hi = Math.max(hi, a, b); }
  const pad = (hi - lo || 1) * 0.08;
  const cur = RB.reg.VR * 100;
  drawPlot(ctx, w, h, { x0: 0, x1: 1, y0: lo - pad, y1: hi + pad, xlabel: 'power factor (load = ' + S.x.toFixed(2) + ' × rated)', ylabel: 'regulation (%)',
    series: [{ xs, ys: lag, color: COL.gap }, { xs, ys: led, color: COL.core }, { xs: [0, 1], ys: [0, 0], color: '#5a6a85', width: 1, dash: [3, 3] }],
    marks: [{ x: S.pf, y: cur, color: COL.ok, label: (RB.lead ? 'leading ' : 'lagging ') + cur.toFixed(2) + '%' }],
    legend: [{ color: COL.gap, text: 'lagging' }, { color: COL.core, text: 'leading' }] });
}
function plotEff(cv) {
  if (!RB || !RB.p.ok) return; const p = RB.p;
  const { ctx, w, h } = cv.prep(); if (w < 40) return;
  const xs = [], ex = [], ap = []; let lo = 100;
  for (let i = 0; i <= 150; i++) { const x = Math.max(0.01, i / 100); xs.push(x); const e = M.effExact(p, x, S.pf, RB.lead) * 100, a = M.effApprox(p, x, S.pf) * 100; ex.push(e); ap.push(a); if (x >= 0.1 && fin(e)) lo = Math.min(lo, e); }
  const y0 = Math.max(0, Math.floor((lo - 1) / 5) * 5), y1 = Math.min(100, Math.ceil(Math.max.apply(null, ex.concat(ap).filter(fin)) / 1 + 0.5));
  drawPlot(ctx, w, h, { x0: 0, x1: 1.5, y0, y1: Math.max(y1, y0 + 1), xlabel: 'load fraction (pf = ' + S.pf.toFixed(2) + (RB.lead ? ' leading)' : ' lagging)'), ylabel: 'efficiency (%)',
    series: [{ xs, ys: ap, color: COL.gap, dash: [6, 4], width: 1.5 }, { xs, ys: ex, color: COL.core }],
    vlines: [{ x: RB.me.x, color: COL.par, label: 'x* = ' + RB.me.x.toFixed(3) }],
    marks: [{ x: S.x, y: RB.fl.eta * 100, color: COL.ok, label: (RB.fl.eta * 100).toFixed(2) + '%' }, { x: RB.me.x, y: RB.me.eta * 100, color: COL.par }],
    legend: [{ color: COL.core, text: 'exact T-model' }, { color: COL.gap, text: 'P_core + x²P_cu approx.' }] });
}
function plotPhasor(cv) {
  if (!RB || !RB.p.ok) return; const p = RB.p, f = RB.fl;
  const { ctx, w, h } = cv.prep(); if (w < 40) return;
  ctx.clearRect(0, 0, w, h); ctx.font = '11px ' + monoFont();
  const rot = f.V2.arg() * -1, rc = z => z.mul(C(Math.cos(rot), Math.sin(rot)));
  const V2 = rc(f.V2), I2 = rc(f.I2), Ie = rc(f.Ie), E = rc(f.Vm), V1 = rc(f.V1), I1 = rc(f.I1);
  const Z1 = C(p.R1, p.X1), Z2 = C(p.R2p, p.X2p);
  const d2 = I2.mul(Z2), d1 = I1.mul(Z1);
  const mag = V2.abs();
  const maxDrop = Math.max(d1.abs(), d2.abs()) / mag, kd = Math.min(40, Math.max(1, 0.18 / (maxDrop || 1)));
  const ke = Ie.abs() > 0 && I2.abs() > 0 ? Math.min(40, Math.max(1, 0.12 * I2.abs() / Ie.abs())) : 1;
  // schematic (exaggerated) points
  const vE = V2.add(d2.scale(kd)), vV1 = vE.add(d1.scale(kd));
  const iE = Ie.scale(ke), i1 = I2.add(iE);
  const Lv = Math.min(0.62 * w, 0.9 * h), sV = Lv / mag, Li = 0.42 * Lv, iMax = Math.max(I2.abs(), i1.abs(), 1e-12), sI = Li / iMax;
  const ox = 14, oy = h / 2 + 4;
  const P = (z, s) => [ox + z.re * s, oy - z.im * s];
  const vec = (z, s, col, label, dash, from, lw, dy) => {
    const a = from ? P(from, s) : [ox, oy], b = P(from ? from.add(z) : z, s);
    ctx.strokeStyle = col; ctx.fillStyle = col; ctx.lineWidth = lw || 2; ctx.setLineDash(dash || []); arrow(ctx, a[0], a[1], b[0], b[1], 8); ctx.setLineDash([]);
    if (label) { ctx.textAlign = 'left'; ctx.textBaseline = 'middle'; const lx = Math.min(w - 4 - ctx.measureText(label).width, b[0] + 5); ctx.fillText(label, lx, b[1] + (dy !== undefined ? dy : (z.im >= 0 ? -8 : 8))); }
  };
  // axes hint
  ctx.strokeStyle = GRID; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(ox, 8); ctx.lineTo(ox, h - 8); ctx.moveTo(ox, oy); ctx.lineTo(w - 6, oy); ctx.stroke();
  // voltages
  vec(V2, sV, COL.core, "V2'", null);
  vec(d2.scale(kd), sV, COL.par, "I2'Z2'", [4, 3], V2);
  vec(vE, sV, COL.gap, "E1=E2'", null, null, 2, 12);
  vec(d1.scale(kd), sV, COL.par, 'I1Z1', [4, 3], vE);
  vec(vV1, sV, COL.fg, 'V1', null);
  // currents
  vec(I2, sI, COL.ok, "I2'", null);
  vec(iE, sI, '#b388ff', 'Ie', null, I2);
  vec(i1, sI, '#ffd166', 'I1', null);
  ctx.fillStyle = TXT; ctx.textAlign = 'left'; ctx.textBaseline = 'top';
  ctx.fillText('drops ×' + kd.toFixed(kd < 10 ? 1 : 0) + ', Ie ×' + ke.toFixed(ke < 10 ? 1 : 0) + ' (exaggerated)', 6, 4);
  ctx.textBaseline = 'bottom'; ctx.fillText('|V1| = ' + fE(f.V1.abs(), 'V', 5) + '   |V2\'| = ' + fE(f.V2.abs(), 'V', 5), 6, h - 3);
}

function workingB() {
  const p = RB.p, f = RB.fl, L = [];
  L.push('Turns ratio a = V1/V2 = ' + fN(S.V1) + '/' + fN(S.V2) + ' = ' + fN(p.a, 6) + '   (impedances measured on the secondary are multiplied by a²)');
  L.push('OC test (' + (S.ocSide === 'sec' ? 'secondary' : 'primary') + ' side): R_c = V²/P = ' + fN(S.ocV) + '²/' + fN(S.ocP) + ';  I_c = P/V = ' + fN(S.ocP / S.ocV) + ' A;  I_m = √(I²−I_c²) = ' + fN(Math.sqrt(Math.max(0, S.ocI * S.ocI - Math.pow(S.ocP / S.ocV, 2)))) + ' A;  X_m = V/I_m');
  L.push('  → referred to primary: R_c = ' + fE(p.Rc, 'Ω', 6) + ',  X_m = ' + fE(p.Xm, 'Ω', 6));
  L.push('SC test (' + (S.scSide === 'sec' ? 'secondary' : 'primary') + ' side): R_eq = P/I² = ' + fN(S.scP) + '/' + fN(S.scI) + '²;  Z_eq = V/I = ' + fN(S.scV / S.scI) + ' Ω;  X_eq = √(Z_eq² − R_eq²)');
  L.push('  → referred to primary: R_eq = ' + fE(p.Req, 'Ω', 6) + ',  X_eq = ' + fE(p.Xeq, 'Ω', 6) + '   (' + fP(p.Rpu, 2) + ' + j' + fP(p.Xpu, 2) + ' per unit)');
  L.push('Split assumed equal: R1 = R2\' = R_eq/2 = ' + fE(p.R1, 'Ω', 5) + ';  X1 = X2\' = X_eq/2 = ' + fE(p.X1, 'Ω', 5));
  L.push('Exact T-circuit: Z_in = Z1 + (Y_m + 1/(Z2\'+Z_L\'))⁻¹, Y_m = 1/R_c − j/X_m;  I1 = V1/Z_in, E1 = V1 − I1 Z1, I2\' = E1/(Z2\'+Z_L\')');
  L.push('Load: |Z_L\'| = V1²/(x·S) = ' + fE(Math.hypot(RB.ZL.re, RB.ZL.im), 'Ω', 5) + ' at ' + (RB.lead ? '−' : '+') + 'acos(' + S.pf + ') = ' + (RB.ZL.arg() * 180 / Math.PI).toFixed(2) + '°');
  L.push('Result: V2 = |V2\'|/a = ' + fE(f.V2.abs() / p.a, 'V', 6) + ';  no-load V2 = ' + fE(RB.reg.V2nl, 'V', 6) + ';  VR = (V2,nl − V2)/V2,nl = ' + fP(RB.reg.VR, 4));
  L.push('Losses: P_cu = |I1|²R1 + |I2\'|²R2\' = ' + fE(f.Pcu, 'W', 5) + ';  P_core = |E1|²/R_c = ' + fE(f.Pcore, 'W', 5) + ';  P_in − P_out − P_loss = ' + (f.Pin - f.Pout - f.Pcu - f.Pcore).toExponential(1) + ' W');
  L.push('Max efficiency at P_cu = P_core: x* = √(P_core/P_cu,rated) = √(' + fN(p.Pcore) + '/' + fN(p.PcuRated) + ') = ' + fN(RB.me.x, 6) + ';  η_max = x*·S·pf/(x*·S·pf + 2 P_core) = ' + fP(RB.me.eta, 4));
  return L;
}

function updateB(U) {
  if (!inited) return;
  computeB();
  if (!RB.p.ok) { U.warn.hidden = false; U.warn.textContent = RB.p.msg; renderHud(U.hud, [['Parameters', '—']]); U.working.set(RB.p.msg); return; }
  U.warn.hidden = true;
  const p = RB.p, f = RB.fl;
  const rows = [
    ['a = V1/V2', fN(p.a, 6)], ['I1 rated / I2 rated', fE(p.I1r, 'A', 4) + ' / ' + fE(p.I2r, 'A', 4)],
    ['R1 = R2\'', fE(p.R1, 'Ω', 5)], ['X1 = X2\'', fE(p.X1, 'Ω', 5)], ['R_c (pri)', fE(p.Rc, 'Ω', 5)], ['X_m (pri)', fE(p.Xm, 'Ω', 5)],
    ['R_eq / X_eq (pri)', fE(p.Req, 'Ω', 4) + ' / ' + fE(p.Xeq, 'Ω', 4)], ['Z_eq per unit', fP(p.Zeq * p.I1r / p.V1r, 3)],
    ['V2 at load', fE(f.V2.abs() / p.a, 'V', 5)], ['No-load V2', fE(RB.reg.V2nl, 'V', 5)],
    ['Regulation', (RB.reg.VR * 100).toFixed(3) + ' %', RB.reg.VR < 0 ? 'warn' : 'ok'],
    ['I1 / I2', fE(f.I1.abs(), 'A', 4) + ' / ' + fE(f.I2.abs() * p.a, 'A', 4)],
    ['P_in / P_out', fE(f.Pin, 'W', 5) + ' / ' + fE(f.Pout, 'W', 5)], ['P_cu / P_core', fE(f.Pcu, 'W', 4) + ' / ' + fE(f.Pcore, 'W', 4)],
    ['Efficiency', fP(f.eta, 3)], ['Max-η load x*', fN(RB.me.x, 5)], ['η_max (this pf)', fP(RB.me.eta, 3)],
  ];
  renderHud(U.hud, rows);
  U.working.set(workingB());
  plotVR(U.vr); plotEff(U.eff); plotPhasor(U.ph);
}

/* ---------- panel construction ---------- */
function touch() { FSP.state.touch(); }
let UA = null, UB = null;
function redrawAll() { if (!inited || !active) return; if (S.sub === 'a') updateA(UA); else updateB(UB); }

function numCtl(parent, key, refresh) {
  const o = CTL[key];
  const sl = FSP.ui.slider(parent, { label: o.label, min: o.min, max: o.max, step: o.step, log: o.log, unit: o.unit, value: S[key], digits: 6, onInput: v => { S[key] = o.int ? Math.round(v) : v; if (o.int) sl.set(S[key], true); touch(); refresh(); } });
  sl.el.title = o.tip || ''; widgets.push(() => sl.set(S[key], true));
  return sl;
}
function selCtl(parent, key, label, options, refresh, after) {
  const s = FSP.ui.select(parent, label, options, S[key], v => { S[key] = v; if (after) after(); touch(); refresh(); });
  widgets.push(() => { s.value = S[key]; }); return s;
}
function segBtns(parent, items, current, onPick) {
  const wrap = FSP.ui.el('div', { class: 'seg', role: 'group' }); const btns = {};
  items.forEach(([id, text]) => { const b = FSP.ui.el('button', { type: 'button', class: 'seg-btn', text, 'aria-pressed': 'false', onclick: () => onPick(id) }); btns[id] = b; wrap.appendChild(b); });
  parent.appendChild(wrap);
  return { el: wrap, set(id) { Object.keys(btns).forEach(k => { btns[k].setAttribute('aria-pressed', k === id ? 'true' : 'false'); btns[k].classList.toggle('active', k === id); }); } };
}

function buildA(root) {
  const U = FSP.ui, el = U.el;
  const layout = el('div', { class: 'layout' }); root.appendChild(layout);
  const ctl = el('div', { class: 'controls' }), stage = el('div', { class: 'stage' }); layout.appendChild(ctl); layout.appendChild(stage);
  const refresh = () => { redrawAll(); };
  const fsC = U.fieldset(ctl, 'Winding & core geometry');
  numCtl(fsC, 'N', refresh); numCtl(fsC, 'lc', refresh); numCtl(fsC, 'Ac', refresh); numCtl(fsC, 'lg', refresh); numCtl(fsC, 'fr', refresh);
  const fsM = U.fieldset(ctl, 'Core material');
  const murRow = el('div'), tabRow = el('div');
  selCtl(fsM, 'mode', 'Material', [['const', 'Constant μr'], ['table', 'B–H table (monotone interp.)']], refresh, () => { vis(); });
  fsM.appendChild(murRow); fsM.appendChild(tabRow);
  numCtl(murRow, 'mur', refresh);
  const ta = el('textarea', { rows: '9', spellcheck: 'false', 'aria-label': 'B-H table, one H B pair per line', style: 'width:100%;background:var(--panel2);color:var(--text);border:1px solid var(--border);border-radius:6px;font:12px var(--mono);padding:4px' });
  ta.value = S.table;
  const tmsg = el('div', { class: 'note', text: '' });
  const applyTable = () => { S.table = ta.value; tableParsed = M.parseTable(S.table); tmsg.textContent = tableParsed.ok ? 'Table OK: ' + (tableParsed.H.length - 1) + ' points, B_max = ' + tableParsed.B[tableParsed.B.length - 1] + ' T. Beyond the last point the slope is μ0.' : 'Invalid: ' + tableParsed.msg; touch(); refresh(); };
  ta.addEventListener('change', applyTable);
  tabRow.appendChild(el('div', { class: 'note', text: 'One "H(A/m) B(T)" pair per line (monotone cubic interpolation; (0,0) is implied).' }));
  tabRow.appendChild(ta); tabRow.appendChild(tmsg);
  const trow = el('div', { class: 'row' }); tabRow.appendChild(trow);
  U.button(trow, 'Reset to silicon steel', () => { ta.value = tableToText(M.DEFAULT_TABLE); applyTable(); });
  widgets.push(() => { ta.value = S.table; tableParsed = M.parseTable(S.table); tmsg.textContent = tableParsed.ok ? '' : 'Invalid: ' + tableParsed.msg; });
  const fsP = U.fieldset(ctl, 'Parallel return legs (optional)');
  const parRow = el('div');
  const parSel = selCtl(fsP, 'par', 'Parallel legs', [['0', 'None (single path)'], ['1', 'Two parallel legs R2‖R3']], refresh, () => vis());
  fsP.appendChild(parRow); ['l2', 'A2', 'l3', 'A3'].forEach(k => numCtl(parRow, k, refresh));
  const fsS = U.fieldset(ctl, 'Operating point');
  const iRow = el('div'), bRow = el('div');
  selCtl(fsS, 'solve', 'Solve', [['I', 'Given I → find φ (nonlinear solve)'], ['B', 'Given B_core → find I']], refresh, () => vis());
  fsS.appendChild(iRow); fsS.appendChild(bRow); numCtl(iRow, 'I', refresh); numCtl(bRow, 'Bc', refresh);
  function vis() { murRow.hidden = S.mode !== 'const'; tabRow.hidden = S.mode !== 'table'; parRow.hidden = S.par !== '1'; iRow.hidden = S.solve !== 'I'; bRow.hidden = S.solve !== 'B'; }
  widgets.push(vis); vis();
  const warn = el('div', { class: 'msg warn', hidden: '' }); stage.appendChild(warn);
  const g = el('div', { class: 'grid2' }); stage.appendChild(g);
  const w1 = el('div', { class: 'canvas-wrap' }), w2 = el('div', { class: 'canvas-wrap' }); g.appendChild(w1); g.appendChild(w2);
  const net = U.canvas(w1, { height: 250 }), bh = U.canvas(w2, { height: 250 });
  net.cv.setAttribute('aria-label', 'Reluctance network'); bh.cv.setAttribute('aria-label', 'B-H curve with operating point');
  const hud = el('div', { class: 'hud' }); stage.appendChild(hud);
  const working = U.working(stage);
  stage.appendChild(el('div', { class: 'note', text: 'Core length is the full mean path (the gap is not subtracted). Energy W = ∫NI dφ is the field energy; coenergy is NIφ − W. Dashed orange line: gap load line H_c = NI/l_c − (l_g A_c/(l_c A_g μ0))·B_c, whose intersection with the B–H curve is the operating point.' }));
  net.onResize(() => drawNetwork(net)); bh.onResize(() => drawBH(bh));
  return { hud, working, warn, net, bh };
}

function buildB(root) {
  const U = FSP.ui, el = U.el;
  const layout = el('div', { class: 'layout' }); root.appendChild(layout);
  const ctl = el('div', { class: 'controls' }), stage = el('div', { class: 'stage' }); layout.appendChild(ctl); layout.appendChild(stage);
  const refresh = () => { redrawAll(); };
  const fsR = U.fieldset(ctl, 'Ratings');
  ['kva', 'V1', 'V2', 'f'].forEach(k => numCtl(fsR, k, refresh));
  const fsO = U.fieldset(ctl, 'Open-circuit test');
  selCtl(fsO, 'ocSide', 'Measured on', [['sec', 'Secondary (LV) side'], ['pri', 'Primary (HV) side']], refresh);
  ['ocV', 'ocI', 'ocP'].forEach(k => numCtl(fsO, k, refresh));
  const fsS = U.fieldset(ctl, 'Short-circuit test');
  selCtl(fsS, 'scSide', 'Measured on', [['pri', 'Primary (HV) side'], ['sec', 'Secondary (LV) side']], refresh);
  ['scV', 'scI', 'scP'].forEach(k => numCtl(fsS, k, refresh));
  const fsL = U.fieldset(ctl, 'Load');
  numCtl(fsL, 'x', refresh); numCtl(fsL, 'pf', refresh);
  selCtl(fsL, 'lead', 'Power factor', [['0', 'Lagging (inductive)'], ['1', 'Leading (capacitive)']], refresh);
  const warn = el('div', { class: 'msg warn', hidden: '' }); stage.appendChild(warn);
  const hud = el('div', { class: 'hud' }); stage.appendChild(hud);
  const g = el('div', { class: 'grid2' }); stage.appendChild(g);
  const mk = (title, h) => { const w = el('div', { class: 'canvas-wrap' }); w.appendChild(el('div', { class: 'note', text: title, style: 'padding:4px 8px 0' })); g.appendChild(w); return U.canvas(w, { height: h }); };
  const vr = mk('Voltage regulation vs power factor', 240), eff = mk('Efficiency vs load', 240), ph = mk('Phasor diagram (referred to primary)', 280);
  vr.cv.setAttribute('aria-label', 'Voltage regulation versus power factor'); eff.cv.setAttribute('aria-label', 'Efficiency versus load'); ph.cv.setAttribute('aria-label', 'Phasor diagram');
  const working = U.working(stage);
  stage.appendChild(el('div', { class: 'note', text: 'Equivalent circuit: R1+jX1 in series, then R_c ∥ jX_m (exact T-model), then R2\'+jX2\' and the load. The tests give only the totals R_eq and X_eq, so the split R1=R2\', X1=X2\' is the usual assumption. Regulation is at constant rated primary voltage.' }));
  vr.onResize(() => plotVR(vr)); eff.onResize(() => plotEff(eff)); ph.onResize(() => plotPhasor(ph));
  return { hud, working, warn, vr, eff, ph };
}

/* ---------- URL state ---------- */
const NUMK = Object.keys(CTL), STRK = Object.keys(STR).filter(k => k !== 'table');
function stateGet() {
  const o = {}; NUMK.forEach(k => { o[k] = S[k]; }); STRK.forEach(k => { o[k] = S[k]; });
  o.table = S.table.trim().split(/\n+/).map(l => l.trim().replace(/\s+/g, ',')).join(';');
  return o;
}
const ENUMS = { sub: ['a', 'b'], mode: ['const', 'table'], solve: ['I', 'B'], par: ['0', '1'], ocSide: ['pri', 'sec'], scSide: ['pri', 'sec'], lead: ['0', '1'] };
function stateSet(o) {
  if (!o) return;
  NUMK.forEach(k => { if (o[k] !== undefined) { const v = parseFloat(o[k]); if (Number.isFinite(v)) { S[k] = Math.min(CTL[k].max, Math.max(CTL[k].min, CTL[k].int ? Math.round(v) : v)); } } });
  STRK.forEach(k => { if (o[k] !== undefined && ENUMS[k].indexOf(String(o[k])) >= 0) S[k] = String(o[k]); });
  if (typeof o.table === 'string' && o.table) { const txt = o.table.split(';').map(r => r.replace(/,/g, ' ')).join('\n'); if (M.parseTable(txt).ok) S.table = txt; }
  tableParsed = M.parseTable(S.table);
  if (inited) { widgets.forEach(fn => fn()); showSub(); redrawAll(); }
}

let subA = null, subB = null, seg = null;
function showSub() { if (!inited) return; subA.hidden = S.sub !== 'a'; subB.hidden = S.sub !== 'b'; seg.set(S.sub); }

FSP.registerTab({
  id: 'magnetics', title: 'Magnetics & Transformer',
  init(panel) {
    rootEl = panel;
    const U = FSP.ui, el = U.el;
    const head = el('div', { class: 'row' }); panel.appendChild(head);
    subA = el('div'); subB = el('div'); panel.appendChild(subA); panel.appendChild(subB);
    seg = segBtns(head, [['a', 'Magnetic circuit'], ['b', 'Single-phase transformer']], S.sub, id => { S.sub = id; showSub(); touch(); redrawAll(); });
    UA = buildA(subA); UB = buildB(subB);
    inited = true; showSub();
    FSP.state.bind('magnetics', { get: stateGet, set: stateSet });
    widgets.forEach(fn => fn());
  },
  activate() { active = true; redrawAll(); if (typeof requestAnimationFrame === 'function') requestAnimationFrame(redrawAll); },
  deactivate() { active = false; },
});

/* =============================================================
   Tests
   ============================================================= */
FSP.registerTests('magnetics', t => {
  const cm = (o) => M.makeCircuit(Object.assign({ N: 500, lc: 0.3, Ac: 1e-3, lg: 1e-3, fringe: 1, mat: M.makeMaterial({ mode: 'const', mur: 2000 }), par: null }, o || {}));
  const tab = () => { const p = M.parseTable(tableToText(M.DEFAULT_TABLE)); return M.makeMaterial({ mode: 'table', H: p.H, B: p.B }); };

  // --- spec example
  const c0 = cm(), a0 = M.analyze(c0, { NI: 500 * 1 });
  t.check('L = 0.2732 H for N=500, l=0.3, A=1e-3, mu_r=2000, 1 mm gap (±0.2%)', t.rel(a0.Lsec, 0.2732, 0.002), 'L=' + a0.Lsec.toFixed(5));
  t.check('L matches closed form N²/(lc/(μ0μrA) + lg/(μ0A)) (1e-12)', t.rel(a0.Lsec, 250000 / (0.3 / (MU0 * 2000 * 1e-3) + 1e-3 / (MU0 * 1e-3)), 1e-12));
  t.check('independent value R_total = 915140.92 A/Wb (python/numpy)', t.rel(a0.Rtot, 915140.9223, 1e-9), 'R=' + a0.Rtot.toFixed(2));
  const cNo = cm({ lg: 0 }), aNo = M.analyze(cNo, { NI: 500 });
  t.check('gap=0 raises L by 1 + lg·μr·Ac/(lc·Ag) = 7.6667', t.rel(aNo.Lsec / a0.Lsec, 1 + (1e-3 * 2000) / 0.3, 1e-10), 'ratio=' + (aNo.Lsec / a0.Lsec).toFixed(5));
  t.check('gap=0: L = N²μ0μrA/lc = 1.0053 H', t.rel(aNo.Lsec, 250000 * MU0 * 2000 * 1e-3 / 0.3, 1e-12));
  t.check('linear solve: φ = NI/R_total (rel 1e-12)', t.rel(a0.phi, 500 / a0.Rtot, 1e-12) && a0.sol.converged);
  t.check('solver reports convergence with bounded iterations', a0.sol.converged && a0.sol.iter > 5 && a0.sol.iter < 200, 'iter=' + a0.sol.iter);
  t.check('linear: incremental L equals secant L', t.rel(a0.Linc, a0.Lsec, 1e-6));
  t.check('linear: energy equals ½LI² (1e-9)', t.rel(a0.W, 0.5 * a0.Lsec * a0.I * a0.I, 1e-9), 'W=' + a0.W.toExponential(5));
  t.check('linear: B_gap = B_core with equal areas; Fg/F = Rg/Rtot', t.near(a0.m.Bg, a0.m.Bc, 1e-15) && t.rel(a0.m.Fg / a0.NI, a0.Rgap / a0.Rtot, 1e-10));
  const cF = cm({ fringe: 1.2 }), aF = M.analyze(cF, { NI: 500 });
  t.check('fringing (A_g = 1.2 A_c) lowers R_gap by 1.2', t.rel(aF.Rgap, a0.Rgap / 1.2, 1e-10));
  t.check('gap energy fraction = R_gap/R_total (linear)', t.rel(a0.Wgap / a0.W, a0.Rgap / a0.Rtot, 1e-6));

  // --- parallel legs
  const cp = cm({ par: { l2: 0.2, A2: 5e-4, l3: 0.4, A3: 8e-4 } }), ap = M.analyze(cp, { NI: 300 });
  t.check('parallel legs (linear): R_par = R2‖R3 and R_total closed form', t.rel(ap.Rtot, M.linearReluctance(cp), 1e-9), 'R=' + ap.Rtot.toFixed(1));
  t.check('parallel legs: flux splits as A/l and sums to φ', t.rel(ap.m.par.phi2 / ap.m.par.phi3, (5e-4 / 0.2) / (8e-4 / 0.4), 1e-9) && t.rel(ap.m.par.phi2 + ap.m.par.phi3, ap.phi, 1e-9));
  t.check('parallel legs: energy = ½LI²', t.rel(ap.W, 0.5 * ap.Lsec * ap.I * ap.I, 1e-8));

  // --- B-H table
  const mat = tab();
  const Hs = M.DEFAULT_TABLE.H, Bs = M.DEFAULT_TABLE.B;
  t.check('table: B(H) reproduces all nodes (1e-10)', Hs.every((h, i) => t.near(mat.Bof(h), Bs[i], 1e-10)));
  let mono = true, prev = -1; for (let i = 0; i <= 4000; i++) { const b = mat.Bof(i * 6); if (b < prev - 1e-12) mono = false; prev = b; }
  t.check('table: interpolation is monotone (no overshoot) on 4001 points', mono);
  t.check('table: H(B(H)) round trip (rel 1e-9) at 6 values', [3, 15, 55, 175, 800, 7000, 30000].every(h => t.rel(mat.Hof(mat.Bof(h)), h, 1e-9)));
  t.check('table: beyond last node slope is μ0', t.rel(mat.Bof(30000) - mat.Bof(20000), MU0 * 10000, 1e-12));
  t.check('table: odd symmetry', t.near(mat.Bof(-100), -mat.Bof(100), 1e-15));
  const cT = cm({ mat }), cTp = cm({ mat, par: { l2: 0.25, A2: 6e-4, l3: 0.35, A3: 4e-4 } });
  const cases = [[cT, 20], [cT, 400], [cT, 5000], [cT, 80000], [cTp, 50], [cTp, 700], [cTp, 30000], [cm({ mat, lg: 0 }), 90], [cm({ mat, lg: 0 }), 30000]];
  const worst = cases.map(([c, NI]) => { const r = M.analyze(c, { NI }); return Math.abs(M.mmf(c, r.phi).NI - NI) / NI; });
  t.check('nonlinear solve: forward re-evaluation reproduces NI within 1e-6 (9 cases incl. saturation, parallel legs, no gap)', Math.max.apply(null, worst) < 1e-6, 'worst=' + Math.max.apply(null, worst).toExponential(2));
  t.check('nonlinear solve: converged flag true on all cases', cases.every(([c, NI]) => M.analyze(c, { NI }).sol.converged));
  const aT1 = M.analyze(cm({ mat, lg: 0 }), { NI: 40 }), aT2 = M.analyze(cm({ mat, lg: 0 }), { NI: 20000 });
  t.check('table, no gap: NI=40 → H=133.3 A/m, B from table', t.near(aT1.m.Hc, 40 / 0.3, 1e-6) && t.rel(aT1.m.Bc, mat.Bof(40 / 0.3), 1e-12));
  const aT0 = M.analyze(cm({ mat, lg: 0 }), { NI: 6 });
  t.check('saturation indicator: linear at low B, knee/linear mid, saturated at high NI', aT0.sat === 'linear' && aT1.sat !== 'saturated' && aT2.sat === 'saturated', aT0.sat + ' / ' + aT1.sat + ' / ' + aT2.sat);
  t.check('saturation: incremental L < secant L; L falls as current rises', aT2.Linc < aT2.Lsec && aT2.Lsec < aT1.Lsec);
  t.check('mmf is monotone in flux (table + parallel)', (() => { let p = -1; for (let i = 0; i < 300; i++) { const v = M.mmf(cTp, i * 1e-5).NI; if (v < p) return false; p = v; } return true; })());
  const aE = M.analyze(cT, { NI: 3000 });
  t.check('nonlinear energy: 0 < W < NI·φ (coenergy positive)', aE.W > 0 && aE.W < aE.NI * aE.phi);
  const trap = (() => { let s = 0, N = 20000; for (let i = 0; i < N; i++) s += 0.5 * (M.mmf(cT, i * aE.phi / N).NI + M.mmf(cT, (i + 1) * aE.phi / N).NI) * aE.phi / N; return s; })();
  t.check('nonlinear energy agrees with fine trapezoid (1e-5)', t.rel(aE.W, trap, 1e-5));
  t.check('given-φ mode inverts given-NI mode', (() => { const r2 = M.analyze(cT, { phi: aE.phi }); return t.rel(r2.NI, 3000, 1e-6); })());
  t.check('parseTable accepts commas/semicolons and rejects non-monotone / duplicate H', M.parseTable('10,0.5;20,0.9;40,1.2').ok && !M.parseTable('10 0.5\n20 0.4\n30 0.9').ok && !M.parseTable('10 0.5\n10 0.6\n30 0.9').ok && !M.parseTable('hello').ok);

  // --- transformer
  const spec = { S: 10000, V1: 2400, V2: 240, oc: { V: 240, I: 1.5, P: 120, side: 'sec' }, sc: { V: 60, I: 10000 / 2400, P: 150, side: 'pri' } };
  const p = M.xfmrFromTests(spec);
  t.check('xfmr params from OC/SC (python): R_c=48000, X_m=16970.56 Ω; R_eq=8.64, X_eq=11.52 Ω', p.ok && t.rel(p.Rc, 48000, 1e-9) && t.rel(p.Xm, 16970.5627, 1e-8) && t.rel(p.Req, 8.64, 1e-3) && t.rel(p.Xeq, 11.52, 1e-3), p.ok ? 'Req=' + p.Req.toFixed(4) + ' Xeq=' + p.Xeq.toFixed(4) : p.msg);
  const spec2 = { S: 10000, V1: 2400, V2: 240, oc: { V: 2400, I: 0.15, P: 120, side: 'pri' }, sc: { V: 6, I: 41.6667, P: 150, side: 'sec' } };
  const p2 = M.xfmrFromTests(spec2);
  t.check('side referral: same transformer measured on the other sides gives same R_c, X_m, R_eq, X_eq (1e-4)', p2.ok && t.rel(p2.Rc, p.Rc, 1e-9) && t.rel(p2.Xm, p.Xm, 1e-8) && t.rel(p2.Req, p.Req, 1e-3) && t.rel(p2.Xeq, p.Xeq, 1e-3));
  t.check('invalid test data rejected (P_oc > V·I, P_sc > V·I, zero, NaN)', !M.xfmrFromTests(Object.assign({}, spec, { oc: { V: 240, I: 0.4, P: 120, side: 'sec' } })).ok && !M.xfmrFromTests(Object.assign({}, spec, { sc: { V: 10, I: 4, P: 150, side: 'pri' } })).ok && !M.xfmrFromTests(Object.assign({}, spec, { S: 0 })).ok && !M.xfmrFromTests(Object.assign({}, spec, { V2: NaN })).ok);
  const nl = M.xfmrForward(p, C(2400, 0), null);
  t.check('no-load secondary reproduces rated V2 within 0.5%', t.rel(nl.V2.abs() / p.a, 240, 0.005), 'V2nl=' + (nl.V2.abs() / p.a).toFixed(4));
  t.check('no-load input reproduces OC test: P≈P_oc (1%), |I|≈I_oc (1%)', t.rel(nl.Pin, 120, 0.01) && t.rel(nl.I1.abs() * p.a, 1.5, 0.01), 'P=' + nl.Pin.toFixed(3) + ' I=' + (nl.I1.abs() * p.a).toFixed(4));
  const sc = M.xfmrForward(p, C(60, 0), C(0, 0));
  t.check('shorted secondary at V_sc reproduces SC test: I1≈I_sc (0.5%), P≈P_sc (0.5%)', t.rel(sc.I1.abs(), 10000 / 2400, 0.005) && t.rel(sc.Pin, 150, 0.005), 'I=' + sc.I1.abs().toFixed(4) + ' P=' + sc.Pin.toFixed(3));
  const worstV = [[1, 0.8, false], [1, 0.8, true], [1, 1, false], [0.5, 0.6, false], [1.2, 0.9, true]].map(([x, pf, ld]) => {
    const ZL = M.loadZ(p, x, pf, ld), back = M.xfmrBack(p, C(2400, 0), ZL), fwd = M.xfmrForward(p, back.V1, ZL);
    return Math.abs(fwd.V2.abs() - 2400) / 2400;
  });
  t.check('load analysis at rated load reproduces rated V2: back-solve V1, forward V2 within 0.5% (5 loads)', Math.max.apply(null, worstV) < 0.005, 'worst=' + Math.max.apply(null, worstV).toExponential(2));
  t.check('rated-load back-solve: V1 exceeds rated by about the regulation (0.8 lag, 2.2–2.8%)', (() => { const b = M.xfmrBack(p, C(2400, 0), M.loadZ(p, 1, 0.8, false)); const d = b.V1.abs() / 2400 - 1; return d > 0.022 && d < 0.028; })());
  const fl = M.xfmrForward(p, C(2400, 0), M.loadZ(p, 1, 0.8, false));
  t.check('power balance P_in = P_out + P_cu + P_core (1e-9 rel)', t.rel(fl.Pin, fl.Pout + fl.Pcu + fl.Pcore, 1e-9));
  t.check('KCL at magnetising branch: I1 = I2\' + Ie (1e-12)', fl.I1.sub(fl.I2).sub(fl.Ie).abs() < 1e-12 * fl.I1.abs() && fl.Ie.sub(fl.Ic).sub(fl.Im).abs() < 1e-12);
  t.check('KVL: V1 = E1 + I1·Z1 and E1 = V2 + I2·Z2\' (1e-12)', fl.V1.sub(fl.Vm).sub(fl.I1.mul(C(p.R1, p.X1))).abs() < 1e-9 && fl.Vm.sub(fl.V2).sub(fl.I2.mul(C(p.R2p, p.X2p))).abs() < 1e-9);
  const reg = (x, pf, ld) => M.regulation(p, x, pf, ld).VR;
  const vrApprox = (x, pf, ld) => x * (p.Rpu * pf + (ld ? -1 : 1) * p.Xpu * Math.sqrt(1 - pf * pf));
  t.check('regulation ≈ x(R_pu cosφ ± X_pu sinφ) within 0.1 percentage points (0.8 lag, 0.8 lead, unity)', [[1, 0.8, false], [1, 0.8, true], [1, 1, false], [0.5, 0.6, false]].every(a => Math.abs(reg.apply(null, a) - vrApprox.apply(null, a)) < 0.001), 'VR(0.8lag)=' + (100 * reg(1, 0.8, false)).toFixed(4) + '%');
  t.check('regulation sign: lagging > 0; leading at low pf < 0; unity small positive', reg(1, 0.8, false) > 0.02 && reg(1, 0.3, true) < 0 && reg(1, 1, false) > 0 && reg(1, 1, false) < reg(1, 0.8, false));
  t.check('regulation independent value (python complex solve): 0.8 lag = 2.345496 %', t.near(100 * reg(1, 0.8, false), 2.345496, 1e-5), (100 * reg(1, 0.8, false)).toFixed(5));
  const me = M.maxEff(p, 0.8);
  t.check('max-efficiency load fraction equals √(P_core/P_cu,rated) (1e-6)', t.near(me.x, Math.sqrt(p.Pcore / p.PcuRated), 1e-6) && t.near(me.x, 0.894427, 1e-6), 'x*=' + me.x.toFixed(6));
  const gss = (() => { let a = 0.1, b = 1.5; const g = (Math.sqrt(5) - 1) / 2; for (let i = 0; i < 200; i++) { const c1 = b - g * (b - a), c2 = a + g * (b - a); if (M.effApprox(p, c1, 0.8) > M.effApprox(p, c2, 0.8)) b = c2; else a = c1; } return 0.5 * (a + b); })();
  t.check('numerical maximum of η(x) (golden section) lands on x* (1e-6) with P_cu = P_core', t.near(gss, me.x, 1e-6) && t.rel(me.Pcu, me.Pcore, 1e-12), 'gss=' + gss.toFixed(7));
  t.check('η_max matches closed form x·S·pf/(x·S·pf+2·P_core) and approx η(x*)', t.rel(me.eta, M.effApprox(p, me.x, 0.8), 1e-12));
  t.check('exact-model peak efficiency within 0.2% of approximation and located within 5% of x*', (() => { let bx = 0, be = 0; for (let x = 0.3; x <= 1.5; x += 0.001) { const e = M.effExact(p, x, 0.8, false); if (e > be) { be = e; bx = x; } } return Math.abs(be - me.eta) < 0.002 && Math.abs(bx / me.x - 1) < 0.05; })());
  t.check('efficiency at rated load 0.8 pf lag: approx 96.735 %, exact T-model 96.6469 % (python)', t.near(100 * M.effApprox(p, 1, 0.8), 100 * 8000 / (8000 + 120 + 150), 1e-9) && t.near(100 * fl.eta, 96.646868, 1e-5), 'exact=' + (100 * fl.eta).toFixed(3));
  t.check('no NaN/Infinity from outputs at extreme loads (x=0.01,1.5, pf=0.05)', [0.01, 1.5].every(x => [0.05, 1].every(pf => [false, true].every(ld => { const r = M.regulation(p, x, pf, ld), e = M.effExact(p, x, pf, ld); return Number.isFinite(r.VR) && Number.isFinite(e); }))));
});
})();
