/* Phase 3: Transmission lines + Smith chart (tab 'tline').
   Conventions: e^{jwt}; d / l = distance from the LOAD toward the generator, in wavelengths unless stated.
   Gamma(d) = Gamma_L e^{-2 gamma d}, Zin = Z0 (1+Gamma)/(1-Gamma) (no tan(), so no singularities).
   Z0 may be complex (voltage-wave reflection coefficient (ZL-Z0)/(ZL+Z0) is used). */
(function () {
  'use strict';
  const TWO_PI = 2 * Math.PI, NP_PER_DB = Math.log(10) / 20;
  const INF = Infinity;
  const mod = (x, m) => ((x % m) + m) % m;

  /* ===================== pure math ===================== */
  const M = {};
  M.gammaL = (Z0, ZL) => {
    if (!Number.isFinite(ZL.re) || !Number.isFinite(ZL.im)) return C(1, 0);   // open circuit
    return ZL.sub(Z0).div(ZL.add(Z0));
  };
  // alpha: attenuation in Np per wavelength (0 = lossless). d in wavelengths.
  M.gammaAt = (GL, d, alpha) => {
    const a = alpha || 0, ph = -2 * TWO_PI * d, m = Math.exp(-2 * a * d);
    return GL.mul(C(m * Math.cos(ph), m * Math.sin(ph)));
  };
  M.zinFromGamma = (Z0, G) => {
    const den = C(1, 0).sub(G);
    if (den.abs() < 1e-15) return C(INF, 0);
    return Z0.mul(C(1, 0).add(G)).div(den);
  };
  M.zin = (Z0, ZL, d, alpha) => M.zinFromGamma(Z0, M.gammaAt(M.gammaL(Z0, ZL), d, alpha));
  // classic tan form, used only as an independent cross-check in the tests
  M.zinTan = (Z0, ZL, d) => { const t = Math.tan(TWO_PI * d), jt = C(0, t); return Z0.mul(ZL.add(Z0.mul(jt))).div(Z0.add(ZL.mul(jt))); };
  M.vswr = G => { const m = G.abs(); return m >= 1 - 1e-12 ? INF : (1 + m) / (1 - m); };
  M.returnLossDb = G => { const m = G.abs(); return m < 1e-15 ? INF : -20 * Math.log10(m); };
  M.mismatchLossDb = G => { const q = 1 - G.abs2(); return q <= 1e-15 ? INF : -10 * Math.log10(q); };
  // Normalised V (with V+ = 1 at the load) and Z0*I along the line.
  M.vi = (GL, d, alpha) => {
    const a = alpha || 0, g = TWO_PI * d;
    const f = C(Math.exp(a * d) * Math.cos(g), Math.exp(a * d) * Math.sin(g));       // forward wave e^{gamma d}
    const r = GL.mul(C(Math.exp(-a * d) * Math.cos(-g), Math.exp(-a * d) * Math.sin(-g))); // reflected
    return { V: f.add(r), ZI: f.sub(r) };
  };
  M.zToGamma = (z) => z.sub(C(1, 0)).div(z.add(C(1, 0)));
  M.gammaToZ = (G) => M.zinFromGamma(C(1, 0), G);
  // Wavelengths-toward-generator reading of a reflection-coefficient angle (0 at Gamma=-1, 0.25 at Gamma=+1).
  M.wtg = G => mod((Math.PI - Math.atan2(G.im, G.re)) / (4 * Math.PI), 0.5);

  const isReal = Z => Math.abs(Z.im) <= 1e-9 * Math.max(1, Math.abs(Z.re));
  // Quarter-wave transformer placed at the first voltage maximum (where Zin is real = Z0*VSWR).
  M.quarterWave = (Z0, ZL) => {
    if (!isReal(Z0) || Z0.re <= 0) return null;
    const G = M.gammaL(Z0, ZL), m = G.abs();
    if (m >= 1 - 1e-9) return null;
    if (m < 1e-12) return { matched: true, d: 0, Rin: Z0.re, Zt: Z0.re };
    const d = mod(Math.atan2(G.im, G.re) / (4 * Math.PI), 0.5);
    const Rin = Z0.re * (1 + m) / (1 - m);
    return { matched: false, d, Rin, Zt: Math.sqrt(Z0.re * Rin) };
  };
  M.quarterWaveCheck = (Z0, ZL, sol) => {   // Gamma seen at the generator side of the transformer
    const Zd = M.zin(Z0, ZL, sol.d, 0);
    const Zq = M.zin(C(sol.Zt, 0), Zd, 0.25, 0);
    return M.gammaL(Z0, Zq);
  };
  // Single shunt stub (same Z0, short- or open-circuited). Returns up to two solutions sorted by d.
  M.shuntStub = (Z0, ZL, type) => {
    if (!isReal(Z0) || Z0.re <= 0) return null;
    const G = M.gammaL(Z0, ZL), m = G.abs();
    if (m >= 1 - 1e-9 || m < 1e-12) return m < 1e-12 ? [] : null;
    const th = Math.atan2(G.im, G.re), out = [];
    [1, -1].forEach(s => {
      const phi = s * Math.acos(-m), d = mod((th - phi) / (4 * Math.PI), 0.5);
      const Gd = M.gammaAt(G, d, 0), y = C(1, 0).sub(Gd).div(C(1, 0).add(Gd)), b = y.im;
      const ls = type === 'open' ? mod(Math.atan(-b) / TWO_PI, 0.5) : mod(Math.atan2(1, b) / TWO_PI, 0.5);
      out.push({ d, ls, b, y, Gd });
    });
    return out.sort((p, q) => p.d - q.d);
  };
  // Normalised admittance of load+line+stub at the stub position (should be 1 + j0), via the tan form.
  M.stubCheck = (Z0, ZL, sol, type) => {
    const Zd = M.zinTan(Z0, ZL, sol.d), y = C(Z0.re, 0).div(Zd);
    const t = Math.tan(TWO_PI * sol.ls);
    const ys = type === 'open' ? C(0, t) : C(0, -1 / t);
    return y.add(ys);
  };
  FSP.math.tline = M;

  /* ===================== tests ===================== */
  FSP.registerTests('tline', t => {
    const Z0 = C(50, 0);
    let G = M.gammaL(Z0, C(50, 0));
    t.check('ZL=Z0 gives Gamma=0 and VSWR=1', G.abs() < 1e-15 && M.vswr(G) === 1, 'VSWR=' + M.vswr(G));
    G = M.gammaL(Z0, C(100, 0));
    t.check('Z0=50, ZL=100: Gamma=1/3', t.near(G.re, 1 / 3, 1e-12) && t.near(G.im, 0, 1e-12));
    t.check('Z0=50, ZL=100: VSWR=2', t.near(M.vswr(G), 2, 1e-12));
    t.check('Z0=50, ZL=100: return loss 9.542 dB', t.near(M.returnLossDb(G), 9.542, 0.001), M.returnLossDb(G).toFixed(4));
    t.check('Z0=50, ZL=100: mismatch loss 0.5115 dB', t.near(M.mismatchLossDb(G), -10 * Math.log10(8 / 9), 1e-9) && t.near(M.mismatchLossDb(G), 0.5115, 1e-4));
    const Gs = M.gammaL(Z0, C(0, 0));
    t.check('short: |Gamma|=1, VSWR infinite', t.near(Gs.abs(), 1, 1e-15) && M.vswr(Gs) === INF);
    const zq = M.zin(Z0, C(0, 0), 0.25, 0);
    t.check('short: Zin at lambda/4 is open (|Zin| > 1e6 Z0)', zq.abs() > 1e6 * 50 || !Number.isFinite(zq.abs()), '|Zin|=' + zq.abs().toExponential(2));
    const zo = M.zin(Z0, C(INF, 0), 0.25, 0);
    t.check('open: Zin at lambda/4 is a short (|Zin| < 1e-6 Z0)', zo.abs() < 1e-6 * 50, '|Zin|=' + zo.abs().toExponential(2));
    const qw = M.quarterWave(Z0, C(100, 0));
    t.check('quarter-wave 50 -> 100 needs Zt=70.71 ohm', t.near(qw.Zt, 70.7107, 1e-3), qw.Zt.toFixed(4));
    t.check('quarter-wave transformer gives Gamma=0 (50 -> 100)', M.quarterWaveCheck(Z0, C(100, 0), qw).abs() < 1e-12);
    const ZLc = C(25, 50), qc = M.quarterWave(Z0, ZLc);
    t.check('quarter-wave on complex load (at voltage max) gives Gamma=0', M.quarterWaveCheck(Z0, ZLc, qc).abs() < 1e-12, '|G|=' + M.quarterWaveCheck(Z0, ZLc, qc).abs().toExponential(2));
    // periodicity
    let worst = 0;
    [[25, 50], [100, -30], [3, 0.5], [500, 800], [50, 0], [12, -90]].forEach(([r, x]) => [0.013, 0.2, 0.377, 1.1].forEach(l => {
      const a = M.zin(Z0, C(r, x), l, 0), b = M.zin(Z0, C(r, x), l + 0.5, 0);
      worst = Math.max(worst, a.sub(b).abs() / Math.max(1, a.abs()));
    }));
    t.check('Zin periodic with lambda/2 (rel error < 1e-9)', worst < 1e-9, worst.toExponential(2));
    // cross-check vs tan form (python: Z0=50, ZL=25+j50, l=0.1 lambda -> 184.7522+70.2557j)
    const zt = M.zin(Z0, ZLc, 0.1, 0);
    t.check('Zin(0.1 lambda) = 184.752+j70.256 (numpy reference)', t.near(zt.re, 184.75223623645647, 1e-9) && t.near(zt.im, 70.25570694848419, 1e-9), zt.re.toFixed(4) + '+j' + zt.im.toFixed(4));
    let dtan = 0;
    [0.01, 0.07, 0.19, 0.31, 0.43].forEach(l => dtan = Math.max(dtan, M.zin(Z0, ZLc, l, 0).sub(M.zinTan(Z0, ZLc, l)).abs()));
    t.check('Gamma form equals tan form away from singularities (<1e-9)', dtan < 1e-9, dtan.toExponential(2));
    const z4 = M.zin(Z0, ZLc, 0.25, 0), zq4 = C(2500, 0).div(ZLc);
    t.check('quarter-wave line inverts: Zin = Z0^2/ZL', z4.sub(zq4).abs() < 1e-9);
    const z8 = M.zin(Z0, C(0, 0), 0.125, 0);
    t.check('shorted line, lambda/8: Zin = jZ0', t.near(z8.re, 0, 1e-9) && t.near(z8.im, 50, 1e-9));
    t.check('Zin at l=0 equals ZL', M.zin(Z0, ZLc, 0, 0).sub(ZLc).abs() < 1e-12);
    // complex Z0
    const Zc = C(50, -3);
    t.check('complex Z0: ZL=Z0 gives Gamma=0 at every length', [0, 0.17, 0.5].every(l => M.zin(Zc, Zc, l, 0).sub(Zc).abs() < 1e-12));
    t.check('complex Z0: Zin(l=0)=ZL, periodic lambda/2', M.zin(Zc, C(80, 20), 0, 0).sub(C(80, 20)).abs() < 1e-12 && M.zin(Zc, C(80, 20), 0.5, 0).sub(C(80, 20)).abs() < 1e-9);
    // stub matching
    const ex = M.shuntStub(Z0, ZLc, 'short');
    t.check('short stub, ZL=25+j50: d = 0.2933 / 0.4369 lambda, L = 0.0898 / 0.4102 (numpy)',
      ex.length === 2 && t.near(ex[0].d, 0.293338613545214, 1e-9) && t.near(ex[1].d, 0.4368696743745032, 1e-9) && t.near(ex[0].ls, 0.08975425899284405, 1e-9) && t.near(ex[1].ls, 0.410245741007156, 1e-9));
    const eo = M.shuntStub(Z0, ZLc, 'open');
    t.check('open stub lengths 0.3398 / 0.1602 (numpy)', t.near(eo[0].ls, 0.339754258992844, 1e-9) && t.near(eo[1].ls, 0.16024574100715602, 1e-9));
    let wst = 0;
    [[25, 50], [100, -80], [10, 5], [400, 0], [60, -200]].forEach(([r, x]) => ['short', 'open'].forEach(ty => M.shuntStub(Z0, C(r, x), ty).forEach(s => {
      wst = Math.max(wst, M.stubCheck(Z0, C(r, x), s, ty).sub(C(1, 0)).abs());
    })));
    t.check('stub solutions give y = 1 + j0 (short and open, 5 loads, tan-form check)', wst < 1e-9, wst.toExponential(2));
    t.check('stub helper: none for matched / lossless-reactive loads', M.shuntStub(Z0, C(50, 0), 'short').length === 0 && M.shuntStub(Z0, C(0, 30), 'short') === null);
    // V/I standing wave
    let vmax = 0, vmin = INF, imax = 0;
    const GLt = M.gammaL(Z0, C(25, 50));
    for (let i = 0; i <= 2000; i++) { const v = M.vi(GLt, i / 2000 * 0.5, 0); const a = v.V.abs(); vmax = Math.max(vmax, a); vmin = Math.min(vmin, a); imax = Math.max(imax, v.ZI.abs()); }
    t.check('max|V|/min|V| over a half-wave = VSWR', t.rel(vmax / vmin, M.vswr(GLt), 1e-4), (vmax / vmin).toFixed(5) + ' vs ' + M.vswr(GLt).toFixed(5));
    t.check('max|V| = 1+|Gamma|, max|Z0 I| = 1+|Gamma|', t.near(vmax, 1 + GLt.abs(), 1e-6) && t.near(imax, 1 + GLt.abs(), 1e-6));
    const v0 = M.vi(GLt, 0, 0), zl0 = v0.V.div(v0.ZI).scale(50);
    t.check('V/I at the load equals ZL', zl0.sub(C(25, 50)).abs() < 1e-9);
    // Smith helpers
    t.check('Smith: z=1 -> Gamma=0, z=0 -> -1, wavelength scale 0 at short / 0.25 at open',
      M.zToGamma(C(1, 0)).abs() < 1e-15 && t.near(M.zToGamma(C(0, 0)).re, -1, 1e-15) && t.near(M.wtg(C(-1, 0)), 0, 1e-12) && t.near(M.wtg(C(1, 0)), 0.25, 1e-12) && t.near(M.wtg(C(0, 1)), 0.125, 1e-12));
    const rt = M.gammaToZ(M.zToGamma(C(0.3, -1.7)));
    t.check('Smith: z -> Gamma -> z round trip', rt.sub(C(0.3, -1.7)).abs() < 1e-12);
    // lossy line: |Gamma| shrinks as exp(-2 alpha d)
    const Gl = M.gammaAt(C(0.5, 0), 0.25, 0.1);
    t.check('lossy line: |Gamma(d)| = |Gamma_L| e^{-2 alpha d}', t.near(Gl.abs(), 0.5 * Math.exp(-2 * 0.1 * 0.25), 1e-12));
  });

  /* ===================== UI ===================== */
  const C0 = 299792458;
  const S = { z0r: 50, z0i: 0, zlr: 25, zli: 50, f: 100e6, vf: 1, len: 0.3, unit: 'lam', loss: 0, match: 'qw' };
  const UNITS = ['lam', 'm'], MATCH = ['qw', 'short', 'open'];
  let ui = null, smith = null, vi = null, geo = null, active = false;

  const fnum = (x, d) => Number.isFinite(x) ? (Math.abs(x) >= 1e5 ? x.toExponential(2) : x.toFixed(d === undefined ? 3 : d)) : '∞';
  const fC = (z, d) => { if (!Number.isFinite(z.re) || !Number.isFinite(z.im) || z.abs() > 1e9) return '∞'; return fnum(z.re, d) + (z.im < 0 ? ' − j' : ' + j') + fnum(Math.abs(z.im), d); };
  const fPolar = (z, d) => Number.isFinite(z.abs()) ? fnum(z.abs(), d === undefined ? 4 : d) + ' ∠ ' + fnum(z.arg() * 180 / Math.PI, 2) + '°' : '—';
  const fLam = x => Number.isFinite(x) ? x.toFixed(4) + ' λ' : '—';
  const Z0c = () => C(S.z0r, S.z0i), ZLc = () => C(S.zlr, S.zli);
  const alphaNp = () => S.loss * NP_PER_DB;                 // loss in dB per wavelength -> Np/λ
  const lambdaM = () => S.vf * C0 / S.f;
  const parseNum = s => { s = String(s).trim().toLowerCase(); if (/^[+]?(inf|infinity|∞|open)$/.test(s)) return 1e15; if (/^-(inf|infinity|∞)$/.test(s)) return -1e15; return FSP.parseSI(s); };
  const clamp = (v, a, b) => Math.min(b, Math.max(a, v));

  function compute() {
    const Z0 = Z0c(), ZL = ZLc(), a = alphaNp();
    const GL = M.gammaL(Z0, ZL), Gin = M.gammaAt(GL, S.len, a), Zin = M.zinFromGamma(Z0, Gin);
    const lam = lambdaM();
    return { Z0, ZL, GL, Gin, Zin, lam, a, m: GL.abs(), valid: Z0.re > 0 && Number.isFinite(GL.re) && Number.isFinite(GL.im) };
  }

  function numField(parent, label, key, unit, extra) {
    const inp = FSP.ui.el('input', { type: 'text', inputmode: 'decimal', 'aria-label': label, size: 9 });
    inp.style.cssText = 'background:var(--panel2);color:var(--text);border:1px solid var(--border);border-radius:6px;padding:3px 5px;font:12px var(--mono);width:100%;min-width:0';
    const row = FSP.ui.el('div', { class: 'ctl' }, FSP.ui.el('label', { text: label }), inp, FSP.ui.el('span', { text: unit || '' }));
    const show = () => { const v = S[key]; inp.value = Math.abs(v) >= 1e12 ? (v > 0 ? '∞' : '−∞') : String(+v.toPrecision(6)); };
    inp.addEventListener('change', () => { const v = parseNum(inp.value.replace('−', '-')); if (Number.isFinite(v)) { S[key] = v; if (extra) extra(v); changed(); } show(); });
    parent.appendChild(row); show();
    return { show, el: inp };
  }

  function init(panel) {
    const U = FSP.ui;
    const layout = U.el('div', { class: 'layout' }), aside = U.el('aside', { class: 'controls' }), stage = U.el('div', { class: 'stage' });
    layout.appendChild(aside); layout.appendChild(stage); panel.appendChild(layout);
    ui = { fields: {} };

    const fl = U.fieldset(aside, 'Line');
    ui.fields.z0r = numField(fl, 'Z0 real', 'z0r', 'Ω');
    ui.fields.z0i = numField(fl, 'Z0 imag', 'z0i', 'Ω');
    ui.f = U.slider(fl, { label: 'Frequency', min: 1e6, max: 1e10, value: S.f, log: true, unit: 'Hz', digits: 4, onInput: v => { S.f = v; changed(); } });
    ui.vf = U.slider(fl, { label: 'Vel. factor', min: 0.3, max: 1, step: 0.01, value: S.vf, digits: 3, onInput: v => { S.vf = v; changed(); } });
    ui.loss = U.slider(fl, { label: 'Loss', min: 0, max: 2, step: 0.01, value: S.loss, unit: 'dB/λ', digits: 3, onInput: v => { S.loss = v; changed(); } });
    // length: slider in wavelengths (0..2) + text box in chosen unit
    const lenSl = U.el('input', { type: 'range', min: 0, max: 2, step: 0.001, 'aria-label': 'Line length slider' });
    const lenBox = U.el('input', { type: 'text', inputmode: 'decimal', 'aria-label': 'Line length', size: 8 });
    lenBox.style.cssText = 'background:var(--panel2);color:var(--text);border:1px solid var(--border);border-radius:6px;padding:3px 5px;font:12px var(--mono);width:78px';
    const unitSel = U.el('select', { 'aria-label': 'Length unit' }, U.el('option', { value: 'lam', text: 'λ' }), U.el('option', { value: 'm', text: 'm' }));
    fl.appendChild(U.el('div', { class: 'ctl' }, U.el('label', { text: 'Length l' }), lenSl, U.el('span', null, lenBox, unitSel)));
    ui.showLen = () => { lenSl.value = clamp(S.len, 0, 2); unitSel.value = S.unit; lenBox.value = S.unit === 'lam' ? String(+S.len.toPrecision(5)) : String(+(S.len * lambdaM()).toPrecision(5)); };
    lenSl.addEventListener('input', () => { S.len = parseFloat(lenSl.value); ui.showLen(); changed(); });
    lenBox.addEventListener('change', () => { const v = FSP.parseSI(lenBox.value); if (Number.isFinite(v) && v >= 0) { S.len = S.unit === 'lam' ? v : v / lambdaM(); changed(); } ui.showLen(); });
    unitSel.addEventListener('change', () => { S.unit = unitSel.value; ui.showLen(); changed(); });
    ui.lamNote = U.el('div', { class: 'note' }); fl.appendChild(ui.lamNote);

    aside.querySelectorAll('input[type=text]').forEach(i => { i.style.cssText = 'background:var(--panel2);color:var(--text);border:1px solid var(--border);border-radius:6px;padding:3px 5px;font:12px var(--mono);width:78px;min-width:0'; });
    const ll = U.fieldset(aside, 'Load ZL (click the Smith chart to set)');
    ui.fields.zlr = numField(ll, 'R', 'zlr', 'Ω');
    ui.fields.zli = numField(ll, 'X', 'zli', 'Ω');
    const seg = U.el('div', { class: 'seg' }); ll.appendChild(seg);
    const presets = [['Match', () => [S.z0r, S.z0i]], ['Short', () => [0, 0]], ['Open', () => [1e15, 0]], ['2·Z0', () => [2 * S.z0r, 0]], ['25 + j50', () => [25, 50]], ['jZ0', () => [0, S.z0r]], ['−jZ0', () => [0, -S.z0r]]];
    presets.forEach(([n, f]) => U.el('span') && seg.appendChild(U.el('button', { type: 'button', class: 'seg-btn', text: n, onclick: () => { const v = f(); S.zlr = v[0]; S.zli = v[1]; syncLoad(); changed(); } })));

    const fm = U.fieldset(aside, 'Matching helper');
    U.select(fm, 'Method', [['qw', 'Quarter-wave transformer'], ['short', 'Shunt stub (short-circuit)'], ['open', 'Shunt stub (open-circuit)']], S.match, v => { S.match = v; changed(); });
    ui.matchOut = U.el('div', { class: 'note mono' }); fm.appendChild(ui.matchOut);

    ui.hud = U.el('div', { class: 'hud', 'aria-live': 'off' });
    const gw = U.el('div', { class: 'grid2' });
    smith = U.canvas(null, { height: 400 }); smith.cv.setAttribute('tabindex', '0'); smith.cv.setAttribute('aria-label', 'Smith chart. Click or drag to set the load; arrow keys nudge it.');
    vi = U.canvas(null, { height: 300 }); vi.cv.setAttribute('aria-label', 'Voltage and current magnitude along the line');
    const w1 = U.el('div', { class: 'canvas-wrap' }); w1.appendChild(smith.cv);
    const w2 = U.el('div', { class: 'canvas-wrap' }); w2.appendChild(vi.cv);
    gw.appendChild(w1); gw.appendChild(w2);
    stage.appendChild(ui.hud); stage.appendChild(gw);
    ui.work = U.working(stage);
    const th = U.el('details', { class: 'theory' }, U.el('summary', { text: 'Theory & conventions' }), U.el('ul', null,
      U.el('li', { text: 'Γ_L = (Z_L − Z0)/(Z_L + Z0); Γ(d) = Γ_L e^(−2γd) with d measured from the load toward the generator, γ = α + jβ, β = 2π/λ, λ = v_f·c/f.' }),
      U.el('li', { text: 'Z_in(d) = Z0 (1 + Γ(d))/(1 − Γ(d)): the Γ form has no tan(βl) singularity, so open and shorted lines at λ/4 behave. Z_in is periodic with λ/2.' }),
      U.el('li', { text: 'VSWR = (1+|Γ|)/(1−|Γ|), RL = −20 log|Γ|, mismatch loss = −10 log(1−|Γ|²). With a complex Z0 the same voltage-wave Γ is used and the chart is normalised to Z0.' }),
      U.el('li', { text: 'Smith chart: moving toward the generator rotates clockwise, one full turn = λ/2. The outer scale reads wavelengths toward generator (0 at the short-circuit point, 0.25 at the open-circuit point).' }),
      U.el('li', { text: 'Stub matching (real Z0): find d where Re y = 1 (cos φ = −|Γ|), then a shunt stub supplies −jb. Short stub: ℓ = atan2(1, b)/2π; open stub: ℓ = atan(−b)/2π (mod ½). Quarter-wave: insert at the first voltage maximum where Z = Z0·VSWR, with Zt = √(Z0·Rmax).' })));
    stage.appendChild(th);

    // Smith interaction
    const setFromEvent = e => {
      if (!geo) return; const r = smith.cv.getBoundingClientRect();
      let gx = (e.clientX - r.left - geo.cx) / geo.R, gy = -(e.clientY - r.top - geo.cy) / geo.R;
      const m = Math.hypot(gx, gy); if (m > 0.9995) { gx *= 0.9995 / m; gy *= 0.9995 / m; }
      setLoadFromGamma(C(gx, gy));
    };
    let drag = false;
    smith.cv.addEventListener('pointerdown', e => { drag = true; try { smith.cv.setPointerCapture(e.pointerId); } catch (x) { /* ignore */ } setFromEvent(e); smith.cv.focus(); });
    smith.cv.addEventListener('pointermove', e => { if (drag) setFromEvent(e); });
    const up = () => { drag = false; };
    smith.cv.addEventListener('pointerup', up); smith.cv.addEventListener('pointercancel', up);
    smith.cv.addEventListener('keydown', e => {
      const k = { ArrowLeft: [-1, 0], ArrowRight: [1, 0], ArrowUp: [0, 1], ArrowDown: [0, -1] }[e.key]; if (!k) return;
      e.preventDefault(); const c = compute(); let G = c.GL; if (!c.valid) return;
      G = C(G.re + 0.02 * k[0], G.im + 0.02 * k[1]); const m = G.abs(); if (m > 0.9995) G = G.scale(0.9995 / m);
      setLoadFromGamma(G);
    });
    smith.onResize(() => { if (active) drawSmith(); }); vi.onResize(() => { if (active) drawVI(); });

    FSP.state.bind('tline', {
      get: () => ({ z0r: S.z0r, z0i: S.z0i, zlr: S.zlr, zli: S.zli, f: S.f, vf: S.vf, len: S.len, unit: S.unit, loss: S.loss, match: S.match }),
      set: p => {
        const num = (k, lo, hi) => { if (p && p[k] !== undefined) { const v = Number(p[k]); if (Number.isFinite(v) && v >= lo && v <= hi) S[k] = v; } };
        num('z0r', 1e-3, 1e6); num('z0i', -1e6, 1e6); num('zlr', -1e15, 1e15); num('zli', -1e15, 1e15); num('f', 1e6, 1e10); num('vf', 0.3, 1); num('len', 0, 1e4); num('loss', 0, 2);
        if (p && UNITS.indexOf(p.unit) >= 0) S.unit = p.unit;
        if (p && MATCH.indexOf(p.match) >= 0) S.match = p.match;
        if (ui) syncAll();
      },
    });
    syncAll(); recalc();
  }

  function setLoadFromGamma(G) {
    const c = compute(); const z = M.gammaToZ(G), ZL = c.Z0.mul(z);
    if (!ZL.isFinite()) return;
    S.zlr = +ZL.re.toPrecision(5); S.zli = +ZL.im.toPrecision(5); syncLoad(); changed();
  }
  function syncLoad() { ui.fields.zlr.show(); ui.fields.zli.show(); }
  function syncAll() {
    ['z0r', 'z0i', 'zlr', 'zli'].forEach(k => ui.fields[k].show());
    ui.f.set(S.f, true); ui.vf.set(S.vf, true); ui.loss.set(S.loss, true); ui.showLen();
    const sel = ui.matchOut.parentNode.querySelector('select'); if (sel) sel.value = S.match;
  }
  function changed() { recalc(); FSP.state.touch(); }
  function recalc() { if (!ui) return; if (ui.showLen) ui.showLen(); updateText(); if (active) { drawSmith(); drawVI(); } }

  /* ----- text outputs ----- */
  function row(a, b) { return FSP.ui.el('div', null, FSP.ui.el('span', { text: a }), FSP.ui.el('span', { text: b })); }
  function updateText() {
    const c = compute(), hud = ui.hud; hud.textContent = '';
    ui.lamNote.textContent = 'λ = ' + FSP.fmtEng(c.lam, 'm') + '   l = ' + fnum(S.len, 4) + ' λ = ' + FSP.fmtEng(S.len * c.lam, 'm') + ' (' + fnum(S.len * 360, 1) + '°)';
    if (!c.valid) {
      hud.appendChild(row('Status', 'Z0 must have Re > 0 and ZL must be finite'));
      ui.matchOut.textContent = '—'; ui.work.set('Invalid inputs.'); return;
    }
    const vs = M.vswr(c.GL), rl = M.returnLossDb(c.GL), ml = M.mismatchLossDb(c.GL), yin = c.Zin.abs() > 1e9 ? null : C(1, 0).div(c.Zin);
    [['Γ_L', fPolar(c.GL)], ['Γ_L (re, im)', fnum(c.GL.re, 4) + (c.GL.im < 0 ? ' − j' : ' + j') + fnum(Math.abs(c.GL.im), 4)],
      ['VSWR', vs === INF ? '∞' : fnum(vs, 3)], ['Return loss', rl === INF ? '∞ dB' : fnum(rl, 3) + ' dB'], ['Mismatch loss', ml === INF ? '∞ dB' : fnum(ml, 4) + ' dB'],
      ['Z_in(l)', fC(c.Zin, 3) + ' Ω'], ['Γ_in(l)', fPolar(c.Gin)], ['Y_in(l)', yin ? fC(yin.scale(1), 5) + ' S' : '0 S'],
      ['Load on scale', fLam(M.wtg(c.GL)) + ' WTG'], ['Input on scale', fLam(M.wtg(c.Gin)) + ' WTG']].forEach(r => hud.appendChild(row(r[0], r[1])));

    // matching helper
    const mo = [];
    if (S.match === 'qw') {
      const q = M.quarterWave(c.Z0, c.ZL);
      if (!q) mo.push('Needs real Z0 and |Γ| < 1.');
      else if (q.matched) mo.push('Already matched.');
      else {
        const chk = M.quarterWaveCheck(c.Z0, c.ZL, q);
        mo.push('Insert at d = ' + fLam(q.d) + ' (' + FSP.fmtEng(q.d * c.lam, 'm') + ') from the load,', 'where Z = ' + fnum(q.Rin, 3) + ' Ω (real).', 'Zt = ' + fnum(q.Zt, 3) + ' Ω, length λ/4 = ' + FSP.fmtEng(c.lam / 4, 'm') + '.', 'Residual |Γ| = ' + chk.abs().toExponential(1));
      }
    } else {
      const sols = M.shuntStub(c.Z0, c.ZL, S.match);
      if (sols === null) mo.push('No solution (needs real Z0 and |Γ| < 1).');
      else if (!sols.length) mo.push('Already matched.');
      else sols.forEach((s, i) => mo.push('Sol. ' + 'AB'[i] + ': d = ' + fLam(s.d) + ', stub ℓ = ' + fLam(s.ls), '      (' + FSP.fmtEng(s.d * c.lam, 'm') + ', ' + FSP.fmtEng(s.ls * c.lam, 'm') + '), b = ' + fnum(s.b, 3)));
    }
    ui.matchOut.textContent = ''; mo.forEach(l => ui.matchOut.appendChild(FSP.ui.el('div', { text: l })));

    // working
    const w = [];
    w.push('λ = v_f c / f = ' + fnum(S.vf, 3) + ' × ' + C0 + ' / ' + FSP.fmtNum(S.f, 5) + ' = ' + FSP.fmtEng(c.lam, 'm'));
    w.push('Γ_L = (ZL − Z0)/(ZL + Z0) = (' + fC(c.ZL, 3) + ' − ' + fC(c.Z0, 3) + ')/(' + fC(c.ZL, 3) + ' + ' + fC(c.Z0, 3) + ')');
    w.push('    = ' + fPolar(c.GL, 5));
    w.push('VSWR = (1+|Γ|)/(1−|Γ|) = ' + (vs === INF ? '∞' : fnum(vs, 4)) + ';  RL = −20 log10|Γ| = ' + (rl === INF ? '∞' : fnum(rl, 4)) + ' dB;  ML = −10 log10(1−|Γ|²) = ' + (ml === INF ? '∞' : fnum(ml, 4)) + ' dB');
    w.push('Γ(l) = Γ_L e^(−2γl), 2βl = ' + fnum(4 * Math.PI * S.len, 4) + ' rad' + (c.a ? ', e^(−2αl) = ' + fnum(Math.exp(-2 * c.a * S.len), 5) : ''));
    w.push('    = ' + fPolar(c.Gin, 5));
    w.push('Z_in = Z0 (1+Γ)/(1−Γ) = ' + fC(c.Zin, 4) + ' Ω');
    w.push('Normalised load z = ' + fC(c.ZL.div(c.Z0), 4) + ';  voltage max at d = ' + fLam(mod(Math.atan2(c.GL.im, c.GL.re) / (4 * Math.PI), 0.5)) + ' (|V| = ' + fnum(1 + c.m, 4) + '), min at +λ/4 (|V| = ' + fnum(1 - c.m, 4) + ')');
    ui.work.set(w);
  }

  /* ----- drawing ----- */
  const COL = { grid: '#2b3a55', gridHi: '#3d5278', text: '#8b98ab', fg: '#e6edf3', acc: '#4cc9f0', org: '#f9a03f', pink: '#f72585', ok: '#2ecc71' };
  const mono = () => (typeof getComputedStyle !== 'undefined' ? getComputedStyle(document.body).getPropertyValue('--mono') : '') || 'monospace';

  function drawSmith() {
    if (!smith) return; const { ctx, w, h } = smith.prep(); const c = compute();
    ctx.clearRect(0, 0, w, h);
    const R = Math.max(40, Math.min(w, h) / 2 - 34), cx = w / 2, cy = h / 2; geo = { cx, cy, R };
    const P = G => [cx + R * G.re, cy - R * G.im];
    ctx.lineWidth = 1; ctx.font = '10px ' + mono(); ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
    // grid (clipped to unit circle)
    ctx.save(); ctx.beginPath(); ctx.arc(cx, cy, R, 0, TWO_PI); ctx.clip();
    ctx.strokeStyle = COL.grid;
    [0.2, 0.5, 1, 2, 5].forEach(r => { ctx.strokeStyle = r === 1 ? COL.gridHi : COL.grid; ctx.beginPath(); ctx.arc(cx + R * r / (1 + r), cy, R / (1 + r), 0, TWO_PI); ctx.stroke(); });
    [0.2, 0.5, 1, 2, 5].forEach(x => [1, -1].forEach(s => { ctx.strokeStyle = x === 1 ? COL.gridHi : COL.grid; ctx.beginPath(); ctx.arc(cx + R, cy - s * R / x, R / x, 0, TWO_PI); ctx.stroke(); }));
    ctx.beginPath(); ctx.moveTo(cx - R, cy); ctx.lineTo(cx + R, cy); ctx.strokeStyle = COL.gridHi; ctx.stroke();
    ctx.restore();
    ctx.strokeStyle = COL.gridHi; ctx.beginPath(); ctx.arc(cx, cy, R, 0, TWO_PI); ctx.stroke();
    ctx.fillStyle = COL.text;
    [0.2, 0.5, 1, 2, 5].forEach(r => { const g = (r - 1) / (r + 1); ctx.fillText(String(r), cx + R * g - (r === 5 ? 0 : 0), cy + 8); });
    [0.5, 1, 2].forEach(x => { [1, -1].forEach(s => { const G = M.zToGamma(C(0, s * x)); const p = P(G); ctx.fillText((s < 0 ? '−' : '') + x, p[0] + (G.re > 0 ? 9 : -9), p[1] + (s > 0 ? -7 : 7)); }); });
    // wavelength scale
    ctx.strokeStyle = COL.text;
    for (let i = 0; i < 100; i++) {
      const wv = i / 200, phi = Math.PI - 4 * Math.PI * wv, big = i % 10 === 0, mid = i % 5 === 0, len = big ? 9 : mid ? 6 : 3;
      ctx.beginPath(); ctx.moveTo(cx + R * Math.cos(phi), cy - R * Math.sin(phi)); ctx.lineTo(cx + (R + len) * Math.cos(phi), cy - (R + len) * Math.sin(phi)); ctx.stroke();
      if (big && R > 70) { ctx.fillStyle = COL.text; ctx.fillText(wv.toFixed(2).replace(/^0/, ''), cx + (R + 19) * Math.cos(phi), cy - (R + 19) * Math.sin(phi)); }
    }
    if (R > 70) { ctx.fillStyle = COL.text; ctx.textAlign = 'left'; ctx.fillText('λ toward generator ↻', 6, 10); ctx.textAlign = 'center'; }
    if (!c.valid) { ctx.fillStyle = COL.text; ctx.fillText('Invalid inputs', cx, cy - 20); return; }
    // constant |Gamma| circle
    if (c.m > 1e-9) { ctx.strokeStyle = COL.org; ctx.setLineDash([5, 4]); ctx.beginPath(); ctx.arc(cx, cy, R * Math.min(1, c.m), 0, TWO_PI); ctx.stroke(); ctx.setLineDash([]); }
    // locus of Zin(l') for l' in [0, l]
    ctx.strokeStyle = COL.acc; ctx.lineWidth = 2; ctx.beginPath();
    const N = Math.max(2, Math.min(4000, Math.ceil(S.len * 360)));
    for (let i = 0; i <= N; i++) { const p = P(M.gammaAt(c.GL, S.len * i / N, c.a)); if (i) ctx.lineTo(p[0], p[1]); else ctx.moveTo(p[0], p[1]); }
    ctx.stroke(); ctx.lineWidth = 1;
    // stub / transformer points (voltage max for qw)
    let sols = [];
    if (S.match === 'qw') { const q = M.quarterWave(c.Z0, c.ZL); if (q && !q.matched) sols = [{ Gd: M.gammaAt(c.GL, q.d, 0), tag: 'Q' }]; }
    else { const s = M.shuntStub(c.Z0, c.ZL, S.match); if (s) sols = s.map((x, i) => ({ Gd: x.Gd, tag: 'AB'[i] })); }
    sols.forEach(s => { const p = P(s.Gd); ctx.strokeStyle = COL.ok; ctx.beginPath(); ctx.moveTo(p[0], p[1] - 6); ctx.lineTo(p[0] + 6, p[1]); ctx.lineTo(p[0], p[1] + 6); ctx.lineTo(p[0] - 6, p[1]); ctx.closePath(); ctx.stroke(); ctx.fillStyle = COL.ok; ctx.fillText(s.tag, p[0] + 12, p[1] - 7); });
    // scale pointers
    [[c.GL, COL.org], [c.Gin, COL.pink]].forEach(([G, col]) => { const phi = Math.atan2(G.im, G.re), q = G.abs() < 1e-9 ? 0 : 1; if (!q) return; ctx.strokeStyle = col; ctx.lineWidth = 2; ctx.beginPath(); ctx.moveTo(cx + R * Math.cos(phi), cy - R * Math.sin(phi)); ctx.lineTo(cx + (R + 12) * Math.cos(phi), cy - (R + 12) * Math.sin(phi)); ctx.stroke(); ctx.lineWidth = 1; });
    // markers
    const pl = P(c.GL), pi = P(c.Gin);
    ctx.fillStyle = COL.org; ctx.beginPath(); ctx.arc(pl[0], pl[1], 5, 0, TWO_PI); ctx.fill();
    ctx.fillStyle = COL.pink; ctx.beginPath(); ctx.arc(pi[0], pi[1], 4, 0, TWO_PI); ctx.fill();
    ctx.textAlign = 'left'; ctx.fillStyle = COL.org; ctx.fillText('ZL', pl[0] + 8, pl[1] - 8); ctx.fillStyle = COL.pink; ctx.fillText('Zin', pi[0] + 8, pi[1] + 9);
    ctx.textAlign = 'center';
  }

  function drawVI() {
    if (!vi) return; const { ctx, w, h } = vi.prep(); const c = compute();
    ctx.clearRect(0, 0, w, h);
    const L = 40, Rm = 12, T = 24, B = 34, pw = w - L - Rm, ph = h - T - B;
    ctx.font = '11px ' + mono(); ctx.textBaseline = 'middle';
    if (!c.valid) { ctx.fillStyle = COL.text; ctx.textAlign = 'center'; ctx.fillText('Invalid inputs', w / 2, h / 2); return; }
    const dmax = Math.max(S.len, 0.5), N = 500, V = [], I = [], YCAP = 1e6;
    // On a lossy line the forward wave (V+ = 1 at the load) grows as e^{alpha d} toward the generator and can overflow
    // for long lines; clamp to a finite cap so the axis loops below are always bounded (an unbounded ymax froze the page).
    let ymax = 1, clipped = false;
    const fin = v => (Number.isFinite(v) && v <= YCAP) ? v : (clipped = true, YCAP);
    for (let i = 0; i <= N; i++) { const o = M.vi(c.GL, dmax * i / N, c.a), v = fin(o.V.abs()), j = fin(o.ZI.abs()); V.push(v); I.push(j); ymax = Math.max(ymax, v, j); }
    ymax = Math.min(YCAP, Math.ceil(ymax * 10 * 1.05) / 10);
    const X = d => L + pw * d / dmax, Y = v => T + ph * (1 - v / ymax);
    const niceStep = raw => { const p = Math.pow(10, Math.floor(Math.log10(raw))), f = raw / p; return (f < 1.5 ? 1 : f < 3.5 ? 2 : f < 7.5 ? 5 : 10) * p; };
    ctx.strokeStyle = COL.grid; ctx.fillStyle = COL.text; ctx.textAlign = 'right';
    const ystep = ymax > 1.2 ? niceStep(ymax / 4) : 0.25;                      // at most ~8 grid lines whatever the range
    for (let v = 0; v <= ymax + 1e-9; v += ystep) { ctx.beginPath(); ctx.moveTo(L, Y(v)); ctx.lineTo(L + pw, Y(v)); ctx.stroke(); ctx.fillText(v >= 1e4 ? v.toExponential(1) : FSP.fmtNum(v, 3), L - 4, Y(v)); }
    ctx.textAlign = 'center';
    const xstep = dmax > 1.6 ? niceStep(dmax / 8) : dmax > 0.8 ? 0.25 : 0.125;
    for (let d = 0; d <= dmax + 1e-9; d += xstep) { ctx.beginPath(); ctx.moveTo(X(d), T); ctx.lineTo(X(d), T + ph); ctx.stroke(); ctx.fillText(FSP.fmtNum(d, 4), X(d), T + ph + 12); }
    if (clipped) { ctx.fillStyle = COL.org; ctx.textAlign = 'right'; ctx.fillText('amplitude clipped at ' + YCAP.toExponential(0) + ' (loss × length too large)', L + pw - 4, T + 10); ctx.fillStyle = COL.text; ctx.textAlign = 'center'; }
    ctx.fillText('distance from load d (λ)', L + pw / 2, h - 8);
    const line = (arr, col) => { ctx.strokeStyle = col; ctx.lineWidth = 2; ctx.beginPath(); arr.forEach((v, i) => { const x = X(dmax * i / N), y = Y(v); if (i) ctx.lineTo(x, y); else ctx.moveTo(x, y); }); ctx.stroke(); ctx.lineWidth = 1; };
    line(V, COL.acc); line(I, COL.org);
    // generator end marker
    if (S.len < dmax - 1e-9) { ctx.strokeStyle = COL.pink; ctx.setLineDash([4, 3]); ctx.beginPath(); ctx.moveTo(X(S.len), T); ctx.lineTo(X(S.len), T + ph); ctx.stroke(); ctx.setLineDash([]); ctx.fillStyle = COL.pink; ctx.textAlign = 'left'; ctx.fillText('input end (l)', X(S.len) + 4, T + 8); }
    ctx.fillStyle = COL.text; ctx.textAlign = 'left'; ctx.fillText('load', L + 3, T + ph - 8);
    ctx.fillStyle = COL.acc; ctx.fillText('|V| / |V+|', L + 4, 10); ctx.fillStyle = COL.org; ctx.fillText('|I|·|Z0| / |V+|', L + 90, 10);
  }

  FSP.registerTab({
    id: 'tline', title: 'Transmission Lines',
    init,
    activate() { active = true; recalc(); },
    deactivate() { active = false; },
  });
})();
