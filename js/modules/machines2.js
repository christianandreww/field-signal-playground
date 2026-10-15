/* rotating UI. Numerical implementation: js/math/machines2.js. */
(function () {
  'use strict';
  const { CD, DC, IM, PI, SC, TP, cabs, cx, fin } = FSP.math.rotating._shared;
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
