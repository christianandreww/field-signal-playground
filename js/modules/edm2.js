/* edm2 UI. Numerical implementation: js/math/edm2.js. */
(function () {
  'use strict';
  const E = FSP.math.edm2;
  const { MU0, NAMES, PI, SH, arrowOf, cap, dropLines, fin, fx } = E._shared;
  /* =====================================================================
     UI
     ===================================================================== */
  const css = n => { try { return getComputedStyle(document.documentElement).getPropertyValue(n).trim() || null; } catch (e) { return null; } };
  function theme() {
    return { bg: css('--panel2') || '#0f141f', panel: css('--panel') || '#121826', text: css('--text') || '#e6edf3', muted: css('--muted') || '#8b98ab', border: css('--border') || '#243047',
      a: css('--accent2') || '#f9a03f', b: css('--ok') || '#2ecc71', c: css('--accent') || '#4cc9f0', pink: css('--pink') || '#f72585', warn: css('--warn') || '#ffb347', bad: css('--bad') || '#ff5c5c', mono: css('--mono') || 'monospace' };
  }
  const f2 = (x, d) => (fin(x) ? x.toFixed(d === undefined ? 2 : d) : '—');
  const fE = (x, u, d) => (fin(x) ? FSP.fmtEng(x, u, d || 4) : '—');
  const tickStr = v => (Math.abs(v) < 1e-12 ? '0' : Math.abs(v) >= 1e5 || Math.abs(v) < 1e-3 ? v.toExponential(1) : String(+v.toPrecision(4)));
  function niceTicks(lo, hi, n) {
    const span = hi - lo; if (!(span > 0)) return [lo];
    const raw = span / n, mag = Math.pow(10, Math.floor(Math.log10(raw))), f = raw / mag, step = (f < 1.5 ? 1 : f < 3 ? 2 : f < 7 ? 5 : 10) * mag, out = [];
    for (let v = Math.ceil(lo / step - 1e-9) * step; v <= hi + step * 1e-9; v += step) out.push(Math.abs(v) < step * 1e-9 ? 0 : +v.toPrecision(12));
    return out;
  }
  function withAlpha(col, a) {
    const m = /^#([0-9a-f]{6})$/i.exec(col || ''); if (!m) return col;
    const n = parseInt(m[1], 16); return 'rgba(' + (n >> 16) + ',' + ((n >> 8) & 255) + ',' + (n & 255) + ',' + a + ')';
  }
  // generic plot. o: {xmin,xmax,ymin,ymax,y2min,y2max,xlabel,ylabel,y2label,series:[{x,y,color,width,dash,axis}],fills:[{x,y,color}],bars:[{x0,x1,y,color}],
  //   points:[{x,y,label,color,below,axis}],vlines:[{x,label,color}],hlines:[{y,label,color}],legend:[{text,color}]}
  function drawPlot(c, o) {
    const g = c.prep(), ctx = g.ctx, w = g.w, h = g.h, T = theme(), y2 = o.y2min !== undefined;
    const L = 52, R = y2 ? 42 : 10, Tp = 10, B = 32, pw = w - L - R, ph = h - Tp - B;
    ctx.clearRect(0, 0, w, h); ctx.fillStyle = T.bg; ctx.fillRect(0, 0, w, h);
    if (pw < 20 || ph < 20 || !(o.xmax > o.xmin) || !(o.ymax > o.ymin)) return null;
    const X = x => L + (x - o.xmin) / (o.xmax - o.xmin) * pw, Y = y => Tp + (1 - (y - o.ymin) / (o.ymax - o.ymin)) * ph;
    const Y2 = y => Tp + (1 - (y - o.y2min) / (o.y2max - o.y2min)) * ph;
    ctx.font = '11px ' + T.mono; ctx.lineWidth = 1; ctx.textBaseline = 'middle'; ctx.textAlign = 'right';
    niceTicks(o.ymin, o.ymax, 5).forEach(v => { const yy = Y(v); ctx.strokeStyle = T.border; ctx.beginPath(); ctx.moveTo(L, yy); ctx.lineTo(L + pw, yy); ctx.stroke(); ctx.fillStyle = T.muted; ctx.fillText(tickStr(v), L - 4, yy); });
    if (y2) { ctx.textAlign = 'left'; niceTicks(o.y2min, o.y2max, 5).forEach(v => { ctx.fillStyle = T.muted; ctx.fillText(tickStr(v), L + pw + 4, Y2(v)); }); }
    ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    niceTicks(o.xmin, o.xmax, Math.max(3, Math.floor(pw / 70))).forEach(v => { const xx = X(v); ctx.strokeStyle = T.border; ctx.beginPath(); ctx.moveTo(xx, Tp); ctx.lineTo(xx, Tp + ph); ctx.stroke(); ctx.fillStyle = T.muted; ctx.fillText(tickStr(v), xx, Tp + ph + 3); });
    ctx.strokeStyle = T.muted; if (o.ymin < 0 && o.ymax > 0) { ctx.beginPath(); ctx.moveTo(L, Y(0)); ctx.lineTo(L + pw, Y(0)); ctx.stroke(); }
    ctx.fillStyle = T.muted; ctx.textAlign = 'center'; ctx.textBaseline = 'bottom'; ctx.fillText(o.xlabel || '', L + pw / 2, h - 1);
    ctx.save(); ctx.translate(10, Tp + ph / 2); ctx.rotate(-PI / 2); ctx.textBaseline = 'top'; ctx.fillText(o.ylabel || '', 0, 0); ctx.restore();
    if (y2) { ctx.save(); ctx.translate(w - 3, Tp + ph / 2); ctx.rotate(-PI / 2); ctx.textBaseline = 'bottom'; ctx.fillText(o.y2label || '', 0, 0); ctx.restore(); }
    ctx.save(); ctx.beginPath(); ctx.rect(L, Tp, pw, ph); ctx.clip();
    (o.bars || []).forEach(b => { if (!fin(b.y)) return; ctx.fillStyle = b.color; const y0 = Y(Math.max(o.ymin, 0)), y1 = Y(b.y); ctx.fillRect(X(b.x0) + 1, Math.min(y0, y1), Math.max(1, X(b.x1) - X(b.x0) - 2), Math.abs(y1 - y0)); });
    (o.fills || []).forEach(f => { ctx.fillStyle = f.color; ctx.beginPath(); f.x.forEach((xv, i) => { const px = X(xv), py = Y(f.y[i]); if (i) ctx.lineTo(px, py); else ctx.moveTo(px, py); }); ctx.closePath(); ctx.fill(); });
    (o.vlines || []).forEach(v => { if (!fin(v.x)) return; ctx.strokeStyle = v.color || T.muted; ctx.setLineDash([4, 4]); ctx.beginPath(); ctx.moveTo(X(v.x), Tp); ctx.lineTo(X(v.x), Tp + ph); ctx.stroke(); ctx.setLineDash([]); });
    (o.hlines || []).forEach(v => { if (!fin(v.y)) return; ctx.strokeStyle = v.color || T.muted; ctx.setLineDash([4, 4]); ctx.beginPath(); ctx.moveTo(L, Y(v.y)); ctx.lineTo(L + pw, Y(v.y)); ctx.stroke(); ctx.setLineDash([]); });
    (o.series || []).forEach(s => {
      const fy = s.axis === 'y2' ? Y2 : Y; ctx.strokeStyle = s.color; ctx.lineWidth = s.width || 2; ctx.setLineDash(s.dash || []); ctx.beginPath(); let pen = false;
      for (let i = 0; i < s.x.length; i++) { const xv = s.x[i], yv = s.y[i]; if (!fin(xv) || !fin(yv)) { pen = false; continue; } const px = X(xv), py = Math.max(-1e4, Math.min(1e4, fy(yv))); if (pen) ctx.lineTo(px, py); else { ctx.moveTo(px, py); pen = true; } }
      ctx.stroke(); ctx.setLineDash([]);
    });
    ctx.restore();
    ctx.font = '11px ' + T.mono;
    (o.vlines || []).forEach((v, i) => { if (!fin(v.x) || !v.label) return; const xx = X(v.x), right = xx < L + pw * 0.6; ctx.fillStyle = v.color || T.muted; ctx.textAlign = right ? 'left' : 'right'; ctx.textBaseline = 'top'; ctx.fillText(v.label, xx + (right ? 4 : -4), Tp + 3 + (v.row || 0) * 13); });
    (o.hlines || []).forEach(v => { if (!fin(v.y) || !v.label) return; ctx.fillStyle = v.color || T.muted; ctx.textAlign = 'left'; ctx.textBaseline = 'bottom'; ctx.fillText(v.label, L + 4, Y(v.y) - 2); });
    (o.points || []).forEach(p => {
      if (!fin(p.x) || !fin(p.y)) return; const px = X(p.x), py = (p.axis === 'y2' ? Y2 : Y)(p.y);
      if (px < L - 1 || px > L + pw + 1 || py < Tp - 1 || py > Tp + ph + 1) return;
      ctx.fillStyle = p.color; ctx.strokeStyle = T.bg; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.arc(px, py, 5, 0, 2 * PI); ctx.fill(); ctx.stroke();
      if (p.label) { const right = px < L + pw * 0.62; ctx.textAlign = right ? 'left' : 'right'; ctx.textBaseline = p.below ? 'top' : 'bottom'; ctx.fillStyle = T.text; ctx.fillText(p.label, px + (right ? 8 : -8), py + (p.below ? 6 : -6)); }
    });
    (o.legend || []).forEach((lg, i) => { ctx.fillStyle = lg.color; ctx.textAlign = 'right'; ctx.textBaseline = 'top'; ctx.fillText(lg.text, L + pw - 4, Tp + 4 + i * 14); });
    return { X, Y, L, Tp, pw, ph };
  }
  function hudSet(el, rows) {
    while (el.firstChild) el.removeChild(el.firstChild);
    rows.forEach(r => el.appendChild(FSP.ui.el('div', null, FSP.ui.el('span', { text: r[0] }), FSP.ui.el('span', { text: r[1] }))));
  }
  const wrapCanvas = (parent, label, height) => { const wrap = FSP.ui.el('div', { class: 'canvas-wrap' }); parent.appendChild(wrap); const c = FSP.ui.canvas(wrap, { height }); c.cv.setAttribute('role', 'img'); c.cv.setAttribute('aria-label', label); return c; };
  function layout(root) { const lay = FSP.ui.el('div', { class: 'layout' }); root.appendChild(lay); const ctl = FSP.ui.el('div', { class: 'controls' }), stage = FSP.ui.el('div', { class: 'stage' }); lay.appendChild(ctl); lay.appendChild(stage); return { ctl, stage }; }
  // defs: [{k,l,min,max,v,u,log,step}]
  function makeSliders(parent, defs, st, onChange) {
    defs.forEach(d => { st[d.k] = FSP.ui.slider(parent, { label: d.l, min: d.min, max: d.max, step: d.step, value: d.v, unit: d.u, log: d.log, onInput: () => { onChange(); FSP.state.touch(); } }); });
  }
  const sel = (parent, label, opts, v, onChange) => FSP.ui.select(parent, label, opts, v, () => { onChange(); FSP.state.touch(); });
  const rowOf = selEl => selEl.parentNode;
  function stateIO(st, sels, extra) {
    return {
      get() { const o = {}; Object.keys(st).forEach(k => { o[k] = st[k].get(); }); Object.keys(sels).forEach(k => { o[k] = sels[k].value; }); if (extra) Object.assign(o, extra.get()); return o; },
      set(o) {
        o = o || {};
        Object.keys(st).forEach(k => { if (o[k] !== undefined) { const v = parseFloat(o[k]); if (fin(v)) st[k].set(v, true); } });
        Object.keys(sels).forEach(k => { if (o[k] !== undefined && Array.prototype.some.call(sels[k].options, op => op.value === o[k])) sels[k].value = o[k]; });
        if (extra) extra.set(o);
      },
    };
  }
  function msgBox(parent) { const m = FSP.ui.el('div', { class: 'msg bad', hidden: '' }); parent.appendChild(m); return { el: m, show(t) { m.textContent = t; m.hidden = !t; } }; }
  function arrowHead(ctx, x0, y0, x1, y1, s) {
    const a = Math.atan2(y1 - y0, x1 - x0); ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1, y1); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(x1, y1); ctx.lineTo(x1 - s * Math.cos(a - 0.4), y1 - s * Math.sin(a - 0.4)); ctx.lineTo(x1 - s * Math.cos(a + 0.4), y1 - s * Math.sin(a + 0.4)); ctx.closePath(); ctx.fill();
  }

  /* ---------- shared B–H table (used by sub-topics 1 and 2) ---------- */
  const shared = { text: E.tableToText(E.DEFAULT_TABLE), parsed: E.parseTable(E.tableToText(E.DEFAULT_TABLE)), listeners: [] };
  shared.tableMat = () => (shared.parsed.ok ? E.makeMaterial({ mode: 'table', H: shared.parsed.H, B: shared.parsed.B }) : { ok: false });
  shared.setText = (t, silent) => { shared.text = String(t); shared.parsed = E.parseTable(shared.text); if (!silent) shared.listeners.forEach(f => f()); };

  /* ---------- 1. multi-limb magnetic circuit ---------- */
  function buildMC(root) {
    const st = {}, ui = {}, sels = {}, { ctl, stage } = layout(root);
    const fsP = FSP.ui.fieldset(ctl, 'Problem');
    sels.mode = sel(fsP, 'Solve for', [['I', 'coil current (given Φ or B)'], ['phi', 'fluxes (given currents)']], 'I', update);
    sels.g = sel(fsP, 'Given in', [['0', 'left limb'], ['1', 'centre limb'], ['2', 'right limb']], '2', update);
    sels.gq = sel(fsP, 'Given', [['B', 'flux density B'], ['phi', 'flux Φ']], 'B', update);
    makeSliders(fsP, [{ k: 'gB', l: 'B given', min: 0, max: 2.5, v: 0.75, step: 0.0001, u: 'T' }, { k: 'gPhi', l: 'Φ given', min: 0, max: 50, v: 0.3, step: 0.0001, u: 'mWb' }], st, update);
    sels.gd = sel(fsP, 'Direction', [['-1', '↓ down'], ['1', '↑ up']], '-1', update);
    sels.x = sel(fsP, 'Unknown coil', [['0', 'on left limb'], ['1', 'on centre limb'], ['2', 'on right limb']], '1', update);
    const fsM = FSP.ui.fieldset(ctl, 'Core material & geometry');
    sels.mat = sel(fsM, 'Material', [['table', 'B–H table (linear interp.)'], ['const', 'constant μr']], 'table', update);
    makeSliders(fsM, [{ k: 'mur', l: 'μr', min: 1, max: 100000, v: 1000, u: '', log: true }, { k: 'd', l: 'depth', min: 0.1, max: 50, v: 2, u: 'cm', log: true }], st, update);
    sels.fr = sel(fsM, 'Fringing', [['0', 'none (A_g = A)'], ['1', 'A_g = (a+l_g)(b+l_g)']], '0', update);
    ui.tblWrap = FSP.ui.el('div'); fsM.appendChild(ui.tblWrap);
    ui.tblWrap.appendChild(FSP.ui.el('div', { class: 'note', text: 'B–H table, one "H B" pair per row (H in A/m, B in T). (0,0) is implied; straight lines between points.' }));
    ui.tbl = FSP.ui.el('textarea', { rows: '8', spellcheck: 'false', 'aria-label': 'B–H table' }); ui.tbl.value = shared.text;
    ui.tbl.style.cssText = 'width:100%;box-sizing:border-box;background:var(--panel2);color:var(--text);border:1px solid var(--border);border-radius:6px;font:12px var(--mono);padding:4px';
    ui.tbl.addEventListener('change', () => { shared.setText(ui.tbl.value); update(); FSP.state.touch(); }); ui.tblWrap.appendChild(ui.tbl);
    const rr = FSP.ui.el('div', { class: 'row' }); ui.tblWrap.appendChild(rr);
    FSP.ui.button(rr, 'Past-paper table', () => { ui.tbl.value = E.tableToText(E.DEFAULT_TABLE); shared.setText(ui.tbl.value); update(); FSP.state.touch(); }).className = 'btn';
    ui.tblMsg = msgBox(ui.tblWrap);
    const limbDefs = [
      { on: '1', l: 0.3, w: 2, g: 0, N: 0, I: 0, s: '1' }, { on: '1', l: 0.1, w: 4, g: 0, N: 500, I: 0.5, s: '1' }, { on: '1', l: 0.3, w: 2, g: 0.2, N: 0, I: 0, s: '1' }];
    ui.limbRows = [];
    limbDefs.forEach((D, k) => {
      const fs = FSP.ui.fieldset(ctl, cap(NAMES[k]) + ' limb'), p = SH[k];
      sels['on' + p] = sel(fs, 'Present', [['1', 'yes'], ['0', 'no (removed)']], D.on, update);
      makeSliders(fs, [
        { k: 'l' + p, l: 'mean l', min: 0.01, max: 2, v: D.l, u: 'm', log: true }, { k: 'w' + p, l: 'width a', min: 0.1, max: 50, v: D.w, u: 'cm', log: true },
        { k: 'g' + p, l: 'air gap', min: 0, max: 10, v: D.g, step: 0.01, u: 'mm' }, { k: 'N' + p, l: 'turns N', min: 0, max: 5000, v: D.N, step: 1, u: '' },
        { k: 'I' + p, l: 'current I', min: 0, max: 100, v: D.I, step: 0.0001, u: 'A' }], st, update);
      sels['s' + p] = sel(fs, 'Coil drives', [['1', 'flux ↑'], ['-1', 'flux ↓']], D.s, update);
      ui.limbRows.push(fs);
    });
    ctl.appendChild(FSP.ui.el('div', { class: 'note', text: 'Mean length l = path along the limb between the two junctions (include its share of the yokes). Iron length = l − l_g. Area A = width a × depth b.' }));
    ui.msg = msgBox(ctl);
    ui.cv = wrapCanvas(stage, 'Three-limb magnetic circuit with flux arrows', 320); ui.cv.onResize(update);
    ui.hud = FSP.ui.el('div', { class: 'hud' }); stage.appendChild(ui.hud);
    ui.warn = FSP.ui.el('div', { class: 'msg warn', hidden: '' }); stage.appendChild(ui.warn);
    ui.work = FSP.ui.working(stage); ui.work.el.open = true;
    shared.listeners.push(() => { if (ui.tbl.value !== shared.text) ui.tbl.value = shared.text; });

    function circuit() {
      const fringe = sels.fr.value === '1', mat = sels.mat.value === 'table' ? shared.tableMat() : E.makeMaterial({ mode: 'const', mur: st.mur.get() });
      const limbs = SH.map(p => E.makeLimb({ on: sels['on' + p].value === '1', l: st['l' + p].get(), w: st['w' + p].get() / 100, d: st.d.get() / 100, lg: st['g' + p].get() / 1000, fringe, N: Math.round(st['N' + p].get()), I: st['I' + p].get(), sense: +sels['s' + p].value }));
      return { mat, limbs };
    }
    let last = null;
    function update() {
      if (!ui.cv) return;
      const modeI = sels.mode.value === 'I', g = +sels.g.value, x = +sels.x.value;
      [sels.g, sels.gq, sels.gd, sels.x].forEach(s => { rowOf(s).hidden = !modeI; });
      st.gB.el.hidden = !modeI || sels.gq.value !== 'B'; st.gPhi.el.hidden = !modeI || sels.gq.value !== 'phi';
      st.mur.el.hidden = sels.mat.value !== 'const'; ui.tblWrap.hidden = sels.mat.value !== 'table';
      ui.tblMsg.show(shared.parsed.ok ? '' : shared.parsed.msg);
      SH.forEach((p, k) => { const on = sels['on' + p].value === '1'; ['l', 'w', 'g', 'N', 'I'].forEach(q => { st[q + p].el.hidden = !on || (q === 'I' && modeI && k === x); }); rowOf(sels['s' + p]).hidden = !on; });
      const c = circuit(); let r;
      if (modeI) {
        const L = c.limbs[g], dir = +sels.gd.value, phiG = sels.gq.value === 'B' ? dir * st.gB.get() * L.A : dir * st.gPhi.get() * 1e-3;
        r = !(Math.abs(phiG) > 0) ? { ok: false, msg: 'Enter a non-zero given flux / flux density.' } : E.mcSolveForCurrent(c, g, phiG, x);
      } else r = E.mcSolveCurrents(c);
      last = { c, r };
      ui.msg.show(r.ok ? '' : r.msg);
      drawCore(ui.cv, c, r.ok ? r : null, modeI ? x : -1);
      if (!r.ok) { hudSet(ui.hud, []); ui.work.set(r.msg); ui.warn.hidden = true; return; }
      const rows = [];
      if (modeI) rows.push(['coil current I (' + NAMES[x] + ')', fx(r.I, 4) + ' A'], ['coil MMF N·I', fx(Math.abs(r.sNI), 5) + ' At']);
      r.P.forEach(k => { const d = r.drops[k]; rows.push(['Φ ' + NAMES[k], fx(Math.abs(d.phi) * 1e3, 4) + ' mWb ' + arrowOf(d.phi)], ['B / H ' + NAMES[k], fx(Math.abs(d.B), 4) + ' T / ' + fx(Math.abs(d.H), 4) + ' A/m']); });
      rows.push(['U (yoke-to-yoke MMF)', fx(r.U, 5) + ' At'], ['ΣΦ check', fx(r.kcl, 3) + ' Wb']);
      if (!modeI) rows.push(['method', r.method], ['iterations', String(r.iter)]);
      hudSet(ui.hud, rows);
      ui.warn.hidden = !r.outside; ui.warn.textContent = r.outside ? 'A flux density exceeds the last B–H table point: H was extrapolated along the last segment, so the answer is outside the given data.' : '';
      if (modeI) ui.work.set(r.steps);
      else {
        const lines = ['Convention: Φ and coil MMF positive ↑ in every limb.  For each limb: s·N·I − D(Φ) = U (common yoke-to-yoke MMF);  flux law ΣΦ = 0.',
          'Coil MMFs: ' + r.P.map(k => SH[k] + ': ' + fx(r.S[k]) + ' At').join(',  ')];
        if (r.oth && r.oth.length === 2) lines.push('Unknowns Φ_L, Φ_R with Φ_C = −(Φ_L + Φ_R). Window equations (Kirchhoff MMF law):', '  left window:  [S_L − D_L(Φ_L)] − [S_C − D_C(Φ_C)] = 0', '  right window: [S_R − D_R(Φ_R)] − [S_C − D_C(Φ_C)] = 0');
        else lines.push('Two limbs: single loop, Φ_a = −Φ_b, S_a − D_a(Φ_a) = S_b − D_b(Φ_b).');
        lines.push('D(Φ) is nonlinear (B–H table) → solved by ' + r.method + (r.hist ? '; residual history (At): ' + r.hist.map(v => fx(v, 3)).join(' → ') : '') + '.');
        lines.push('Solution check (forward chain with the solved fluxes):');
        r.P.forEach(k => dropLines(c, k, r.phi[k]).forEach(s => lines.push(s)));
        lines.push('U per limb: ' + r.P.map(k => SH[k] + ': ' + fx(r.S[k]) + ' − (' + fx(r.drops[k].F) + ') = ' + fx(r.S[k] - r.drops[k].F)).join(';  ') + ' At');
        lines.push('ΣΦ = ' + fx(r.kcl, 3) + ' Wb (flux law satisfied)');
        ui.work.set(lines);
      }
    }
    // core drawing
    function drawCore(cv, c, r, xUnknown) {
      const gg = cv.prep(), ctx = gg.ctx, w = gg.w, h = gg.h, T = theme();
      ctx.clearRect(0, 0, w, h); ctx.fillStyle = T.bg; ctx.fillRect(0, 0, w, h); if (w < 120) return;
      const x0 = 12, x1 = w - 12, y0 = 46, y1 = h - 46, wmax = Math.max.apply(null, c.limbs.map(L => L.w || 0)) || 1;
      const lw = c.limbs.map(L => Math.max(14, Math.min(Math.min(48, w * 0.1), (L.w / wmax) * Math.min(48, w * 0.1))));
      const yk = Math.max(16, Math.min(30, Math.max.apply(null, lw) * 0.7));
      const cxs = [x0 + lw[0] / 2, (x0 + x1) / 2, x1 - lw[2] / 2];
      ctx.fillStyle = withAlpha(T.muted, 0.35); ctx.strokeStyle = T.muted; ctx.lineWidth = 1;
      ctx.fillRect(x0, y0, x1 - x0, yk); ctx.fillRect(x0, y1 - yk, x1 - x0, yk);
      const gapPx = 7, yGap = y0 + yk + (y1 - y0 - 2 * yk) * 0.3;
      c.limbs.forEach((L, k) => {
        const xl = cxs[k] - lw[k] / 2;
        if (!L.on) { ctx.setLineDash([4, 4]); ctx.strokeRect(xl, y0 + yk, lw[k], y1 - y0 - 2 * yk); ctx.setLineDash([]); return; }
        ctx.fillStyle = withAlpha(T.muted, 0.35); ctx.fillRect(xl, y0 + yk, lw[k], y1 - y0 - 2 * yk);
        if (L.lg > 0) { ctx.fillStyle = T.bg; ctx.fillRect(xl - 1, yGap - gapPx / 2, lw[k] + 2, gapPx); ctx.fillStyle = T.a; ctx.font = '10px ' + T.mono; ctx.textAlign = k === 2 ? 'right' : 'left'; ctx.textBaseline = 'middle'; ctx.fillText('g ' + fx(L.lg * 1e3, 3) + 'mm', k === 2 ? xl - 4 : xl + lw[k] + 4, yGap); }
      });
      ctx.strokeStyle = T.muted; ctx.strokeRect(x0, y0, x1 - x0, y1 - y0);
      const yc0 = y0 + yk + (y1 - y0 - 2 * yk) * 0.52, yc1 = y0 + yk + (y1 - y0 - 2 * yk) * 0.86;
      ctx.font = '10px ' + T.mono;
      c.limbs.forEach((L, k) => {
        if (!L.on) return; const cxk = cxs[k];
        // coil
        if (L.N > 0) {
          const cw = lw[k] + 12; ctx.strokeStyle = k === xUnknown ? T.pink : T.a; ctx.lineWidth = 1.6;
          for (let yy = yc0; yy <= yc1; yy += 5) { ctx.beginPath(); ctx.moveTo(cxk - cw / 2, yy); ctx.lineTo(cxk + cw / 2, yy + 2.5); ctx.stroke(); }
          ctx.fillStyle = k === xUnknown ? T.pink : T.a; const ax = k === 2 ? cxk - cw / 2 - 7 : cxk + cw / 2 + 7, up = L.sense > 0;
          arrowHead(ctx, ax, up ? yc1 : yc0, ax, up ? yc0 : yc1, 6);
          ctx.textAlign = 'center'; ctx.textBaseline = 'top';
          const cur = r && k === xUnknown ? r.I : L.I;
          ctx.fillText('N=' + L.N, cxk, y1 + 4); ctx.fillText((k === xUnknown ? 'I=' + fx(cur, 3) : 'I=' + fx(cur, 3)) + 'A', cxk, y1 + 17);
        }
        // flux arrow and label
        if (r) {
          const d = r.drops[k], f = d.phi; if (Math.abs(f) > 0) {
            ctx.strokeStyle = T.c; ctx.fillStyle = T.c; ctx.lineWidth = 2.4; const ya = y0 + yk + 6, yb = yc0 - 8;
            arrowHead(ctx, cxk, f > 0 ? yb : ya, cxk, f > 0 ? ya : yb, 8);
          }
          ctx.fillStyle = T.text; ctx.textAlign = k === 0 ? 'left' : k === 2 ? 'right' : 'center'; ctx.textBaseline = 'bottom';
          const tx = k === 0 ? x0 : k === 2 ? x1 : cxk;
          ctx.fillText('Φ' + SH[k] + '=' + fx(Math.abs(f) * 1e3, 3) + 'mWb' + arrowOf(f), tx, y0 - 15);
          ctx.fillStyle = d.outside ? T.bad : T.muted; ctx.fillText('B=' + fx(Math.abs(d.B), 3) + 'T', tx, y0 - 3);
        }
      });
    }
    const io = stateIO(st, sels, { get: () => ({ tbl: shared.text.split(/\r?\n/).join(';') }), set(o) { if (typeof o.tbl === 'string' && o.tbl.length < 4000) { shared.setText(o.tbl.split(';').join('\n'), false); ui.tbl.value = shared.text; } } });
    return { update, get: io.get, set: io.set };
  }

  /* ---------- 2. inductance, energy, force, losses ---------- */
  function buildIE(root) {
    const st = {}, sels = {}, ui = {}, { ctl, stage } = layout(root);
    const fs0 = FSP.ui.fieldset(ctl, 'View');
    sels.view = sel(fs0, 'Topic', [['core', 'gapped core: L & energy'], ['relay', 'relay / plunger force'], ['rot', 'rotary actuator torque'], ['sine', 'sinusoidal excitation & core loss']], 'core', update);
    sels.mat = sel(fs0, 'Material', [['const', 'constant μr'], ['table', 'B–H table (sub-topic 1)']], 'const', update);
    makeSliders(fs0, [{ k: 'mur', l: 'μr', min: 10, max: 100000, v: 2000, u: '', log: true }], st, update);
    ui.fsCore = FSP.ui.fieldset(ctl, 'Gapped core');
    makeSliders(ui.fsCore, [{ k: 'N', l: 'turns N', min: 1, max: 5000, v: 500, step: 1, u: '' }, { k: 'lc', l: 'l iron', min: 0.01, max: 2, v: 0.3, u: 'm', log: true },
      { k: 'A', l: 'area A', min: 0.1, max: 500, v: 10, u: 'cm²', log: true }, { k: 'lg', l: 'air gap', min: 0, max: 10, v: 1, step: 0.01, u: 'mm' },
      { k: 'i', l: 'current i', min: 0.001, max: 50, v: 1, u: 'A', log: true }], st, update);
    ui.fsRel = FSP.ui.fieldset(ctl, 'Relay (U-core + armature)');
    makeSliders(ui.fsRel, [{ k: 'rN', l: 'turns N', min: 1, max: 10000, v: 1000, step: 1, u: '' }, { k: 'rlc', l: 'l iron', min: 0.01, max: 2, v: 0.2, u: 'm', log: true },
      { k: 'rA', l: 'pole A', min: 0.1, max: 100, v: 4, u: 'cm²', log: true }, { k: 'rx', l: 'gap x', min: 0.01, max: 10, v: 1, u: 'mm', log: true },
      { k: 'ri', l: 'current i', min: 0.001, max: 20, v: 0.5, u: 'A', log: true }, { k: 'rFs', l: 'spring F', min: 0.01, max: 1000, v: 20, u: 'N', log: true },
      { k: 'rx0', l: 'open gap x0', min: 0.05, max: 10, v: 2, u: 'mm', log: true }], st, update);
    sels.ng = sel(ui.fsRel, 'Gaps', [['2', '2 (armature across both poles)'], ['1', '1 (single gap)']], '2', update);
    ui.fsRot = FSP.ui.fieldset(ctl, 'Rotary reluctance actuator');
    makeSliders(ui.fsRot, [{ k: 'Ld', l: 'L aligned', min: 0.001, max: 5, v: 0.2, u: 'H', log: true }, { k: 'Lq', l: 'L unalign.', min: 0.0005, max: 5, v: 0.05, u: 'H', log: true },
      { k: 'ti', l: 'current i', min: 0, max: 50, v: 2, step: 0.01, u: 'A' }, { k: 'th', l: 'angle θ', min: 0, max: 180, v: 30, step: 0.5, u: '°' }], st, update);
    ui.fsSine = FSP.ui.fieldset(ctl, 'Sinusoidal excitation & losses');
    makeSliders(ui.fsSine, [{ k: 'V', l: 'V rms', min: 1, max: 50000, v: 230, u: 'V', log: true }, { k: 'f', l: 'f', min: 1, max: 2000, v: 50, u: 'Hz', log: true },
      { k: 'sN', l: 'turns N', min: 1, max: 10000, v: 200, step: 1, u: '' }, { k: 'sA', l: 'area A', min: 0.1, max: 1000, v: 25, u: 'cm²', log: true },
      { k: 'vol', l: 'core vol.', min: 1, max: 1e6, v: 3000, u: 'cm³', log: true }, { k: 'n', l: 'Steinmetz n', min: 1, max: 3, v: 1.6, step: 0.01, u: '' },
      { k: 'f1', l: 'test f1', min: 1, max: 2000, v: 50, u: 'Hz', log: true }, { k: 'P1', l: 'test P1', min: 0.001, max: 1e6, v: 100, u: 'W', log: true },
      { k: 'f2', l: 'test f2', min: 1, max: 2000, v: 60, u: 'Hz', log: true }, { k: 'P2', l: 'test P2', min: 0.001, max: 1e6, v: 130, u: 'W', log: true }], st, update);
    ui.fsSine.appendChild(FSP.ui.el('div', { class: 'note', text: 'Both test readings must be at the same Bmax (V/f kept constant). Bmax used = value from V, f, N, A above.' }));
    ui.msg = msgBox(ctl);
    ui.cv = wrapCanvas(stage, 'Inductance, force or loss plot', 300); ui.cv.onResize(update);
    ui.hud = FSP.ui.el('div', { class: 'hud' }); stage.appendChild(ui.hud);
    ui.warn = FSP.ui.el('div', { class: 'msg warn', hidden: '' }); stage.appendChild(ui.warn);
    ui.work = FSP.ui.working(stage);
    shared.listeners.push(() => update());

    function material() { return sels.mat.value === 'table' ? shared.tableMat() : E.makeMaterial({ mode: 'const', mur: st.mur.get() }); }
    function update() {
      if (!ui.cv) return;
      const v = sels.view.value, T = theme();
      ui.fsCore.hidden = v !== 'core'; ui.fsRel.hidden = v !== 'relay'; ui.fsRot.hidden = v !== 'rot'; ui.fsSine.hidden = v !== 'sine';
      const needMat = v === 'core' || v === 'relay'; rowOf(sels.mat).hidden = !needMat; st.mur.el.hidden = !needMat || sels.mat.value !== 'const';
      ui.warn.hidden = true; ui.msg.show('');
      const mat = material();
      if (needMat && !mat.ok) { ui.msg.show('B–H table invalid: ' + shared.parsed.msg + ' (edit it in sub-topic 1).'); drawPlot(ui.cv, {}); hudSet(ui.hud, []); ui.work.set('—'); return; }
      if (v === 'core') return core(mat, T);
      if (v === 'relay') return relay(mat, T);
      if (v === 'rot') return rot(T);
      return sine(T);
    }
    function core(mat, T) {
      const k = { N: Math.round(st.N.get()), lc: st.lc.get(), A: st.A.get() * 1e-4, lg: st.lg.get() * 1e-3, mat }; k.Ag = k.A;
      const i = st.i.get(), a = E.coreAnalyze(k, i);
      const n = 120, is = [], ls = []; for (let j = 0; j <= n; j++) { const ii = 1.5 * i * j / n; is.push(ii); ls.push(k.N * E.coreFlux(k, ii)); }
      const ii0 = is.filter(z => z <= i + 1e-15), ll0 = ls.slice(0, ii0.length);
      const lmax = Math.max.apply(null, ls) * 1.1 || 1;
      drawPlot(ui.cv, {
        xmin: 0, xmax: 1.5 * i, ymin: 0, ymax: lmax, xlabel: 'current i (A)', ylabel: 'flux linkage λ (Wb-t)',
        fills: [{ x: ii0.concat([i]), y: ll0.concat([0]), color: withAlpha(T.b, 0.25) }, { x: [0].concat(ii0), y: [a.lam].concat(ll0), color: withAlpha(T.a, 0.25) }],
        series: [{ x: is, y: ls, color: T.c, width: 2.2 }], points: [{ x: i, y: a.lam, label: 'operating point', color: T.text }],
        legend: [{ text: 'energy W (above curve)', color: T.a }, { text: "co-energy W' (below)", color: T.b }],
      });
      const lin = mat.mode === 'const';
      hudSet(ui.hud, [['flux Φ', fx(a.phi, 4) + ' Wb'], ['B core / B gap', fx(a.B, 4) + ' / ' + fx(a.Bg, 4) + ' T'], ['λ = NΦ', fx(a.lam, 4) + ' Wb-t'],
        ['L = λ/i (secant)', fE(a.Lsec, 'H')], ['L incremental dλ/di', fE(a.Linc, 'H')], ['R core / R gap', fx(a.Rc, 4) + ' / ' + fx(a.Rg, 4) + ' A/Wb'],
        ['energy W', fE(a.W, 'J')], ['  in core / in gap', fE(a.Wcore, 'J') + ' / ' + fE(a.Wgap, 'J')], ["co-energy W'", fE(a.Wco, 'J')], ['gap share of W', f2(100 * a.Wgap / a.W, 1) + ' %']]);
      if (a.outside) { ui.warn.hidden = false; ui.warn.textContent = 'Core B is beyond the B–H table: extrapolated.'; }
      ui.work.set([
        'R_core = ' + (lin ? 'l_c/(μ0 μr A) = ' + fx(k.lc) + '/(4π×10⁻⁷ × ' + fx(mat.mur) + ' × ' + fx(k.A) + ')' : 'H_c·l_c/Φ (secant, nonlinear)') + ' = ' + fx(a.Rc) + ' A/Wb',
        'R_gap = l_g/(μ0 A) = ' + fx(k.lg) + '/(4π×10⁻⁷ × ' + fx(k.A) + ') = ' + fx(a.Rg) + ' A/Wb',
        lin ? 'L = N²/(R_core + R_gap) = ' + k.N + '²/' + fx(a.Rc + a.Rg) + ' = ' + fE(a.Lsec, 'H') : 'Nonlinear: Φ solved from N·i = H(Φ/A)·l_c + Φ·R_gap (bisection);  L = λ/i = ' + fE(a.Lsec, 'H'),
        'Φ = ' + fx(a.phi) + ' Wb,  B = Φ/A = ' + fx(a.B) + ' T,  λ = NΦ = ' + fx(a.lam) + ' Wb-t',
        'W_gap = (B²/2μ0)·A·l_g = ' + fx(a.Bg) + '²/(2 × 4π×10⁻⁷) × ' + fx(k.A) + ' × ' + fx(k.lg) + ' = ' + fE(a.Wgap, 'J'),
        'W_core = ' + (lin ? '(B²/2μ0μr)·A·l_c' : '(∫₀^B H dB)·A·l_c (exact trapezoids on the table)') + ' = ' + fE(a.Wcore, 'J'),
        'W = ∫ i dλ = ' + fE(a.W, 'J') + (lin ? ' = ½ L i² = ½ × ' + fx(a.Lsec) + ' × ' + fx(i) + '² ✓' : ''),
        "Co-energy W' = λ·i − W = " + fE(a.Wco, 'J') + (lin ? ' (= W for a linear core)' : ' (≠ W: saturating core)'),
        'Ratio W_gap/W_core = ' + fx(a.Wgap / a.Wcore) + (lin ? ' = R_gap/R_core: most energy is stored in the gap.' : ''),
      ]);
    }
    function relay(mat, T) {
      const r = { N: Math.round(st.rN.get()), lc: st.rlc.get(), A: st.rA.get() * 1e-4, ng: +sels.ng.value, mat };
      const i = st.ri.get(), x = st.rx.get() * 1e-3, Fs = st.rFs.get(), x0 = st.rx0.get() * 1e-3, lin = mat.mode === 'const';
      const F = E.relayForce(r, i, x), ip = E.relayPullIn(r, x0, Fs);
      const xm = Math.max(x, x0) * 1.5, n = 100, xs = [], fs = [], fp = [];
      for (let j = 1; j <= n; j++) { const xx = xm * j / n; xs.push(xx * 1e3); fs.push(Math.abs(E.relayForce(r, i, xx))); if (fin(ip)) fp.push(Math.abs(E.relayForce(r, ip, xx))); }
      const ymax = Math.min(Math.max(Fs * 3, Math.abs(F) * 1.3), Math.max.apply(null, fs.concat(fp).concat([Fs])) * 1.05) || 1;
      drawPlot(ui.cv, {
        xmin: 0, xmax: xm * 1e3, ymin: 0, ymax, xlabel: 'gap x (mm)', ylabel: 'attractive force |F| (N)',
        series: [{ x: xs, y: fs, color: T.c, width: 2.2 }].concat(fp.length ? [{ x: xs, y: fp, color: T.pink, width: 1.5, dash: [5, 3] }] : []),
        hlines: [{ y: Fs, label: 'spring/load F_s', color: T.warn }], vlines: [{ x: x0 * 1e3, label: 'x0 (open)', color: T.muted }],
        points: [{ x: x * 1e3, y: Math.abs(F), label: f2(Math.abs(F), 2) + ' N', color: T.text }],
        legend: [{ text: 'F at i = ' + fx(i, 3) + ' A', color: T.c }].concat(fp.length ? [{ text: 'F at pull-in ' + fx(ip, 3) + ' A', color: T.pink }] : []),
      });
      const Lx = lin ? E.relayL(r, x) : NaN, dL = lin ? E.relaydLdx(r, x) : NaN, B = E.coreFlux(E.relayCore(r, x), i) / r.A;
      hudSet(ui.hud, [['force F (+ = opening)', fx(F, 4) + ' N'], ['attractive |F|', fx(Math.abs(F), 4) + ' N'], ['gap B', fx(B, 4) + ' T'],
        ['L(x)', lin ? fE(Lx, 'H') : fE(r.N * E.coreFlux(E.relayCore(r, x), i) / i, 'H') + ' (secant)'], ['dL/dx', lin ? fx(dL, 4) + ' H/m' : '— (nonlinear)'],
        ['pull-in current at x0', fin(ip) ? fx(ip, 4) + ' A' : '—'], ['closes at present i?', fin(ip) ? (i >= ip ? 'yes (i ≥ i_pull)' : 'no (i < i_pull)') : '—']]);
      if (mat.outside(B)) { ui.warn.hidden = false; ui.warn.textContent = 'Iron B exceeds the B–H table: extrapolated.'; }
      const Ri = lin ? r.lc / (MU0 * mat.mur * r.A) : NaN, Rg = r.ng * x / (MU0 * r.A);
      ui.work.set(lin ? [
        'Reluctance: R(x) = l_c/(μ0 μr A) + n·x/(μ0 A) = ' + fx(Ri) + ' + ' + r.ng + '×' + fx(x) + '/(4π×10⁻⁷ × ' + fx(r.A) + ') = ' + fx(Ri + Rg) + ' A/Wb',
        'L(x) = N²/R(x) = ' + r.N + '²/' + fx(Ri + Rg) + ' = ' + fE(Lx, 'H'),
        'dL/dx = −N²·(n/(μ0 A))/R(x)² = ' + fx(dL) + ' H/m',
        'F = ½ i² dL/dx = ½ × ' + fx(i) + '² × (' + fx(dL) + ') = ' + fx(F) + ' N  (negative → pulls the armature in, closing the gap)',
        'Check (iron ideal): F ≈ n·B²A/(2μ0) with B = μ0 N i/(n x) — equal when R_iron ≪ R_gap.',
        'Pull-in at x0 = ' + fx(x0 * 1e3) + ' mm: ½ i² |dL/dx(x0)| = F_s → i = √(2F_s/|dL/dx|) = √(2 × ' + fx(Fs) + '/' + fx(Math.abs(E.relaydLdx(r, x0))) + ') = ' + fx(ip) + ' A',
        'Force ∝ i² and ∝ 1/R(x)²: it rises steeply as the gap closes, so once pulled in the relay stays closed (hold current < pull-in current).',
      ] : [
        'Nonlinear core: λ(i, x) from N·i = H(λ/(NA))·l_c + (λ/N)·n x/(μ0 A) (bisection).',
        "Co-energy W'(i, x) = ∫₀^i λ(i', x) di' (Simpson, 120 steps);  F = ∂W'/∂x at constant i (central difference).",
        "F = " + fx(F) + ' N.  Equivalent route: F = −∂W/∂x at constant λ (same value; see tests).',
        'Pull-in current found by bisection on |F(i, x0)| = F_s: i = ' + fx(ip) + ' A',
      ]);
    }
    function rot(T) {
      const Ld = st.Ld.get(), Lq = st.Lq.get(), i = st.ti.get(), th = st.th.get() * PI / 180;
      if (!(Ld > Lq)) { ui.msg.show('Need L_aligned > L_unaligned for a reluctance torque.'); }
      const ths = [], Ts = [], Lsr = []; for (let j = 0; j <= 180; j++) { ths.push(j); Ts.push(E.rotTorque(Ld, Lq, i, j * PI / 180)); Lsr.push(E.rotL(Ld, Lq, j * PI / 180)); }
      const Tm = 0.5 * i * i * Math.abs(Ld - Lq) || 1, Tn = E.rotTorque(Ld, Lq, i, th);
      drawPlot(ui.cv, { xmin: 0, xmax: 180, ymin: -1.2 * Tm, ymax: 1.2 * Tm, y2min: 0, y2max: Math.max(Ld, Lq) * 1.1, xlabel: 'rotor angle θ from aligned axis (deg)', ylabel: 'torque T (N·m)', y2label: 'L(θ) (H)',
        series: [{ x: ths, y: Ts, color: T.c, width: 2.2 }, { x: ths, y: Lsr, color: T.a, width: 1.5, dash: [5, 3], axis: 'y2' }],
        points: [{ x: st.th.get(), y: Tn, label: fx(Tn, 3) + ' N·m', color: T.text }], legend: [{ text: 'T(θ)', color: T.c }, { text: 'L(θ)', color: T.a }] });
      hudSet(ui.hud, [['L(θ)', fE(E.rotL(Ld, Lq, th), 'H')], ['dL/dθ', fx(-(Ld - Lq) * Math.sin(2 * th), 4) + ' H/rad'], ['torque T', fx(Tn, 4) + ' N·m'], ['max |T| (θ = 45°)', fx(0.5 * i * i * (Ld - Lq), 4) + ' N·m']]);
      ui.work.set([
        'L(θ) = (L_d + L_q)/2 + (L_d − L_q)/2 · cos 2θ  (θ measured from the aligned position)',
        'T = ½ i² dL/dθ = −½ i² (L_d − L_q) sin 2θ = −½ × ' + fx(i) + '² × ' + fx(Ld - Lq) + ' × sin(' + fx(2 * st.th.get()) + '°) = ' + fx(Tn) + ' N·m',
        'Negative T for 0 < θ < 90°: the torque pulls the rotor back toward alignment (minimum reluctance). Peak at θ = 45°.',
        'Linear magnetics assumed (L independent of i), so co-energy W\' = ½ L(θ) i² and T = ∂W\'/∂θ at constant i.',
      ]);
    }
    function sine(T) {
      const V = st.V.get(), f = st.f.get(), N = Math.round(st.sN.get()), A = st.sA.get() * 1e-4, vol = st.vol.get() * 1e-6, n = st.n.get();
      const sf = E.sineFlux(V, f, N, A), sp = E.separateLosses({ f1: st.f1.get(), P1: st.P1.get(), f2: st.f2.get(), P2: st.P2.get(), Bm: sf.Bm, n, vol });
      if (!sp.ok) ui.msg.show(sp.msg);
      const fmax = Math.max(st.f1.get(), st.f2.get(), f) * 1.4, fsr = [], yy = [], ph = [], pe = [];
      for (let j = 0; j <= 100; j++) { const ff = fmax * j / 100; fsr.push(ff); yy.push(sp.a + sp.b * ff); }
      const ymax = Math.max(sp.y1, sp.y2, sp.a + sp.b * fmax) * 1.15 || 1;
      drawPlot(ui.cv, { xmin: 0, xmax: fmax, ymin: 0, ymax: fin(ymax) && ymax > 0 ? ymax : 1, xlabel: 'frequency f (Hz)', ylabel: 'P/f (J per cycle)',
        series: sp.ok ? [{ x: fsr, y: yy, color: T.c, width: 2.2 }] : [],
        points: [{ x: st.f1.get(), y: sp.y1, label: 'test 1', color: T.a }, { x: st.f2.get(), y: sp.y2, label: 'test 2', color: T.a, below: true }].concat(sp.ok ? [{ x: 0, y: sp.a, label: 'a = hysteresis/cycle', color: T.b }] : []),
        legend: [{ text: 'P/f = a + b·f', color: T.c }] });
      const atF = sp.ok ? { Ph: sp.a * f, Pe: sp.b * f * f } : { Ph: NaN, Pe: NaN };
      hudSet(ui.hud, [['Φmax', fx(sf.phim * 1e3, 4) + ' mWb'], ['Bmax', fx(sf.Bm, 4) + ' T'], ['a (intercept)', fx(sp.a, 4) + ' J'], ['b (slope)', fx(sp.b, 4) + ' J·s'],
        ['k_h', fx(sp.kh, 4) + ' J/(m³·T^n)'], ['k_e', fx(sp.ke, 4) + ' J·s/(m³·T²)'], ['P_h / P_e at f1', fx(sp.Ph1, 4) + ' / ' + fx(sp.Pe1, 4) + ' W'], ['P_h / P_e at f2', fx(sp.Ph2, 4) + ' / ' + fx(sp.Pe2, 4) + ' W'],
        ['at f = ' + fx(f) + ' Hz (same Bmax)', fx(atF.Ph, 4) + ' + ' + fx(atF.Pe, 4) + ' W']]);
      if (sf.Bm > 2) { ui.warn.hidden = false; ui.warn.textContent = 'Bmax > 2 T: a real silicon-steel core would be deeply saturated — increase N or A.'; }
      ui.work.set([
        'e = N dΦ/dt with Φ = Φmax sin ωt → E_rms = ωNΦmax/√2 = (2π/√2) f N Φmax = 4.443 f N Φmax (≈ 4.44)',
        'Φmax = V/(4.443 f N) = ' + fx(V) + '/(4.443 × ' + fx(f) + ' × ' + N + ') = ' + fx(sf.phim) + ' Wb;  Bmax = Φmax/A = ' + fx(sf.Bm) + ' T',
        'Steinmetz: P_h = k_h f Bmax^n V,   eddy: P_e = k_e f² Bmax² V   (n = ' + fx(n) + ')',
        'At constant Bmax:  P/f = a + b f  with a = k_h Bmax^n V and b = k_e Bmax² V',
        'Test 1: P1/f1 = ' + fx(st.P1.get()) + '/' + fx(st.f1.get()) + ' = ' + fx(sp.y1) + ';   Test 2: P2/f2 = ' + fx(st.P2.get()) + '/' + fx(st.f2.get()) + ' = ' + fx(sp.y2),
        'b = (' + fx(sp.y2) + ' − ' + fx(sp.y1) + ')/(' + fx(st.f2.get()) + ' − ' + fx(st.f1.get()) + ') = ' + fx(sp.b) + ';   a = ' + fx(sp.y1) + ' − b × ' + fx(st.f1.get()) + ' = ' + fx(sp.a),
        'P_h(f1) = a f1 = ' + fx(sp.Ph1) + ' W,  P_e(f1) = b f1² = ' + fx(sp.Pe1) + ' W  (sum = ' + fx(sp.Ph1 + sp.Pe1) + ' W ✓)',
        'k_h = a/(Bmax^n V) = ' + fx(sp.a) + '/(' + fx(sf.Bm) + '^' + fx(n) + ' × ' + fx(vol) + ') = ' + fx(sp.kh) + ';   k_e = b/(Bmax² V) = ' + fx(sp.ke),
      ]);
    }
    const io = stateIO(st, sels);
    return { update, get: io.get, set: io.set };
  }

  /* ---------- 3. autotransformer ---------- */
  function buildAuto(root) {
    const st = {}, sels = {}, ui = {}, { ctl, stage } = layout(root);
    const fs1 = FSP.ui.fieldset(ctl, 'Two-winding transformer');
    makeSliders(fs1, [{ k: 'V1', l: 'V1 winding', min: 1, max: 500000, v: 220, u: 'V', log: true }, { k: 'V2', l: 'V2 winding', min: 1, max: 500000, v: 110, u: 'V', log: true },
      { k: 'S', l: 'S rating', min: 0.01, max: 100000, v: 5, u: 'kVA', log: true }, { k: 'Pc', l: 'P core', min: 0, max: 100000, v: 50, step: 0.1, u: 'W' }, { k: 'Pcu', l: 'P cu (FL)', min: 0, max: 100000, v: 100, step: 0.1, u: 'W' }], st, update);
    const fs2 = FSP.ui.fieldset(ctl, 'Auto connection & load');
    sels.common = sel(fs2, 'Common wdg', [['V1', 'V1 winding'], ['V2', 'V2 winding']], 'V1', update);
    sels.dir = sel(fs2, 'Use as', [['down', 'step-down (H → L)'], ['up', 'step-up (L → H)']], 'down', update);
    makeSliders(fs2, [{ k: 'x', l: 'load', min: 0, max: 1.5, v: 1, step: 0.01, u: '× rated' }, { k: 'pf', l: 'pf', min: 0.05, max: 1, v: 0.8, step: 0.01, u: '' }], st, update);
    fs2.appendChild(FSP.ui.el('div', { class: 'note', text: 'Windings connected series-aiding: V_H = V_common + V_series, V_L = V_common.' }));
    ui.msg = msgBox(ctl);
    ui.cv = wrapCanvas(stage, 'Autotransformer connection with currents', 300); ui.cv.onResize(update);
    ui.hud = FSP.ui.el('div', { class: 'hud' }); stage.appendChild(ui.hud);
    ui.work = FSP.ui.working(stage);
    function update() {
      if (!ui.cv) return;
      const r = E.auto({ V1: st.V1.get(), V2: st.V2.get(), S: st.S.get() * 1e3, common: sels.common.value, dir: sels.dir.value, x: st.x.get(), pf: st.pf.get(), Pcore: st.Pc.get(), Pcu: st.Pcu.get() });
      ui.msg.show(r.ok ? '' : r.msg); if (!r.ok) { hudSet(ui.hud, []); ui.work.set(r.msg); return; }
      const x = st.x.get(), down = sels.dir.value === 'down';
      draw(r, x, down);
      hudSet(ui.hud, [['voltage ratio', fx(r.VH, 5) + ' / ' + fx(r.VL, 5) + ' V  (' + fx(r.ratio, 4) + ':1)'], ['S auto rating', fx(r.Sa / 1e3, 5) + ' kVA (' + fx(r.gain, 4) + '× S)'],
        ['I_H (= series wdg)', fx(r.IH, 5) + ' A'], ['I_L', fx(r.IL, 5) + ' A'], ['common wdg I_L − I_H', fx(r.Icom, 5) + ' A'],
        ['transformed S', fx(r.Strans / 1e3, 5) + ' kVA'], ['conducted S', fx(r.Scond / 1e3, 5) + ' kVA'],
        ['η two-winding @ x', fx(100 * r.eta2w, 5) + ' %'], ['η auto @ x', fx(100 * r.etaAuto, 5) + ' %']]);
      ui.work.set([
        'Common winding: ' + fx(r.Vc) + ' V (rated ' + fx(r.IcomRated) + ' A); series winding: ' + fx(r.Vs) + ' V (rated ' + fx(r.IserRated) + ' A)',
        'V_H = V_common + V_series = ' + fx(r.Vc) + ' + ' + fx(r.Vs) + ' = ' + fx(r.VH) + ' V;   V_L = ' + fx(r.VL) + ' V;   ratio a = ' + fx(r.ratio),
        'The series winding carries the HV-side current: I_H = I_series,rated = S/V_series = ' + fx(r.S || r.Sw) + '/' + fx(r.Vs) + ' = ' + fx(r.IH) + ' A',
        'S_auto = V_H·I_H = V_H·S/(V_H − V_L) = ' + fx(r.VH) + ' × ' + fx(r.Sw) + '/' + fx(r.VH - r.VL) + ' = ' + fx(r.Sa) + ' VA',
        'I_L = S_auto/V_L = ' + fx(r.IL) + ' A;  common winding I_L − I_H = ' + fx(r.Icom) + ' A = S/V_common ✓ (both windings exactly at rating)',
        'Transformed (magnetically) = (V_H − V_L)·I_H = ' + fx(r.Strans) + ' VA = S (two-winding rating);  conducted = V_L·I_H = ' + fx(r.Scond) + ' VA',
        'Losses unchanged (same winding currents and flux): P = P_core + x²P_cu = ' + fx(st.Pc.get()) + ' + ' + fx(x) + '² × ' + fx(st.Pcu.get()) + ' = ' + fx(r.loss) + ' W',
        'η_2w = xS·pf/(xS·pf + P) = ' + fx(r.P2w) + '/' + fx(r.P2w + r.loss) + ' = ' + fx(100 * r.eta2w) + ' %;   η_auto = ' + fx(r.Pauto) + '/' + fx(r.Pauto + r.loss) + ' = ' + fx(100 * r.etaAuto) + ' %',
        'Trade-off: no electrical isolation between HV and LV, and a higher fault level (lower per-unit impedance on the auto base).',
      ]);
    }
    function draw(r, x, down) {
      const g = ui.cv.prep(), ctx = g.ctx, w = g.w, h = g.h, T = theme(); ctx.clearRect(0, 0, w, h); ctx.fillStyle = T.bg; ctx.fillRect(0, 0, w, h); if (w < 160) return;
      const cx = w / 2, yH = 40, yT = h * 0.47, yN = h - 40, xl = Math.max(30, cx - Math.min(170, w * 0.4)), xr = Math.min(w - 30, cx + Math.min(170, w * 0.4));
      ctx.lineWidth = 2; ctx.strokeStyle = T.muted;
      const coil = (ya, yb, col) => { const n = 6, step = (yb - ya) / n; ctx.strokeStyle = col; ctx.lineWidth = 2; ctx.beginPath(); for (let k = 0; k < n; k++) ctx.arc(cx, ya + step * (k + 0.5), step / 2, -PI / 2, PI / 2, false); ctx.stroke(); };
      ctx.strokeStyle = T.muted; ctx.beginPath(); ctx.moveTo(xl, yH); ctx.lineTo(cx, yH); ctx.lineTo(cx, yH + 12); ctx.moveTo(cx, yT - 12); ctx.lineTo(cx, yT + 12); ctx.moveTo(cx, yN - 12); ctx.lineTo(cx, yN); ctx.lineTo(xl, yN); ctx.moveTo(cx, yN); ctx.lineTo(xr, yN); ctx.moveTo(cx, yT); ctx.lineTo(xr, yT); ctx.stroke();
      coil(yH + 12, yT - 12, T.a); coil(yT + 12, yN - 12, T.c);
      // core bar
      ctx.strokeStyle = T.muted; ctx.lineWidth = 3; ctx.beginPath(); ctx.moveTo(cx + 26, yH + 10); ctx.lineTo(cx + 26, yN - 10); ctx.stroke(); ctx.lineWidth = 1.5;
      const T0 = (t, xx, yy, al, col) => { ctx.fillStyle = col || T.text; ctx.font = '11px ' + T.mono; ctx.textAlign = al || 'left'; ctx.textBaseline = 'middle'; ctx.fillText(t, xx, yy); };
      [[xl, yH], [xl, yN], [xr, yT], [xr, yN]].forEach(p => { ctx.fillStyle = T.text; ctx.beginPath(); ctx.arc(p[0], p[1], 3.5, 0, 2 * PI); ctx.fill(); });
      T0('H', xl - 4, yH - 12, 'left'); T0('L', xr + 4, yT - 12, 'right');
      T0('V_H = ' + fx(r.VH, 4) + ' V', xl, (yH + yN) / 2, 'left', T.muted); T0('V_L = ' + fx(r.VL, 4) + ' V', xr, (yT + yN) / 2, 'right', T.muted);
      T0('series ' + fx(r.Vs, 4) + ' V', cx - 14, (yH + yT) / 2, 'right', T.a); T0('common ' + fx(r.Vc, 4) + ' V', cx - 14, (yT + yN) / 2, 'right', T.c);
      T0(down ? 'source' : 'load', xl, yH + 16, 'left', T.muted); T0(down ? 'load' : 'source', xr, yT + 16, 'right', T.muted);
      // currents
      ctx.strokeStyle = T.pink; ctx.fillStyle = T.pink; ctx.lineWidth = 2;
      const mxH = (xl + cx) / 2, mxL = (cx + xr) / 2, IH = x * r.IH, IL = x * r.IL, Ic = x * r.Icom;
      if (down) arrowHead(ctx, mxH - 18, yH, mxH + 18, yH, 7); else arrowHead(ctx, mxH + 18, yH, mxH - 18, yH, 7);
      if (down) arrowHead(ctx, mxL - 18, yT, mxL + 18, yT, 7); else arrowHead(ctx, mxL + 18, yT, mxL - 18, yT, 7);
      const ya = (yT + yN) / 2; if (down) arrowHead(ctx, cx + 12, ya + 16, cx + 12, ya - 16, 7); else arrowHead(ctx, cx + 12, ya - 16, cx + 12, ya + 16, 7);
      const yb = (yH + yT) / 2; if (down) arrowHead(ctx, cx + 12, yb - 16, cx + 12, yb + 16, 7); else arrowHead(ctx, cx + 12, yb + 16, cx + 12, yb - 16, 7);
      ctx.font = '11px ' + T.mono; ctx.textBaseline = 'bottom'; ctx.textAlign = 'center';
      ctx.fillText('I_H=' + fx(IH, 4) + 'A', mxH, yH - 6); ctx.fillText('I_L=' + fx(IL, 4) + 'A', mxL, yT - 6);
      ctx.textAlign = 'left'; ctx.textBaseline = 'middle'; ctx.fillText(fx(Ic, 4) + 'A', cx + 32, ya); ctx.fillText(fx(IH, 4) + 'A', cx + 32, yb);
      T0('S_auto = ' + fx(r.Sa / 1e3, 4) + ' kVA  (load ' + fx(x, 3) + '×)', 8, 12, 'left', T.text);
    }
    const io = stateIO(st, sels);
    return { update, get: io.get, set: io.set };
  }

  /* ---------- 4. three-phase transformers ---------- */
  function buildTP(root) {
    const st = {}, sels = {}, ui = {}, { ctl, stage } = layout(root);
    const fs1 = FSP.ui.fieldset(ctl, 'Bank of three single-phase units');
    sels.conn = sel(fs1, 'HV-LV', [['YY', 'Y – Y'], ['YD', 'Y – Δ'], ['DY', 'Δ – Y'], ['DD', 'Δ – Δ'], ['VV', 'open-Δ (V – V)']], 'YD', update);
    makeSliders(fs1, [{ k: 'S', l: 'S bank (3 units)', min: 0.1, max: 1e6, v: 150, u: 'kVA', log: true }, { k: 'VH', l: 'V line HV', min: 1, max: 1e6, v: 11000, u: 'V', log: true },
      { k: 'VL', l: 'V line LV', min: 1, max: 1e6, v: 415, u: 'V', log: true }], st, update);
    sels.conv = sel(fs1, 'Phase shift', [['ANSI', 'ANSI: LV lags HV by 30° (Yd1/Dy1)'], ['IEC11', 'IEC clock 11: LV leads 30° (Dyn11/YNd11)']], 'ANSI', update);
    const fs2 = FSP.ui.fieldset(ctl, 'Each unit: series impedance');
    makeSliders(fs2, [{ k: 'R', l: 'R_eq (HV wdg)', min: 0, max: 10000, v: 12, step: 0.01, u: 'Ω' }, { k: 'X', l: 'X_eq (HV wdg)', min: 0, max: 10000, v: 30, step: 0.01, u: 'Ω' }], st, update);
    fs2.appendChild(FSP.ui.el('div', { class: 'note', text: 'Impedance of one single-phase unit referred to its own HV winding.' }));
    ui.msg = msgBox(ctl);
    ui.cv = wrapCanvas(stage, 'HV and LV line-voltage phasors', 300); ui.cv.onResize(update);
    ui.hud = FSP.ui.el('div', { class: 'hud' }); stage.appendChild(ui.hud);
    ui.work = FSP.ui.working(stage);
    const nm = { Y: 'Y', D: 'Δ' };
    function update() {
      if (!ui.cv) return;
      const conn = sels.conn.value, r = E.threePhase({ conn, S: st.S.get() * 1e3, VHL: st.VH.get(), VLL: st.VL.get(), conv: sels.conv.value, Req: st.R.get(), Xeq: st.X.get() });
      ui.msg.show(r.ok ? '' : r.msg); if (!r.ok) { hudSet(ui.hud, []); ui.work.set(r.msg); return; }
      rowOf(sels.conv).hidden = !(conn === 'YD' || conn === 'DY');
      draw(r);
      const hn = conn === 'VV' ? 'Δ (open)' : nm[conn[0]], ln = conn === 'VV' ? 'Δ (open)' : nm[conn[1]];
      hudSet(ui.hud, [['HV: line / winding V', fx(st.VH.get(), 5) + ' / ' + fx(r.VHw, 5) + ' V'], ['LV: line / winding V', fx(st.VL.get(), 5) + ' / ' + fx(r.VLw, 5) + ' V'],
        ['winding ratio a', fx(r.a, 5)], ['line-voltage ratio', fx(r.lineRatio, 5)], ['unit rating', fx(r.Sunit / 1e3, 5) + ' kVA'], ['bank capacity', fx(r.Sbank / 1e3, 5) + ' kVA'],
        ['HV line / winding I', fx(r.IHL, 5) + ' / ' + fx(r.IHw, 5) + ' A'], ['LV line / winding I', fx(r.ILL, 5) + ' / ' + fx(r.ILw, 5) + ' A'],
        ['LV vs HV line voltage', (r.shift > 0 ? '+' : '') + r.shift + '°'], ['Z per phase (Y-eq) HV', fx(r.ZYH.re, 4) + ' + j' + fx(r.ZYH.im, 4) + ' Ω'],
        ['Z per phase (Y-eq) LV', fx(r.ZYL.re, 4) + ' + j' + fx(r.ZYL.im, 4) + ' Ω'], ['Z per unit', fx(r.Zpu.re, 4) + ' + j' + fx(r.Zpu.im, 4) + ' pu']]);
      const L = [
        'Connection ' + hn + ' (HV) – ' + ln + ' (LV). Y: V_winding = V_line/√3, I_winding = I_line.  Δ: V_winding = V_line, I_winding = I_line/√3.',
        'HV winding V = ' + fx(r.VHw) + ' V;  LV winding V = ' + fx(r.VLw) + ' V;  turns ratio a = ' + fx(r.a),
        'Line-voltage ratio V_HL/V_LL = ' + fx(r.lineRatio) + (conn === 'YD' ? ' = √3·a' : conn === 'DY' ? ' = a/√3' : ' = a'),
      ];
      if (conn === 'VV') L.push('Open-Δ: two units supply the full 3-phase voltage set; each carries the line current. Capacity = √3·V_L·I_rated = √3·S_unit = ' + fx(r.Sbank / 1e3) + ' kVA = ' + fx(100 * r.vvFraction) + ' % of the Δ-Δ bank (3 S_unit); unit utilisation √3/2 = ' + fx(100 * r.vvUtil) + ' %.');
      else L.push('Unit rating = S/3 = ' + fx(r.Sunit / 1e3) + ' kVA.');
      L.push('I_line = S/(√3 V_line):  HV ' + fx(r.Sbank) + '/(√3 × ' + fx(st.VH.get()) + ') = ' + fx(r.IHL) + ' A;  LV = ' + fx(r.ILL) + ' A;  winding currents ' + fx(r.IHw) + ' A / ' + fx(r.ILw) + ' A');
      if (conn === 'YD' || conn === 'DY') L.push('Phase shift: ' + (sels.conv.value === 'ANSI' ? 'ANSI C57.12 standard — the HV line voltage leads the LV line voltage by 30° (V_ab lags V_AB by 30°); IEC equivalent Yd1 / Dy1.' : 'IEC clock number 11 (e.g. Dyn11): LV line voltage leads HV by 30° (LV phasor at "11 o\'clock" when HV is at 12). Achieved by reversing the Δ connection relative to ANSI.'));
      else L.push('Y-Y / Δ-Δ (and V-V): no phase shift between HV and LV line voltages (clock 0).');
      L.push('Per-phase (Y-equivalent, line-to-neutral) model referred to HV: Z = ' + (r.hY ? 'Z_winding (Y on HV)' : 'Z_winding/3 (Δ→Y)') + ' = ' + fx(r.ZYH.re) + ' + j' + fx(r.ZYH.im) + ' Ω;  referred to LV × (V_LL/V_HL)² = ' + fx(r.ZYL.re) + ' + j' + fx(r.ZYL.im) + ' Ω');
      L.push('Per unit: Z_base(HV) = V_HL²/S = ' + fx(r.Zbase) + ' Ω → Z = ' + fx(r.Zpu.re) + ' + j' + fx(r.Zpu.im) + ' pu (same on either side). Per-phase analysis then uses V_line/√3 and ignores the 30° shift (add it back at the end).');
      ui.work.set(L);
    }
    function draw(r) {
      const g = ui.cv.prep(), ctx = g.ctx, w = g.w, h = g.h, T = theme(); ctx.clearRect(0, 0, w, h); ctx.fillStyle = T.bg; ctx.fillRect(0, 0, w, h); if (w < 120) return;
      const cx = w / 2, cy = h / 2 + 8, R = Math.min(w, h) / 2 - 34;
      ctx.strokeStyle = T.border; ctx.lineWidth = 1; ctx.beginPath(); ctx.arc(cx, cy, R, 0, 2 * PI); ctx.stroke();
      const ph = (ang, len, col, lab, dash) => {
        const a = ang * PI / 180, x1 = cx + len * Math.cos(a), y1 = cy - len * Math.sin(a); ctx.strokeStyle = col; ctx.fillStyle = col; ctx.lineWidth = 2.2; ctx.setLineDash(dash || []);
        arrowHead(ctx, cx, cy, x1, y1, 9); ctx.setLineDash([]); ctx.font = '11px ' + T.mono; ctx.textAlign = 'center'; ctx.textBaseline = 'middle'; ctx.fillText(lab, cx + (len + 16) * Math.cos(a), cy - (len + 14) * Math.sin(a));
      };
      const base = 90; // V_AB drawn at 12 o'clock
      ['AB', 'BC', 'CA'].forEach((n, k) => ph(base - 120 * k, R, T.c, 'V' + n));
      ['ab', 'bc', 'ca'].forEach((n, k) => ph(base - 120 * k + r.shift, R * 0.72, T.a, 'V' + n, [6, 3]));
      ctx.fillStyle = T.text; ctx.font = '11px ' + T.mono; ctx.textAlign = 'left'; ctx.textBaseline = 'top';
      ctx.fillText('HV (solid) at 12 o\'clock; LV (dashed) ' + (r.shift === 0 ? 'in phase' : r.shift < 0 ? 'lags 30°' : 'leads 30°'), 6, 5);
      ctx.fillStyle = T.muted; ctx.fillText('lengths normalised; ABC sequence, CCW rotation', 6, 19);
      if (r.shift !== 0) { ctx.strokeStyle = T.pink; ctx.lineWidth = 1.2; ctx.beginPath(); const a0 = -base * PI / 180, a1 = -(base + r.shift) * PI / 180; ctx.arc(cx, cy, R * 0.45, Math.min(a0, a1), Math.max(a0, a1)); ctx.stroke(); }
    }
    const io = stateIO(st, sels);
    return { update, get: io.get, set: io.set };
  }

  /* ---------- 5. all-day efficiency ---------- */
  function buildAD(root) {
    const st = {}, sels = {}, ui = {}, { ctl, stage } = layout(root);
    ctl.appendChild(FSP.ui.el('div', { class: 'note', text: 'Regulation, OC/SC tests and η(load) at one operating point are in the Magnetics & Transformer tab. This adds the all-day (energy) efficiency over a 24-h load cycle.' }));
    const fs1 = FSP.ui.fieldset(ctl, 'Transformer');
    makeSliders(fs1, [{ k: 'S', l: 'S rated', min: 0.1, max: 1e5, v: 100, u: 'kVA', log: true }, { k: 'Pc', l: 'P core', min: 0, max: 1000, v: 1, step: 0.001, u: 'kW' },
      { k: 'Pcu', l: 'P cu (FL)', min: 0.001, max: 1000, v: 1.5, u: 'kW', log: true }, { k: 'pf', l: 'pf for η(x)', min: 0.05, max: 1, v: 0.8, step: 0.01, u: '' }], st, update);
    const fs2 = FSP.ui.fieldset(ctl, '24-h load profile');
    fs2.appendChild(FSP.ui.el('div', { class: 'note', text: 'One row per period: hours  kW  pf. Remaining hours are no-load (core loss only — the transformer stays energised).' }));
    ui.prof = FSP.ui.el('textarea', { rows: '6', spellcheck: 'false', 'aria-label': 'load profile' }); ui.prof.value = '6 80 0.8\n10 40 1\n8 0 1';
    ui.prof.style.cssText = 'width:100%;box-sizing:border-box;background:var(--panel2);color:var(--text);border:1px solid var(--border);border-radius:6px;font:12px var(--mono);padding:4px';
    ui.prof.addEventListener('change', () => { update(); FSP.state.touch(); }); fs2.appendChild(ui.prof);
    ui.msg = msgBox(ctl);
    ui.cv1 = wrapCanvas(stage, '24-hour load profile', 220); ui.cv1.onResize(update);
    ui.cv2 = wrapCanvas(stage, 'Efficiency versus load fraction', 240); ui.cv2.onResize(update);
    ui.hud = FSP.ui.el('div', { class: 'hud' }); stage.appendChild(ui.hud);
    ui.warn = FSP.ui.el('div', { class: 'msg warn', hidden: '' }); stage.appendChild(ui.warn);
    ui.work = FSP.ui.working(stage);
    function update() {
      if (!ui.cv1) return;
      const T = theme(), S = st.S.get(), Pc = st.Pc.get(), Pcu = st.Pcu.get(), pf = st.pf.get(), pr = E.parseProfile(ui.prof.value);
      const ad = pr.ok ? E.allDay({ S, Pcore: Pc, Pcu, rows: pr.rows }) : pr;
      ui.msg.show(ad.ok ? '' : ad.msg);
      const me = E.maxEffPoint(S, Pc, Pcu, pf);
      // efficiency curve
      const xs = [], es = []; for (let j = 1; j <= 150; j++) { const x = j / 100; xs.push(x); es.push(100 * E.etaAt(S, Pc, Pcu, x, pf)); }
      const emin = Math.max(0, Math.min.apply(null, es.filter(fin)) - 1);
      drawPlot(ui.cv2, { xmin: 0, xmax: 1.5, ymin: emin, ymax: 100, xlabel: 'load fraction x = S_load/S_rated', ylabel: 'efficiency η (%) at pf ' + fx(pf, 3),
        series: [{ x: xs, y: es, color: T.c, width: 2.2 }], vlines: [{ x: me.x, label: 'x* = ' + fx(me.x, 4), color: T.warn }, { x: 1, color: T.muted }],
        points: [{ x: me.x, y: 100 * me.eta, label: 'η max ' + fx(100 * me.eta, 4) + '%', color: T.warn, below: true }].concat(ad.ok ? ad.rows.filter(r => r.kW > 0).map(r => ({ x: r.x, y: 100 * E.etaAt(S, Pc, Pcu, r.x, r.pf), color: T.a })) : []) });
      if (!ad.ok) { drawPlot(ui.cv1, {}); hudSet(ui.hud, [['max-efficiency load x*', fx(me.x, 4)]]); ui.work.set(ad.msg); ui.warn.hidden = true; return; }
      // profile bars (kVA fraction) over 24 h
      const bars = []; let t0 = 0; ad.rows.forEach(r => { bars.push({ x0: t0, x1: t0 + r.h, y: r.x, color: withAlpha(T.c, 0.7) }); t0 += r.h; });
      if (t0 < 24) bars.push({ x0: t0, x1: 24, y: 0, color: T.muted });
      const ymax = Math.max(1.1, Math.max.apply(null, ad.rows.map(r => r.x)) * 1.15);
      drawPlot(ui.cv1, { xmin: 0, xmax: 24, ymin: 0, ymax, xlabel: 'hour of day', ylabel: 'load fraction x', bars,
        hlines: [{ y: me.x, label: 'x* (max η)', color: T.warn }, { y: 1, label: 'rated', color: T.muted }] });
      ui.warn.hidden = !ad.overload; ui.warn.textContent = ad.overload ? 'A period exceeds rated kVA (x > 1).' : '';
      hudSet(ui.hud, [['energy out', fx(ad.Eout, 5) + ' kWh'], ['core loss energy (24 h)', fx(ad.Ecore, 5) + ' kWh'], ['copper loss energy', fx(ad.Ecu, 5) + ' kWh'],
        ['all-day efficiency', fx(100 * ad.eta, 5) + ' %'], ['max-η load fraction x*', fx(me.x, 4) + '  (' + fx(me.Skva, 4) + ' kVA)'], ['η max at pf ' + fx(pf, 3), fx(100 * me.eta, 5) + ' %'],
        ['idle hours', fx(ad.idle, 4) + ' h']]);
      const L = ['All-day (energy) efficiency η_AD = E_out / (E_out + E_core + E_cu) over 24 h.',
        'Core loss runs all 24 h: E_core = P_core × 24 = ' + fx(Pc) + ' × 24 = ' + fx(ad.Ecore) + ' kWh'];
      ad.rows.forEach((r, i) => L.push('  period ' + (i + 1) + ': ' + fx(r.h) + ' h at ' + fx(r.kW) + ' kW, pf ' + fx(r.pf) + ' → S = ' + fx(r.kW / r.pf) + ' kVA, x = ' + fx(r.x) + ';  E_out = ' + fx(r.Eout) + ' kWh;  E_cu = x²·P_cu·h = ' + fx(r.x) + '² × ' + fx(Pcu) + ' × ' + fx(r.h) + ' = ' + fx(r.Ecu) + ' kWh'));
      if (ad.idle > 0) L.push('  remaining ' + fx(ad.idle) + ' h: no load (E_out = 0, E_cu = 0)');
      L.push('E_out = ' + fx(ad.Eout) + ' kWh;  E_cu = ' + fx(ad.Ecu) + ' kWh;  η_AD = ' + fx(ad.Eout) + '/(' + fx(ad.Eout) + ' + ' + fx(ad.Ecore) + ' + ' + fx(ad.Ecu) + ') = ' + fx(100 * ad.eta) + ' %');
      L.push('Maximum efficiency where copper loss = core loss: x*² P_cu = P_core → x* = √(P_core/P_cu,FL) = √(' + fx(Pc) + '/' + fx(Pcu) + ') = ' + fx(me.x) + '  (' + fx(me.Skva) + ' kVA)');
      L.push('η_max = x*·S·pf / (x*·S·pf + 2 P_core) = ' + fx(100 * me.eta) + ' % at pf ' + fx(pf) + '.  Distribution transformers (lightly loaded most of the day) are designed with low core loss, i.e. small x*.');
      ui.work.set(L);
    }
    const io = stateIO(st, sels, { get: () => ({ prof: ui.prof.value.split(/\r?\n/).join(';') }), set(o) { if (typeof o.prof === 'string' && o.prof.length < 4000) ui.prof.value = o.prof.split(';').join('\n'); } });
    return { update, get: io.get, set: io.set };
  }

  /* ---------- tab ---------- */
  const host = { active: false, sub: 'mc' };
  let subs = null; const panels = {}, btns = {};
  function showSub(id, silent) {
    host.sub = id;
    Object.keys(panels).forEach(k => { panels[k].hidden = k !== id; btns[k].setAttribute('aria-pressed', String(k === id)); btns[k].classList.toggle('active', k === id); });
    if (host.active) subs[id].update();
    if (!silent) FSP.state.touch();
  }
  FSP.registerTab({
    id: 'edm2', title: 'Transformers & Actuators',
    init(panel) {
      const bar = FSP.ui.el('div', { class: 'seg row', role: 'group', 'aria-label': 'Transformers and actuators topic' }); panel.appendChild(bar);
      const defs = [['mc', 'Multi-limb magnetic circuit'], ['ie', 'Inductance, force & losses'], ['auto', 'Autotransformer'], ['tp', 'Three-phase banks'], ['ad', 'All-day efficiency']];
      defs.forEach(d => { btns[d[0]] = FSP.ui.el('button', { type: 'button', class: 'seg-btn', 'aria-pressed': 'false', text: d[1], onclick: () => showSub(d[0]) }); bar.appendChild(btns[d[0]]); });
      defs.forEach(d => { panels[d[0]] = FSP.ui.el('div', { hidden: '' }); panel.appendChild(panels[d[0]]); });
      subs = { mc: buildMC(panels.mc), ie: buildIE(panels.ie), auto: buildAuto(panels.auto), tp: buildTP(panels.tp), ad: buildAD(panels.ad) };
      const pre = { mc: 'mc_', ie: 'ie_', auto: 'at_', tp: 'tp_', ad: 'ad_' };
      FSP.state.bind('edm2', {
        get() { const o = { sub: host.sub }; Object.keys(subs).forEach(s => { const g = subs[s].get(); Object.keys(g).forEach(k => { o[pre[s] + k] = g[k]; }); }); return o; },
        set(o) {
          o = o || {};
          Object.keys(subs).forEach(s => { const r = {}; Object.keys(o).forEach(k => { if (k.indexOf(pre[s]) === 0) r[k.slice(pre[s].length)] = o[k]; }); try { subs[s].set(r); } catch (e) { /* ignore bad state */ } });
          showSub(panels[o.sub] ? o.sub : 'mc', true);
        },
      });
      showSub(host.sub, true);
    },
    activate() { host.active = true; if (subs) subs[host.sub].update(); },
    deactivate() { host.active = false; },
  });
})();
