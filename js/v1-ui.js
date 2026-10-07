/* ===================== UI ===================== */
'use strict';
const $ = id => document.getElementById(id);
function prepCanvas(cv) {
  const rect = cv.getBoundingClientRect(), dpr = window.devicePixelRatio || 1;
  const w = Math.max(1, Math.round(rect.width)), h = Math.max(1, Math.round(rect.height));
  const pw = Math.round(w * dpr), ph = Math.round(h * dpr);
  if (cv.width !== pw || cv.height !== ph) { cv.width = pw; cv.height = ph; }
  const ctx = cv.getContext('2d'); ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  return { ctx, w, h, dpr };
}
function clamp(v, a, b) { return Math.min(b, Math.max(a, v)); }
function arrow(ctx, x0, y0, x1, y1, head) {
  ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1, y1); ctx.stroke();
  const a = Math.atan2(y1 - y0, x1 - x0); head = head || 6;
  ctx.beginPath(); ctx.moveTo(x1, y1); ctx.lineTo(x1 - head * Math.cos(a - 0.45), y1 - head * Math.sin(a - 0.45)); ctx.lineTo(x1 - head * Math.cos(a + 0.45), y1 - head * Math.sin(a + 0.45)); ctx.closePath(); ctx.fill();
}
function labelBox(ctx, text, x, y, align) {
  ctx.font = '12px ' + getComputedStyle(document.body).getPropertyValue('--mono');
  const wt = ctx.measureText(text).width + 8; const xx = align === 'right' ? x - wt : align === 'center' ? x - wt / 2 : x;
  ctx.fillStyle = 'rgba(11,15,23,0.72)'; ctx.fillRect(xx, y - 13, wt, 17);
  ctx.fillStyle = '#e6edf3'; ctx.textAlign = 'left'; ctx.fillText(text, xx + 4, y);
}

/* ---------- tabs & animation loops ---------- */
let activeModule = 'em', emRaf = 0, dspRaf = 0;
FSP.mountTabs(document.querySelector('.tabs'), document.querySelector('main'));
const tabs = Array.from(document.querySelectorAll('[role=tab]'));
function activateTab(id, focus) {
  const next = id.replace(/^tab-/, '');
  tabs.forEach(t => { const on = t.id === id; t.setAttribute('aria-selected', on ? 'true' : 'false'); t.tabIndex = on ? 0 : -1; $(t.getAttribute('aria-controls')).hidden = !on; if (on && focus) t.focus(); });
  if (next !== activeModule) {
    const prev = activeModule; activeModule = next;
    stopLoops();
    const pd = FSP.tabs.find(d => d.id === prev);
    if (pd && pd.deactivate) { try { pd.deactivate(); } catch (e) { console.error(e); } }
    if (prev === 'dsp') stopAudio('tab switch');
    if (next === 'em') { em.dirty = true; }
    else if (next === 'dsp') { dspRedrawAll(); }
    else { const nd = FSP.tabs.find(d => d.id === next); if (nd && nd.activate) { try { nd.activate(); } catch (e) { console.error(e); } } }
    startLoops();
    FSP.state.setActive(next);
  }
}
tabs.forEach((t, i) => {
  t.addEventListener('click', () => activateTab(t.id, false));
  t.addEventListener('keydown', e => {
    let j = null;
    if (e.key === 'ArrowRight') j = (i + 1) % tabs.length; else if (e.key === 'ArrowLeft') j = (i - 1 + tabs.length) % tabs.length;
    else if (e.key === 'Home') j = 0; else if (e.key === 'End') j = tabs.length - 1;
    if (j !== null) { e.preventDefault(); activateTab(tabs[j].id, true); }
  });
});
function startLoops() {
  if (document.hidden) return;
  if (activeModule === 'em') { if (!emRaf) emRaf = requestAnimationFrame(emFrame); }
  else if (activeModule === 'dsp' && !dspRaf) dspRaf = requestAnimationFrame(dspFrame);
}
function stopLoops() { if (emRaf) cancelAnimationFrame(emRaf); if (dspRaf) cancelAnimationFrame(dspRaf); emRaf = 0; dspRaf = 0; em.lastTs = 0; }
document.addEventListener('visibilitychange', () => { if (document.hidden) stopLoops(); else startLoops(); });

/* ===================== MODULE 1: EM ===================== */
const em = { theta: 45, pol: 'TE', logf: 9, er1: 1, mr1: 1, s1: 0, er2: 2.25, mr2: 1, s2: 0, span: 10, speed: 1, playing: true,
  quality: 240, arrows: 'inst', field: 'A', phase: 0, lastTs: 0, dirty: true, cache: null, an: null, m1: null, m2: null, oc: null, img: null, hudT: 0, rtDirty: true };
const EM_LUT = (() => { // diverging colormap: blue (−1) → dark (0) → orange (+1)
  const stops = [[-1, [60, 140, 255]], [-0.5, [34, 70, 150]], [0, [14, 18, 30]], [0.5, [160, 60, 40]], [1, [255, 150, 70]]];
  const lut = new Uint8ClampedArray(256 * 3);
  for (let i = 0; i < 256; i++) {
    const v = i / 255 * 2 - 1; let k = 0; while (k < stops.length - 2 && v > stops[k + 1][0]) k++;
    const [v0, c0] = stops[k], [v1, c1] = stops[k + 1], f = (v - v0) / (v1 - v0);
    for (let c = 0; c < 3; c++) lut[i * 3 + c] = c0[c] + (c1[c] - c0[c]) * f;
  }
  return lut;
})();
const EM_PRESETS = {
  airglass: { er1: 1, mr1: 1, s1: 0, er2: 2.25, mr2: 1, s2: 0, theta: 45 },
  glassair: { er1: 2.25, mr1: 1, s1: 0, er2: 1, mr2: 1, s2: 0, theta: 60 },
  seawater: { er1: 1, mr1: 1, s1: 0, er2: 72, mr2: 1, s2: 4, theta: 45, logf: 9 },
  copper: { er1: 1, mr1: 1, s1: 0, er2: 1, mr2: 1, s2: 5.8e7, theta: 45, logf: 9 },
  magnetic: { er1: 1, mr1: 1, s1: 0, er2: 2, mr2: 3, s2: 0, theta: 37.8 }
};
function emReadInputs() {
  em.er1 = Math.max(1e-9, +$('em-er1').value || 1); em.mr1 = Math.max(1e-9, +$('em-mr1').value || 1); em.s1 = Math.max(0, +$('em-s1').value || 0);
  em.er2 = Math.max(1e-9, +$('em-er2').value || 1); em.mr2 = Math.max(1e-9, +$('em-mr2').value || 1); em.s2 = Math.max(0, +$('em-s2').value || 0);
  em.logf = +$('em-freq').value; em.span = +$('em-span').value;
  em.dirty = true; em.rtDirty = true;
}
function emSetTheta(v) { v = clamp(Math.round(v * 10) / 10, 0, 89); em.theta = v; $('em-theta').value = v; $('em-theta-n').value = v; em.dirty = true; }
function emSetPol(p) { em.pol = p; $('em-pol-te').setAttribute('aria-pressed', p === 'TE'); $('em-pol-tm').setAttribute('aria-pressed', p === 'TM'); $('em-field').disabled = p !== 'TM'; if (p !== 'TM') { em.field = 'A'; $('em-field').value = 'A'; } em.dirty = true; em.rtDirty = true; }
function emApplyPreset(p) {
  for (const k of ['er1', 'mr1', 's1', 'er2', 'mr2', 's2']) $('em-' + k).value = p[k];
  if (p.logf !== undefined) $('em-freq').value = p.logf;
  emReadInputs(); emSetTheta(p.theta);
}
['em-er1', 'em-mr1', 'em-s1', 'em-er2', 'em-mr2', 'em-s2', 'em-freq', 'em-span'].forEach(id => $(id).addEventListener('input', emReadInputs));
$('em-theta').addEventListener('input', e => emSetTheta(+e.target.value));
$('em-theta-n').addEventListener('change', e => emSetTheta(+e.target.value));
$('em-pol-te').addEventListener('click', () => emSetPol('TE'));
$('em-pol-tm').addEventListener('click', () => emSetPol('TM'));
$('em-presets').addEventListener('click', e => { const b = e.target.closest('[data-preset]'); if (b) emApplyPreset(EM_PRESETS[b.dataset.preset]); });
$('em-speed').addEventListener('input', e => { em.speed = +e.target.value; });
$('em-play').addEventListener('click', () => { em.playing = !em.playing; $('em-play').textContent = em.playing ? 'Pause' : 'Play'; $('em-play').setAttribute('aria-pressed', em.playing); });
$('em-quality').addEventListener('change', e => { em.quality = +e.target.value; em.dirty = true; });
$('em-arrows').addEventListener('change', e => { em.arrows = e.target.value; });
$('em-field').addEventListener('change', e => { em.field = e.target.value; });
$('em-field').disabled = true;
$('em-canvas').addEventListener('pointerdown', e => { e.currentTarget.setPointerCapture(e.pointerId); em.drag = true; emPointerTheta(e); });
$('em-canvas').addEventListener('pointermove', e => { if (em.drag) emPointerTheta(e); });
$('em-canvas').addEventListener('pointerup', () => { em.drag = false; });
$('em-canvas').addEventListener('pointercancel', () => { em.drag = false; });
function emPointerTheta(e) {
  const r = e.currentTarget.getBoundingClientRect(), x = e.clientX - r.left - r.width / 2, z = e.clientY - r.top - r.height / 2;
  if (z >= -2) return; // only above the interface (medium 1)
  emSetTheta(Math.atan2(Math.abs(x), -z) * 180 / Math.PI);
}
new ResizeObserver(() => { em.dirty = true; em.rtDirty = true; }).observe($('em-canvas'));

function emRecompute() {
  const f = Math.pow(10, em.logf);
  em.m1 = makeMedium(em.er1, em.mr1, em.s1, f); em.m2 = makeMedium(em.er2, em.mr2, em.s2, f);
  em.an = emAnalyze(em.theta, em.pol, em.m1, em.m2);
  const { w, h } = prepCanvas($('em-canvas'));
  const iw = em.quality, ih = Math.max(2, Math.round(iw * h / w));
  em.cache = buildFieldCache(em.an.sol, iw, ih, em.span, h / w);
  if (!em.oc || em.oc.width !== iw || em.oc.height !== ih) { em.oc = document.createElement('canvas'); em.oc.width = iw; em.oc.height = ih; em.octx = em.oc.getContext('2d'); em.img = em.octx.createImageData(iw, ih); }
  em.dirty = false;
  emUpdateHud();
  if (em.rtDirty) { emDrawRT(); em.rtDirty = false; }
}
function emUpdateHud() {
  const an = em.an, m1 = em.m1, m2 = em.m2, f = m1.freq;
  $('em-freq-o').value = fmtSI(f, 'Hz'); $('em-lambda0').textContent = fmtSI(an.lambda0, 'm'); $('em-span-o').value = em.span.toFixed(1);
  $('em-n1').textContent = m1.nReal.toFixed(4) + (m1.lossy ? '  (n_c = ' + fmtComplexN(m1.nc) + ')' : '');
  $('em-n2').textContent = m2.nReal.toFixed(4) + (m2.lossy ? '  (n_c = ' + fmtComplexN(m2.nc) + ')' : '');
  const badge = an.conserv === 'pass' ? '<span class="badge ok">PASS</span>' : an.conserv === 'fail' ? '<span class="badge bad">FAIL</span>' : '<span class="badge warn">n/a (σ1&gt;0: interference flux)</span>';
  const b = an.brewster;
  const rows = [
    ['Polarization', em.pol + (em.pol === 'TE' ? ' (E ∥ y)' : ' (H ∥ y)')],
    ['θi / θr', fmtDeg(em.theta, 1) + ' / ' + fmtDeg(em.theta, 1)],
    ['θt', an.evanescent ? 'evanescent' : fmtDeg(an.thetaT)],
    ['r = Γ', fmtPolar(an.r)],
    ['t = τ', fmtPolar(an.t)],
    ['refl. phase shift', fmtDeg(an.phaseDeg)],
    ['R = |r|²', fmtNum(an.R, 6)],
    ['T (Poynting flux)', fmtNum(an.T, 6)],
    ['R + T', fmtNum(an.sum, 10) + ' ' + badge],
    [b.kind === 'exact' ? 'Brewster θB (TM)' : 'pseudo-Brewster (min |r_TM|)', b.tm === null || !Number.isFinite(b.tm) ? 'N/A' : fmtDeg(b.tm, 3) + (b.kind === 'pseudo' ? ' (|r|=' + b.rmin.toFixed(4) + ')' : '')],
    ['Brewster θB (TE)', b.kind === 'exact' && b.te !== null ? fmtDeg(b.te, 3) : 'N/A'],
    ['critical θc', an.thetaC === null ? 'N/A' : fmtDeg(an.thetaC, 3)],
    [an.evanescent ? 'evanescent decay δ' : 'decay length δ', Number.isFinite(an.delta) ? fmtNum(an.delta, 4) + ' λ0 = ' + fmtSI(an.deltaPhys, 'm') : '— (propagating)'],
    ['λ0 / λ1 / λ2', fmtSI(an.lambda0, 'm') + ' / ' + fmtSI(an.lambda0 / Math.max(m1.nc.re, 1e-300), 'm') + ' / ' + fmtSI(an.lambda0 / Math.max(m2.nc.re, 1e-300), 'm')],
    ...(() => { const sd = skinDepth(m2); return sd ? [['skin depth δ (medium 2)', fmtSI(sd.exact, 'm') + '  (good-conductor ' + fmtSI(sd.approx, 'm') + ', ratio ' + fmtNum(sd.ratio, 6) + ')'], ['surface impedance Zs', fmtComplexN(sd.Zs) + ' Ω  (= (1+j)/(σδ); |η2| = ' + fmtNum(sd.eta.abs(), 4) + ' Ω)']] : []; })(),
    ['η1 / η2', fmtPolar(m1.eta).replace(/^([\d.]+)/, (s, a) => (+a).toFixed(1)) + ' Ω / ' + fmtPolar(m2.eta).replace(/^([\d.]+)/, (s, a) => (+a).toFixed(1)) + ' Ω'],
    ['ωt', (em.phase * 180 / Math.PI).toFixed(0) + '°' + (em.playing ? '' : ' (paused)')]
  ];
  $('em-hud').innerHTML = rows.map(r => '<div><span>' + r[0] + '</span><span>' + r[1] + '</span></div>').join('');
}
function emFrame(ts) {
  emRaf = 0;
  if (activeModule !== 'em' || document.hidden) return;
  if (em.lastTs && em.playing) em.phase = (em.phase + em.speed * Math.min(0.1, (ts - em.lastTs) / 1000) * Math.PI) % (2 * Math.PI);
  em.lastTs = ts;
  if (em.dirty) emRecompute();
  emRender();
  if (ts - em.hudT > 100) { em.hudT = ts; emUpdateHud(); emDrawRTMarker(); }
  emRaf = requestAnimationFrame(emFrame);
}
function emRender() {
  const cv = $('em-canvas'), { ctx, w, h } = prepCanvas(cv), cache = em.cache, an = em.an, sol = an.sol;
  const cphi = Math.cos(em.phase), sphi = Math.sin(em.phase), n = cache.iw * cache.ih, data = em.img.data;
  const emag = em.pol === 'TM' && em.field === 'E';
  const { Are, Aim, Bre, Bim, Dre, Dim } = cache;
  if (emag) {
    const inv = 127 / cache.vmaxE;
    for (let k = 0, p = 0; k < n; k++, p += 4) {
      const ex = Bre[k] * cphi - Bim[k] * sphi, ez = Dre[k] * cphi - Dim[k] * sphi;
      const idx = 128 + Math.min(127, Math.sqrt(ex * ex + ez * ez) * inv) | 0;
      data[p] = EM_LUT[idx * 3]; data[p + 1] = EM_LUT[idx * 3 + 1]; data[p + 2] = EM_LUT[idx * 3 + 2]; data[p + 3] = 255;
    }
  } else {
    const inv = 127.5 / cache.vmaxA;
    for (let k = 0, p = 0; k < n; k++, p += 4) {
      let v = (Are[k] * cphi - Aim[k] * sphi) * inv; if (v > 127.5) v = 127.5; else if (v < -127.5) v = -127.5;
      const idx = (v + 127.5) | 0;
      data[p] = EM_LUT[idx * 3]; data[p + 1] = EM_LUT[idx * 3 + 1]; data[p + 2] = EM_LUT[idx * 3 + 2]; data[p + 3] = 255;
    }
  }
  em.octx.putImageData(em.img, 0, 0);
  ctx.imageSmoothingEnabled = true; ctx.drawImage(em.oc, 0, 0, w, h);

  const scale = w / em.span, ox = w / 2, oy = h / 2;
  // Poynting arrows
  if (em.arrows !== 'off') {
    const g = 40, pts = [];
    let smax = 0;
    for (let gy = g / 2; gy < h; gy += g) for (let gx = g / 2; gx < w; gx += g) {
      const i = Math.min(cache.iw - 1, (gx / w * cache.iw) | 0), j = Math.min(cache.ih - 1, (gy / h * cache.ih) | 0), k = j * cache.iw + i;
      const F = { A: C(Are[k], Aim[k]), B: C(Bre[k], Bim[k]), D: C(Dre[k], Dim[k]) };
      const S = em.arrows === 'avg' ? poyntingAvg(em.pol, F) : poyntingInst(em.pol, F, em.phase);
      const m = Math.hypot(S.x, S.z); if (Number.isFinite(m)) { pts.push([gx, gy, S.x, S.z, m]); if (m > smax) smax = m; }
    }
    ctx.strokeStyle = 'rgba(255,255,255,0.75)'; ctx.fillStyle = 'rgba(255,255,255,0.75)'; ctx.lineWidth = 1.2;
    const cap = g * 0.85;
    for (const [gx, gy, sx, sz, m] of pts) {
      if (!(m > 1e-9 * smax) || smax === 0) continue;
      const L = Math.min(cap, cap * m / smax); if (L < 2) continue;
      arrow(ctx, gx - sx / m * L / 2, gy - sz / m * L / 2, gx + sx / m * L / 2, gy + sz / m * L / 2, 5);
    }
  }
  // interface, normal
  ctx.lineWidth = 2; ctx.strokeStyle = '#e6edf3'; ctx.beginPath(); ctx.moveTo(0, oy); ctx.lineTo(w, oy); ctx.stroke();
  ctx.setLineDash([6, 5]); ctx.lineWidth = 1; ctx.strokeStyle = 'rgba(230,237,243,0.6)'; ctx.beginPath(); ctx.moveTo(ox, 0); ctx.lineTo(ox, h); ctx.stroke(); ctx.setLineDash([]);
  ctx.fillStyle = 'rgba(230,237,243,0.8)'; arrow(ctx, ox, oy, ox, oy - 34, 7);
  labelBox(ctx, 'n̂', ox + 8, oy - 36);
  // rays
  const th = em.theta * Math.PI / 180, L = Math.min(w, h) * 0.46, si = Math.sin(th), ci = Math.cos(th);
  const drawRay = (x0, y0, x1, y1, color, label, lx, ly) => { ctx.strokeStyle = color; ctx.fillStyle = color; ctx.lineWidth = 2; arrow(ctx, x0, y0, x1, y1, 9); if (label && w >= 720) labelBox(ctx, label, clamp(lx, 40, w - 40), clamp(ly, 50, h - 24), 'center'); };
  const wavefronts = (dirx, diry, kvec, fromX, fromY, len) => { // ticks perpendicular to the ray, spaced by λ along it
    const kmag = Math.hypot(kvec.re, kvec.im2); if (!(kmag > 1e-9)) return;
    const lam = 2 * Math.PI / kmag * scale; if (lam < 4) return;
    ctx.strokeStyle = 'rgba(255,255,255,0.35)'; ctx.lineWidth = 1;
    for (let d = lam; d < len; d += lam) { const x = fromX + dirx * d, y = fromY + diry * d; ctx.beginPath(); ctx.moveTo(x - diry * 12, y + dirx * 12); ctx.lineTo(x + diry * 12, y - dirx * 12); ctx.stroke(); }
  };
  // incident (from upper-left to origin), reflected (origin to upper-right)
  drawRay(ox - L * si, oy - L * ci, ox, oy, '#4cc9f0', 'incident', ox - L * si - 10, oy - L * ci - 6);
  wavefronts(si, ci, { re: sol.kx.re, im2: sol.kz1.re }, ox - L * si, oy - L * ci, L);
  drawRay(ox, oy, ox + L * si, oy - L * ci, '#f9a03f', 'reflected', ox + L * si + 10, oy - L * ci - 6);
  wavefronts(si, -ci, { re: sol.kx.re, im2: sol.kz1.re }, ox, oy, L);
  if (!an.evanescent && Number.isFinite(an.thetaT)) {
    const tt = an.thetaT * Math.PI / 180, st = Math.sin(tt), ct = Math.cos(tt);
    drawRay(ox, oy, ox + L * st, oy + L * ct, '#f72585', 'transmitted', ox + L * st + 14, oy + L * ct + 14);
    wavefronts(st, ct, { re: sol.kx.re, im2: sol.kz2.re }, ox, oy, L);
    ctx.strokeStyle = '#f72585'; ctx.lineWidth = 1.2; ctx.beginPath(); ctx.arc(ox, oy, 48, Math.PI / 2, Math.PI / 2 - tt, true); ctx.stroke();
    labelBox(ctx, 'θt=' + fmtDeg(an.thetaT, 1), ox + 62 * Math.sin(tt / 2) + 4, oy + 62 * Math.cos(tt / 2) + 6);
  } else {
    ctx.strokeStyle = '#f72585'; ctx.fillStyle = '#f72585'; ctx.lineWidth = 2; ctx.setLineDash([8, 6]); arrow(ctx, ox, oy + 8, ox + L * 0.8, oy + 8, 9); ctx.setLineDash([]);
    labelBox(ctx, 'evanescent (no propagating θt)', ox + 10, oy + 28);
  }
  // angle arcs θi, θr
  ctx.lineWidth = 1.2; ctx.strokeStyle = '#4cc9f0'; ctx.beginPath(); ctx.arc(ox, oy, 48, -Math.PI / 2, -Math.PI / 2 - th, true); ctx.stroke();
  labelBox(ctx, 'θi=' + fmtDeg(em.theta, 1), ox - 64 * Math.sin(th / 2) - 4, oy - 64 * Math.cos(th / 2) + 4, 'right');
  ctx.strokeStyle = '#f9a03f'; ctx.beginPath(); ctx.arc(ox, oy, 48, -Math.PI / 2, -Math.PI / 2 + th, false); ctx.stroke();
  labelBox(ctx, 'θr=' + fmtDeg(em.theta, 1), ox + 64 * Math.sin(th / 2) + 4, oy - 64 * Math.cos(th / 2) + 4);
  // decay envelope
  if (Number.isFinite(an.delta)) {
    const dpx = an.delta * scale;
    ctx.setLineDash([4, 4]); ctx.strokeStyle = 'rgba(255,255,255,0.5)'; ctx.lineWidth = 1;
    for (let m = 1; m <= 3; m++) { const y = oy + m * dpx; if (y > h - 4) break; ctx.beginPath(); ctx.moveTo(0, y); ctx.lineTo(w, y); ctx.stroke(); labelBox(ctx, 'z = ' + m + 'δ', 8, y - 3); }
    ctx.setLineDash([]);
    ctx.strokeStyle = '#ffd166'; ctx.lineWidth = 1.5; ctx.beginPath();
    for (let y = oy; y <= h; y += 2) { const x = w - 12 - 56 * Math.exp(-(y - oy) / Math.max(dpx, 1e-6)); if (y === oy) ctx.moveTo(x, y); else ctx.lineTo(x, y); }
    ctx.stroke(); labelBox(ctx, 'e^(−z/δ)', w - 8, oy + 16, 'right');
  }
  // medium labels & shown field
  const m1 = em.m1, m2 = em.m2;
  labelBox(ctx, 'Medium 1  εr=' + m1.er + '  μr=' + m1.mr + '  σ=' + m1.sigma + ' S/m  n=' + m1.nReal.toFixed(3), 8, 18);
  labelBox(ctx, 'Medium 2  εr=' + m2.er + '  μr=' + m2.mr + '  σ=' + m2.sigma + ' S/m  n=' + m2.nReal.toFixed(3), 8, h - 8);
  const showTxt = 'showing ' + (emag ? '|E| (Ex,Ez)' : em.pol === 'TE' ? 'Ey' : 'Hy') + ' · ' + em.pol + ' · arrows: ' + (em.arrows === 'avg' ? '⟨S⟩' : em.arrows === 'inst' ? 'S(t)' : 'off');
  if (w < 720) labelBox(ctx, showTxt, 8, 36); else labelBox(ctx, showTxt, w - 8, 18, 'right');
}
function emPower(th, pol, m1, m2) { // R and T only (cheap, for the sweep plot)
  const sol = solveInterface(th, pol, m1, m2);
  const SzInc = poyntingAvg(pol, fieldPhasors(sol, 0, 0, 'inc')).z, SzTr = poyntingAvg(pol, fieldPhasors(sol, 0, 0, 'tr')).z;
  return { R: sol.r.isFinite() ? sol.r.abs2() : NaN, T: Math.abs(SzInc) > 1e-300 ? SzTr / SzInc : NaN };
}
function emDrawRT() { // curves rendered once into an offscreen canvas; the marker is drawn per frame
  const cv = $('em-rt'), rect = cv.getBoundingClientRect(), dpr = window.devicePixelRatio || 1;
  const w = Math.max(1, Math.round(rect.width)), h = Math.max(1, Math.round(rect.height));
  if (!em.rtOc) em.rtOc = document.createElement('canvas');
  em.rtOc.width = Math.round(w * dpr); em.rtOc.height = Math.round(h * dpr);
  const ctx = em.rtOc.getContext('2d'); ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  const L = 40, R = 12, T = 10, B = 22, pw = w - L - R, ph = h - T - B;
  ctx.fillStyle = '#121826'; ctx.fillRect(0, 0, w, h);
  ctx.strokeStyle = '#243047'; ctx.lineWidth = 1; ctx.font = '11px ' + getComputedStyle(document.body).getPropertyValue('--mono'); ctx.fillStyle = '#8b98ab';
  for (let v = 0; v <= 1; v += 0.25) { const y = T + ph * (1 - v); ctx.beginPath(); ctx.moveTo(L, y); ctx.lineTo(w - R, y); ctx.stroke(); ctx.textAlign = 'right'; ctx.fillText(v.toFixed(2), L - 4, y + 4); }
  for (let d = 0; d <= 90; d += 15) { const x = L + pw * d / 90; ctx.beginPath(); ctx.moveTo(x, T); ctx.lineTo(x, T + ph); ctx.stroke(); ctx.textAlign = 'center'; ctx.fillText(d + '°', x, h - 6); }
  const series = [['TE', 'R', '#4cc9f0'], ['TE', 'T', '#2ecc71'], ['TM', 'R', '#f9a03f'], ['TM', 'T', '#f72585']];
  const N = 181;
  for (const [pol, key, col] of series) {
    ctx.strokeStyle = col; ctx.lineWidth = 1.5; ctx.beginPath(); let started = false;
    for (let i = 0; i < N; i++) { const th = 89.9 * i / (N - 1); const an = emPower(th, pol, em.m1, em.m2); const v = an[key]; if (!Number.isFinite(v)) { started = false; continue; } const x = L + pw * th / 90, y = T + ph * (1 - clamp(v, 0, 1)); if (started) ctx.lineTo(x, y); else { ctx.moveTo(x, y); started = true; } }
    ctx.stroke();
  }
  ctx.textAlign = 'left'; let lx = L + 6; for (const [pol, key, col] of series) { ctx.fillStyle = col; ctx.fillText(key + '_' + pol, lx, T + 12); lx += 46; }
  em.rtCanvasInfo = { L, pw, T, ph, w, h };
}
function emDrawRTMarker() {
  const cv = $('em-rt'), { ctx, w, h } = prepCanvas(cv), i = em.rtCanvasInfo; if (!i || !em.rtOc) return;
  ctx.drawImage(em.rtOc, 0, 0, w, h);
  const x = i.L + i.pw * em.theta / 90;
  ctx.strokeStyle = '#ffffff'; ctx.lineWidth = 1; ctx.setLineDash([3, 3]); ctx.beginPath(); ctx.moveTo(x, i.T); ctx.lineTo(x, i.T + i.ph); ctx.stroke(); ctx.setLineDash([]);
  const an = em.an; if (!an) return;
  ctx.fillStyle = em.pol === 'TE' ? '#4cc9f0' : '#f9a03f'; ctx.beginPath(); ctx.arc(x, i.T + i.ph * (1 - clamp(an.R, 0, 1)), 4, 0, 7); ctx.fill();
  ctx.fillStyle = em.pol === 'TE' ? '#2ecc71' : '#f72585'; ctx.beginPath(); ctx.arc(x, i.T + i.ph * (1 - clamp(an.T, 0, 1)), 4, 0, 7); ctx.fill();
}
/* ===================== MODULE 2: DSP ===================== */
const dsp = { type: 'lp', fc: 1000, fs: 48000, Q: 0.7071, gainDb: 0, N: 128, custom: false, normalize: true, autoRange: false,
  coefs: null, roots: null, refLevel: 1, resp: null, h: null, hLen: 128, stable: true, lastStable: null, sel: null, drag: null, cursorF: null };
const MONO = () => getComputedStyle(document.body).getPropertyValue('--mono');
function dspMaxLogFc() { return Math.log10(0.49 * dsp.fs); }
function dspSetFc(f) {
  f = clamp(f, 10, 0.49 * dsp.fs); dsp.fc = f;
  $('dsp-fc').max = dspMaxLogFc().toFixed(4); $('dsp-fc').value = Math.log10(f); $('dsp-fc-n').value = +f.toPrecision(6); $('dsp-fc-n').max = 0.49 * dsp.fs;
}
function dspSetFs(fs, fromAudio) {
  const sel = $('dsp-fs');
  if (!Array.from(sel.options).some(o => +o.value === fs)) { const o = document.createElement('option'); o.value = fs; o.textContent = fs; sel.appendChild(o); }
  sel.value = fs; dsp.fs = fs; dspSetFc(dsp.fc);
  $('dsp-fs-note').textContent = fromAudio ? 'locked to AudioContext rate' : '';
  sel.disabled = !!fromAudio;
}
function dspApplyCookbook() {
  dsp.custom = false;
  dsp.coefs = rbjCoefs(dsp.type, dsp.fc, dsp.fs, dsp.Q, dsp.gainDb);
  dsp.roots = rootsFromCoefs(dsp.coefs);
  dsp.refLevel = freqResponseAt(dsp.coefs, 2 * Math.PI * refFreqFor(dsp.type, dsp.fc, dsp.fs) / dsp.fs).abs();
  dspUpdateAll();
}
function dspApplyRoots() {
  dsp.custom = true;
  const wref = 2 * Math.PI * refFreqFor(dsp.type, dsp.fc, dsp.fs) / dsp.fs;
  if (dsp.normalize) {
    const unit = freqResponseAt(coefsFromRoots({ poles: dsp.roots.poles, zeros: dsp.roots.zeros, K: 1 }), wref).abs();
    if (unit > 1e-12 && dsp.refLevel > 1e-12) dsp.roots.K = dsp.refLevel / unit;
  }
  if (!(Math.abs(dsp.roots.K) > 1e-12)) dsp.roots.K = 1e-12;
  dsp.coefs = coefsFromRoots(dsp.roots);
  dspUpdateAll();
}
function dspUpdateAll() {
  const c = dsp.coefs;
  dsp.stable = isStable(c);
  if (dsp.stable) { dsp.lastStable = c; audioPushCoefs(c); }
  dsp.resp = freqResponse(c, dsp.fs, 1024, 20);
  dsp.h = impulseResponse(c, dsp.N); dsp.hLen = Math.min(dsp.N, trimLength(dsp.h, 1e-6));
  const f = x => (x >= 0 ? ' ' : '') + x.toFixed(9);
  $('dsp-coefs').textContent = 'b0 =' + f(c.b0) + '\nb1 =' + f(c.b1) + '\nb2 =' + f(c.b2) + '\na0 =' + f(1) + '\na1 =' + f(c.a1) + '\na2 =' + f(c.a2) + '\n(pre-normalization a0 = ' + c.a0raw.toFixed(9) + ')';
  const rad = poleRadii(c);
  $('dsp-stable').textContent = dsp.stable ? 'STABLE' : 'UNSTABLE — audio keeps last stable set';
  $('dsp-stable').className = 'badge ' + (dsp.stable ? 'ok' : 'bad');
  $('dsp-poles-note').textContent = '|p| = ' + rad.map(x => x.toFixed(5)).join(', ') + (dsp.roots.zeroNote ? ' · ' + dsp.roots.zeroNote : '');
  $('dsp-custom').hidden = !dsp.custom;
  document.querySelectorAll('#dsp-types .seg-btn').forEach(b => b.classList.toggle('active', !dsp.custom && b.dataset.type === dsp.type));
  dspRedrawAll();
}
function dspRedrawAll() { if (!dsp.resp) return; dspDrawMag(); dspDrawPhase(); dspDrawImpulse(); dspDrawZplane(); }
// controls
document.querySelectorAll('#dsp-types .seg-btn').forEach(b => b.addEventListener('click', () => { dsp.type = b.dataset.type; dspApplyCookbook(); }));
$('dsp-fc').addEventListener('input', e => { dspSetFc(Math.pow(10, +e.target.value)); dspApplyCookbook(); });
$('dsp-fc-n').addEventListener('change', e => { dspSetFc(+e.target.value || 1000); dspApplyCookbook(); });
$('dsp-fs').addEventListener('change', e => { dspSetFs(+e.target.value, false); dspApplyCookbook(); });
$('dsp-q').addEventListener('input', e => { dsp.Q = +e.target.value; $('dsp-q-n').value = dsp.Q; dspApplyCookbook(); });
$('dsp-q-n').addEventListener('change', e => { dsp.Q = clamp(+e.target.value || 0.7071, 0.1, 20); $('dsp-q').value = dsp.Q; $('dsp-q-n').value = dsp.Q; dspApplyCookbook(); });
$('dsp-gain').addEventListener('input', e => { dsp.gainDb = +e.target.value; $('dsp-gain-n').value = dsp.gainDb; dspApplyCookbook(); });
$('dsp-gain-n').addEventListener('change', e => { dsp.gainDb = clamp(+e.target.value || 0, -24, 24); $('dsp-gain').value = dsp.gainDb; $('dsp-gain-n').value = dsp.gainDb; dspApplyCookbook(); });
$('dsp-n').addEventListener('input', e => { dsp.N = +e.target.value; $('dsp-n-o').value = dsp.N; dsp.h = impulseResponse(dsp.coefs, dsp.N); dsp.hLen = Math.min(dsp.N, trimLength(dsp.h, 1e-6)); dspDrawImpulse(); });
$('dsp-autorange').addEventListener('change', e => { dsp.autoRange = e.target.checked; dspDrawMag(); });
$('dsp-normgain').addEventListener('change', e => { dsp.normalize = e.target.checked; });

/* ---------- plots ---------- */
function logX(fmin, fmax, x0, x1) { const k = (x1 - x0) / Math.log10(fmax / fmin); return f => x0 + Math.log10(f / fmin) * k; }
function invLogX(fmin, fmax, x0, x1) { return x => fmin * Math.pow(fmax / fmin, (x - x0) / (x1 - x0)); }
const F_TICKS = [10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000];
const fmtHz = f => f >= 1000 ? (f / 1000).toFixed(f >= 10000 ? 0 : 1).replace(/\.0$/, '') + 'k' : f.toFixed(0);
function plotFrame(cv, L, R, T, B) {
  const { ctx, w, h } = prepCanvas(cv);
  ctx.fillStyle = '#121826'; ctx.fillRect(0, 0, w, h); ctx.font = '11px ' + MONO(); ctx.lineWidth = 1;
  return { ctx, w, h, L, R, T, B, pw: w - L - R, ph: h - T - B };
}
function drawLogGrid(P, fmin, fmax, nice) {
  const { ctx } = P, xf = logX(fmin, fmax, P.L, P.L + P.pw);
  ctx.strokeStyle = '#243047'; ctx.fillStyle = '#8b98ab'; ctx.textAlign = 'center';
  for (const f of F_TICKS) { if (f < fmin || f > fmax) continue; const x = xf(f); ctx.beginPath(); ctx.moveTo(x, P.T); ctx.lineTo(x, P.T + P.ph); ctx.stroke(); ctx.fillText(fmtHz(f), x, P.h - 6); }
  if (nice) { for (let d = 1; d < 10; d++) for (const e of [10, 100, 1000, 10000]) { const f = d * e; if (f > fmin && f < fmax && !F_TICKS.includes(f)) { const x = xf(f); ctx.strokeStyle = '#1a2233'; ctx.beginPath(); ctx.moveTo(x, P.T); ctx.lineTo(x, P.T + P.ph); ctx.stroke(); } } }
  return xf;
}
function drawYGrid(P, y0, y1, step, unit) {
  const { ctx } = P, yf = v => P.T + P.ph * (1 - (v - y0) / (y1 - y0));
  ctx.strokeStyle = '#243047'; ctx.fillStyle = '#8b98ab'; ctx.textAlign = 'right';
  for (let v = Math.ceil(y0 / step) * step; v <= y1 + 1e-9; v += step) { const y = yf(v); ctx.beginPath(); ctx.moveTo(P.L, y); ctx.lineTo(P.L + P.pw, y); ctx.stroke(); ctx.fillText(v.toFixed(0) + unit, P.L - 4, y + 4); }
  return yf;
}
function fcMarker(P, xf) { const x = xf(dsp.fc); P.ctx.strokeStyle = 'rgba(249,160,63,0.8)'; P.ctx.setLineDash([4, 4]); P.ctx.beginPath(); P.ctx.moveTo(x, P.T); P.ctx.lineTo(x, P.T + P.ph); P.ctx.stroke(); P.ctx.setLineDash([]); P.ctx.fillStyle = '#f9a03f'; P.ctx.textAlign = 'left'; P.ctx.fillText('fc', x + 3, P.T + 11); }
function cursorMarker(P, xf) { if (dsp.cursorF === null) return; const x = xf(dsp.cursorF); P.ctx.strokeStyle = 'rgba(255,255,255,0.7)'; P.ctx.beginPath(); P.ctx.moveTo(x, P.T); P.ctx.lineTo(x, P.T + P.ph); P.ctx.stroke(); }
function dspDrawMag() {
  const P = plotFrame($('dsp-mag'), 46, 12, 10, 20), r = dsp.resp;
  let y0 = -60, y1 = 30;
  if (dsp.autoRange) { let mn = Infinity, mx = -Infinity; for (const v of r.mag) { if (v < mn) mn = v; if (v > mx) mx = v; } y0 = Math.floor((mn - 5) / 10) * 10; y1 = Math.ceil((mx + 5) / 10) * 10; if (y1 - y0 < 20) y1 = y0 + 20; }
  const xf = drawLogGrid(P, r.fmin, r.fmax, true), yf = drawYGrid(P, y0, y1, (y1 - y0) > 100 ? 20 : 10, '');
  P.ctx.fillStyle = '#8b98ab'; P.ctx.textAlign = 'left'; P.ctx.fillText('|H| dB', P.L + 4, P.T + 11);
  fcMarker(P, xf);
  P.ctx.strokeStyle = '#4cc9f0'; P.ctx.lineWidth = 2; P.ctx.beginPath();
  for (let i = 0; i < r.freqs.length; i++) { const x = xf(r.freqs[i]), y = clamp(yf(r.mag[i]), P.T - 2, P.T + P.ph + 2); if (i) P.ctx.lineTo(x, y); else P.ctx.moveTo(x, y); }
  P.ctx.stroke(); cursorMarker(P, xf);
}
function dspDrawPhase() {
  const P = plotFrame($('dsp-phase'), 46, 12, 10, 20), r = dsp.resp;
  let mn = Infinity, mx = -Infinity; for (const v of r.phase) { if (v < mn) mn = v; if (v > mx) mx = v; }
  const y0 = Math.floor(mn / 45) * 45, y1 = Math.max(y0 + 45, Math.ceil(mx / 45) * 45);
  const xf = drawLogGrid(P, r.fmin, r.fmax, false), yf = drawYGrid(P, y0, y1, (y1 - y0) > 360 ? 90 : 45, '°');
  P.ctx.fillStyle = '#8b98ab'; P.ctx.textAlign = 'left'; P.ctx.fillText('∠H (unwrapped)', P.L + 4, P.T + 11);
  fcMarker(P, xf);
  P.ctx.strokeStyle = '#f72585'; P.ctx.lineWidth = 2; P.ctx.beginPath();
  for (let i = 0; i < r.freqs.length; i++) { const x = xf(r.freqs[i]), y = yf(r.phase[i]); if (i) P.ctx.lineTo(x, y); else P.ctx.moveTo(x, y); }
  P.ctx.stroke(); cursorMarker(P, xf);
}
function dspDrawImpulse() {
  const P = plotFrame($('dsp-impulse'), 46, 12, 10, 20), h = dsp.h, n = dsp.hLen;
  let mx = 1e-12; for (let i = 0; i < n; i++) mx = Math.max(mx, Math.abs(h[i]));
  const yf = v => P.T + P.ph * (1 - (v + mx) / (2 * mx)), xf = i => P.L + P.pw * (i + 0.5) / n;
  const { ctx } = P; ctx.strokeStyle = '#243047'; ctx.fillStyle = '#8b98ab';
  for (const v of [-1, -0.5, 0, 0.5, 1]) { const y = yf(v * mx); ctx.beginPath(); ctx.moveTo(P.L, y); ctx.lineTo(P.L + P.pw, y); ctx.stroke(); ctx.textAlign = 'right'; ctx.fillText((v * mx).toPrecision(2), P.L - 4, y + 4); }
  ctx.textAlign = 'center'; const step = n <= 32 ? 4 : n <= 128 ? 16 : n <= 512 ? 64 : 128;
  for (let i = 0; i < n; i += step) ctx.fillText(String(i), xf(i), P.h - 6);
  labelBox(ctx, 'h[n] via DF-II-T · showing ' + n + ' of ' + dsp.N + (n < dsp.N ? ' (auto-trimmed: |h| < 1e-6 beyond)' : ''), P.L + 4, P.T + 13);
  ctx.strokeStyle = '#2ecc71'; ctx.fillStyle = '#2ecc71'; ctx.lineWidth = n > 300 ? 1 : 1.5;
  const y0 = yf(0), rad = n > 300 ? 1 : n > 128 ? 1.5 : 2.5;
  for (let i = 0; i < n; i++) { const x = xf(i), y = yf(h[i]); ctx.beginPath(); ctx.moveTo(x, y0); ctx.lineTo(x, y); ctx.stroke(); ctx.beginPath(); ctx.arc(x, y, rad, 0, 7); ctx.fill(); }
}
function bodeCursor(e) {
  const cv = e.currentTarget, rect = cv.getBoundingClientRect(), r = dsp.resp; if (!r) return;
  const x = e.clientX - rect.left, L = 46, pw = rect.width - 46 - 12;
  if (x < L || x > L + pw) { dsp.cursorF = null; } else { dsp.cursorF = invLogX(r.fmin, r.fmax, L, L + pw)(x); }
  if (dsp.cursorF !== null) {
    let i = Math.round(Math.log10(dsp.cursorF / r.fmin) / Math.log10(r.fmax / r.fmin) * (r.freqs.length - 1)); i = clamp(i, 0, r.freqs.length - 1);
    $('dsp-cursor').textContent = 'f = ' + fmtSI(r.freqs[i], 'Hz', 4) + '   |H| = ' + r.mag[i].toFixed(2) + ' dB   ∠H = ' + r.phase[i].toFixed(1) + '°';
  } else $('dsp-cursor').textContent = 'Hover or touch the Bode plots for a readout.';
  dspDrawMag(); dspDrawPhase();
}
for (const id of ['dsp-mag', 'dsp-phase']) { const cv = $(id); cv.addEventListener('pointermove', bodeCursor); cv.addEventListener('pointerdown', bodeCursor); cv.addEventListener('pointerleave', () => { dsp.cursorF = null; $('dsp-cursor').textContent = 'Hover or touch the Bode plots for a readout.'; dspDrawMag(); dspDrawPhase(); }); }
/* ---------- z-plane editor ---------- */
function zGeom() {
  const cv = $('dsp-zplane'), rect = cv.getBoundingClientRect(), w = rect.width, h = rect.height;
  let maxZ = 1; for (const z of dsp.roots.zeros) maxZ = Math.max(maxZ, z.abs());
  const view = Math.max(1.5, maxZ + 0.3), R = Math.max(1, (Math.min(w, h) / 2 - 14) / view);
  return { w, h, cx: w / 2, cy: h / 2, R, view };
}
function dspDrawZplane() {
  if ($('dsp-zplane').getBoundingClientRect().width < 40) return; // hidden panel
  const { ctx, w, h } = prepCanvas($('dsp-zplane')), g = zGeom(), { cx, cy, R } = g;
  ctx.fillStyle = '#121826'; ctx.fillRect(0, 0, w, h); ctx.font = '11px ' + MONO();
  ctx.strokeStyle = '#1e293b'; ctx.lineWidth = 1;
  for (let v = -Math.floor(g.view * 2) / 2; v <= g.view; v += 0.5) { const x = cx + v * R, y = cy - v * R; ctx.beginPath(); ctx.moveTo(x, 0); ctx.lineTo(x, h); ctx.stroke(); ctx.beginPath(); ctx.moveTo(0, y); ctx.lineTo(w, y); ctx.stroke(); }
  ctx.strokeStyle = '#3b4a63'; ctx.beginPath(); ctx.moveTo(0, cy); ctx.lineTo(w, cy); ctx.stroke(); ctx.beginPath(); ctx.moveTo(cx, 0); ctx.lineTo(cx, h); ctx.stroke();
  ctx.strokeStyle = '#8b98ab'; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.arc(cx, cy, R, 0, 7); ctx.stroke();
  ctx.fillStyle = '#8b98ab'; ctx.textAlign = 'left'; ctx.fillText('Re', w - 20, cy - 4); ctx.fillText('Im', cx + 4, 12); ctx.fillText('|z|=1', cx + R * 0.72, cy - R * 0.72);
  ctx.fillText('z-plane · ' + (dsp.custom ? 'custom roots' : FILTER_NAMES[dsp.type]), 6, h - 6);
  const drawSet = (arr, kind) => {
    const dbl = arr[0].sub(arr[1]).abs() < 1e-7;
    arr.forEach((z, i) => {
      if (dbl && i === 1) return;
      const x = cx + z.re * R, y = cy - z.im * R, sel = dsp.sel && dsp.sel.kind === kind && (dsp.sel.idx === i || (dbl && dsp.sel.idx === 1));
      ctx.lineWidth = sel ? 3 : 2; ctx.strokeStyle = kind === 'pole' ? '#ff5c5c' : '#4cc9f0';
      if (sel) { ctx.fillStyle = 'rgba(255,255,255,0.12)'; ctx.beginPath(); ctx.arc(x, y, 13, 0, 7); ctx.fill(); }
      if (kind === 'pole') { ctx.beginPath(); ctx.moveTo(x - 6, y - 6); ctx.lineTo(x + 6, y + 6); ctx.moveTo(x + 6, y - 6); ctx.lineTo(x - 6, y + 6); ctx.stroke(); }
      else { ctx.beginPath(); ctx.arc(x, y, 6, 0, 7); ctx.stroke(); }
      if (dbl) { ctx.fillStyle = ctx.strokeStyle; ctx.fillText('×2', x + 8, y - 8); }
    });
  };
  drawSet(dsp.roots.zeros, 'zero'); drawSet(dsp.roots.poles, 'pole');
}
function zHit(px, py) {
  const g = zGeom(); let best = null, bd = 14;
  for (const kind of ['pole', 'zero']) (kind === 'pole' ? dsp.roots.poles : dsp.roots.zeros).forEach((z, i) => { const d = Math.hypot(g.cx + z.re * g.R - px, g.cy - z.im * g.R - py); if (d < bd) { bd = d; best = { kind, idx: i }; } });
  return best;
}
function zMoveSelected(x, y, snapTol) {
  applyRootMove(dsp.roots, dsp.sel.kind, dsp.sel.idx, x, y, snapTol);
  dspApplyRoots();
  const z = (dsp.sel.kind === 'pole' ? dsp.roots.poles : dsp.roots.zeros)[dsp.sel.idx];
  $('dsp-sel-note').textContent = 'Selected: ' + dsp.sel.kind + ' ' + (dsp.sel.idx + 1) + ' at ' + fmtComplexN(z) + '  |z| = ' + z.abs().toFixed(4);
}
{
  const cv = $('dsp-zplane');
  cv.addEventListener('pointerdown', e => {
    const rect = cv.getBoundingClientRect(), px = e.clientX - rect.left, py = e.clientY - rect.top;
    const hit = zHit(px, py); dsp.sel = hit; cv.focus();
    if (hit) { dsp.drag = true; cv.setPointerCapture(e.pointerId); zMoveSelected(...zFromPx(px, py), 8 / zGeom().R); }
    else { $('dsp-sel-note').textContent = 'Selected: none'; dspDrawZplane(); }
  });
  cv.addEventListener('pointermove', e => { if (!dsp.drag || !dsp.sel) return; const rect = cv.getBoundingClientRect(); zMoveSelected(...zFromPx(e.clientX - rect.left, e.clientY - rect.top), 8 / zGeom().R); });
  const end = () => { dsp.drag = false; };
  cv.addEventListener('pointerup', end); cv.addEventListener('pointercancel', end);
  cv.addEventListener('keydown', e => {
    if (!dsp.sel) return;
    const step = e.shiftKey ? 0.02 : 0.005; let dx = 0, dy = 0;
    if (e.key === 'ArrowLeft') dx = -step; else if (e.key === 'ArrowRight') dx = step; else if (e.key === 'ArrowUp') dy = step; else if (e.key === 'ArrowDown') dy = -step; else return;
    e.preventDefault();
    const z = (dsp.sel.kind === 'pole' ? dsp.roots.poles : dsp.roots.zeros)[dsp.sel.idx];
    let ny = z.im + dy; if (z.im !== 0 && Math.sign(ny) !== Math.sign(z.im)) ny = 0; // crossing the axis merges
    zMoveSelected(z.re + dx, ny, 1e-12);
  });
  function zFromPx(px, py) { const g = zGeom(); return [(px - g.cx) / g.R, (g.cy - py) / g.R]; }
}
new ResizeObserver(() => dspRedrawAll()).observe($('panel-dsp'));

/* ---------- audio ---------- */
const audio = { ctx: null, mode: null, node: null, noise: null, iir: { node: null, gain: null }, analyserIn: null, analyserOut: null, wet: null, dry: null, master: null, shaper: null,
  src: null, srcKind: null, micStream: null, bypass: false, inData: null, outData: null, buffers: {}, iirTimer: 0, lastCoefs: null };
function tanhCurve(points, ceiling) { // deterministic soft clip: y = ceiling · tanh(x / ceiling)
  const c = new Float32Array(points);
  for (let i = 0; i < points; i++) { const x = (i / (points - 1)) * 2 - 1; c[i] = ceiling * Math.tanh(x / ceiling); }
  return c;
}
function auMsg(text, cls) { const d = document.createElement('div'); d.className = 'msg ' + (cls || ''); d.textContent = text; $('au-msgs').appendChild(d); while ($('au-msgs').children.length > 4) $('au-msgs').firstChild.remove(); }
function auClearMsgs() { $('au-msgs').innerHTML = ''; }
function auState() { const b = $('au-state'); b.textContent = audio.ctx ? audio.ctx.state + ' @ ' + audio.ctx.sampleRate + ' Hz' + (audio.mode ? ' · ' + audio.mode : '') : 'no context'; b.className = 'badge ' + (audio.ctx && audio.ctx.state === 'running' ? 'ok' : audio.ctx ? 'warn' : ''); $('au-stop').disabled = !audio.ctx; $('au-start').textContent = audio.ctx && audio.ctx.state === 'suspended' ? 'Resume audio' : audio.ctx ? 'Running' : 'Start audio'; }
async function startAudio() {
  auClearMsgs();
  if (audio.ctx) { try { await audio.ctx.resume(); } catch (e) { auMsg('Resume failed: ' + e.message, 'err'); } auState(); return; }
  const AC = window.AudioContext || window.webkitAudioContext;
  if (!AC) { auMsg('Web Audio API is not available in this browser.', 'err'); return; }
  let ctx; try { ctx = new AC({ sampleRate: dsp.fs }); } catch (e) { ctx = new AC(); }
  audio.ctx = ctx; ctx.addEventListener('statechange', auState);
  if (ctx.sampleRate !== dsp.fs) auMsg('AudioContext runs at ' + ctx.sampleRate + ' Hz (selected fs was ' + dsp.fs + ' Hz). Coefficient math now uses the context rate.', 'warn');
  dspSetFs(ctx.sampleRate, true);
  audio.analyserIn = ctx.createAnalyser(); audio.analyserOut = ctx.createAnalyser();
  for (const a of [audio.analyserIn, audio.analyserOut]) { a.fftSize = 4096; a.smoothingTimeConstant = 0.8; a.minDecibels = -100; a.maxDecibels = 0; }
  audio.inData = new Float32Array(audio.analyserIn.frequencyBinCount); audio.outData = new Float32Array(audio.analyserOut.frequencyBinCount);
  audio.wet = ctx.createGain(); audio.dry = ctx.createGain(); audio.master = ctx.createGain(); audio.shaper = ctx.createWaveShaper();
  audio.wet.gain.value = audio.bypass ? 0 : 1; audio.dry.gain.value = audio.bypass ? 1 : 0; audio.master.gain.value = +$('au-master').value;
  audio.shaper.curve = tanhCurve(4096, 0.9); audio.shaper.oversample = '2x'; // output ceiling 0.9, no makeup gain
  try {
    if (!ctx.audioWorklet) throw new Error('AudioWorklet unavailable');
    const url = URL.createObjectURL(new Blob([WORKLET_SRC], { type: 'application/javascript' }));
    try { await ctx.audioWorklet.addModule(url); } finally { URL.revokeObjectURL(url); }
    audio.node = new AudioWorkletNode(ctx, 'biquad-df2t', { numberOfInputs: 1, numberOfOutputs: 1, outputChannelCount: [2] });
    audio.mode = 'AudioWorklet';
    audio.analyserIn.connect(audio.node); audio.node.connect(audio.wet);
  } catch (e) {
    audio.mode = 'IIRFilterNode fallback';
    auMsg('AudioWorklet unavailable (' + e.message + '). Using IIRFilterNode recreated on change with a 30 ms crossfade; noise sources loop a 10 s buffer with a crossfaded seam.', 'warn');
  }
  audio.analyserIn.connect(audio.dry); audio.dry.connect(audio.analyserOut); audio.wet.connect(audio.analyserOut);
  audio.analyserOut.connect(audio.master); audio.master.connect(audio.shaper); audio.shaper.connect(ctx.destination);
  dspApplyCookbookIfNotCustom();
  audioPushCoefs(dsp.lastStable || dsp.coefs, true);
  try { await ctx.resume(); } catch (e) { auMsg('Context could not resume: ' + e.message, 'err'); }
  await setSource($('au-source').value);
  auState();
}
function dspApplyCookbookIfNotCustom() { if (!dsp.custom) dspApplyCookbook(); else dspApplyRoots(); }
function audioPushCoefs(c, force) {
  if (!audio.ctx || !c) return;
  if (!force && audio.lastCoefs && coefArray(c).every((v, i) => v === audio.lastCoefs[i])) return;
  audio.lastCoefs = coefArray(c);
  if (audio.node) audio.node.port.postMessage({ type: 'coefs', c: audio.lastCoefs });
  else if (audio.mode) { clearTimeout(audio.iirTimer); audio.iirTimer = setTimeout(() => iirSwap(c), 40); }
}
function iirSwap(c) {
  const ctx = audio.ctx; if (!ctx) return;
  let node; try { node = ctx.createIIRFilter([c.b0, c.b1, c.b2], [1, c.a1, c.a2]); } catch (e) { auMsg('IIRFilterNode rejected coefficients: ' + e.message, 'err'); return; }
  const g = ctx.createGain(), now = ctx.currentTime; g.gain.setValueAtTime(0, now); g.gain.linearRampToValueAtTime(1, now + 0.03);
  audio.analyserIn.connect(node); node.connect(g); g.connect(audio.wet);
  const old = audio.iir;
  if (old.node) { old.gain.gain.setValueAtTime(old.gain.gain.value, now); old.gain.gain.linearRampToValueAtTime(0, now + 0.03); setTimeout(() => { try { old.node.disconnect(); old.gain.disconnect(); } catch (e) {} }, 80); }
  audio.iir = { node, gain: g };
}
function makeNoiseFallbackBuffer(kind, ctx) { // no-AudioWorklet fallback: 10 s loop with an equal-power crossfade at the seam
  const fs = ctx.sampleRate, n = fs * 10, F = Math.round(0.05 * fs), buf = ctx.createBuffer(1, n, fs), d = buf.getChannelData(0);
  const g = new Float32Array(n + F), pink = NOISE.makePinkFilter();
  for (let i = 0; i < n + F; i++) g[i] = kind === 'pink' ? pink(Math.random() * 2 - 1) : NOISE.whiteSample(Math.random);
  for (let i = 0; i < n; i++) {
    if (i < F) { const w = (i / F) * Math.PI / 2; d[i] = g[i] * Math.sin(w) + g[n + i] * Math.cos(w); } // tail fades out, head fades in
    else d[i] = g[i];
  }
  return buf;
}
function makeSweepBuffer(ctx, dur) {
  const fs = ctx.sampleRate, n = Math.round(fs * dur), buf = ctx.createBuffer(1, n, fs), d = buf.getChannelData(0);
  const f0 = 20, f1 = Math.min(20000, 0.45 * fs), k = Math.log(f1 / f0), fade = Math.round(0.01 * fs);
  for (let i = 0; i < n; i++) { const t = i / fs, ph = 2 * Math.PI * f0 * dur / k * (Math.exp(t * k / dur) - 1); const env = Math.min(1, i / fade, (n - 1 - i) / fade); d[i] = 0.5 * env * Math.sin(ph); }
  if (f1 < 20000) auMsg('Sweep capped at ' + f1.toFixed(0) + ' Hz for this sample rate.', 'warn');
  return buf;
}
function stopSource() {
  if (audio.src) { try { audio.src.stop(); } catch (e) {} try { audio.src.disconnect(); } catch (e) {} audio.src = null; }
  if (audio.noise) { try { audio.noise.disconnect(); } catch (e) {} audio.noise.port.postMessage({ type: 'mode', mode: 'off' }); }
  if (audio.micStream) { audio.micStream.getTracks().forEach(t => t.stop()); audio.micStream = null; }
  audio.srcKind = null;
}
async function setSource(kind) {
  const ctx = audio.ctx; if (!ctx) return;
  stopSource();
  const warn = document.querySelector('#au-msgs .msg.feedback'); if (warn) warn.remove();
  if (kind === 'mic') {
    if (!window.isSecureContext || !navigator.mediaDevices || !navigator.mediaDevices.getUserMedia) { auMsg('Microphone needs a secure context (https:// or localhost) and getUserMedia support.', 'err'); return; }
    try {
      audio.micStream = await navigator.mediaDevices.getUserMedia({ audio: { echoCancellation: false, noiseSuppression: false, autoGainControl: false }, video: false });
      if (!audio.ctx) { audio.micStream.getTracks().forEach(t => t.stop()); audio.micStream = null; return; }
      audio.src = ctx.createMediaStreamSource(audio.micStream);
      const w = document.createElement('div'); w.className = 'msg warn feedback'; w.textContent = '⚠ Microphone → speakers can feed back. Use headphones or keep the master gain low.'; $('au-msgs').appendChild(w);
    } catch (e) {
      const m = e.name === 'NotAllowedError' || e.name === 'SecurityError' ? 'Microphone permission denied.' : e.name === 'NotFoundError' || e.name === 'OverconstrainedError' ? 'No microphone device found.' : 'Microphone error: ' + e.message;
      auMsg(m, 'err'); return;
    }
  } else {
    const s = ctx.createBufferSource(); s.loop = true;
    if (kind === 'sweep') { const dur = +$('au-sweep').value; if (!audio.buffers.sweep || audio.buffers.sweepDur !== dur || audio.buffers.sweepFs !== ctx.sampleRate) { audio.buffers.sweep = makeSweepBuffer(ctx, dur); audio.buffers.sweepDur = dur; audio.buffers.sweepFs = ctx.sampleRate; } s.buffer = audio.buffers.sweep; }
    else if (audio.node) { // AudioWorklet available: never-repeating generator
      if (!audio.noise) audio.noise = new AudioWorkletNode(ctx, 'noise-gen', { numberOfInputs: 0, numberOfOutputs: 1, outputChannelCount: [1] });
      audio.noise.port.postMessage({ type: 'mode', mode: kind });
      audio.noise.connect(audio.analyserIn); audio.srcKind = kind; return;
    } else {
      auMsg('AudioWorklet unavailable: ' + kind + ' noise loops a 10 s buffer with a crossfaded seam.', 'warn');
      const key = kind + ctx.sampleRate; if (!audio.buffers[key]) audio.buffers[key] = makeNoiseFallbackBuffer(kind, ctx); s.buffer = audio.buffers[key];
    }
    s.start(); audio.src = s;
  }
  audio.srcKind = kind;
  audio.src.connect(audio.analyserIn);
}
function stopAudio(reason) {
  if (!audio.ctx) return;
  stopSource(); clearTimeout(audio.iirTimer);
  for (const n of [audio.node, audio.noise, audio.iir.node, audio.iir.gain, audio.analyserIn, audio.analyserOut, audio.wet, audio.dry, audio.master, audio.shaper]) { try { n && n.disconnect(); } catch (e) {} }
  const ctx = audio.ctx;
  audio.ctx = null; audio.node = null; audio.noise = null; audio.iir = { node: null, gain: null }; audio.mode = null; audio.lastCoefs = null; audio.buffers = {};
  try { ctx.close(); } catch (e) {}
  dspSetFs(dsp.fs, false);
  auState();
  if (reason) auMsg('Audio stopped (' + reason + ').');
}
$('au-start').addEventListener('click', () => { startAudio().catch(e => { auMsg('Audio start failed: ' + e.message, 'err'); auState(); }); });
$('au-stop').addEventListener('click', () => stopAudio('stop'));
$('au-source').addEventListener('change', e => { if (audio.ctx) setSource(e.target.value).catch(err => auMsg(err.message, 'err')); });
$('au-sweep').addEventListener('input', e => { $('au-sweep-o').value = (+e.target.value).toFixed(1) + ' s'; if (audio.ctx && audio.srcKind === 'sweep') setSource('sweep'); });
$('au-master').addEventListener('input', e => { $('au-master-o').value = (+e.target.value).toFixed(2); if (audio.master) audio.master.gain.setTargetAtTime(+e.target.value, audio.ctx.currentTime, 0.02); });
$('au-bypass').addEventListener('click', () => {
  audio.bypass = !audio.bypass; const b = $('au-bypass'); b.setAttribute('aria-pressed', audio.bypass); b.textContent = audio.bypass ? 'Bypass: ON (dry signal)' : 'Bypass: off (filter active)';
  if (audio.ctx) { const t = audio.ctx.currentTime; audio.wet.gain.setTargetAtTime(audio.bypass ? 0 : 1, t, 0.01); audio.dry.gain.setTargetAtTime(audio.bypass ? 1 : 0, t, 0.01); }
});
window.addEventListener('pagehide', () => stopAudio());
/* ---------- spectrum analyzer ---------- */
function dspDrawSpectrum() {
  const P = plotFrame($('dsp-spectrum'), 46, 12, 10, 20), { ctx } = P;
  const fs = dsp.fs, fmin = 20, fmax = fs / 2;
  const xf = drawLogGrid(P, fmin, fmax, true), yf = drawYGrid(P, -100, 0, 20, '');
  ctx.fillStyle = '#8b98ab'; ctx.textAlign = 'left';
  if (!audio.ctx || audio.ctx.state !== 'running') { ctx.fillText('Spectrum analyzer (fftSize 4096) — start audio to see input (dim) and filtered output (bright).', P.L + 4, P.T + 11); fcMarker(P, xf); return; }
  audio.analyserIn.getFloatFrequencyData(audio.inData); audio.analyserOut.getFloatFrequencyData(audio.outData);
  const binHz = audio.ctx.sampleRate / audio.analyserIn.fftSize, nb = audio.inData.length, cols = Math.max(2, Math.floor(P.pw));
  const inv = invLogX(fmin, fmax, 0, cols), sample = (arr, f) => { const b = f / binHz, i = Math.floor(b); if (i >= nb - 1) return arr[nb - 1]; const t = b - i; return arr[i] * (1 - t) + arr[i + 1] * t; };
  const inCol = new Float32Array(cols), outCol = new Float32Array(cols), fcol = new Float32Array(cols);
  for (let i = 0; i < cols; i++) { const f = inv(i + 0.5); fcol[i] = f; inCol[i] = clamp(sample(audio.inData, f), -100, 0); outCol[i] = clamp(sample(audio.outData, f), -100, 0); }
  const line = (arr, color, lw) => { ctx.strokeStyle = color; ctx.lineWidth = lw; ctx.beginPath(); for (let i = 0; i < cols; i++) { const x = P.L + i + 0.5, y = yf(arr[i]); if (i) ctx.lineTo(x, y); else ctx.moveTo(x, y); } ctx.stroke(); };
  line(inCol, 'rgba(139,152,171,0.55)', 1); line(outCol, '#4cc9f0', 1.6);
  if ($('au-showh').checked && dsp.lastStable) {
    const sorted = Array.from(inCol).filter(v => v > -99).sort((a, b) => a - b), off = sorted.length ? sorted[sorted.length >> 1] : -60;
    const hc = new Float32Array(cols); for (let i = 0; i < cols; i++) hc[i] = clamp(off + magDb(freqResponseAt(dsp.lastStable, 2 * Math.PI * fcol[i] / fs)), -100, 0);
    ctx.setLineDash([5, 4]); line(hc, '#f9a03f', 1.2); ctx.setLineDash([]);
    ctx.fillStyle = '#f9a03f'; ctx.fillText('|H| overlay (offset to input median ' + off.toFixed(0) + ' dB)', P.L + 4, P.T + 24);
  }
  ctx.fillStyle = '#8b98ab'; ctx.fillText('in (dim) / out (bright) · ' + audio.srcKind + (audio.bypass ? ' · BYPASS' : '') + ' · dBFS', P.L + 4, P.T + 11);
  fcMarker(P, xf);
}
function dspFrame() {
  dspRaf = 0;
  if (activeModule !== 'dsp' || document.hidden) return;
  dspDrawSpectrum();
  dspRaf = requestAnimationFrame(dspFrame);
}

/* ---------- init ---------- */
dspSetFs(48000, false); dspSetFc(1000); dspApplyCookbook();
emReadInputs(); emSetTheta(45); emSetPol('TE');
startLoops();
/* ---------- shared state (URL hash) + self-tests ---------- */
FSP.state.bind('em', {
  get: () => ({ th: em.theta, pol: em.pol, f: em.logf, e1: em.er1, m1: em.mr1, s1: em.s1, e2: em.er2, m2: em.mr2, s2: em.s2, span: em.span }),
  set: p => {
    const n = (k, id) => { const v = parseFloat(p[k]); if (Number.isFinite(v)) $(id).value = v; };
    n('e1', 'em-er1'); n('m1', 'em-mr1'); n('s1', 'em-s1'); n('e2', 'em-er2'); n('m2', 'em-mr2'); n('s2', 'em-s2'); n('f', 'em-freq'); n('span', 'em-span');
    emReadInputs(); const th = parseFloat(p.th); if (Number.isFinite(th)) emSetTheta(th); if (p.pol === 'TE' || p.pol === 'TM') emSetPol(p.pol);
  }
});
FSP.state.bind('dsp', {
  get: () => ({ type: dsp.type, fc: dsp.fc, fs: dsp.fs, q: dsp.Q, g: dsp.gainDb }),
  set: p => {
    if (FILTER_NAMES[p.type]) dsp.type = p.type;
    const fs = parseFloat(p.fs); if (Number.isFinite(fs) && fs >= 3000 && fs <= 192000 && !audio.ctx) dspSetFs(fs, false);
    const fc = parseFloat(p.fc); if (Number.isFinite(fc)) dspSetFc(fc);
    const q = parseFloat(p.q); if (Number.isFinite(q)) { dsp.Q = clamp(q, 0.1, 20); $('dsp-q').value = dsp.Q; $('dsp-q-n').value = dsp.Q; }
    const g = parseFloat(p.g); if (Number.isFinite(g)) { dsp.gainDb = clamp(g, -24, 24); $('dsp-gain').value = dsp.gainDb; $('dsp-gain-n').value = dsp.gainDb; }
    dspApplyCookbook();
  }
});
['input', 'change', 'click'].forEach(ev => document.querySelector('main').addEventListener(ev, () => FSP.state.touch()));
{
  const hdr = document.querySelector('header'), b = document.createElement('button');
  b.type = 'button'; b.textContent = 'Copy link'; b.className = 'copy-link'; b.title = 'Copy a link that restores the current tab and parameters';
  b.addEventListener('click', () => { const url = FSP.state.link(); const done = () => { b.textContent = 'Copied'; setTimeout(() => { b.textContent = 'Copy link'; }, 1200); }; if (navigator.clipboard && navigator.clipboard.writeText) navigator.clipboard.writeText(url).then(done, () => window.prompt('Copy this link', url)); else window.prompt('Copy this link', url); });
  hdr.appendChild(b);
}
{ const want = FSP.state.restore(); if (want && want !== 'em' && $('tab-' + want)) activateTab('tab-' + want, false); else FSP.state.setActive('em'); }
if (/selftest/.test(location.search)) { const res = FSP.runSelfTests(); const pre = $('selftest-out'); pre.hidden = false; pre.textContent = res.lines.join('\n'); }
window.runSelfTests = FSP.runSelfTests;
