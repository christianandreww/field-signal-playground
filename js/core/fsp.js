/* FSP core: namespace, tab registry, test registry, URL-hash state, SI parsing, small UI helpers.
   Loads in browsers and under Node (tests/run.js). No DOM access at load time. */
'use strict';
var FSP = (typeof globalThis !== 'undefined' ? globalThis : this).FSP = (typeof FSP !== 'undefined' && FSP) || {};
FSP.math = FSP.math || {};
FSP.tabs = FSP.tabs || [];       // registered tab descriptors, in order
FSP.testSets = FSP.testSets || []; // [{name, fn}]

/* ---------- registries ---------- */
// desc: {id, title, init(panelEl), activate(), deactivate()}; init is called once at mount time.
FSP.registerTab = function (desc) { if (!desc || !desc.id) throw new Error('registerTab: id required'); FSP.tabs.push(desc); return desc; };
// fn(t): t.check(name, ok, detail), t.near(a,b,tol) -> bool, t.rel(a,b,relTol) -> bool
FSP.registerTests = function (name, fn) { FSP.testSets.push({ name: name, fn: fn }); };

FSP.runSelfTests = function () {
  const lines = [], results = [];
  const record = (ok, name, detail) => {
    ok = !!ok; results.push({ name, ok });
    const l = (ok ? 'PASS' : 'FAIL') + ' — ' + name + (detail ? '  [' + detail + ']' : '');
    lines.push(l);
  };
  if (typeof runV1Tests === 'function') {
    const r = runV1Tests();
    r.results.forEach((x, i) => { results.push(x); lines.push(r.lines[i]); });
  }
  FSP.testSets.forEach(set => {
    const t = {
      check: (name, ok, detail) => record(ok, set.name + ': ' + name, detail),
      near: (a, b, tol) => Number.isFinite(a) && Number.isFinite(b) && Math.abs(a - b) <= tol,
      rel: (a, b, tol) => Number.isFinite(a) && Number.isFinite(b) && Math.abs(a - b) <= tol * Math.max(Math.abs(b), 1e-300),
    };
    try { set.fn(t); } catch (e) { record(false, set.name + ': threw ' + (e && e.message || e)); }
  });
  const passed = results.filter(r => r.ok).length;
  lines.push(passed + '/' + results.length + ' checks passed');
  return { passed, failed: results.length - passed, results, lines };
};

/* ---------- URL-hash state ---------- */
// Format: #t=<activeTab>&<tabId>.<key>=<value>&...   (values are URI-encoded strings)
FSP.state = (function () {
  const binds = {}; let timer = 0; let activeTab = null; let restoring = false;
  function encode(active, params) {
    const out = []; if (active) out.push('t=' + encodeURIComponent(active));
    Object.keys(params || {}).forEach(tab => {
      const p = params[tab] || {};
      Object.keys(p).forEach(k => { const v = p[k]; if (v === undefined || v === null || (typeof v === 'number' && !Number.isFinite(v))) return; out.push(encodeURIComponent(tab) + '.' + encodeURIComponent(k) + '=' + encodeURIComponent(String(v))); });
    });
    return out.join('&');
  }
  function decode(str) {
    const res = { active: null, params: {} };
    if (typeof str !== 'string') return res;
    str = str.replace(/^#/, '');
    str.split('&').forEach(pair => {
      if (!pair) return; const i = pair.indexOf('='); if (i < 0) return;
      let k, v; try { k = decodeURIComponent(pair.slice(0, i)); v = decodeURIComponent(pair.slice(i + 1)); } catch (e) { return; }
      if (k === 't') { res.active = v; return; }
      const d = k.indexOf('.'); if (d <= 0 || d === k.length - 1) return;
      const tab = k.slice(0, d), key = k.slice(d + 1);
      (res.params[tab] = res.params[tab] || {})[key] = v;
    });
    return res;
  }
  function collect() { const params = {}; Object.keys(binds).forEach(tab => { try { params[tab] = binds[tab].get(); } catch (e) { /* ignore */ } }); return params; }
  function write() {
    timer = 0; if (typeof location === 'undefined' || typeof history === 'undefined') return;
    try { history.replaceState(null, '', '#' + encode(activeTab, collect())); } catch (e) { /* file:// or sandboxed */ }
  }
  return {
    encode, decode,
    // bind(tabId, {get: () => flatObject, set: (flatObjectOfStrings) => void}); set must tolerate missing/invalid keys.
    bind(tabId, io) { binds[tabId] = io; },
    // modules call touch() after any user parameter change
    touch() { if (restoring || typeof window === 'undefined') return; clearTimeout(timer); timer = setTimeout(write, 200); },
    setActive(id) { activeTab = id; this.touch(); },
    restore(hash) {
      const d = decode(hash === undefined ? (typeof location !== 'undefined' ? location.hash : '') : hash);
      restoring = true;
      try { Object.keys(d.params).forEach(tab => { if (binds[tab]) { try { binds[tab].set(d.params[tab]); } catch (e) { /* ignore bad state */ } } }); } finally { restoring = false; }
      return d.active;
    },
    link() { return (typeof location !== 'undefined' ? location.href.split('#')[0] : '') + '#' + encode(activeTab, collect()); },
  };
})();

/* ---------- parsing / formatting helpers ---------- */
// "4.7k" -> 4700, "10u" -> 1e-5, "2.2Meg"/"2.2M" -> 2.2e6 (case: m = milli, M = mega); returns NaN on junk.
FSP.parseSI = function (s) {
  if (typeof s === 'number') return s;
  if (typeof s !== 'string') return NaN;
  const m = /^\s*([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)\s*(Meg|meg|MEG|[TGMkKmuµnpf])?\s*[A-Za-zΩ°%]*\s*$/.exec(s);
  if (!m) return NaN;
  const mult = { T: 1e12, G: 1e9, Meg: 1e6, meg: 1e6, MEG: 1e6, M: 1e6, k: 1e3, K: 1e3, m: 1e-3, u: 1e-6, 'µ': 1e-6, n: 1e-9, p: 1e-12, f: 1e-15 };
  return parseFloat(m[1]) * (m[2] ? mult[m[2]] : 1);
};
FSP.fmt = function (x, d) { return Number.isFinite(x) ? x.toFixed(d === undefined ? 4 : d) : '—'; };
FSP.fmtEng = function (x, unit, digits) {
  if (!Number.isFinite(x)) return '—'; unit = unit || ''; digits = digits === undefined ? 4 : digits;
  if (x === 0) return '0 ' + unit;
  const ax = Math.abs(x), pre = [[1e12, 'T'], [1e9, 'G'], [1e6, 'M'], [1e3, 'k'], [1, ''], [1e-3, 'm'], [1e-6, 'µ'], [1e-9, 'n'], [1e-12, 'p']];
  for (const [s, p] of pre) if (ax >= s * 0.99995) return (x / s).toPrecision(digits).replace(/\.?0+$/, m => m.includes('.') ? '' : m) + ' ' + p + unit;
  return x.toExponential(2) + ' ' + unit;
};

/* ---------- DOM helpers (browser only; call from init/activate, never at load) ---------- */
FSP.ui = {
  el(tag, props) {
    const e = document.createElement(tag); props = props || {};
    Object.keys(props).forEach(k => { if (k === 'class') e.className = props[k]; else if (k === 'text') e.textContent = props[k]; else if (k === 'html') e.innerHTML = props[k]; else if (k.slice(0, 2) === 'on') e.addEventListener(k.slice(2), props[k]); else e.setAttribute(k, props[k]); });
    for (let i = 2; i < arguments.length; i++) { const c = arguments[i]; if (c) e.appendChild(typeof c === 'string' ? document.createTextNode(c) : c); }
    return e;
  },
  // Slider + numeric box. opts: {label, min, max, step, value, log, unit, digits, onInput(v)} -> {el, get(), set(v, silent)}
  slider(parent, opts) {
    const U = FSP.ui, log = !!opts.log, N = 1000; let val = opts.value;
    const toPos = v => log ? Math.round(N * Math.log(v / opts.min) / Math.log(opts.max / opts.min)) : v;
    const fromPos = p => log ? opts.min * Math.pow(opts.max / opts.min, p / N) : p;
    const rng = U.el('input', { type: 'range', min: log ? 0 : opts.min, max: log ? N : opts.max, step: log ? 1 : (opts.step || 'any') });
    const box = U.el('input', { type: 'text', inputmode: 'decimal', 'aria-label': opts.label, size: 9 });
    const unit = U.el('span', { text: opts.unit || '' });
    const row = U.el('div', { class: 'ctl' }, U.el('label', { text: opts.label }), rng, U.el('span', null, box, unit));
    function show() { rng.value = toPos(val); box.value = FSP.fmtNum ? FSP.fmtNum(val, opts.digits) : String(+val.toPrecision(opts.digits || 6)); }
    function set(v, silent) { if (!Number.isFinite(v)) return; v = Math.min(opts.max, Math.max(opts.min, v)); if (opts.step && !log) v = Math.round(v / opts.step) * opts.step; val = v; show(); if (!silent && opts.onInput) opts.onInput(val); }
    rng.addEventListener('input', () => { let v = fromPos(parseFloat(rng.value)); if (opts.step && !log) v = Math.round(v / opts.step) * opts.step; val = v; show(); if (opts.onInput) opts.onInput(val); });
    box.addEventListener('change', () => { const v = FSP.parseSI(box.value); if (Number.isFinite(v)) set(v); else show(); });
    show(); if (parent) parent.appendChild(row);
    return { el: row, get: () => val, set };
  },
  // HiDPI canvas that re-renders on resize. -> {cv, prep() -> {ctx,w,h,dpr}, onResize(fn)}
  canvas(parent, opts) {
    opts = opts || {}; const cv = FSP.ui.el('canvas', { class: 'plot' }); cv.style.width = '100%'; cv.style.height = (opts.height || 260) + 'px'; cv.style.display = 'block';
    if (parent) parent.appendChild(cv); let cb = null;
    if (typeof ResizeObserver !== 'undefined') new ResizeObserver(() => { if (cb) cb(); }).observe(cv);
    return {
      cv, onResize(fn) { cb = fn; },
      prep() { const r = cv.getBoundingClientRect(), dpr = window.devicePixelRatio || 1, w = Math.max(1, Math.round(r.width)), h = Math.max(1, Math.round(r.height)); const pw = Math.round(w * dpr), ph = Math.round(h * dpr); if (cv.width !== pw || cv.height !== ph) { cv.width = pw; cv.height = ph; } const ctx = cv.getContext('2d'); ctx.setTransform(dpr, 0, 0, dpr, 0, 0); return { ctx, w, h, dpr }; },
    };
  },
  // Collapsible "Show working" panel. -> {el, set(linesArrayOrString)}
  working(parent) {
    const pre = FSP.ui.el('pre', { class: 'working' }); const d = FSP.ui.el('details', { class: 'working-box' }, FSP.ui.el('summary', { text: 'Show working' }), pre);
    if (parent) parent.appendChild(d);
    return { el: d, set(x) { pre.textContent = Array.isArray(x) ? x.join('\n') : String(x); } };
  },
  fieldset(parent, title) { const f = FSP.ui.el('fieldset', null, FSP.ui.el('legend', { text: title })); if (parent) parent.appendChild(f); return f; },
  readout(parent, label) { const v = FSP.ui.el('output', { text: '—' }); parent.appendChild(FSP.ui.el('div', { class: 'ctl' }, FSP.ui.el('span', { text: label }), v, FSP.ui.el('span'))); return v; },
  select(parent, label, options, value, onChange) { const s = FSP.ui.el('select', { 'aria-label': label }); options.forEach(o => { const v = Array.isArray(o) ? o[0] : o, t = Array.isArray(o) ? o[1] : o; s.appendChild(FSP.ui.el('option', { value: v, text: t })); }); s.value = value; s.addEventListener('change', () => onChange(s.value)); parent.appendChild(FSP.ui.el('div', { class: 'ctl' }, FSP.ui.el('label', { text: label }), s, FSP.ui.el('span'))); return s; },
  button(parent, text, onClick) { const b = FSP.ui.el('button', { type: 'button', text, onclick: onClick }); parent.appendChild(b); return b; },
};
FSP.fmtNum = function (v, d) { if (!Number.isFinite(v)) return '—'; const a = Math.abs(v); if (a !== 0 && (a >= 1e6 || a < 1e-3)) return v.toExponential(3); return String(+v.toPrecision(d || 6)); };

/* ---------- tab mounting (called once by js/v1-ui.js) ---------- */
FSP.mountTabs = function (tablistEl, mainEl) {
  const created = [];
  FSP.tabs.forEach(d => {
    const btn = FSP.ui.el('button', { role: 'tab', id: 'tab-' + d.id, 'aria-selected': 'false', 'aria-controls': 'panel-' + d.id, tabindex: '-1', text: d.title });
    const sec = FSP.ui.el('section', { id: 'panel-' + d.id, role: 'tabpanel', 'aria-labelledby': 'tab-' + d.id, hidden: '' });
    tablistEl.appendChild(btn); mainEl.insertBefore(sec, document.getElementById('selftest-out'));
    try { if (d.init) d.init(sec); } catch (e) { sec.appendChild(FSP.ui.el('div', { class: 'msg bad', text: 'Failed to initialise ' + d.title + ': ' + e.message })); if (typeof console !== 'undefined') console.error(e); }
    created.push(d);
  });
  return created;
};
if (typeof module !== 'undefined' && module.exports) module.exports = FSP;
