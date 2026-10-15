#!/usr/bin/env node
/* Real navigation, saved-link, and audio lifecycle regression checks. */
'use strict';
const path = require('node:path'), assert = require('node:assert/strict');
const { pathToFileURL } = require('node:url');
const http = require('node:http'), fs = require('node:fs/promises');
// Missing browser tooling must fail locally and in CI, never silently skip.
const pw = require(process.env.PW_MODULE || 'playwright');
const url = pathToFileURL(path.join(__dirname, '..', 'index.html')).href;
const root = path.resolve(__dirname, '..');
const server = http.createServer(async (req, res) => {
  try {
    const pathname = decodeURIComponent(new URL(req.url, 'http://localhost').pathname);
    const file = path.resolve(root, '.' + pathname);
    const relative = path.relative(root, file);
    if (relative.startsWith('..') || path.isAbsolute(relative)) { res.writeHead(403).end(); return; }
    const body = await fs.readFile(file);
    res.setHeader('Content-Type', ({ '.js': 'text/javascript', '.css': 'text/css', '.html': 'text/html' })[path.extname(file)] || 'application/octet-stream');
    res.end(body);
  } catch (e) { res.writeHead(404).end(); }
});
(async () => {
  const browser = await pw.chromium.launch({ ...(process.env.CHROMIUM ? { executablePath: process.env.CHROMIUM } : {}), args: ['--mute-audio'] });
  try {
  const page = await browser.newPage({ viewport: { width: 1280, height: 900 } });
  page.setDefaultTimeout(15000);
  const errs = []; page.on('pageerror', e => errs.push('pageerror: ' + e.message)); page.on('console', m => { if (m.type() === 'error') errs.push('console: ' + m.text()); });
  // Permission timing is mocked; the audio contexts and graphs are real Web Audio.
  await page.addInitScript(() => {
    window.__mic = { requests: [], attached: 0, stopped: 0 };
    const devices = navigator.mediaDevices || {};
    Object.defineProperty(navigator, 'mediaDevices', { value: devices, configurable: true });
    devices.getUserMedia = () => new Promise((resolve, reject) => {
      __mic.requests.push({ resolve: () => resolve({ getTracks: () => [{ stop: () => __mic.stopped++ }] }), reject });
    });
    AudioContext.prototype.createMediaStreamSource = function () { __mic.attached++; return this.createGain(); };
    window.__worklets = { hold: false, requests: [] };
    const addModule = AudioWorklet.prototype.addModule;
    AudioWorklet.prototype.addModule = function (...args) {
      const loaded = addModule.apply(this, args);
      return __worklets.hold ? loaded.then(() => new Promise(resolve => __worklets.requests.push(resolve))) : loaded;
    };
  });
  await page.goto(url + '?selftest');
  await page.waitForFunction(() => /\d+\/\d+ checks passed/.test(document.getElementById('selftest-out').textContent));
  assert(!/FAIL/.test(await page.locator('#selftest-out').textContent()));
  const groups = await page.locator('.tab-group').evaluateAll(els => els.map(el => [...el.querySelectorAll('[role=tab]')].map(b => b.id)));
  const ids = groups.flat(); assert.equal(ids.length, 13);
  async function select(id) {
    await page.locator('.group-btn').nth(groups.findIndex(g => g.includes(id))).click();
    await page.locator('#' + id).click();
    await page.locator('#panel-' + id.slice(4)).waitFor({ state: 'visible' });
    assert.equal(await page.locator('#' + id).getAttribute('aria-selected'), 'true');
  }
  for (const id of ids) {
    await select(id);
    const topics = page.locator('#panel-' + id.slice(4) + ' > [role=group] > button');
    for (let i = 0; i < await topics.count(); i++) {
      await topics.nth(i).click();
      assert.equal(await topics.nth(i).getAttribute('aria-pressed'), 'true');
    }
  }
  console.log('PASS: visible desktop navigation for 13 tools and their top-level subtopics');
  // DOM order is the navigation order, even though the subject groups move buttons.
  await select(ids[0]); await page.locator('#' + ids[0]).focus();
  for (let i = 1; i <= ids.length; i++) {
    await page.keyboard.press('ArrowRight');
    assert.equal(await page.evaluate(() => document.activeElement.id), ids[i % ids.length]);
  }
  await page.keyboard.press('End'); assert.equal(await page.evaluate(() => document.activeElement.id), ids.at(-1));
  await page.keyboard.press('Home'); assert.equal(await page.evaluate(() => document.activeElement.id), ids[0]);
  await page.keyboard.press('ArrowLeft'); assert.equal(await page.evaluate(() => document.activeElement.id), ids.at(-1));
  await page.setViewportSize({ width: 390, height: 844 });
  for (const id of ids) { await page.locator('#tab-picker').selectOption(id); await page.locator('#panel-' + id.slice(4)).waitFor({ state: 'visible' }); }
  await page.setViewportSize({ width: 1280, height: 900 });
  console.log('PASS: cross-subject keyboard focus, wraparound, Home/End, and mobile picker');

  await select('tab-dsp');
  const pole = await page.evaluate(() => { const g = zGeom(), p = dsp.roots.poles[0]; return { x: g.cx + p.re * g.R, y: g.cy - p.im * g.R }; });
  await page.locator('#dsp-zplane').click({ position: pole }); await page.keyboard.press('Shift+ArrowLeft');
  await page.locator('#dsp-custom').waitFor({ state: 'visible' });
  await page.locator('#dsp-normgain').uncheck(); await page.locator('#dsp-autorange').check();
  const before = await page.evaluate(() => ({ coefs: coefArray(dsp.coefs), response: [0, 0.1, 1, 2, Math.PI].map(w => freqResponseAt(dsp.coefs, w).abs()), link: FSP.state.link() }));
  await page.goto(before.link);
  await page.locator('#dsp-custom').waitFor({ state: 'visible' });
  const after = await page.evaluate(() => ({ coefs: coefArray(dsp.coefs), response: [0, 0.1, 1, 2, Math.PI].map(w => freqResponseAt(dsp.coefs, w).abs()), custom: dsp.custom, normalize: dsp.normalize, auto: dsp.autoRange }));
  assert.deepEqual(after.coefs, before.coefs); assert.deepEqual(after.response, before.response);
  assert(after.custom && !after.normalize && after.auto);
  for (const cs of ['not-json', '[1,0,0,-3,2]', '[1,0,0]', '[null,0,0,0,0]', '[1,100,0,0,0]']) {
    await page.goto(url + '?case=' + encodeURIComponent(cs) + '#t=dsp&dsp.custom=1&dsp.cs=' + encodeURIComponent(cs));
    assert.equal(await page.evaluate(() => dsp.custom), false); assert(await page.evaluate(() => isStable(dsp.coefs)));
  }
  await page.goto(url + '?case=hostile#t=dsp&__proto__.fspPolluted=yes&constructor.x=no');
  assert.equal(await page.evaluate(() => ({}).fspPolluted), undefined);
  console.log('PASS: custom response round trip, invalid filter fallback, hostile links');

  await select('tab-audio');
  const panel = page.locator('#panel-audio'), input = panel.getByRole('combobox', { name: 'Input', exact: true });
  await input.selectOption('mic');
  async function requestMic() {
    const n = await page.evaluate(() => __mic.requests.length);
    await panel.getByRole('button', { name: 'Start', exact: true }).click();
    await page.waitForFunction(i => __mic.requests.length === i + 1, n); return n;
  }
  async function resolveMic(i) { await page.evaluate(async i => { __mic.requests[i].resolve(); await new Promise(r => setTimeout(r, 0)); }, i); }
  let request = await requestMic();
  await panel.getByRole('button', { name: 'Cancel', exact: true }).click(); await resolveMic(request);
  assert.deepEqual(await page.evaluate(() => [__mic.attached, __mic.stopped]), [0, 1]);
  request = await requestMic(); await select('tab-em'); await resolveMic(request);
  assert.deepEqual(await page.evaluate(() => [__mic.attached, __mic.stopped]), [0, 2]);
  await select('tab-audio'); request = await requestMic(); await input.selectOption('tone'); await resolveMic(request);
  assert.deepEqual(await page.evaluate(() => [__mic.attached, __mic.stopped]), [0, 3]);
  await input.selectOption('mic'); const old = await requestMic();
  await panel.getByRole('button', { name: 'Cancel', exact: true }).click(); const newer = await requestMic();
  await resolveMic(newer); await resolveMic(old);
  assert.deepEqual(await page.evaluate(() => [__mic.attached, __mic.stopped]), [1, 4]);
  await panel.getByRole('button', { name: 'Stop', exact: true }).click();
  assert.equal(await page.evaluate(() => __mic.stopped), 5);
  console.log('PASS: mic cancellation, tab/source changes, and out-of-order permissions');

  await select('tab-dsp'); await page.locator('#au-source').selectOption('mic'); await page.locator('#au-start').click();
  await page.waitForFunction(() => __mic.requests.length === 6);
  await page.locator('#au-source').selectOption('white'); await resolveMic(5);
  assert.deepEqual(await page.evaluate(() => [__mic.attached, __mic.stopped]), [1, 6]);
  await page.locator('#au-stop').click();
  console.log('PASS: DSP rejects stale microphone sources');

  // file:// deliberately exercises the IIR fallback; real worklets also need
  // a localhost run because Chromium may reject opaque-origin blob modules.
  await new Promise(resolve => server.listen(0, '127.0.0.1', resolve));
  await page.goto('http://127.0.0.1:' + server.address().port + '/index.html');
  for (const [tab, panelId, start, stop] of [
    ['tab-dsp', '#panel-dsp', '#au-start', '#au-stop'],
    ['tab-filters', '#panel-filters', 'button:text-is("Start audio")', 'button:text-is("Stop")'],
  ]) {
    await select(tab); await page.evaluate(() => { __worklets.hold = true; __worklets.requests = []; });
    const p = page.locator(panelId);
    console.log('Checking delayed worklet:', tab);
    await p.locator(start).click();
    await page.waitForFunction(() => __worklets.requests.length === 1).catch(async error => { console.error(await p.locator('.msg, .badge').allTextContents()); throw error; });
    await p.locator(stop).click(); await p.locator(start).click();
    await page.waitForFunction(() => __worklets.requests.length === 2);
    await page.evaluate(async () => { __worklets.requests[0](); await new Promise(r => setTimeout(r, 0)); });
    await page.evaluate(async () => { __worklets.requests[1](); await new Promise(r => setTimeout(r, 50)); });
    const states = await p.locator('.badge').allTextContents();
    assert(states.some(s => /running.*AudioWorklet/.test(s)), states.join('; '));
    await p.locator(stop).click(); await page.evaluate(() => { __worklets.hold = false; });
  }
  console.log('PASS: worklet loading survives stop/restart in both filter tools');
  assert.deepEqual(errs, [], 'browser console/page errors');
  console.log('Browser regression suite passed.');
  } finally { await browser.close(); if (server.listening) await new Promise(resolve => server.close(resolve)); }
})().catch(error => { console.error(error); process.exitCode = 1; });
