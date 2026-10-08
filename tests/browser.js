#!/usr/bin/env node
/* Headless smoke test (needs Playwright + Chromium): loads index.html over file://, visits every tab, fails on console/page errors.
   Usage: node tests/browser.js   (set PW_MODULE to the playwright module path if not globally resolvable) */
'use strict';
const path = require('path');
let pw; try { pw = require(process.env.PW_MODULE || 'playwright'); } catch (e) { console.log('SKIP: playwright not available'); process.exit(0); }
(async () => {
  const exe = process.env.CHROMIUM || '/opt/pw-browsers/chromium';
  const browser = await pw.chromium.launch({ executablePath: require('fs').existsSync(exe) ? exe : undefined, args: ['--no-sandbox', '--autoplay-policy=no-user-gesture-required'] });
  const page = await browser.newPage({ viewport: { width: 1280, height: 900 } });
  const errs = []; page.on('pageerror', e => errs.push('pageerror: ' + e.message)); page.on('console', m => { if (m.type() === 'error') errs.push('console: ' + m.text()); });
  await page.goto('file://' + path.join(__dirname, '..', 'index.html') + '?selftest');
  await page.waitForTimeout(500);
  const ids = await page.$$eval('[role=tab]', els => els.map(e => e.id));
  for (const id of ids) { await page.$eval('#' + id, e => e.click()); await page.waitForTimeout(400); const vis = await page.$eval('#' + id, e => document.getElementById(e.getAttribute('aria-controls')).hidden); if (vis) errs.push('panel hidden for ' + id); }
  const st = await page.$eval('#selftest-out', e => e.textContent);
  console.log('tabs: ' + ids.join(', ')); console.log(st.split('\n').slice(-1)[0]);
  if (/FAIL/.test(st)) errs.push('self-test failures in page');
  await browser.close();
  if (errs.length) { console.log(errs.join('\n')); process.exit(1); } console.log('browser smoke OK');
})();
