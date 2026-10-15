#!/usr/bin/env node
/* Shared UI behavior: themes, real results, saved designs, and responsive reachability. */
'use strict';
const assert = require('node:assert/strict');
const path = require('node:path');
const { pathToFileURL } = require('node:url');
const { chromium } = require(process.env.PW_MODULE || 'playwright');
const url = pathToFileURL(path.join(__dirname, '..', 'index.html')).href;

(async () => {
  const browser = await chromium.launch(process.env.CHROMIUM ? { executablePath: process.env.CHROMIUM } : {});
  try {
    const page = await browser.newPage({ viewport: { width: 1440, height: 1000 } });
    const errors = []; page.on('pageerror', e => errors.push(e.message));
    await page.goto(url);
    const groups = await page.locator('.tab-group').evaluateAll(els => els.map(el => [...el.querySelectorAll('[role=tab]')].map(tab => tab.id)));
    async function select(id) {
      if (page.viewportSize().width <= 900) await page.locator('#tab-picker').selectOption(id);
      else { await page.locator('.group-btn').nth(groups.findIndex(g => g.includes(id))).click(); await page.locator('#' + id).click(); }
      await page.locator('#panel-' + id.slice(4)).waitFor({ state: 'visible' });
    }
    // All tools remain reachable at desktop, tablet, and narrow phone widths.
    const courses = ['EE3101', 'IE3014', 'IE3017', 'EE2005'];
    for (const width of [1440, 1024, 900, 768, 390, 320]) {
      await page.setViewportSize({ width, height: 1000 });
      for (let group = 0; group < groups.length; group++) for (const id of groups[group]) {
        await select(id);
        assert.equal(await page.evaluate(() => document.documentElement.dataset.subject), courses[group]);
        assert((await page.locator('#workspace-course').textContent()).includes(courses[group]));
        const size = await page.evaluate(() => [document.documentElement.scrollWidth, innerWidth]);
        assert(size[0] <= size[1] + 2, 'horizontal overflow: ' + id + ' at ' + width);
      }
    }
    console.log('PASS: all subject identities and 13 tools fit desktop, tablet, and phone widths');
    await page.setViewportSize({ width: 1440, height: 1000 }); await select('tab-em');
    const materials = page.locator('#panel-em .control-details').filter({ has: page.locator('#em-er2') });
    assert.equal(await materials.getAttribute('open'), null);
    assert.equal(await page.locator('#panel-em .result-details').getAttribute('open'), null);
    await page.locator('#em-theta-n').fill('45'); await page.locator('#em-theta-n').dispatchEvent('change');
    await page.waitForFunction(() => document.getElementById('em-refracted').textContent === '28.1°');
    assert.equal(await page.locator('#em-reflected').textContent(), '9.2%');
    assert.equal(await page.locator('#em-transmitted').textContent(), '90.8%');
    await page.locator('[data-preset=glassair]').click();
    await page.waitForFunction(() => document.getElementById('em-refracted').textContent === 'Evanescent');
    assert.equal(await page.locator('#em-transmitted').textContent(), '0.0%');
    await materials.locator('summary').click();
    assert(await page.locator('#em-er2').isVisible());
    await page.locator('#panel-em .result-details > summary').click();
    await page.waitForFunction(() => document.getElementById('em-rt').width > 100);
    assert(await page.locator('#em-hud').isVisible());
    console.log('PASS: simplified wave results match the real solver; advanced controls and sweep remain usable');

    await select('tab-dsp');
    const before = await page.evaluate(() => JSON.stringify(coefArray(dsp.coefs)));
    await page.locator('#appearance-picker').selectOption('dark');
    assert.equal(await page.evaluate(() => document.documentElement.dataset.appearance), 'dark');
    assert.equal(await page.evaluate(() => JSON.stringify(coefArray(dsp.coefs))), before);
    assert.equal(await page.evaluate(() => audio.ctx), null);
    const darkPixel = await page.locator('#dsp-mag').evaluate(cv => [...cv.getContext('2d').getImageData(2, 2, 1, 1).data].slice(0, 3));
    await page.locator('#appearance-picker').selectOption('light');
    const lightPixel = await page.locator('#dsp-mag').evaluate(cv => [...cv.getContext('2d').getImageData(2, 2, 1, 1).data].slice(0, 3));
    assert.notDeepEqual(darkPixel, lightPixel, 'canvas must redraw when appearance changes');
    await page.locator('#appearance-picker').selectOption('dark');
    await page.waitForFunction(() => FSP.state.decode(location.hash).active === 'dsp');
    await page.reload();
    assert.equal(await page.locator('#appearance-picker').inputValue(), 'dark');
    assert.equal(await page.evaluate(() => document.documentElement.dataset.subject), 'IE3014');
    await page.emulateMedia({ colorScheme: 'light' });
    await page.locator('#appearance-picker').selectOption('system');
    assert.equal(await page.evaluate(() => document.documentElement.dataset.appearance), 'light');
    await page.emulateMedia({ colorScheme: 'dark' });
    await page.waitForFunction(() => document.documentElement.dataset.appearance === 'dark');
    console.log('PASS: appearance redraws plots, persists independently, follows the system, and never starts audio');

    await page.setViewportSize({ width: 390, height: 844 }); await select('tab-em');
    const order = await page.evaluate(() => ({ controls: document.querySelector('#panel-em .controls').getBoundingClientRect().top, theory: document.querySelector('#panel-em .em-deep-dive').getBoundingClientRect().top }));
    assert(order.controls < order.theory, 'mobile controls must precede detailed theory');
    const isolated = await browser.newContext({ reducedMotion: 'reduce' });
    const reduced = await isolated.newPage();
    await reduced.addInitScript(() => {
      Object.defineProperty(window, 'localStorage', { get() { throw new DOMException('Storage disabled', 'SecurityError'); } });
    });
    await reduced.goto(url);
    assert.equal(await reduced.locator('#em-play').textContent(), 'Play');
    await reduced.locator('#appearance-picker').selectOption('dark');
    assert.equal(await reduced.evaluate(() => document.documentElement.dataset.appearance), 'dark');
    await isolated.close();
    assert.deepEqual(errors, []);
    console.log('PASS: mobile control order, reduced motion, and disabled-storage fallback');
  } finally { await browser.close(); }
})().catch(error => { console.error(error); process.exitCode = 1; });
