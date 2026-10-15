'use strict';
const assert = require('node:assert/strict');
const { chromium } = require('playwright');
const { pathToFileURL } = require('node:url');
const path = require('node:path');
(async () => {
  const browser = await chromium.launch();
  try {
    const page = await browser.newPage({ viewport: { width: 390, height: 900 } });
    const errors = []; page.on('pageerror', error => errors.push(error.message));
    await page.goto(pathToFileURL(path.join(__dirname, '..', 'index.html')).href);
    async function open(id) {
      await page.locator('#tab-picker').selectOption('tab-' + id);
      const lesson = page.locator('#panel-' + id + ' .foundation');
      await lesson.locator('summary').click(); return lesson;
    }
    const field = await open('fields');
    const flux = () => field.locator('output').nth(1).textContent();
    const before = await flux();
    await field.getByLabel('Gaussian sphere radius (m)').fill('2');
    assert((await flux()).includes('2.247 V/m'));
    assert.equal(before.split('Closed-surface flux = ')[1].split('\n')[0], (await flux()).split('Closed-surface flux = ')[1].split('\n')[0]);
    await field.getByLabel('Field pattern').selectOption('rotation');
    assert((await field.locator('output').first().textContent()).includes('Divergence: 0 · Curl z: 2'));
    const fft = await open('dspdesign');
    await fft.getByLabel('FFT stage (0 = bit-reversed input)').fill('3');
    for (let i = 0; i < 8; i++) assert((await fft.locator('output').textContent()).includes(`${i}: 1.000 + j0.000`));
    await fft.getByLabel('Eight-sample input').selectOption('constant');
    assert((await fft.locator('output').textContent()).includes('0: 8.000 + j0.000'));
    const net = await open('net');
    await net.getByLabel('Journey step').fill('1');
    assert((await net.locator('output').textContent()).includes('Who has 192.168.1.1'));
    await net.getByLabel('Destination').selectOption('local');
    assert((await net.locator('output').textContent()).includes('Who has 192.168.1.20'));
    const motor = await open('rotating');
    assert((await motor.locator('output').textContent()).includes('1500 rpm'));
    await motor.getByLabel('Internal EMF lead angle δ (degrees)').fill('100');
    assert((await motor.locator('output').textContent()).includes('unstable branch'));
    for (const id of ['fields', 'dspdesign', 'net', 'rotating']) {
      await page.locator('#tab-picker').selectOption('tab-' + id);
      assert(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth + 2), id + ' overflow');
    }
    assert.deepEqual(errors, []);
    console.log('PASS: theory results, FFT reference cases, next-hop ARP, synchronous speed/stability, and mobile layout');
  } finally { await browser.close(); }
})().catch(error => { console.error(error); process.exit(1); });
