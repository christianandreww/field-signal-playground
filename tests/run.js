#!/usr/bin/env node
/* Node test harness: loads the same classic scripts the browser loads (no DOM), then runs FSP.runSelfTests().
   Usage: node tests/run.js */
'use strict';
const fs = require('fs'), path = require('path'), vm = require('vm');
const root = path.join(__dirname, '..');
const html = fs.readFileSync(path.join(root, 'index.html'), 'utf8');
// Script order comes from index.html; DOM-only scripts are marked data-dom and skipped here.
const files = [...html.matchAll(/<script\s+([^>]*)src="([^"]+)"([^>]*)><\/script>/g)].filter(m => !/data-dom/.test(m[1] + m[3])).map(m => m[2]);
const ctx = vm.createContext({ console, Math, Number, Float64Array, Float32Array, Uint8Array, Int32Array, Uint32Array, TextEncoder, setTimeout, clearTimeout });
ctx.globalThis = ctx;
for (const f of files) {
  const code = fs.readFileSync(path.join(root, f), 'utf8');
  try { vm.runInContext(code, ctx, { filename: f }); } catch (e) { console.error('FAILED to load ' + f + ': ' + (e && e.stack || e)); process.exit(2); }
}
const res = vm.runInContext('FSP.runSelfTests()', ctx);
console.log(res.lines.join('\n'));
process.exit(res.failed ? 1 : 0);
