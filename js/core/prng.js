/* Seeded PRNG + Gaussian. Pure; no DOM. */
'use strict';
FSP.prng = (function () {
  function mulberry32(seed) { let a = seed >>> 0; return function () { a = (a + 0x6D2B79F5) >>> 0; let t = a; t = Math.imul(t ^ (t >>> 15), t | 1); t ^= t + Math.imul(t ^ (t >>> 7), t | 61); return ((t ^ (t >>> 14)) >>> 0) / 4294967296; }; }
  // Box–Muller with cached spare; returns function() -> N(0,1)
  function gaussian(rng) { let spare = null; return function () { if (spare !== null) { const s = spare; spare = null; return s; } let u = 0; while (u <= 1e-300) u = rng(); const v = rng(), r = Math.sqrt(-2 * Math.log(u)), th = 2 * Math.PI * v; spare = r * Math.sin(th); return r * Math.cos(th); }; }
  return { mulberry32, gaussian };
})();
