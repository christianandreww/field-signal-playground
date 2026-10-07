/* FFT / DFT / convolution. Pure; no DOM. Arrays are Float64Array (or Array) of equal length for re/im. */
'use strict';
FSP.fft = (function () {
  function isPow2(n) { return n > 0 && (n & (n - 1)) === 0; }
  function nextPow2(n) { let p = 1; while (p < n) p <<= 1; return p; }
  // In-place iterative radix-2. inverse=true applies conjugate twiddles and 1/N scaling.
  function fft(re, im, inverse) {
    const n = re.length; if (!isPow2(n) || im.length !== n) throw new Error('fft: length must be a power of two');
    for (let i = 1, j = 0; i < n; i++) { let bit = n >> 1; for (; j & bit; bit >>= 1) j ^= bit; j ^= bit; if (i < j) { let t = re[i]; re[i] = re[j]; re[j] = t; t = im[i]; im[i] = im[j]; im[j] = t; } }
    for (let len = 2; len <= n; len <<= 1) {
      const ang = (inverse ? 2 : -2) * Math.PI / len, half = len >> 1;
      for (let k = 0; k < half; k++) {
        const wr = Math.cos(ang * k), wi = Math.sin(ang * k);
        for (let i = k; i < n; i += len) { const j = i + half, xr = re[j] * wr - im[j] * wi, xi = re[j] * wi + im[j] * wr; re[j] = re[i] - xr; im[j] = im[i] - xi; re[i] += xr; im[i] += xi; }
      }
    }
    if (inverse) for (let i = 0; i < n; i++) { re[i] /= n; im[i] /= n; }
  }
  // Direct O(N²) DFT (tests only). Returns {re, im}.
  function dft(xr, xi) {
    const n = xr.length, re = new Float64Array(n), im = new Float64Array(n);
    for (let k = 0; k < n; k++) { let sr = 0, si = 0; for (let t = 0; t < n; t++) { const a = -2 * Math.PI * ((k * t) % n) / n, c = Math.cos(a), s = Math.sin(a), vi = xi ? xi[t] : 0; sr += xr[t] * c - vi * s; si += xr[t] * s + vi * c; } re[k] = sr; im[k] = si; }
    return { re, im };
  }
  function convolveDirect(a, b) { const out = new Float64Array(a.length + b.length - 1); for (let i = 0; i < a.length; i++) for (let j = 0; j < b.length; j++) out[i + j] += a[i] * b[j]; return out; }
  function convolveFFT(a, b) {
    const L = a.length + b.length - 1, N = nextPow2(L), ar = new Float64Array(N), ai = new Float64Array(N), br = new Float64Array(N), bi = new Float64Array(N);
    ar.set(a); br.set(b); fft(ar, ai, false); fft(br, bi, false);
    for (let k = 0; k < N; k++) { const r = ar[k] * br[k] - ai[k] * bi[k], i = ar[k] * bi[k] + ai[k] * br[k]; ar[k] = r; ai[k] = i; }
    fft(ar, ai, true); return ar.slice(0, L);
  }
  function convolveCircular(a, b) { const N = Math.max(a.length, b.length), out = new Float64Array(N); for (let n = 0; n < N; n++) { let s = 0; for (let k = 0; k < N; k++) s += (a[k] || 0) * (b[((n - k) % N + N) % N] || 0); out[n] = s; } return out; }
  return { isPow2, nextPow2, fft, dft, convolveDirect, convolveFFT, convolveCircular };
})();
