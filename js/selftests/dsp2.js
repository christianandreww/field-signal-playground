/* Numerical regression checks for dspdesign; also loaded by ?selftest. */
(function () {
  'use strict';
  const { C, PI, TAU, cabs, cmul, csub } = FSP.math.dspdesign._shared;
  FSP.registerTests('dspdesign', t => {
    const D = FSP.math.dspdesign;
    const hasRoot = (rs, re, im, tol) => rs.some(o => Math.abs(o.r.re - re) < tol && Math.abs(o.r.im - im) < tol);
    // roots: z^2 + 1 -> ±j
    { const r = D.rootsMult([1, 0, 1]); t.check('roots of z²+1 are ±j', r.length === 2 && hasRoot(r, 0, 1, 1e-14) && hasRoot(r, 0, -1, 1e-14)); }
    // (z−1)²(z+2) = z³ − 3z + 2: double root 1, simple −2
    { const r = D.rootsMult([1, 0, -3, 2]), d = r.find(o => o.m === 2);
      t.check('z³−3z+2: double root at 1 (multiplicity detected), simple root −2', r.length === 2 && d && Math.abs(d.r.re - 1) < 1e-12 && Math.abs(d.r.im) === 0 && hasRoot(r, -2, 0, 1e-13), JSON.stringify(r.map(o => [o.r.re, o.r.im, o.m]))); }
    // z^10 − 1: ten roots of unity
    { const r = D.polyRoots([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1]); let ok = r.length === 10;
      for (let k = 0; k < 10; k++) { const th = TAU * k / 10; if (!r.some(z => Math.hypot(z.re - Math.cos(th), z.im - Math.sin(th)) < 1e-13)) ok = false; }
      t.check('z¹⁰ − 1: all ten 10th roots of unity to 1e-13', ok); }
    // (z − 0.5)^2 (z² + 0.81): coefficients by hand: (z² − z + 0.25)(z² + 0.81) = z⁴ − z³ + 1.06 z² − 0.81 z + 0.2025
    { const r = D.rootsMult([1, -1, 1.06, -0.81, 0.2025]);
      t.check('(z−0.5)²(z²+0.81): double 0.5 and ±0.9j', r.length === 3 && r.some(o => o.m === 2 && Math.abs(o.r.re - 0.5) < 1e-12) && hasRoot(r, 0, 0.9, 1e-12) && hasRoot(r, 0, -0.9, 1e-12)); }
    // complex-coefficient: (z − (1+2j))(z − (−3+0.5j)) = z² − (−2+2.5j) z + (1+2j)(−3+0.5j) = z² + (2 − 2.5j) z + (−3 − 1 + (0.5 − 6)j) = z² + (2−2.5j)z + (−4 − 5.5j)
    { const r = D.rootsMult([C(1, 0), C(2, -2.5), C(-4, -5.5)]); t.check('complex-coefficient quadratic: roots 1+2j, −3+0.5j', hasRoot(r, 1, 2, 1e-12) && hasRoot(r, -3, 0.5, 1e-12)); }
    // order-10 Wilkinson-type with known roots 0.1..1.0 rebuilt from roots
    { const want = []; for (let k = 1; k <= 8; k++) want.push(C(0.95 * Math.cos(k * 0.7), 0.95 * Math.sin(k * 0.7))); want.push(C(-0.3, 0), C(0.6, 0));
      const p = D.polyFromRoots(want), r = D.polyRoots(p); const ok = want.every(w => r.some(z => cabs(csub(z, w)) < 1e-10));
      t.check('order-10 complex polynomial rebuilt from known roots: all recovered to 1e-10', ok && r.length === 10); }
    // triple root (z+0.5)^3 = z³ + 1.5 z² + 0.75 z + 0.125
    { const r = D.rootsMult([1, 1.5, 0.75, 0.125]); t.check('(z+0.5)³: one root −0.5 of multiplicity 3', r.length === 1 && r[0].m === 3 && Math.abs(r[0].r.re + 0.5) < 1e-10); }
    // impulse response of y[n] = 0.5 y[n−1] + x[n]  (a1 = −0.5)
    { const h = D.impulse([1], [1, -0.5], 30); let e = 0; for (let n = 0; n < 30; n++) e = Math.max(e, Math.abs(h[n] - Math.pow(0.5, n)));
      t.check('h[n] of y[n] = 0.5y[n−1] + x[n] is 0.5^n', e < 1e-15, e.toExponential(2)); }
    // step response of same: Σ 0.5^k = 2(1 − 0.5^{n+1})
    { const s = D.stepResp([1], [1, -0.5], 20); t.check('step response = 2(1 − 0.5^(n+1))', s.every((v, n) => Math.abs(v - 2 * (1 - Math.pow(0.5, n + 1))) < 1e-14)); }
    // partial fractions: 1/((1−0.5z^-1)(1−0.25z^-1)): a = [1, −0.75, 0.125]; A1 = 1/(1 − 0.25/0.5) = 2, A2 = 1/(1 − 0.5/0.25) = −1
    { const pf = D.partialFractions([1], [1, -0.75, 0.125]);
      const c5 = pf.terms.find(x => Math.abs(x.d.re - 0.5) < 1e-12), c25 = pf.terms.find(x => Math.abs(x.d.re - 0.25) < 1e-12);
      t.check('PF of 1/((1−0.5z⁻¹)(1−0.25z⁻¹)): residues 2 and −1', c5 && c25 && Math.abs(c5.C.re - 2) < 1e-12 && Math.abs(c25.C.re + 1) < 1e-12);
      let e = 0; for (let n = 0; n < 25; n++) e = Math.max(e, Math.abs(D.inverseZAt(pf, 0.5, n) - (2 * Math.pow(0.5, n) - Math.pow(0.25, n))));
      t.check('causal inverse = 2(0.5)ⁿ − (0.25)ⁿ', e < 1e-14, e.toExponential(2));
      // anticausal (|z| < 0.25): h[n] = −2(0.5)^n + (0.25)^n for n ≤ −1 ; at n = −1: −4 + 4 = 0; n = −2: −8 + 16 = 8
      t.check('anticausal ROC: h[−1] = 0, h[−2] = 8, h[0] = 0', Math.abs(D.inverseZAt(pf, 0, -1)) < 1e-12 && Math.abs(D.inverseZAt(pf, 0, -2) - 8) < 1e-12 && D.inverseZAt(pf, 0, 0) === 0);
      // two-sided 0.25<|z|<0.5: −2(0.5)^n u[−n−1] − (0.25)^n u[n]; check it satisfies the difference equation everywhere
      const res = D.diffEqResidual([1], [1, -0.75, 0.125], n => D.inverseZAt(pf, 0.25, n), -15, 15);
      t.check('two-sided ROC closed form satisfies the difference equation (residual < 1e-9)', res < 1e-9, res.toExponential(2)); }
    // double pole: 1/(1 − 0.5 z^-1)² = Σ (n+1) 0.5^n
    { const pf = D.partialFractions([1], [1, -1, 0.25]); let e = 0; for (let n = 0; n < 30; n++) e = Math.max(e, Math.abs(D.inverseZAt(pf, 0.5, n) - (n + 1) * Math.pow(0.5, n)));
      t.check('double pole 1/(1−0.5z⁻¹)²: h[n] = (n+1)(0.5)ⁿ', e < 1e-12, e.toExponential(2)); }
    // improper: (1 + 2z^-1 + z^-2)/(1 − 0.5z^-1): long division in w: Q = −2w − 8? check via recursion instead (independent algorithm)
    { const b = [1, 2, 1, 0.5], a = [1, -0.5, 0.06], pf = D.partialFractions(b, a), h = D.impulse(b, a, 40); let e = 0; for (let n = 0; n < 40; n++) e = Math.max(e, Math.abs(D.inverseZAt(pf, 0.5, n) - h[n]));
      t.check('improper H (M ≥ N): direct terms + PF matches the recursion', e < 1e-12 && pf.direct.length === 2, e.toExponential(2)); }
    // Kaiser β: A = 60 -> 0.1102·51.3 = 5.65326; A = 40 -> 0.5842·19^0.4 + 0.07886·19 = 1.89692 + 1.49834 = 3.39526; A = 20 -> 0
    t.check('Kaiser β(A=60) = 5.6533', t.near(D.kaiserBeta(60), 5.65326, 1e-5));
    t.check('Kaiser β(A=40) = 3.3953', t.near(D.kaiserBeta(40), 3.39526, 1e-4));
    t.check('Kaiser β(A=20) = 0', D.kaiserBeta(20) === 0);
    // O&S Example: ωp = 0.4π, ωs = 0.6π, δ = 0.001 -> A = 60, M = 52/(2.285·0.2π) = 36.22 -> 37
    t.check('Kaiser order (O&S example): M = 37', D.estimateOrder('kaiser', 0.2 * PI, 60) === 37);
    // Hamming-windowed lowpass: symmetric, DC gain ≈ 1
    { const r = D.designWindow({ resp: 'lp', edges: [0.3 * PI, 0.5 * PI], win: 'hamming', d1: 0.01, As: 50, orderMode: 'auto' });
      const dc = r.h.reduce((s, v) => s + v, 0); let sym = true; for (let n = 0; n <= r.M; n++) if (Math.abs(r.h[n] - r.h[r.M - n]) > 1e-15) sym = false;
      // Hamming main lobe 8π/M ≤ 0.2π -> M ≥ 40
      t.check('Hamming LP: M = 8π/Δω = 40, symmetric, DC gain ≈ 1 (±0.005)', r.M === 40 && sym && Math.abs(dc - 1) < 0.005, 'M=' + r.M + ' dc=' + dc.toFixed(5));
      t.check('Hamming LP meets ≈ −53 dB class stopband (A_s > 50 dB)', r.meas.As > 50, r.meas.As.toFixed(1)); }
    // highpass auto order is forced even (Type I)
    { const r = D.designWindow({ resp: 'hp', edges: [0.4 * PI, 0.6 * PI], win: 'kaiser', d1: 0.001, As: 60, orderMode: 'auto' });
      t.check('Kaiser HP: M bumped 37 -> 38 (Type I)', r.M === 38 && r.type.type === 1); }
    // FIR type classification
    t.check('type I: [1,2,1]', D.firType([1, 2, 1]).type === 1);
    t.check('type II: [1,1] (zero at ω=π)', D.firType([1, 1]).type === 2 && Math.abs(cabs(D.freqAt([1, 1], [1], PI))) < 1e-15);
    t.check('type III: [1,0,−1] (zeros at 0 and π)', D.firType([1, 0, -1]).type === 3);
    t.check('type IV: [1,−1] (zero at ω=0)', D.firType([1, -1]).type === 4);
    t.check('not linear phase: [1,2,3]', D.firType([1, 2, 3]).type === 0);
    // frequency sampling hits its samples exactly: |H(e^{jω_k})| = A_k
    { const r = D.designFreqSamp({ resp: 'lp', edges: [0.3 * PI, 0.4 * PI], M: 32, T: 0.4 }); let e = 0;
      r.Ak.forEach(s => { e = Math.max(e, Math.abs(cabs(D.freqAt(r.h, [1], s.w)) - s.A)); });
      t.check('frequency sampling: |H(e^{jω_k})| = A_k at all samples', e < 1e-12 && r.type.type === 1, e.toExponential(2)); }
    // structures vs recursion: 4th order with poles 0.8e^{±j0.5}, 0.6e^{±j1.8}, b = [1, 2, 3, 2, 1] and a 5th-order with real poles
    { const mkA = ps => D.polyFromRoots(ps).map(c => c.re);
      const cases = [
        [[1, 2, 3, 2, 1], mkA([C(0.8 * Math.cos(0.5), 0.8 * Math.sin(0.5)), C(0.8 * Math.cos(0.5), -0.8 * Math.sin(0.5)), C(0.6 * Math.cos(1.8), 0.6 * Math.sin(1.8)), C(0.6 * Math.cos(1.8), -0.6 * Math.sin(1.8))])],
        [[0.2, -0.1, 0.4, 0.3, 0, 0.1, 0.05], mkA([C(0.9, 0), C(0.8 * Math.cos(0.7), 0.8 * Math.sin(0.7)), C(0.8 * Math.cos(0.7), -0.8 * Math.sin(0.7)), C(-0.5, 0), C(0.3, 0)])],
        [[0, 0, 1, 0.5], [1, -0.9]]];
      let worst = 0, okAll = true;
      cases.forEach(([b, a]) => {
        const x = new Float64Array(80); x[0] = 1; const ref = D.filter(b, a, x), q = D.normBA(b, a);
        const ys = [D.simDF1(q.b, q.a, x), D.simDF2(q.b, q.a, x), D.simTDF2(q.b, q.a, x), D.simCascade(D.cascadeSOS(b, a), x), D.simParallel(D.parallelForm(b, a), x)];
        ys.forEach(y => { const e = D.maxAbsDiff(y, ref); worst = Math.max(worst, e); if (!(e < 1e-12)) okAll = false; });
      });
      t.check('DF-I, DF-II, TDF-II, cascade SOS, parallel all match the recursion to 1e-12', okAll, worst.toExponential(2)); }
    { const h = [0.1, 0.2, 0.4, 0.2, 0.1], x = new Float64Array(20); x[0] = 1; x[3] = -2; x[7] = 0.5;
      t.check('folded linear-phase FIR equals direct convolution', D.maxAbsDiff(D.simFolded(h, x), D.filter(h, [1], x)) < 1e-15);
      t.check('folded FIR (M=4, type I): 3 multipliers, 4 delays', D.foldedCounts(h).mult === 3 && D.foldedCounts(h).delay === 4); }
    { const cas = D.cascadeSOS([1, 2, 3, 2, 1], [1, -1.2, 0.9, -0.4, 0.1]);
      t.check('cascade: 4th order -> 2 SOS, every a-section has a₀ = 1', cas.secs.length === 2 && cas.secs.every(s => s.a[0] === 1)); }
    // DFT of [1,2,3,4] = [10, −2+2j, −2, −2−2j]
    { const X = D.dft([1, 2, 3, 4], 4), want = [[10, 0], [-2, 2], [-2, 0], [-2, -2]];
      t.check('DFT{1,2,3,4} = [10, −2+2j, −2, −2−2j]', want.every((w, k) => Math.abs(X.re[k] - w[0]) < 1e-12 && Math.abs(X.im[k] - w[1]) < 1e-12)); }
    // circular conv of [1,2,3,4] with [1,1]: y[n] = x[n] + x[(n−1) mod 4] = [5,3,5,7]; N = 5 equals linear [1,3,5,7,4]
    { const y4 = D.circConv([1, 2, 3, 4], [1, 1], 4), y5 = D.circConv([1, 2, 3, 4], [1, 1], 5);
      t.check('4-point circular conv [1,2,3,4]⊛[1,1] = [5,3,5,7]', [5, 3, 5, 7].every((v, i) => Math.abs(y4[i] - v) < 1e-12));
      t.check('N = 5 ≥ L+M−1: circular = linear [1,3,5,7,4]', [1, 3, 5, 7, 4].every((v, i) => Math.abs(y5[i] - v) < 1e-12)); }
    { const x = [1, 2, 3, 4], X = D.dft(x, 4), Xs = D.dft(D.circShift(x, 1, 4), 4); let e = 0;
      for (let k = 0; k < 4; k++) { const w = cmul(C(X.re[k], X.im[k]), C(Math.cos(-TAU * k / 4), Math.sin(-TAU * k / 4))); e = Math.max(e, Math.abs(w.re - Xs.re[k]), Math.abs(w.im - Xs.im[k])); }
      t.check('circular shift by 1 ↔ X[k]·W₄^k', e < 1e-12);
      const r = D.timeRev(x, 4); t.check('time reversal mod 4: [1,4,3,2]', [1, 4, 3, 2].every((v, i) => r[i] === v)); }
    // aliasing 7 kHz at 10 kHz -> 3 kHz (phase flips)
    { const a = D.alias(7000, 10000); t.check('7 kHz sampled at 10 kHz aliases to 3 kHz (phase negated)', a.Fa === 3000 && a.sign === -1 && a.aliased); }
    t.check('12 kHz at Fs = 10 kHz aliases to 2 kHz', D.alias(12000, 10000).Fa === 2000 && D.alias(12000, 10000).sign === 1);
    // periodicity: ω = 3π/8 -> f = 3/16 -> N = 16
    { const p = D.periodicityOmega(D.parseOmega('3pi/8')); t.check('ω = 3π/8 is periodic with N = 16 (k = 3)', p.periodic && p.N === 16 && p.k === 3); }
    t.check('ω = 1 rad is not periodic', D.periodicityOmega(D.parseOmega('1')).periodic === false);
    t.check('ω = π/4 → N = 8; ω = 0.75π → N = 8', D.periodicityOmega(D.parseOmega('π/4')).N === 8 && D.periodicityOmega(D.parseOmega('0.75pi')).N === 8);
    t.check('F = 3 kHz, Fs = 8 kHz -> f = 3/8, N = 8', D.periodicityF(3000, 8000).N === 8);
    t.check('parseOmega rejects junk', D.parseOmega('abc') === null && D.parseOmega('3pi/0') === null);
    // quantisation SNR ≈ 6.02B + 1.76 dB
    t.check('8-bit quantiser SNR ≈ 49.9 dB (±0.3)', t.near(D.quantSNR(8), 6.02 * 8 + 1.76, 0.3), D.quantSNR(8).toFixed(2));
    t.check('12-bit quantiser SNR ≈ 74.0 dB (±0.3)', t.near(D.quantSNR(12), 6.02 * 12 + 1.76, 0.3), D.quantSNR(12).toFixed(2));
    // invalid input
    t.check('normBA rejects a₀ = 0 and empty lists', !!D.normBA([1], [0, 1]).err && !!D.normBA([], [1]).err);
    t.check('parseList: fractions ok, junk -> null', JSON.stringify(D.parseList('1, -1/2 0.25')) === '[1,-0.5,0.25]' && D.parseList('1, x') === null);
    t.check('designWindow rejects unordered edges', !!D.designWindow({ resp: 'lp', edges: [0.5 * PI, 0.3 * PI], win: 'hann', d1: 0.01, As: 40, orderMode: 'auto' }).err);
    t.check('stability: pole at 1.1 unstable, at 0.9 stable', D.stability([{ r: C(1.1, 0), m: 1 }]).verdict === 'unstable' && D.stability([{ r: C(0.9, 0), m: 1 }]).verdict === 'stable');
  });


})();
