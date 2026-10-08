/* Communications lab: constellation + Monte Carlo BER, eye diagram, Shannon calculator, CIDR calculator.
   Pure math in FSP.math.comms (no DOM). UI only inside init/activate. */
(function () {
  'use strict';

  /* ===================================================================== math */
  const SQRT2 = Math.SQRT2, SQRTPI = Math.sqrt(Math.PI);

  // erf by the all-positive series  erf(x) = 2/sqrt(pi) * exp(-x^2) * sum_{n>=0} 2^n x^(2n+1) / (1*3*5*...*(2n+1))
  function erfSeries(x) {
    let term = x, sum = x;
    for (let n = 1; n < 500; n++) { term *= 2 * x * x / (2 * n + 1); sum += term; if (term < sum * 1e-17) break; }
    return 2 / SQRTPI * Math.exp(-x * x) * sum;
  }
  // erfc for x >= 2 by the Laplace continued fraction, evaluated bottom-up:
  // erfc(x) = exp(-x^2)/sqrt(pi) * 1/(x + (1/2)/(x + (2/2)/(x + (3/2)/(x + ...))))
  function erfcCF(x) {
    let f = x;
    for (let k = 120; k >= 1; k--) f = x + (k / 2) / f;
    return Math.exp(-x * x) / SQRTPI / f;
  }
  function erfc(x) {
    if (Number.isNaN(x)) return NaN;
    if (x < 0) return 2 - erfc(-x);
    if (x === 0) return 1;
    if (x < 2) return 1 - erfSeries(x);   // erfc >= 4.7e-3 here, so the subtraction costs <3 digits of 16
    if (x > 27) return 0;
    return erfcCF(x);
  }
  function erf(x) { return x < 0 ? -erf(-x) : (x < 2 ? erfSeries(x) : 1 - erfcCF(x)); }
  function Q(x) { return 0.5 * erfc(x / SQRT2); }

  const MODS = {
    bpsk: { name: 'BPSK', k: 1 },
    qpsk: { name: 'QPSK', k: 2 },
    qam16: { name: '16-QAM', k: 4 },
  };
  function ebn0Lin(db) { return Math.pow(10, db / 10); }

  // Theoretical BER at Eb/N0 (dB). BPSK/QPSK exact (Gray); 16-QAM nearest-neighbour approximation.
  function theoryBER(mod, ebn0dB) {
    const g = ebn0Lin(ebn0dB);
    if (mod === 'qam16') return 3 / 8 * erfc(Math.sqrt(0.4 * g));
    return Q(Math.sqrt(2 * g));
  }

  // Per-axis Gray levels (integer amplitude units). m=1: bit 0->-1, 1->+1.  m=2: 00->-3, 01->-1, 11->+1, 10->+3.
  const AX2 = [-3, -1, 3, 1];   // index = b1*2+b0 -> level
  function axisLevel(m, v) { return m === 1 ? (v ? 1 : -1) : AX2[v]; }
  function axisBits(m, x) {      // hard decision on x (amplitude units) -> integer of the bits
    if (m === 1) return x >= 0 ? 1 : 0;
    return x < -2 ? 0 : x < 0 ? 1 : x < 2 ? 3 : 2;   // 00,01,11,10
  }
  function modInfo(mod) {
    // bits per axis, number of quadrature axes, amplitude unit so that mean symbol energy Es = 1
    if (mod === 'bpsk') return { m: 1, axes: 1, unit: 1, k: 1 };
    if (mod === 'qpsk') return { m: 1, axes: 2, unit: 1 / SQRT2, k: 2 };
    if (mod === 'qam16') return { m: 2, axes: 2, unit: 1 / Math.sqrt(10), k: 4 };
    throw new Error('unknown modulation ' + mod);
  }
  // Constellation: [{bits:int, bitStr, re, im}] normalised to Es = 1. bits packed MSB = I-axis bits first.
  function constellation(mod) {
    const mi = modInfo(mod), pts = [], nA = 1 << mi.m;
    for (let a = 0; a < (mi.axes === 2 ? nA : 1); a++) for (let b = 0; b < nA; b++) {
      // symbol bits: I-axis bits (value b), then Q-axis bits (value a) for 2-axis
      const re = axisLevel(mi.m, b) * mi.unit, im = mi.axes === 2 ? axisLevel(mi.m, a) * mi.unit : 0;
      const bits = mi.axes === 2 ? (b << mi.m) | a : b;
      let s = bits.toString(2); while (s.length < mi.k) s = '0' + s;
      pts.push({ bits, bitStr: s, re, im });
    }
    return pts;
  }

  // Wilson 95% interval for errors/n
  function wilson(err, n, z) {
    z = z || 1.959964; if (!(n > 0)) return [NaN, NaN];
    const p = err / n, z2 = z * z, d = 1 + z2 / n, c = p + z2 / (2 * n), h = z * Math.sqrt(p * (1 - p) / n + z2 / (4 * n * n));
    return [Math.max(0, (c - h) / d), Math.min(1, (c + h) / d)];
  }

  // Monte Carlo BER over AWGN. Seeded; deterministic. opts: {mod, ebn0dB, seed, nMax (bits), minErrors, keep (scatter symbols)}
  function simulateBER(opts) {
    const mod = opts.mod, mi = modInfo(mod), k = mi.k;
    const nMax = Math.max(k, Math.floor(opts.nMax || 200000)), minErr = opts.minErrors === undefined ? 100 : opts.minErrors;
    const keep = opts.keep === undefined ? 0 : opts.keep;
    const seed = (opts.seed === undefined ? 1 : opts.seed) >>> 0;
    const bitRng = FSP.prng.mulberry32(seed), noise = FSP.prng.gaussian(FSP.prng.mulberry32((seed ^ 0x9E3779B9) >>> 0));
    // Es = 1, Eb = Es/k, N0 = Eb/(Eb/N0), noise variance per real dimension = N0/2
    const Eb = 1 / k, N0 = Eb / ebn0Lin(opts.ebn0dB), sigma = Math.sqrt(N0 / 2);
    const m = mi.m, two = mi.axes === 2, u = mi.unit, inv = 1 / u;
    const rx = new Float64Array(2 * keep);
    let bits = 0, errors = 0, nsym = 0, kept = 0;
    while (bits + k <= nMax) {
      let bi = 0, bq = 0;
      for (let j = 0; j < m; j++) bi = (bi << 1) | (bitRng() < 0.5 ? 1 : 0);
      if (two) for (let j = 0; j < m; j++) bq = (bq << 1) | (bitRng() < 0.5 ? 1 : 0);
      const ri = axisLevel(m, bi) * u + sigma * noise();
      const rq = two ? axisLevel(m, bq) * u + sigma * noise() : 0;
      if (kept < keep) { rx[2 * kept] = ri; rx[2 * kept + 1] = rq; kept++; }
      let x = axisBits(m, ri * inv) ^ bi, e = 0; while (x) { e += x & 1; x >>= 1; }
      if (two) { x = axisBits(m, rq * inv) ^ bq; while (x) { e += x & 1; x >>= 1; } }
      errors += e; bits += k; nsym++;
      if (errors >= minErr) break;
    }
    const ber = bits > 0 ? errors / bits : NaN, ci = wilson(errors, bits);
    return { mod, ebn0dB: opts.ebn0dB, bits, errors, symbols: nsym, ber, ci, sigma, theory: theoryBER(mod, opts.ebn0dB),
      stopped: errors >= minErr ? 'errors' : 'nmax', rx: rx.subarray(0, 2 * kept) };
  }

  /* ---- raised cosine (T = 1) ---- */
  function sinc(x) { return Math.abs(x) < 1e-12 ? 1 : Math.sin(Math.PI * x) / (Math.PI * x); }
  function raisedCosine(t, beta) {
    if (Math.abs(t) < 1e-12) return 1;
    if (beta <= 0) return sinc(t);
    const d = 1 - 4 * beta * beta * t * t;
    if (Math.abs(d) < 1e-9) return Math.PI / 4 * sinc(1 / (2 * beta));   // removable singularity at |t| = 1/(2 beta)
    return sinc(t) * Math.cos(Math.PI * beta * t) / d;
  }

  // Eye traces for binary PAM (+/-1) with RC pulse; each trace spans t in [-1, 1] symbol periods around a sampling instant.
  // jitter: rms timing error (fraction of T), noise: rms amplitude noise per sample. Returns {t, traces, centre}.
  function eyeTraces(o) {
    const beta = o.beta, n = o.n || 200, sps = o.sps || 20, S = 8, seed = (o.seed === undefined ? 7 : o.seed) >>> 0;
    const rng = FSP.prng.mulberry32(seed), g = FSP.prng.gaussian(FSP.prng.mulberry32((seed ^ 0xA5A5A5A5) >>> 0));
    const a = new Float64Array(n + 2 * S + 2); for (let i = 0; i < a.length; i++) a[i] = rng() < 0.5 ? -1 : 1;
    const nt = 2 * sps + 1, t = new Float64Array(nt); for (let j = 0; j < nt; j++) t[j] = -1 + j / sps;
    const traces = [], centre = [], sym = [];
    for (let i = 0; i < n; i++) {
      const c = i + S + 1, tau = (o.jitter || 0) * g(), y = new Float64Array(nt);
      for (let j = 0; j < nt; j++) {
        let v = 0; for (let d = -S; d <= S; d++) v += a[c + d] * raisedCosine(t[j] + tau - d, beta);
        y[j] = v + (o.noise || 0) * g();
      }
      traces.push(y); centre.push(y[sps]); sym.push(a[c]);
    }
    return { t, traces, centre, sym, sps };
  }
  // Vertical eye opening at the nominal sampling instant (t = 0): min(+1 branch) - max(-1 branch). Negative = closed.
  function eyeOpening(e) {
    let lo = Infinity, hi = -Infinity;
    for (let i = 0; i < e.centre.length; i++) { if (e.sym[i] > 0) lo = Math.min(lo, e.centre[i]); else hi = Math.max(hi, e.centre[i]); }
    return Number.isFinite(lo) && Number.isFinite(hi) ? lo - hi : NaN;
  }

  /* ---- Shannon ---- */
  // inputs: B (Hz), snrDb, C (bit/s); `solve` = 'B' | 'snr' | 'C' selects the unknown.
  function shannon(inp, solve) {
    const bad = m => ({ ok: false, error: m, B: NaN, snrDb: NaN, C: NaN });
    let B = inp.B, snrDb = inp.snrDb, C = inp.C;
    if (solve !== 'C') { if (!Number.isFinite(C) || C <= 0) return bad('Capacity C must be a positive number.'); }
    if (solve !== 'B') { if (!Number.isFinite(B) || B <= 0) return bad('Bandwidth B must be a positive number.'); }
    if (solve !== 'snr') { if (!Number.isFinite(snrDb)) return bad('SNR (dB) must be a number.'); }
    if (solve === 'C') C = B * Math.log2(1 + Math.pow(10, snrDb / 10));
    else if (solve === 'B') B = C / Math.log2(1 + Math.pow(10, snrDb / 10));
    else {
      const eta = C / B; if (eta > 1000) return bad('C/B is too large (> 1000 bit/s/Hz); SNR out of range.');
      snrDb = 10 * Math.log10(Math.pow(2, eta) - 1);
    }
    if (![B, snrDb, C].every(Number.isFinite)) return bad('Result is out of numeric range.');
    const eta = C / B, ebn0min = (Math.pow(2, eta) - 1) / eta;
    return { ok: true, B, snrDb, C, snrLin: Math.pow(10, snrDb / 10), eta, ebn0minDb: 10 * Math.log10(ebn0min) };
  }

  /* ---- IPv4 / CIDR ---- */
  function parseIPv4(s) {
    if (typeof s !== 'string') return { ok: false, error: 'Address is required.' };
    s = s.trim(); if (!s) return { ok: false, error: 'Address is required.' };
    const p = s.split('.');
    if (p.length !== 4) return { ok: false, error: 'IPv4 address needs exactly 4 dot-separated octets.' };
    let v = 0;
    for (let i = 0; i < 4; i++) {
      if (!/^\d{1,3}$/.test(p[i])) return { ok: false, error: 'Octet ' + (i + 1) + ' ("' + p[i] + '") is not a decimal number 0-255.' };
      if (p[i].length > 1 && p[i][0] === '0') return { ok: false, error: 'Octet ' + (i + 1) + ' has a leading zero (ambiguous octal); write it without.' };
      const o = parseInt(p[i], 10); if (o > 255) return { ok: false, error: 'Octet ' + (i + 1) + ' = ' + o + ' is greater than 255.' };
      v = v * 256 + o;
    }
    return { ok: true, value: v };
  }
  function parsePrefix(s) {
    s = String(s).trim().replace(/^\//, '');
    if (!/^\d{1,3}$/.test(s)) return { ok: false, error: 'Prefix length must be an integer 0-32.' };
    const n = parseInt(s, 10); if (n > 32) return { ok: false, error: 'Prefix length ' + n + ' is greater than 32.' };
    return { ok: true, value: n };
  }
  function toDotted(v) { return [Math.floor(v / 16777216) % 256, Math.floor(v / 65536) % 256, Math.floor(v / 256) % 256, v % 256].join('.'); }
  function maskOf(p) { return p === 0 ? 0 : (0xFFFFFFFF - (Math.pow(2, 32 - p) - 1)); }
  function subnetInfo(ip, p) {
    const size = Math.pow(2, 32 - p), mask = maskOf(p), net = Math.floor(ip / size) * size, bc = net + size - 1;
    let first, last, hosts, note = '';
    if (p === 32) { first = last = net; hosts = 1; note = '/32 host route: a single address.'; }
    else if (p === 31) { first = net; last = bc; hosts = 2; note = '/31 point-to-point link (RFC 3021): both addresses usable, no network/broadcast.'; }
    else { first = net + 1; last = bc - 1; hosts = size - 2; }
    return { ip, prefix: p, network: net, broadcast: bc, mask, wildcard: 0xFFFFFFFF - mask, first, last, hosts, total: size, note };
  }
  // "a.b.c.d/n" -> {ok, info} | {ok:false, error}
  function parseCIDR(str) {
    if (typeof str !== 'string' || !str.trim()) return { ok: false, error: 'Enter an address in CIDR form, e.g. 192.168.1.130/26.' };
    const parts = str.trim().split('/');
    if (parts.length !== 2) return { ok: false, error: 'Use CIDR form address/prefix, e.g. 192.168.1.130/26.' };
    const ip = parseIPv4(parts[0]); if (!ip.ok) return ip;
    const pf = parsePrefix(parts[1]); if (!pf.ok) return pf;
    return { ok: true, info: subnetInfo(ip.value, pf.value) };
  }
  // Split a block into N equal subnets (N rounded up to a power of two). limit caps the returned list length.
  function splitSubnet(info, N, limit) {
    if (!Number.isInteger(N) || N < 1) return { ok: false, error: 'Number of subnets must be a whole number >= 1.' };
    const bits = N === 1 ? 0 : Math.ceil(Math.log2(N)), np = info.prefix + bits;
    if (np > 32) return { ok: false, error: 'Cannot split a /' + info.prefix + ' into ' + N + ' subnets (would need /' + np + ').' };
    const count = Math.pow(2, bits), size = Math.pow(2, 32 - np), shown = Math.min(count, limit || 256), list = [];
    for (let i = 0; i < shown; i++) list.push(subnetInfo(info.network + i * size, np));
    return { ok: true, count, newPrefix: np, list, truncated: shown < count, rounded: count !== N };
  }

  // "100, 50 25" -> {ok, hosts:[100,50,25]}
  function parseHostList(str) {
    if (typeof str !== 'string' || !str.trim()) return { ok: false, error: 'Enter the required host counts, e.g. 100, 50, 25.' };
    const toks = str.split(/[\s,;]+/).filter(Boolean), hosts = [];
    for (const t of toks) {
      if (!/^\d+$/.test(t)) return { ok: false, error: '"' + t + '" is not a whole number of hosts.' };
      const n = parseInt(t, 10); if (n < 1) return { ok: false, error: 'Each subnet needs at least 1 host.' };
      if (n > 4294967294) return { ok: false, error: n + ' hosts cannot fit in IPv4.' };
      hosts.push(n);
    }
    if (hosts.length > 64) return { ok: false, error: 'At most 64 subnets.' };
    return { ok: true, hosts };
  }
  // Largest aligned CIDR blocks covering [start, end] (inclusive)
  function rangeToBlocks(start, end) {
    const out = [];
    while (start <= end) {
      let size = 1;
      while (size * 2 <= 4294967296 && start % (size * 2) === 0 && start + size * 2 - 1 <= end) size *= 2;
      out.push(subnetInfo(start, 32 - Math.round(Math.log2(size)))); start += size;
    }
    return out;
  }
  // VLSM: allocate subnets largest-first from the network of `info`. Each needs hosts + 2 addresses (network + broadcast),
  // rounded up to a power of two (classic VLSM; /31 point-to-point is not used, so 1-2 hosts -> /30).
  function vlsm(info, hosts) {
    if (!Array.isArray(hosts) || !hosts.length) return { ok: false, error: 'No host counts given.' };
    const base = info.network, end = info.network + info.total - 1;
    const order = hosts.map((h, i) => ({ h, i })).sort((a, b) => b.h - a.h || a.i - b.i);
    let cur = base; const list = [];
    for (const o of order) {
      const bits = Math.ceil(Math.log2(o.h + 2)), size = Math.pow(2, bits), prefix = 32 - bits;
      if (prefix < info.prefix) return { ok: false, error: 'Subnet ' + (o.i + 1) + ' needs ' + o.h + ' hosts → /' + prefix + ' (' + size + ' addresses), larger than the whole /' + info.prefix + ' block (' + info.total + ').' };
      cur = Math.ceil((cur - base) / size) * size + base;
      if (cur + size - 1 > end) {
        const need = order.reduce((s, x) => s + Math.pow(2, Math.ceil(Math.log2(x.h + 2))), 0);
        return { ok: false, error: 'Does not fit: the subnets need ' + need.toLocaleString() + ' addresses in total but /' + info.prefix + ' has ' + info.total.toLocaleString() + ' (ran out at subnet ' + (o.i + 1) + ', ' + o.h + ' hosts).' };
      }
      const s = subnetInfo(cur, prefix);
      list.push(Object.assign(s, { index: o.i, need: o.h, size, wasted: s.hosts - o.h }));
      cur += size;
    }
    const used = cur - base, leftover = cur <= end ? rangeToBlocks(cur, end) : [];
    return { ok: true, list, used, free: info.total - used, leftover, base: info };
  }

  FSP.math.comms = { erf, erfc, Q, theoryBER, constellation, modInfo, simulateBER, wilson, raisedCosine, eyeTraces, eyeOpening, shannon,
    parseIPv4, parsePrefix, parseCIDR, subnetInfo, splitSubnet, toDotted, MODS, parseHostList, vlsm, rangeToBlocks };

  /* ===================================================================== tests */
  FSP.registerTests('comms', function (t) {
    const M = FSP.math.comms;
    // erfc / Q against scipy
    [[0.1, 0.8875370839817152], [0.5, 0.4795001221869535], [1, 0.15729920705028516], [1.999, null], [2, 0.004677734981047266],
      [2.5, 0.00040695201744495886], [4, 1.541725790028002e-08], [6, 2.1519736712498913e-17], [10, 2.0884875837625446e-45]].forEach(([x, ref]) => {
      if (ref !== null) t.check('erfc(' + x + ') rel err < 1e-12', t.rel(M.erfc(x), ref, 1e-12), M.erfc(x) + ' vs ' + ref);
    });
    t.check('erfc either side of the x=2 seam (scipy)', t.rel(M.erfc(1.9999), 0.00467980209297061, 1e-12) && t.rel(M.erfc(2.0001), 0.0046756686958033394, 1e-12) && t.rel(M.erfc(1.5), 0.03389485352468927, 1e-12), '');
    t.check('erfc(-x) = 2 - erfc(x), erfc(0)=1', t.near(M.erfc(-0.7), 2 - M.erfc(0.7), 1e-15) && M.erfc(0) === 1, '');
    t.check('Q(0)=0.5, Q(1)=0.15865525393145707', t.near(M.Q(0), 0.5, 1e-15) && t.rel(M.Q(1), 0.15865525393145707, 1e-12), '');
    // theory BER (spec values)
    [[0, 0.07865], [4, 0.012501], [8, 1.909e-4]].forEach(([d, ref]) => t.check('BPSK theory BER at ' + d + ' dB', t.rel(M.theoryBER('bpsk', d), ref, 1e-3), M.theoryBER('bpsk', d)));
    t.check('QPSK theory = BPSK theory (equal Eb/N0)', M.theoryBER('qpsk', 5) === M.theoryBER('bpsk', 5), '');
    t.check('16-QAM theory at 12 dB = (3/8)erfc(sqrt(0.4*10^1.2))', t.rel(M.theoryBER('qam16', 12), 0.00013865868881261898, 1e-10), M.theoryBER('qam16', 12));
    // Gray mapping + unit energy
    ['bpsk', 'qpsk', 'qam16'].forEach(mod => {
      const c = M.constellation(mod), mi = M.modInfo(mod);
      const es = c.reduce((s, p) => s + p.re * p.re + p.im * p.im, 0) / c.length;
      let dmin = Infinity; for (let i = 0; i < c.length; i++) for (let j = i + 1; j < c.length; j++) dmin = Math.min(dmin, Math.hypot(c[i].re - c[j].re, c[i].im - c[j].im));
      let ok = c.length === (1 << mi.k), cnt = 0, uniq = new Set(c.map(p => p.bits)).size === c.length;
      for (let i = 0; i < c.length; i++) for (let j = 0; j < c.length; j++) if (i !== j && Math.hypot(c[i].re - c[j].re, c[i].im - c[j].im) < dmin * 1.0000001) {
        cnt++; let x = c[i].bits ^ c[j].bits, pop = 0; while (x) { pop += x & 1; x >>= 1; } if (pop !== 1) ok = false;
      }
      t.check(mod + ': adjacent points differ by exactly 1 bit (' + cnt + ' neighbour pairs)', ok && uniq && cnt > 0, '');
      t.check(mod + ': mean symbol energy = 1', t.near(es, 1, 1e-12), es);
    });
    // Monte Carlo
    const a = M.simulateBER({ mod: 'bpsk', ebn0dB: 4, seed: 12345, nMax: 200000, minErrors: Infinity });
    const b = M.simulateBER({ mod: 'bpsk', ebn0dB: 4, seed: 12345, nMax: 200000, minErrors: Infinity });
    t.check('MC BPSK 4 dB, 2e5 bits, within 15% of theory', t.rel(a.ber, M.theoryBER('bpsk', 4), 0.15), 'BER ' + a.ber + ' vs ' + a.theory);
    t.check('MC BPSK bit-identical across two runs', a.errors === b.errors && a.bits === b.bits && a.ber === b.ber, a.errors + ' vs ' + b.errors);
    const c2 = M.simulateBER({ mod: 'bpsk', ebn0dB: 4, seed: 999, nMax: 200000, minErrors: Infinity });
    t.check('MC different seed gives different error count', c2.errors !== a.errors, '');
    const q = M.simulateBER({ mod: 'qpsk', ebn0dB: 4, seed: 9001, nMax: 400000, minErrors: Infinity });
    const bp = M.simulateBER({ mod: 'bpsk', ebn0dB: 4, seed: 4242, nMax: 400000, minErrors: Infinity });
    t.check('MC QPSK BER == BPSK BER at equal Eb/N0 (within 8%)', t.rel(q.ber, bp.ber, 0.08) && t.rel(q.ber, M.theoryBER('qpsk', 4), 0.08), q.ber + ' vs ' + bp.ber);
    const q16 = M.simulateBER({ mod: 'qam16', ebn0dB: 8, seed: 77, nMax: 400000, minErrors: Infinity });
    t.check('MC 16-QAM at 8 dB within 10% of (3/8)erfc approx', t.rel(q16.ber, M.theoryBER('qam16', 8), 0.10), q16.ber + ' vs ' + q16.theory);
    const st = M.simulateBER({ mod: 'bpsk', ebn0dB: 0, seed: 3, nMax: 1e6, minErrors: 100 });
    t.check('MC stops once >= 100 errors (reason "errors")', st.stopped === 'errors' && st.errors >= 100 && st.errors < 101 + 1 && st.bits < 1e6, st.errors + ' errors in ' + st.bits);
    const nm = M.simulateBER({ mod: 'bpsk', ebn0dB: 12, seed: 3, nMax: 5000, minErrors: 100 });
    t.check('MC stops at Nmax with few errors (reason "nmax")', nm.stopped === 'nmax' && nm.bits === 5000, nm.errors);
    t.check('Wilson CI brackets the estimate and theory', a.ci[0] < a.ber && a.ber < a.ci[1] && a.ci[0] < a.theory * 1.2 && a.ci[1] > a.theory * 0.8, a.ci.join(','));
    // Raised cosine
    [0.25, 0.5, 1].forEach(beta => {
      let m = 0; for (let n = 1; n <= 12; n++) m = Math.max(m, Math.abs(M.raisedCosine(n, beta)), Math.abs(M.raisedCosine(-n, beta)));
      t.check('RC pulse zeros at integer t, beta=' + beta + ' (max |p| = ' + m.toExponential(1) + ')', m < 1e-9 && M.raisedCosine(0, beta) === 1, '');
      const ts = 1 / (2 * beta), pS = M.raisedCosine(ts, beta), pN = M.raisedCosine(ts + 1e-5, beta), pM = M.raisedCosine(ts - 1e-5, beta);
      t.check('RC singularity t=1/(2 beta) finite and continuous, beta=' + beta, Number.isFinite(pS) && Math.abs(pS - pN) < 1e-4 && Math.abs(pS - pM) < 1e-4,
        pS + ' vs ' + pN);
    });
    t.check('RC value at t=0.5, beta=1 equals pi/4 * sinc(0.5)', t.near(M.raisedCosine(0.5, 1), Math.PI / 4 * (2 / Math.PI), 1e-12), '');
    t.check('RC beta=0.5, t=0.3 vs numpy', t.near(M.raisedCosine(0.3, 0.5), 0.840477339891049, 1e-12), '');
    // Eye
    const e0 = M.eyeTraces({ beta: 0.35, n: 300, jitter: 0, noise: 0, seed: 5 });
    t.check('Eye (no jitter/noise) opening = 2 (zero ISI at sampling instant)', t.near(M.eyeOpening(e0), 2, 1e-8), M.eyeOpening(e0));
    const e1 = M.eyeTraces({ beta: 0.35, n: 300, jitter: 0.1, noise: 0.05, seed: 5 });
    t.check('Eye opening shrinks with jitter + noise', M.eyeOpening(e1) < 1.9 && M.eyeOpening(e1) > -2, M.eyeOpening(e1));
    const e2 = M.eyeTraces({ beta: 0.35, n: 300, jitter: 0.1, noise: 0.05, seed: 5 });
    t.check('Eye traces deterministic for a fixed seed', e1.traces[17][9] === e2.traces[17][9], '');
    // Shannon
    const s1 = M.shannon({ B: 1e6, snrDb: 30 }, 'C');
    t.check('Shannon B=1 MHz, SNR=30 dB -> 9.967 Mbit/s', t.rel(s1.C, 9.967e6, 1e-4) && t.rel(s1.C, 9967226.258835994, 1e-9), s1.C);
    const s2 = M.shannon({ C: s1.C, snrDb: 30 }, 'B'), s3 = M.shannon({ B: 1e6, C: s1.C }, 'snr');
    t.check('Shannon solve B and SNR round-trip', t.rel(s2.B, 1e6, 1e-12) && t.near(s3.snrDb, 30, 1e-9), s2.B + ', ' + s3.snrDb);
    t.check('Shannon C=B => SNR = 0 dB (1 bit/s/Hz); eta=1 gives Eb/N0 min 0 dB', t.near(M.shannon({ B: 5, C: 5 }, 'snr').snrDb, 0, 1e-12) && t.near(M.shannon({ B: 5, C: 5 }, 'snr').ebn0minDb, 0, 1e-12), '');
    t.check('Shannon Eb/N0 min -> -1.59 dB as eta->0', t.near(M.shannon({ B: 1e6, C: 1 }, 'snr').ebn0minDb, 10 * Math.log10(Math.LN2), 1e-4), '');
    t.check('Shannon rejects bad input', !M.shannon({ B: -1, snrDb: 10 }, 'C').ok && !M.shannon({ B: 1, C: 5000 }, 'snr').ok && !M.shannon({ B: NaN, snrDb: 10 }, 'C').ok, '');
    // CIDR
    const c1 = M.parseCIDR('192.168.1.130/26').info;
    t.check('192.168.1.130/26 -> net .128, bcast .191, 62 hosts, mask 255.255.255.192, wildcard 0.0.0.63',
      M.toDotted(c1.network) === '192.168.1.128' && M.toDotted(c1.broadcast) === '192.168.1.191' && c1.hosts === 62 && M.toDotted(c1.mask) === '255.255.255.192' && M.toDotted(c1.wildcard) === '0.0.0.63'
      && M.toDotted(c1.first) === '192.168.1.129' && M.toDotted(c1.last) === '192.168.1.190', '');
    const c31 = M.parseCIDR('10.0.0.4/31').info, c32 = M.parseCIDR('10.0.0.7/32').info, c8 = M.parseCIDR('10.0.0.0/8').info, c0 = M.parseCIDR('1.2.3.4/0').info;
    t.check('/31 -> 2 usable (RFC 3021), first=network, last=broadcast', c31.hosts === 2 && M.toDotted(c31.first) === '10.0.0.4' && M.toDotted(c31.last) === '10.0.0.5', '');
    t.check('/32 -> 1 host', c32.hosts === 1 && M.toDotted(c32.first) === '10.0.0.7' && M.toDotted(c32.last) === '10.0.0.7' && M.toDotted(c32.mask) === '255.255.255.255', '');
    t.check('10.0.0.0/8 -> 16,777,214 hosts', c8.hosts === 16777214 && M.toDotted(c8.broadcast) === '10.255.255.255', c8.hosts);
    t.check('/0 -> mask 0.0.0.0, 4294967294 hosts, high-bit addresses handled', c0.hosts === 4294967294 && M.toDotted(c0.mask) === '0.0.0.0' && M.toDotted(M.parseCIDR('255.255.255.255/1').info.network) === '128.0.0.0', '');
    t.check('Invalid: octet 256, prefix 33, junk rejected with messages',
      !M.parseCIDR('192.168.1.256/24').ok && !M.parseCIDR('192.168.1.1/33').ok && !M.parseCIDR('192.168.1/24').ok && !M.parseCIDR('192.168.1.1').ok
      && !M.parseCIDR('a.b.c.d/8').ok && !M.parseCIDR('1.2.3.4/-1').ok && !M.parseCIDR('1.2.3.4/').ok && !M.parseCIDR('').ok && /256|255/.test(M.parseCIDR('192.168.1.256/24').error)
      && /33|32/.test(M.parseCIDR('1.1.1.1/33').error), '');
    const sp = M.splitSubnet(M.parseCIDR('192.168.1.0/24').info, 4);
    t.check('Split /24 into 4 -> four /26 blocks .0 .64 .128 .192', sp.ok && sp.newPrefix === 26 && sp.list.map(x => M.toDotted(x.network)).join() === '192.168.1.0,192.168.1.64,192.168.1.128,192.168.1.192', '');
    const sp3 = M.splitSubnet(M.parseCIDR('10.0.0.0/24').info, 3);
    t.check('Split into 3 rounds up to 4 and says so; /31 -> 2 x /32; over-split rejected', sp3.ok && sp3.count === 4 && sp3.rounded && M.splitSubnet(M.parseCIDR('1.1.1.0/31').info, 2).newPrefix === 32
      && !M.splitSubnet(M.parseCIDR('1.1.1.1/31').info, 4).ok && !M.splitSubnet(M.parseCIDR('1.1.1.1/24').info, 0).ok && !M.splitSubnet(M.parseCIDR('1.1.1.1/24').info, 2.5).ok, '');
    // VLSM 192.168.1.0/24 for 100, 50, 25 hosts (given out of order):
    // 100+2 -> 128 = /25: .0-.127 (usable .1-.126, 126 hosts); 50+2 -> 64 = /26: .128-.191; 25+2 -> 32 = /27: .192-.223; free .224/27
    const v1 = M.vlsm(M.parseCIDR('192.168.1.0/24').info, M.parseHostList('25, 100 50').hosts), D = M.toDotted;
    t.check('VLSM /24 -> 100, 50, 25 hosts: /25 .0, /26 .128, /27 .192, leftover .224/27', v1.ok
      && v1.list.map(s => D(s.network) + '/' + s.prefix).join() === '192.168.1.0/25,192.168.1.128/26,192.168.1.192/27'
      && v1.list.map(s => s.need).join() === '100,50,25' && v1.list[0].index === 1 && D(v1.list[1].broadcast) === '192.168.1.191' && D(v1.list[2].first) === '192.168.1.193' && D(v1.list[2].last) === '192.168.1.222'
      && D(v1.list[1].mask) === '255.255.255.192' && v1.free === 32 && v1.leftover.length === 1 && D(v1.leftover[0].network) + '/' + v1.leftover[0].prefix === '192.168.1.224/27', '');
    // 200 -> /24 (256) and 100 -> /25 (128): 384 > 256 -> error. 2 hosts -> /30. Leftover of a /24 after one /26 = /26 + /25.
    const v2 = M.vlsm(M.parseCIDR('10.0.0.0/24').info, [200, 100]), v3 = M.vlsm(M.parseCIDR('10.0.0.0/24').info, [2]), v4 = M.vlsm(M.parseCIDR('10.0.0.0/24').info, [60]);
    t.check('VLSM: overflow rejected with message; 2 hosts -> /30; leftover split into aligned blocks', !v2.ok && /384/.test(v2.error) && !M.vlsm(M.parseCIDR('10.0.0.0/24').info, [300]).ok
      && v3.ok && v3.list[0].prefix === 30 && v4.leftover.map(b => D(b.network) + '/' + b.prefix).join() === '10.0.0.64/26,10.0.0.128/25'
      && !M.parseHostList('10, x').ok && !M.parseHostList('0').ok && !M.parseHostList('').ok, v2.error);
  });

  /* ===================================================================== UI */
  let built = false, active = false, ui = {}, redrawTimer = 0;
  const S = { mod: 'qpsk', ebn0: 6, seed: 1, nmax: 200000, beta: 0.35, jitter: 0.03, noise: 0.05, solve: 'C', sB: '1M', sSnr: '30', sC: '10M', cidr: '192.168.1.130/26', split: 4, cmode: 'split', vlsm: '100, 50, 25' };
  let sim = null, sweep = [], eye = null;

  function css(name, fb) { try { const v = getComputedStyle(document.documentElement).getPropertyValue(name).trim(); return v || fb; } catch (e) { return fb; } }
  function fmtBER(x) { return Number.isFinite(x) ? (x === 0 ? '0' : x.toExponential(3)) : '—'; }

  function scheduleRun(now) { clearTimeout(redrawTimer); redrawTimer = setTimeout(runSim, now ? 0 : 120); FSP.state.touch(); }

  function runSim() {
    redrawTimer = 0; if (!built) return;
    sim = FSP.math.comms.simulateBER({ mod: S.mod, ebn0dB: S.ebn0, seed: S.seed, nMax: S.nmax, minErrors: 100, keep: 1500 });
    renderSim();
  }
  function renderSim() {
    const r = sim, M = FSP.math.comms; if (!r) return;
    const ratio = r.theory > 0 && r.errors > 0 ? r.ber / r.theory : NaN;
    ui.hud.innerHTML = '';
    [['Modulation', M.MODS[r.mod].name + ' (' + M.MODS[r.mod].k + ' bit/sym)'], ['Bits simulated', r.bits.toLocaleString()], ['Errors', String(r.errors)],
      ['Measured BER', fmtBER(r.ber)], ['95% CI (Wilson)', fmtBER(r.ci[0]) + ' … ' + fmtBER(r.ci[1])], ['Theory BER' + (r.mod === 'qam16' ? ' (approx)' : ''), fmtBER(r.theory)],
      ['Measured / theory', Number.isFinite(ratio) ? ratio.toFixed(3) : '—'],
      ['Stopped', r.stopped === 'errors' ? 'at ≥100 errors' : r.errors < 100 ? 'at Nmax (only ' + r.errors + ' errors: wide CI)' : 'at Nmax']]
      .forEach(([a, b]) => ui.hud.appendChild(FSP.ui.el('div', null, FSP.ui.el('span', { text: a }), FSP.ui.el('span', { text: b }))));
    ui.simWork.set(['Es = 1 (constellation normalised), k = ' + M.MODS[r.mod].k + ' bit/symbol  =>  Eb = Es/k = ' + FSP.fmt(1 / M.MODS[r.mod].k, 4),
      'Eb/N0 = ' + S.ebn0 + ' dB = ' + FSP.fmt(Math.pow(10, S.ebn0 / 10), 4) + ' (linear)  =>  N0 = Eb/(Eb/N0) = ' + FSP.fmt((1 / M.MODS[r.mod].k) / Math.pow(10, S.ebn0 / 10), 5),
      'Noise std per real dimension  sigma = sqrt(N0/2) = ' + FSP.fmt(r.sigma, 5),
      'Gray mapping per axis; hard decision at the axis thresholds (0 for BPSK/QPSK; 0, ±2/sqrt(10) for 16-QAM).',
      r.mod === 'qam16' ? 'Theory: (3/8) erfc( sqrt(0.4 Eb/N0) ) = ' + fmtBER(r.theory) + '  (nearest-neighbour approximation)' : 'Theory: Q( sqrt(2 Eb/N0) ) = ' + fmtBER(r.theory),
      'Seed ' + S.seed + ' (bits and noise use separate mulberry32 streams); same seed + parameters give bit-identical results.',
      'Wilson 95% interval: p ± z·sqrt(p(1-p)/n + z²/4n²) over (1+z²/n), z = 1.96.']);
    drawConst(); drawBER();
  }

  function drawConst() {
    if (!built) return; const { ctx, w, h } = ui.cvC.prep(); const M = FSP.math.comms;
    ctx.clearRect(0, 0, w, h); const R = 1.7, side = Math.min(w, h), cx = w / 2, cy = h / 2, sc = side / (2 * R);
    const X = x => cx + x * sc, Y = y => cy - y * sc;
    ctx.strokeStyle = css('--border', '#243047'); ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(cx - side / 2, cy); ctx.lineTo(cx + side / 2, cy); ctx.moveTo(cx, cy - side / 2); ctx.lineTo(cx, cy + side / 2); ctx.stroke();
    ctx.fillStyle = css('--muted', '#8b98ab'); ctx.font = '11px ' + css('--mono', 'monospace'); ctx.textAlign = 'left'; ctx.fillText('I', cx + side / 2 - 10, cy - 4); ctx.fillText('Q', cx + 4, cy - side / 2 + 12);
    if (sim) {
      ctx.fillStyle = 'rgba(76,201,240,0.35)'; const rx = sim.rx;
      for (let i = 0; i < rx.length; i += 2) { const x = X(rx[i]), y = Y(rx[i + 1]); ctx.fillRect(x - 1.2, y - 1.2, 2.4, 2.4); }
    }
    const pts = M.constellation(S.mod); ctx.lineWidth = 1.5; ctx.textAlign = 'center';
    pts.forEach(p => {
      ctx.strokeStyle = css('--accent2', '#f9a03f'); ctx.beginPath(); ctx.arc(X(p.re), Y(p.im), 5, 0, 2 * Math.PI); ctx.stroke();
      if (pts.length <= 16) { ctx.fillStyle = css('--text', '#e6edf3'); ctx.fillText(p.bitStr, X(p.re), Y(p.im) - 9); }
    });
  }

  function drawBER() {
    if (!built) return; const { ctx, w, h } = ui.cvB.prep(); const M = FSP.math.comms; ctx.clearRect(0, 0, w, h);
    const L = 44, Rr = 10, T = 10, B = 28, pw = w - L - Rr, ph = h - T - B, x0 = 0, x1 = 14, y0 = -6, y1 = 0;
    const X = d => L + (d - x0) / (x1 - x0) * pw, Y = lg => T + (y1 - lg) / (y1 - y0) * ph;
    ctx.font = '11px ' + css('--mono', 'monospace'); ctx.lineWidth = 1; ctx.textAlign = 'right'; ctx.fillStyle = css('--muted', '#8b98ab');
    for (let e = y0; e <= y1; e++) { ctx.strokeStyle = css('--border', '#243047'); ctx.beginPath(); ctx.moveTo(L, Y(e)); ctx.lineTo(L + pw, Y(e)); ctx.stroke(); ctx.fillText('1e' + e, L - 4, Y(e) + 4); }
    ctx.textAlign = 'center'; for (let d = 0; d <= x1; d += 2) { ctx.strokeStyle = css('--border', '#243047'); ctx.beginPath(); ctx.moveTo(X(d), T); ctx.lineTo(X(d), T + ph); ctx.stroke(); ctx.fillText(String(d), X(d), h - 12); }
    ctx.fillText('Eb/N0 (dB)', L + pw / 2, h - 1);
    const cols = { bpsk: css('--accent', '#4cc9f0'), qpsk: css('--ok', '#2ecc71'), qam16: css('--pink', '#f72585') };
    ctx.save(); ctx.beginPath(); ctx.rect(L, T, pw, ph); ctx.clip();
    ['bpsk', 'qam16'].forEach(m => {
      ctx.strokeStyle = cols[m]; ctx.lineWidth = m === S.mod || (m === 'bpsk' && S.mod === 'qpsk') ? 2.2 : 1; ctx.beginPath(); let first = true;
      for (let d = x0; d <= x1; d += 0.1) { const b = M.theoryBER(m, d); if (!(b > 0)) continue; const x = X(d), y = Y(Math.log10(b)); if (first) { ctx.moveTo(x, y); first = false; } else ctx.lineTo(x, y); }
      ctx.stroke();
    });
    sweep.forEach(p => { if (!(p.errors > 0)) return; ctx.strokeStyle = ctx.fillStyle = cols[p.mod]; const x = X(p.ebn0dB); const lo = Math.max(p.ci[0], 1e-9), hi = p.ci[1];
      ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(x, Y(Math.log10(lo))); ctx.lineTo(x, Y(Math.log10(hi))); ctx.stroke(); ctx.beginPath(); ctx.arc(x, Y(Math.log10(p.ber)), 3, 0, 2 * Math.PI); ctx.fill(); });
    if (sim && sim.errors > 0) {
      const x = X(sim.ebn0dB), y = Y(Math.log10(sim.ber)); ctx.strokeStyle = css('--accent2', '#f9a03f'); ctx.fillStyle = css('--accent2', '#f9a03f'); ctx.lineWidth = 1.5;
      ctx.beginPath(); ctx.moveTo(x, Y(Math.log10(Math.max(sim.ci[0], 1e-9)))); ctx.lineTo(x, Y(Math.log10(sim.ci[1]))); ctx.stroke();
      ctx.beginPath(); ctx.moveTo(x, y - 6); ctx.lineTo(x + 6, y); ctx.lineTo(x, y + 6); ctx.lineTo(x - 6, y); ctx.closePath(); ctx.fill();
    }
    ctx.restore();
    ctx.textAlign = 'left'; ctx.lineWidth = 2; let lx = L + 8;
    [['bpsk', 'BPSK/QPSK theory'], ['qam16', '16-QAM approx']].forEach(([m, lab]) => { ctx.strokeStyle = cols[m]; ctx.beginPath(); ctx.moveTo(lx, T + 10); ctx.lineTo(lx + 16, T + 10); ctx.stroke(); ctx.fillStyle = css('--text', '#e6edf3'); ctx.fillText(lab, lx + 20, T + 14); lx += 20 + ctx.measureText(lab).width + 12; });
  }

  function runSweep() {
    sweep = []; const mods = ['bpsk', 'qpsk', 'qam16'];
    mods.forEach(m => { for (let d = 0; d <= 12; d += 2) sweep.push(FSP.math.comms.simulateBER({ mod: m, ebn0dB: d, seed: S.seed + d, nMax: S.nmax, minErrors: 100 })); });
    drawBER();
  }

  /* ---- eye ---- */
  function runEye() {
    if (!built) return; eye = FSP.math.comms.eyeTraces({ beta: S.beta, jitter: S.jitter, noise: S.noise, n: 250, seed: S.seed }); drawEye();
    const op = FSP.math.comms.eyeOpening(eye);
    ui.eyeOut.value = Number.isFinite(op) ? op.toFixed(3) + (op <= 0 ? '  (closed)' : '  of 2.000 ideal') : '—';
    ui.eyeWork.set(['Binary PAM a_k = ±1, pulse p(t) = sinc(t)·cos(πβt)/(1 − (2βt)²), T = 1, β = ' + S.beta,
      'Singularity at |t| = 1/(2βT): p = (π/4)·sinc(1/(2β)) (L\'Hôpital limit).',
      'y(t) = Σ_d a_(k+d) p(t + τ − d) + n; τ ~ N(0, ' + S.jitter + '² T²) per trace, n ~ N(0, ' + S.noise + '²) per sample.',
      'Zeros of p at every nonzero integer t make the opening exactly 2 at t = 0 when jitter = noise = 0.',
      'Eye opening = min(y(0) | a=+1) − max(y(0) | a=−1) over ' + eye.traces.length + ' traces (nominal sampling instant).']);
  }
  function drawEye() {
    if (!built || !eye) return; const { ctx, w, h } = ui.cvE.prep(); ctx.clearRect(0, 0, w, h);
    const L = 34, Rr = 8, T = 8, B = 22, pw = w - L - Rr, ph = h - T - B, X = t => L + (t + 1) / 2 * pw, Y = v => T + (1 - (v + 2) / 4) * ph;
    ctx.font = '11px ' + css('--mono', 'monospace'); ctx.fillStyle = css('--muted', '#8b98ab'); ctx.strokeStyle = css('--border', '#243047'); ctx.lineWidth = 1;
    [-1, 0, 1].forEach(v => { ctx.beginPath(); ctx.moveTo(L, Y(v)); ctx.lineTo(L + pw, Y(v)); ctx.stroke(); ctx.textAlign = 'right'; ctx.fillText(String(v), L - 4, Y(v) + 4); });
    [-1, -0.5, 0, 0.5, 1].forEach(t => { ctx.beginPath(); ctx.moveTo(X(t), T); ctx.lineTo(X(t), T + ph); ctx.stroke(); ctx.textAlign = 'center'; ctx.fillText(String(t), X(t), h - 8); });
    ctx.textAlign = 'right'; ctx.fillText('t/T', L + pw - 4, T + 12);
    ctx.strokeStyle = 'rgba(76,201,240,0.22)'; ctx.lineWidth = 1;
    eye.traces.forEach(y => { ctx.beginPath(); for (let j = 0; j < y.length; j++) { const px = X(eye.t[j]), py = Y(y[j]); if (j) ctx.lineTo(px, py); else ctx.moveTo(px, py); } ctx.stroke(); });
    ctx.strokeStyle = css('--accent2', '#f9a03f'); ctx.setLineDash([4, 3]); ctx.beginPath(); ctx.moveTo(X(0), T); ctx.lineTo(X(0), T + ph); ctx.stroke(); ctx.setLineDash([]);
  }

  /* ---- shannon ---- */
  function updShannon() {
    const M = FSP.math.comms, num = el => FSP.parseSI(el.value), sv = S.solve;
    const r = M.shannon({ B: num(ui.sh.B), snrDb: num(ui.sh.snr), C: num(ui.sh.C) }, sv);
    ['B', 'snr', 'C'].forEach(k => { ui.sh[k].disabled = k === sv; ui.sh[k].setAttribute('aria-readonly', k === sv ? 'true' : 'false'); });
    ui.shMsg.className = 'msg err'; ui.shMsg.hidden = r.ok;
    if (!r.ok) { ui.shMsg.textContent = r.error; ['snrLin', 'eta', 'eb'].forEach(k => ui.shOut[k].value = '—'); ui.shWork.set('Fix the inputs above.'); const el = ui.sh[sv]; el.value = '—'; return; }
    if (sv === 'B') ui.sh.B.value = FSP.fmtNum(r.B, 7); else if (sv === 'snr') ui.sh.snr.value = FSP.fmtNum(r.snrDb, 7); else ui.sh.C.value = FSP.fmtNum(r.C, 7);
    ui.shOut.snrLin.value = FSP.fmtNum(r.snrLin, 6); ui.shOut.eta.value = FSP.fmtNum(r.eta, 6) + ' bit/s/Hz'; ui.shOut.eb.value = FSP.fmt(r.ebn0minDb, 3) + ' dB';
    ui.shWork.set(['C = B · log2(1 + SNR)', 'SNR = 10^(' + FSP.fmt(r.snrDb, 4) + '/10) = ' + FSP.fmtNum(r.snrLin, 6),
      'log2(1 + SNR) = ' + FSP.fmt(Math.log2(1 + r.snrLin), 6) + ' bit/s/Hz',
      'C = ' + FSP.fmtEng(r.B, 'Hz') + ' × ' + FSP.fmt(r.eta, 6) + ' = ' + FSP.fmtEng(r.C, 'bit/s', 7),
      'Minimum Eb/N0 at this spectral efficiency: (2^η − 1)/η = ' + FSP.fmt(Math.pow(10, r.ebn0minDb / 10), 5) + ' = ' + FSP.fmt(r.ebn0minDb, 3) + ' dB (limit −1.592 dB as η → 0)']);
    S.sB = ui.sh.B.value; S.sSnr = ui.sh.snr.value; S.sC = ui.sh.C.value;
  }

  /* ---- cidr ---- */
  function updCidr() {
    const M = FSP.math.comms, p = M.parseCIDR(S.cidr), D = M.toDotted;
    ui.cidrOut.innerHTML = ''; ui.splitOut.textContent = '';
    ui.cidrIn.setAttribute('aria-invalid', p.ok ? 'false' : 'true');
    if (!p.ok) { ui.cidrMsg.hidden = false; ui.cidrMsg.textContent = p.error; ui.cidrWork.set('Fix the input above.'); ui.splitMsg.hidden = true; ui.vlsmMsg.hidden = true; ui.vlsmOut.textContent = ''; return; }
    ui.cidrMsg.hidden = true; const i = p.info, bin = v => D(v).split('.').map(o => ('00000000' + (+o).toString(2)).slice(-8)).join('.');
    [['Address', D(i.ip)], ['Network', D(i.network) + '/' + i.prefix], ['Broadcast', i.prefix >= 31 ? D(i.broadcast) + (i.prefix === 32 ? ' (n/a)' : ' (n/a, RFC 3021)') : D(i.broadcast)],
      ['Subnet mask', D(i.mask)], ['Wildcard', D(i.wildcard)], ['Host range', D(i.first) + (i.first === i.last ? '' : ' – ' + D(i.last))],
      ['Usable hosts', i.hosts.toLocaleString()], ['Total addresses', i.total.toLocaleString()]]
      .forEach(([a, b]) => ui.cidrOut.appendChild(FSP.ui.el('div', null, FSP.ui.el('span', { text: a }), FSP.ui.el('span', { text: b }))));
    if (i.note) ui.cidrOut.appendChild(FSP.ui.el('div', { class: 'note', text: i.note }));
    ui.cidrWorkBase = (['Address  ' + bin(i.ip), 'Mask     ' + bin(i.mask) + '  (/' + i.prefix + ')', 'Network  ' + bin(i.network) + '  = address AND mask', 'Bcast    ' + bin(i.broadcast) + '  = network OR wildcard',
      'Total = 2^(32 − ' + i.prefix + ') = ' + i.total.toLocaleString() + '; usable = ' + (i.prefix >= 31 ? (i.prefix === 31 ? '2 (RFC 3021 point-to-point)' : '1 (host route)') : 'total − 2 = ' + i.hosts.toLocaleString())]);
    ui.cidrWork.set(ui.cidrWorkBase);
    ui.splitBox.hidden = S.cmode !== 'split'; ui.vlsmBox.hidden = S.cmode !== 'vlsm';
    if (S.cmode === 'vlsm') { updVlsm(i); return; }
    const sp = M.splitSubnet(i, S.split, 64);
    if (!sp.ok) { ui.splitMsg.hidden = false; ui.splitMsg.textContent = sp.error; return; }
    ui.splitMsg.hidden = !sp.rounded; ui.splitMsg.className = 'msg warn'; ui.splitMsg.textContent = 'Rounded up to ' + sp.count + ' subnets (a power of two) so that all blocks are equal.';
    ui.splitOut.textContent = sp.count + ' × /' + sp.newPrefix + '  (' + (sp.list[0].hosts).toLocaleString() + ' usable hosts each)\n' +
      sp.list.map((s, n) => String(n + 1).padStart(2) + '. ' + D(s.network) + '/' + s.prefix + '  ' + D(s.first) + ' – ' + D(s.last) + (s.prefix >= 31 ? '' : '  bc ' + D(s.broadcast))).join('\n') +
      (sp.truncated ? '\n… ' + (sp.count - sp.list.length) + ' more' : '');
  }

  function updVlsm(base) {
    const M = FSP.math.comms, D = M.toDotted, ph = M.parseHostList(S.vlsm);
    ui.vlsmOut.textContent = '';
    const r = ph.ok ? M.vlsm(base, ph.hosts) : ph;
    if (!r.ok) { ui.vlsmMsg.hidden = false; ui.vlsmMsg.textContent = r.error; ui.cidrWork.set(ui.cidrWorkBase.concat(['', 'VLSM: ' + r.error])); return; }
    ui.vlsmMsg.hidden = true;
    const pad = (x, n) => String(x).padEnd(n);
    ui.vlsmOut.textContent = 'Base ' + D(base.network) + '/' + base.prefix + ' (' + base.total.toLocaleString() + ' addresses), allocated largest first\n' +
      pad('#', 3) + pad('need', 7) + pad('subnet', 20) + pad('mask', 17) + pad('usable range', 33) + 'broadcast\n' +
      r.list.map(x => pad(x.index + 1, 3) + pad(x.need, 7) + pad(D(x.network) + '/' + x.prefix, 20) + pad(D(x.mask), 17) + pad(D(x.first) + ' – ' + D(x.last), 33) + D(x.broadcast) + '   (' + x.hosts.toLocaleString() + ' usable, ' + x.wasted.toLocaleString() + ' spare)').join('\n') +
      '\nUsed ' + r.used.toLocaleString() + ' of ' + base.total.toLocaleString() + ' addresses; free ' + r.free.toLocaleString() + (r.leftover.length ? ': ' + r.leftover.slice(0, 16).map(b => D(b.network) + '/' + b.prefix).join(', ') + (r.leftover.length > 16 ? ', …' : '') : '');
    const W = ui.cidrWorkBase.concat(['', 'VLSM (classic: each subnet reserves network + broadcast, so it needs hosts + 2 addresses rounded up to a power of two):',
      '1. Sort requests largest first: ' + r.list.map(x => x.need).join(', ') + '.', '2. For each, smallest 2^k ≥ hosts + 2 gives prefix /(32 − k); place it at the next free address (largest-first keeps every block aligned).']);
    r.list.forEach(x => W.push('   ' + x.need + ' hosts: ' + x.need + ' + 2 = ' + (x.need + 2) + ' ≤ 2^' + (32 - x.prefix) + ' = ' + x.size + ' → /' + x.prefix + ' at ' + D(x.network) + ', broadcast ' + D(x.broadcast) + ', next free ' + D(x.network + x.size)));
    W.push('3. Leftover ' + r.free.toLocaleString() + ' addresses' + (r.leftover.length ? ' = ' + r.leftover.map(b => D(b.network) + '/' + b.prefix).slice(0, 16).join(' + ') : '') + '.');
    ui.cidrWork.set(W);
  }

  /* ---- build ---- */
  function numInput(parent, label, value, onChange, attrs) {
    const inp = FSP.ui.el('input', Object.assign({ type: 'text', inputmode: 'decimal', 'aria-label': label, size: 10 }, attrs || {})); inp.value = value;
    inp.addEventListener('input', () => onChange(inp)); parent.appendChild(FSP.ui.el('div', { class: 'ctl' }, FSP.ui.el('label', { text: label }), inp, FSP.ui.el('span')));
    return inp;
  }

  function build(panel) {
    const U = FSP.ui, M = FSP.math.comms;
    // 1. BER
    const f1 = U.fieldset(panel, '1 · Constellation & Monte Carlo BER (Gray, AWGN)'), g1 = U.el('div', { class: 'grid2' }), ctl1 = U.el('div'), plots1 = U.el('div');
    const wrapC = U.el('div', { class: 'canvas-wrap' }), wrapB = U.el('div', { class: 'canvas-wrap' }); wrapB.style.marginTop = '12px';
    ui.cvC = U.canvas(wrapC, { height: 300 }); ui.cvB = U.canvas(wrapB, { height: 260 }); plots1.appendChild(wrapC); plots1.appendChild(wrapB);
    ui.modSel = U.select(ctl1, 'Modulation', [['bpsk', 'BPSK'], ['qpsk', 'QPSK'], ['qam16', '16-QAM']], S.mod, v => { S.mod = v; scheduleRun(true); });
    ui.ebn0 = U.slider(ctl1, { label: 'Eb/N0', min: 0, max: 14, step: 0.1, value: S.ebn0, unit: 'dB', onInput: v => { S.ebn0 = v; scheduleRun(); } });
    ui.nmaxSel = U.select(ctl1, 'Nmax (bits)', [['20000', '2·10⁴'], ['200000', '2·10⁵'], ['1000000', '10⁶'], ['5000000', '5·10⁶']], String(S.nmax), v => { S.nmax = +v; scheduleRun(true); });
    ui.seed = numInput(ctl1, 'Seed', String(S.seed), el => { const v = parseInt(el.value, 10); if (Number.isFinite(v) && v >= 0) { S.seed = v >>> 0; scheduleRun(); } });
    const rowB = U.el('div', { class: 'row' }); U.button(rowB, 'Re-run', () => runSim()); U.button(rowB, 'Sweep 0–12 dB (all 3)', () => runSweep()); U.button(rowB, 'Clear sweep', () => { sweep = []; drawBER(); }); ctl1.appendChild(rowB);
    ctl1.appendChild(U.el('div', { class: 'note', text: 'Stops at ≥100 errors or Nmax bits. Orange diamond = this run with 95% CI; dots = sweep. QPSK sits on the BPSK curve at equal Eb/N0.' }));
    ui.hud = U.el('div', { class: 'hud' }); ctl1.appendChild(ui.hud);
    g1.appendChild(ctl1); g1.appendChild(plots1); f1.appendChild(g1); ui.simWork = U.working(f1);
    ui.cvC.onResize(drawConst); ui.cvB.onResize(drawBER);

    // 2. Eye
    const f2 = U.fieldset(panel, '2 · Eye diagram (raised-cosine pulse, binary PAM)'), g2 = U.el('div', { class: 'grid2' }), ctl2 = U.el('div'), wrapE = U.el('div', { class: 'canvas-wrap' });
    ui.cvE = U.canvas(wrapE, { height: 260 }); ui.cvE.onResize(drawEye);
    ui.beta = U.slider(ctl2, { label: 'Roll-off β', min: 0, max: 1, step: 0.01, value: S.beta, onInput: v => { S.beta = v; runEye(); FSP.state.touch(); } });
    ui.jit = U.slider(ctl2, { label: 'Jitter rms', min: 0, max: 0.2, step: 0.005, value: S.jitter, unit: 'T', onInput: v => { S.jitter = v; runEye(); FSP.state.touch(); } });
    ui.nz = U.slider(ctl2, { label: 'Noise σ', min: 0, max: 0.5, step: 0.005, value: S.noise, onInput: v => { S.noise = v; runEye(); FSP.state.touch(); } });
    const eo = U.el('div', { class: 'ctl' }, U.el('span', { text: 'Eye opening' })); ui.eyeOut = U.el('output', { text: '—' }); eo.appendChild(ui.eyeOut); eo.appendChild(U.el('span')); ctl2.appendChild(eo);
    ctl2.appendChild(U.el('div', { class: 'note', text: 'Traces span 2 symbol periods; the dashed line is the sampling instant. Low β: narrow bandwidth but large timing sensitivity; β = 1: widest eye horizontally.' }));
    g2.appendChild(ctl2); g2.appendChild(wrapE); f2.appendChild(g2); ui.eyeWork = U.working(f2);

    // 3. Shannon
    const f3 = U.fieldset(panel, '3 · Shannon capacity  C = B·log₂(1+SNR)'), ctl3 = U.el('div'); ui.sh = {};
    U.select(ctl3, 'Solve for', [['C', 'Capacity C'], ['B', 'Bandwidth B'], ['snr', 'SNR']], S.solve, v => { S.solve = v; updShannon(); FSP.state.touch(); });
    ui.sh.B = numInput(ctl3, 'B (Hz)', S.sB, () => { updShannon(); FSP.state.touch(); }); ui.sh.snr = numInput(ctl3, 'SNR (dB)', S.sSnr, () => { updShannon(); FSP.state.touch(); }); ui.sh.C = numInput(ctl3, 'C (bit/s)', S.sC, () => { updShannon(); FSP.state.touch(); });
    ui.shSel = ctl3.querySelector('select');
    ui.shMsg = U.el('div', { class: 'msg err', hidden: '' }); ctl3.appendChild(ui.shMsg);
    ui.shOut = { snrLin: U.readout(ctl3, 'SNR (linear)'), eta: U.readout(ctl3, 'Spectral eff.'), eb: U.readout(ctl3, 'Min Eb/N0') };
    ctl3.appendChild(U.el('div', { class: 'note', text: 'SI suffixes accepted (1M, 56k). The solved-for box is read-only.' }));
    f3.appendChild(ctl3); ui.shWork = U.working(f3);

    // 4. CIDR
    const f4 = U.fieldset(panel, '4 · IPv4 CIDR / subnet calculator'), ctl4 = U.el('div'); ui.cidrIn = numInput(ctl4, 'CIDR', S.cidr, el => { S.cidr = el.value; updCidr(); FSP.state.touch(); }, { inputmode: 'text', size: 18, spellcheck: 'false' });
    ui.cidrMsg = U.el('div', { class: 'msg err', hidden: '', role: 'alert' }); ctl4.appendChild(ui.cidrMsg);
    ui.cidrOut = U.el('div', { class: 'hud' }); ctl4.appendChild(ui.cidrOut);
    ui.cmodeSel = U.select(ctl4, 'Subnetting', [['split', 'Equal split'], ['vlsm', 'VLSM']], S.cmode, v => { S.cmode = v; updCidr(); FSP.state.touch(); });
    ui.splitBox = U.el('div'); ctl4.appendChild(ui.splitBox);
    ui.splitIn = numInput(ui.splitBox, 'Split into N', String(S.split), el => { const v = Number(el.value); S.split = el.value.trim() === '' ? NaN : v; updCidr(); FSP.state.touch(); }, { inputmode: 'numeric' });
    ui.splitMsg = U.el('div', { class: 'msg err', hidden: '' }); ui.splitBox.appendChild(ui.splitMsg);
    ui.splitOut = U.el('pre', { class: 'mono note' }); ui.splitOut.style.overflowX = 'auto'; ui.splitBox.appendChild(ui.splitOut);
    ui.vlsmBox = U.el('div', { hidden: '' }); ctl4.appendChild(ui.vlsmBox);
    ui.vlsmIn = numInput(ui.vlsmBox, 'Hosts per subnet', S.vlsm, el => { S.vlsm = el.value; updCidr(); FSP.state.touch(); }, { inputmode: 'text', size: 18, spellcheck: 'false' });
    ui.vlsmBox.appendChild(U.el('div', { class: 'note', text: 'Comma-separated host counts; the base block is the network of the CIDR above.' }));
    ui.vlsmMsg = U.el('div', { class: 'msg err bad', hidden: '', role: 'alert' }); ui.vlsmBox.appendChild(ui.vlsmMsg);
    ui.vlsmOut = U.el('pre', { class: 'mono note' }); ui.vlsmOut.style.overflowX = 'auto'; ui.vlsmOut.style.whiteSpace = 'pre'; ui.vlsmBox.appendChild(ui.vlsmOut);
    f4.appendChild(ctl4); ui.cidrWork = U.working(f4);
    built = true;
  }
  function syncControls() {
    if (!built) return;
    ui.modSel.value = S.mod; ui.ebn0.set(S.ebn0, true); ui.nmaxSel.value = String(S.nmax); ui.seed.value = String(S.seed);
    ui.beta.set(S.beta, true); ui.jit.set(S.jitter, true); ui.nz.set(S.noise, true);
    ui.shSel.value = S.solve; ui.sh.B.value = S.sB; ui.sh.snr.value = S.sSnr; ui.sh.C.value = S.sC; ui.cidrIn.value = S.cidr; ui.splitIn.value = String(S.split); ui.cmodeSel.value = S.cmode; ui.vlsmIn.value = S.vlsm;
  }
  function refreshAll() { runSim(); runEye(); updShannon(); updCidr(); }

  FSP.state.bind('comms', {
    get: () => ({ mod: S.mod, ebn0: S.ebn0, seed: S.seed, nmax: S.nmax, beta: S.beta, jit: S.jitter, nz: S.noise, solve: S.solve, B: S.sB, snr: S.sSnr, C: S.sC, cidr: S.cidr, split: S.split, cmode: S.cmode, vlsm: S.vlsm }),
    set(p) {
      const num = (k, lo, hi, d) => { const v = parseFloat(p[k]); return Number.isFinite(v) && v >= lo && v <= hi ? v : d; };
      if (p.mod && FSP.math.comms.MODS[p.mod]) S.mod = p.mod;
      S.ebn0 = num('ebn0', 0, 14, S.ebn0); S.seed = Math.floor(num('seed', 0, 4294967295, S.seed)); S.nmax = [20000, 200000, 1000000, 5000000].indexOf(+p.nmax) >= 0 ? +p.nmax : S.nmax;
      S.beta = num('beta', 0, 1, S.beta); S.jitter = num('jit', 0, 0.2, S.jitter); S.noise = num('nz', 0, 0.5, S.noise);
      if (p.solve === 'B' || p.solve === 'C' || p.solve === 'snr') S.solve = p.solve;
      ['B:sB', 'snr:sSnr', 'C:sC'].forEach(s => { const [k, f] = s.split(':'); if (typeof p[k] === 'string' && p[k].length < 40) S[f] = p[k]; });
      if (typeof p.cidr === 'string' && p.cidr.length < 60) S.cidr = p.cidr;
      const sp = parseInt(p.split, 10); if (Number.isFinite(sp) && sp >= 1 && sp <= 1048576) S.split = sp;
      if (p.cmode === 'split' || p.cmode === 'vlsm') S.cmode = p.cmode;
      if (typeof p.vlsm === 'string' && p.vlsm.length < 400) S.vlsm = p.vlsm;
      if (built) { syncControls(); if (active) refreshAll(); }
    },
  });

  FSP.registerTab({
    id: 'comms', title: 'Comms & Networks',
    init(panel) { build(panel); syncControls(); },
    activate() { active = true; refreshAll(); },
    deactivate() { active = false; clearTimeout(redrawTimer); redrawTimer = 0; },
  });
})();
