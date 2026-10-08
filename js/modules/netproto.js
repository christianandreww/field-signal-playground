/* Networks & Protocols (NTU IE3017 Computer Communications): line coding, delays & multiplexing,
   multiple access, ARQ & TCP congestion control, error detection.
   Pure math lives in FSP.math.net (no DOM). All DOM work happens inside init/activate. */
(function () {
  'use strict';
  const fin = Number.isFinite, PI = Math.PI;

  /* ============================================================ 1. line coding */
  const SCHEMES = [['unrz', 'Unipolar NRZ'], ['pnrz', 'Polar NRZ'], ['nrzi', 'NRZ-I'], ['ami', 'Bipolar AMI'],
    ['man', 'Manchester'], ['dman', 'Diff. Manchester'], ['b8zs', 'B8ZS'], ['hdb3', 'HDB3']];

  // "0100 1110" -> {ok, bits:[0,1,0,0,1,1,1,0]}
  function parseBits(s, maxLen) {
    maxLen = maxLen || 64;
    if (typeof s !== 'string') return { ok: false, error: 'Enter a bit string such as 01001110.' };
    const t = s.replace(/[\s_,]/g, '');
    if (!t) return { ok: false, error: 'Enter at least one bit (0 or 1).' };
    if (!/^[01]+$/.test(t)) return { ok: false, error: 'Only the characters 0 and 1 are allowed (spaces are ignored).' };
    if (t.length > maxLen) return { ok: false, error: 'At most ' + maxLen + ' bits (got ' + t.length + ').' };
    return { ok: true, bits: t.split('').map(Number) };
  }

  // Encode to half-bit levels (+1, 0, -1; unipolar uses 1/0). o: {manchester:'course'|'ieee', prior:+1|-1 (line level
  // before bit 0, used by NRZ-I and diff. Manchester), lastPulse:+1|-1 (polarity of the last mark before bit 0, AMI family),
  // hdb3Odd: true if an odd number of pulses has been sent since the last HDB3 substitution}.
  // -> {half:[2n], marks:[{i, label:'V'|'B'}], subs:[{start,len,pattern}], prior}
  function lineEncode(bits, scheme, o) {
    o = o || {};
    const n = bits.length, half = new Array(2 * n).fill(0), marks = [], subs = [];
    const prior = o.prior === 1 ? 1 : -1, lastPulse = o.lastPulse === 1 ? 1 : -1, ieee = o.manchester === 'ieee';
    const put = (i, a, b) => { half[2 * i] = a; half[2 * i + 1] = b === undefined ? a : b; };
    let L = prior, last = lastPulse;
    switch (scheme) {
      case 'unrz': bits.forEach((b, i) => put(i, b ? 1 : 0)); break;
      case 'pnrz': bits.forEach((b, i) => put(i, b ? 1 : -1)); break;
      case 'nrzi': bits.forEach((b, i) => { if (b) L = -L; put(i, L); }); break;
      case 'man': bits.forEach((b, i) => { const hl = ieee ? !b : !!b; put(i, hl ? 1 : -1, hl ? -1 : 1); }); break;
      case 'dman': bits.forEach((b, i) => { const first = b ? L : -L; put(i, first, -first); L = -first; }); break;
      case 'ami': bits.forEach((b, i) => { if (b) { last = -last; put(i, last); } else put(i, 0); }); break;
      case 'b8zs': {
        let run = 0;
        bits.forEach((b, i) => {
          if (b) { last = -last; put(i, last); run = 0; return; }
          put(i, 0); run++;
          if (run === 8) {          // 000VB0VB, V = same polarity as the preceding pulse
            const P = last;
            put(i - 4, P); put(i - 3, -P); put(i - 1, -P); put(i, P);
            marks.push({ i: i - 4, label: 'V' }, { i: i - 3, label: 'B' }, { i: i - 1, label: 'V' }, { i: i, label: 'B' });
            subs.push({ start: i - 7, len: 8, pattern: '000VB0VB' });
            run = 0;                // last pulse polarity is P again
          }
        });
        break;
      }
      case 'hdb3': {
        let run = 0, count = o.hdb3Odd ? 1 : 0;
        bits.forEach((b, i) => {
          if (b) { last = -last; put(i, last); run = 0; count++; return; }
          put(i, 0); run++;
          if (run === 4) {
            if (count % 2 === 1) {  // odd: 000V, V same polarity as the preceding pulse
              put(i, last); marks.push({ i, label: 'V' }); subs.push({ start: i - 3, len: 4, pattern: '000V' });
            } else {                // even: B00V, B obeys AMI (opposite of the preceding pulse), V = same as B
              const B = -last; put(i - 3, B); put(i, B); last = B;
              marks.push({ i: i - 3, label: 'B' }, { i, label: 'V' }); subs.push({ start: i - 3, len: 4, pattern: 'B00V' });
            }
            run = 0; count = 0;
          }
        });
        break;
      }
      default: return null;
    }
    return { scheme, half, marks, subs, prior: (scheme === 'nrzi' || scheme === 'dman') ? prior : undefined };
  }
  // DC (mean level in units of V), number of level transitions (including the one at t = 0 against the prior level for
  // NRZ-I / diff. Manchester), longest constant stretch in bit times.
  function lineStats(enc) {
    const h = enc.half; let sum = 0, tr = 0, run = 1, best = h.length ? 1 : 0;
    for (let i = 0; i < h.length; i++) {
      sum += h[i];
      if (i > 0) { if (h[i] !== h[i - 1]) { tr++; run = 1; } else { run++; if (run > best) best = run; } }
    }
    if (enc.prior !== undefined && h.length && h[0] !== enc.prior) tr++;
    return { dc: h.length ? sum / h.length : NaN, transitions: tr, longestFlat: best / 2 };
  }

  /* ============================================================ 2. delays & multiplexing */
  function nodal(p) { // {L, R, d, v, tproc, tqueue}
    const tt = p.L / p.R, tp = p.d / p.v;
    return { tTrans: tt, tProp: tp, tProc: p.tproc, tQueue: p.tqueue, total: tt + tp + p.tproc + p.tqueue, a: tp / tt, bdp: p.R * tp };
  }
  // Message M bits split into N = ceil(M/P) packets of payload P (last one may be shorter), each with header h,
  // over S store-and-forward switches (S+1 identical links of rate R, propagation tprop, per-hop processing tproc).
  function switching(p) {
    const { M, P, h, S, R, tprop, tproc } = p;
    if (!(M > 0 && P > 0 && h >= 0 && R > 0 && tprop >= 0 && tproc >= 0) || !Number.isInteger(S) || S < 0)
      return { ok: false, error: 'Need M > 0, P > 0, h ≥ 0, R > 0, S a whole number ≥ 0.' };
    const N = Math.ceil(M / P - 1e-12), Plast = M - (N - 1) * P, tPkt = (P + h) / R, tLast = (Plast + h) / R;
    const totalEq = (N + S) * tPkt + (S + 1) * (tprop + tproc);                      // all N packets full-size
    const total = (N + S - 1) * tPkt + tLast + (S + 1) * (tprop + tproc);             // exact with a short last packet
    const totalMsg = (S + 1) * ((M + h) / R + tprop + tproc);
    const Popt = S > 0 ? Math.sqrt(M * h / S) : Infinity;
    return { ok: true, N, Plast, tPkt, tLast, total, totalEq, totalMsg, Popt, overhead: N * h, equal: Math.abs(Plast - P) < 1e-9 };
  }
  // Event-by-event store-and-forward pipeline (independent of the closed form). sizes in bits.
  function simulatePackets(sizes, S, R, tprop, tproc) {
    const links = S + 1, free = new Array(links).fill(0), hops = [];
    let total = 0;
    sizes.forEach(Lb => {
      let ready = 0; const row = [];
      for (let j = 0; j < links; j++) {
        const start = Math.max(ready, free[j]), end = start + Lb / R;
        row.push({ start, end }); free[j] = end; ready = end + tprop + tproc;
      }
      hops.push(row); total = Math.max(total, ready);
    });
    return { hops, total };
  }
  // Binomial pmf by logs (no overflow for N up to ~1e5)
  function binomPmf(N, p) {
    const out = new Array(N + 1).fill(0);
    if (p <= 0) { out[0] = 1; return out; }
    if (p >= 1) { out[N] = 1; return out; }
    const lp = Math.log(p), lq = Math.log1p(-p); let lc = 0;
    for (let n = 0; n <= N; n++) { out[n] = Math.exp(lc + n * lp + (N - n) * lq); lc += Math.log(N - n) - Math.log(n + 1); }
    return out;
  }
  // N users, each active w.p. p at rate r, link capacity C. Circuit switching supports M = floor(C/r) users.
  function statMux(p) {
    const { C, r, N, p: q } = p;
    if (!(C > 0 && r > 0) || !Number.isInteger(N) || N < 1 || !(q >= 0 && q <= 1)) return { ok: false, error: 'Need C > 0, r > 0, N a whole number ≥ 1 and 0 ≤ p ≤ 1.' };
    const M = Math.floor(C / r + 1e-9), pmf = binomPmf(N, q);
    let over = 0; for (let n = N; n > M; n--) over += pmf[n];
    return { ok: true, M, pmf, over, mean: N * q, sd: Math.sqrt(N * q * (1 - q)), circuitUsers: M };
  }
  function nyquistShannon(B, M, snrDb) {
    const snr = Math.pow(10, snrDb / 10);
    return { snr, Rnyq: 2 * B * Math.log2(M), C: B * Math.log2(1 + snr), Mmax: Math.sqrt(1 + snr) };
  }

  /* ============================================================ 3. multiple access */
  const alohaPure = G => G * Math.exp(-2 * G), alohaSlotted = G => G * Math.exp(-G);
  // per-node slotted ALOHA, independent transmit probabilities p_i; k = slot index for "first success in slot k"
  function slottedNodes(ps, k) {
    const idle = ps.reduce((m, p) => m * (1 - p), 1);
    const nodes = ps.map((p, i) => {
      let others = 1; ps.forEach((q, j) => { if (j !== i) others *= 1 - q; });
      const s = p * others;
      return { p, s, firstAtK: Math.pow(1 - s, k - 1) * s, expSlots: s > 0 ? 1 / s : Infinity };
    });
    const eff = nodes.reduce((m, x) => m + x.s, 0);
    return { nodes, eff, idle, collision: Math.max(0, 1 - eff - idle), anyFirstAtK: Math.pow(1 - eff, k - 1) * eff };
  }
  const equalEff = (N, p) => N * p * Math.pow(1 - p, N - 1);
  const equalMax = N => (N === 1 ? 1 : Math.pow(1 - 1 / N, N - 1));
  function csmacd(p) { // {R, d, v, L}
    const tprop = p.d / p.v, ttrans = p.L / p.R, a = tprop / ttrans;
    return { tprop, ttrans, a, Lmin: 2 * p.R * tprop, slot: 2 * tprop, eff: 1 / (1 + 5 * a), effExact: 1 / (1 + 2 * Math.E * a) };
  }
  // Ethernet binary exponential backoff after the k-th collision: K uniform in {0..2^min(k,10) - 1}, wait K*512 bit times.
  function backoff(k, R, slotBits) {
    slotBits = slotBits || 512;
    if (!Number.isInteger(k) || k < 1) return { ok: false, error: 'Collision count k must be a whole number ≥ 1.' };
    if (k >= 16) return { ok: true, abort: true, k };
    const m = Math.min(k, 10), Kmax = Math.pow(2, m) - 1, slot = slotBits / R;
    return { ok: true, abort: false, k, m, Kmax, choices: Kmax + 1, slot, maxWait: Kmax * slot, meanWait: Kmax / 2 * slot };
  }

  /* ============================================================ 4. ARQ & TCP */
  // Stallings (Data & Computer Communications, ARQ performance): a = t_prop / t_frame, P = frame error probability
  function arqU(a, W, P) {
    const w = 2 * a + 1, big = W >= w;
    const sw = (1 - P) / (1 + 2 * a);
    const sr = big ? 1 - P : W * (1 - P) / w;
    const gbn = big ? (1 - P) / (1 + 2 * a * P) : W * (1 - P) / (w * (1 - P + W * P));
    return { sw, gbn, sr, big, w };
  }
  // Stepped ARQ simulation in units of t_frame. o: {proto:'sw'|'gbn'|'sr', n, W, a, lost (frame lost on its first send, -1 none),
  // timeout (measured from end of a frame's transmission; default 2a+1)}. ACKs have zero length.
  // SW/GBN: ACK k = "next frame expected is k" (cumulative), single timer on the oldest outstanding frame.
  // SR: ACK k acknowledges frame k only, one timer per frame, receiver buffers out-of-order frames.
  function arqSim(o) {
    const proto = o.proto, n = o.n, W = proto === 'sw' ? 1 : o.W, a = o.a, lost = o.lost;
    const TO = o.timeout !== undefined ? o.timeout : 2 * a + 1;
    const ev = []; let ek = 0;
    const push = (t, type, data, pri) => ev.push({ t, type, data, pri, k: ek++ });
    const sends = [], acks = [], timeouts = [], delivered = [];
    const firstTx = new Array(n).fill(true), acked = new Array(n).fill(false), retxQ = [], gen = {};
    let base = 0, next = 0, txFree = 0, now = 0, expected = 0, tEnd = NaN;
    const rbuf = new Array(n).fill(false);
    const startTimer = (key, t) => { gen[key] = (gen[key] || 0) + 1; push(t + TO, 'to', { key, g: gen[key] }, 2); };
    const stopTimer = key => { gen[key] = (gen[key] || 0) + 1; };
    function trySend() {
      if (txFree > now + 1e-9) return;
      let f = -1;
      if (proto === 'sr' && retxQ.length) f = retxQ.shift();
      else if (next < n && next < base + W) f = next++;
      if (f < 0) return;
      const t0 = now, t1 = now + 1, isLost = f === lost && firstTx[f], re = !firstTx[f]; firstTx[f] = false;
      sends.push({ f, t0, t1, lost: isLost, re });
      txFree = t1; push(t1, 'free', null, 3);
      if (!isLost) push(t1 + a, 'arr', { f }, 1);
      if (proto === 'sr') startTimer(f, t1); else if (f === base) startTimer('g', t1);
    }
    function sendAck(num, f) { acks.push({ n: num, f, t0: now, t1: now + a }); push(now + a, 'ack', { n: num }, 0); }
    trySend();
    let guard = 0;
    while (ev.length && guard++ < 20000) {
      let bi = 0;
      for (let i = 1; i < ev.length; i++) { const x = ev[i], b = ev[bi]; if (x.t < b.t - 1e-12 || (Math.abs(x.t - b.t) <= 1e-12 && (x.pri < b.pri || (x.pri === b.pri && x.k < b.k)))) bi = i; }
      const e = ev.splice(bi, 1)[0]; now = e.t;
      if (e.type === 'arr') {
        const f = e.data.f;
        if (proto === 'sr') {
          if (f >= expected && f < expected + W) rbuf[f] = true;
          while (expected < n && rbuf[expected]) { delivered.push({ f: expected, t: now }); expected++; }
          if (f < expected + W) sendAck(f, f);
        } else {
          if (f === expected) { delivered.push({ f, t: now }); expected++; }
          sendAck(expected, f);
        }
      } else if (e.type === 'ack') {
        const k = e.data.n;
        if (proto === 'sr') {
          if (!acked[k]) { acked[k] = true; stopTimer(k); const qi = retxQ.indexOf(k); if (qi >= 0) retxQ.splice(qi, 1); }
          while (base < n && acked[base]) base++;
        } else if (k > base) {
          base = k; if (next < base) next = base;
          if (base < next) startTimer('g', now); else stopTimer('g');
        }
        if (base >= n && !fin(tEnd)) tEnd = now;
      } else if (e.type === 'to') {
        const key = e.data.key;
        if (gen[key] !== e.data.g) continue;
        if (proto === 'sr') { if (!acked[key] && retxQ.indexOf(key) < 0) { retxQ.push(key); timeouts.push({ t: now, f: key }); } }
        else if (base < n) { timeouts.push({ t: now, f: base }); next = base; }
      }
      trySend();
      if (fin(tEnd)) break;
    }
    return { sends, acks, timeouts, delivered, tEnd, transmissions: sends.length, W, TO };
  }
  // TCP congestion window per RTT round. o: {rounds, ssthresh0, cwnd0 (default 1), losses: {round: 'to'|'3dup'},
  // variant:'tahoe'|'reno', plus3:bool}. Loss in round t takes effect in round t+1. Slow start: cwnd <- min(2cwnd, ssthresh);
  // congestion avoidance (cwnd >= ssthresh): cwnd + 1. Loss: ssthresh <- cwnd/2; timeout (and any loss in Tahoe): cwnd <- 1;
  // Reno triple-dup-ACK: cwnd <- ssthresh (or ssthresh + 3 with plus3, as in Kurose's Fig. 3.53).
  function tcpTrace(o) {
    const out = []; let c = o.cwnd0 || 1, s = o.ssthresh0;
    for (let t = 1; t <= o.rounds; t++) {
      const ev = o.losses[t], phase = c < s ? 'SS' : 'CA';
      out.push({ round: t, cwnd: c, ssthresh: s, phase, event: ev || '' });
      if (ev) {
        s = c / 2;
        if (ev === 'to' || o.variant === 'tahoe') c = 1; else c = s + (o.plus3 ? 3 : 0);
      } else if (c < s) c = Math.min(2 * c, s); else c = c + 1;
    }
    return out;
  }
  function parseLosses(str, maxRound) {
    const res = {}; if (typeof str !== 'string' || !str.trim()) return { ok: true, losses: res };
    const parts = str.split(/[,;\s]+/).filter(Boolean);
    for (const p of parts) {
      const m = /^(\d+)\s*[:=]?\s*(to|timeout|t|3dup|dup|3d|d)$/i.exec(p);
      if (!m) return { ok: false, error: 'Bad loss event "' + p + '". Use round:type, e.g. 8:3dup, 16:to.' };
      const r = parseInt(m[1], 10); if (r < 1 || r > maxRound) return { ok: false, error: 'Loss round ' + r + ' is outside 1..' + maxRound + '.' };
      res[r] = /^t/i.test(m[2]) ? 'to' : '3dup';
    }
    return { ok: true, losses: res };
  }

  /* ============================================================ 5. error detection */
  function parseHexWords(str) {
    if (typeof str !== 'string' || !str.trim()) return { ok: false, error: 'Enter 16-bit words in hex, e.g. 0001 f203 f4f5 f6f7.' };
    let toks = str.trim().replace(/0x/gi, '').split(/[\s,;]+/).filter(Boolean);
    if (toks.length === 1 && toks[0].length > 4) {  // one long hex string: split into 16-bit words, pad an odd tail with zeros
      const s = toks[0]; toks = []; for (let i = 0; i < s.length; i += 4) toks.push((s.slice(i, i + 4) + '000').slice(0, 4));
    }
    const words = [];
    for (const t of toks) {
      if (!/^[0-9a-fA-F]{1,4}$/.test(t)) return { ok: false, error: '"' + t + '" is not a 16-bit hex word (1–4 hex digits).' };
      words.push(parseInt(t, 16));
    }
    if (words.length > 64) return { ok: false, error: 'At most 64 words.' };
    return { ok: true, words };
  }
  // One's-complement sum with end-around carry; steps for the working
  function onesSum(words) {
    let sum = 0; const steps = [];
    words.forEach((w, i) => {
      if (i === 0) { sum = w; steps.push({ w, raw: w, carry: 0, res: w }); return; }
      const raw = sum + w, carry = raw > 0xFFFF ? 1 : 0, res = (raw & 0xFFFF) + carry;
      steps.push({ w, raw, carry, res }); sum = res;
    });
    return { sum, steps };
  }
  function inetChecksum(words) { const s = onesSum(words); return { sum: s.sum, steps: s.steps, checksum: (~s.sum) & 0xFFFF }; }
  function inetVerify(words) { const s = onesSum(words); return { sum: s.sum, steps: s.steps, ok: s.sum === 0xFFFF }; }

  // Generator as bits "10011" or polynomial "x^4+x+1"
  function parsePoly(str) {
    if (typeof str !== 'string' || !str.trim()) return { ok: false, error: 'Enter a generator, e.g. 10011 or x^4+x+1.' };
    const s = str.replace(/\s/g, '');
    let bits;
    if (/^[01]+$/.test(s)) bits = s.replace(/^0+/, '').split('').map(Number);
    else if (/^[x0-9^+]+$/i.test(s)) {
      const exps = [];
      for (const term of s.toLowerCase().split('+')) {
        let e;
        if (term === '1') e = 0; else if (term === 'x') e = 1; else { const m = /^x\^(\d{1,2})$/.exec(term); if (!m) return { ok: false, error: 'Cannot read polynomial term "' + term + '".' }; e = +m[1]; }
        exps.push(e);
      }
      const deg = Math.max.apply(null, exps); bits = new Array(deg + 1).fill(0);
      exps.forEach(e => { bits[deg - e] ^= 1; });
      while (bits.length && bits[0] === 0) bits.shift();
    } else return { ok: false, error: 'Generator must be bits (10011) or a polynomial (x^4+x+1).' };
    if (bits.length < 2) return { ok: false, error: 'Generator degree must be at least 1.' };
    if (bits.length > 33) return { ok: false, error: 'Generator degree at most 32.' };
    return { ok: true, bits };
  }
  // Mod-2 long division of dividend bits by g; steps record each XOR
  function mod2Div(dividend, g) {
    const r = g.length - 1, reg = dividend.slice(), q = [], steps = [];
    for (let i = 0; i + r < reg.length; i++) {
      if (reg[i]) { steps.push({ pos: i, before: reg.slice(i, i + g.length) }); for (let j = 0; j < g.length; j++) reg[i + j] ^= g[j]; q.push(1); } else q.push(0);
    }
    return { quotient: q, remainder: reg.slice(reg.length - r), reg, steps };
  }
  function crcEncode(msg, g) {
    const r = g.length - 1, d = mod2Div(msg.concat(new Array(r).fill(0)), g);
    return { r, dividend: msg.concat(new Array(r).fill(0)), quotient: d.quotient, remainder: d.remainder, steps: d.steps, frame: msg.concat(d.remainder) };
  }
  function crcCheck(frame, g) { const d = mod2Div(frame, g); return { remainder: d.remainder, steps: d.steps, quotient: d.quotient, ok: d.remainder.every(b => b === 0) }; }
  const parityBit = (bits, odd) => (bits.reduce((s, b) => s ^ b, 0) ^ (odd ? 1 : 0));
  // rows: array of equal-length bit arrays. Returns full (R+1)x(C+1) matrix with parity column/row.
  function parity2D(rows, odd) {
    const R = rows.length, C = rows[0].length, full = rows.map(r => r.concat([parityBit(r, odd)]));
    const last = []; for (let j = 0; j <= C; j++) last.push(parityBit(full.map(r => r[j]), odd));
    full.push(last); return { full, R, C };
  }
  function check2D(full, odd) {
    const want = odd ? 1 : 0, badRows = [], badCols = [];
    full.forEach((r, i) => { if (r.reduce((s, b) => s ^ b, 0) !== want) badRows.push(i); });
    for (let j = 0; j < full[0].length; j++) if (full.map(r => r[j]).reduce((s, b) => s ^ b, 0) !== want) badCols.push(j);
    return { badRows, badCols, ok: !badRows.length && !badCols.length, locatable: badRows.length === 1 && badCols.length === 1 };
  }

  FSP.math.net = {
    SCHEMES, parseBits, lineEncode, lineStats,
    nodal, switching, simulatePackets, binomPmf, statMux, nyquistShannon,
    alohaPure, alohaSlotted, slottedNodes, equalEff, equalMax, csmacd, backoff,
    arqU, arqSim, tcpTrace, parseLosses,
    parseHexWords, onesSum, inetChecksum, inetVerify, parsePoly, mod2Div, crcEncode, crcCheck, parityBit, parity2D, check2D,
  };

  /* ============================================================ tests */
  FSP.registerTests('net', function (t) {
    const M = FSP.math.net, sgn = a => a.map(v => (v > 0 ? '+' : v < 0 ? '-' : '0')).join('');
    const enc = (s, sc, o) => sgn(M.lineEncode(M.parseBits(s).bits, sc, o).half);
    // Manchester (course: 1 = high->low, 0 = low->high), bits 01001110, two half-bits per bit:
    // 0:-+ 1:+- 0:-+ 0:-+ 1:+- 1:+- 1:+- 0:-+
    t.check('Manchester (course) 01001110', enc('01001110', 'man') === '-++--+-++-+-+--+');
    // IEEE 802.3 is the mirror image: 1 = low->high
    t.check('Manchester (IEEE) 01001110 is the inverse', enc('01001110', 'man', { manchester: 'ieee' }) === '+--++-+--+-+-++-');
    // Diff. Manchester, prior level low (-). 0 = transition at start, 1 = none; always a mid-bit transition.
    // b0=0: start flips - -> + ; mid -> -  => +-   end -
    // b1=1: stays -, mid +                 => -+   end +
    // b2=0: flips to -, mid +              => -+   end +
    // b3=0: -+ ; b4=1: stays +, mid - => +- ; b5=1: stays -, => -+ ; b6=1: +- ; b7=0: flips to +, => +-
    t.check('Diff. Manchester 01001110, prior low', enc('01001110', 'dman', { prior: -1 }) === '+--+-+-++--++-+-');
    // NRZ-I, prior low: 0 1 0 0 1 1 1 0 -> levels - + + + - + - -
    t.check('NRZ-I 01001110, prior low', enc('01001110', 'nrzi', { prior: -1 }) === '--++++++--++----');
    // AMI, last pulse negative: 0 1 0 0 1 1 1 0 -> 0 + 0 0 - + - 0
    t.check('AMI 01001110', enc('01001110', 'ami', { lastPulse: -1 }) === '00++0000--++--00');
    // B8ZS 1 00000000 1, last pulse -: first 1 -> +; eight zeros with preceding + -> 000+-0-+ ; final 1 -> - (last pulse still +)
    const b8 = M.lineEncode(M.parseBits('1000000001').bits, 'b8zs', { lastPulse: -1 });
    t.check('B8ZS 1000000001 -> + 0 0 0 + - 0 - + -', sgn(b8.half.filter((_, i) => i % 2 === 0)) === '+000+-0-+-', sgn(b8.half));
    t.check('B8ZS marks V,B,V,B at bits 4,5,7,8', b8.marks.map(m => m.label + m.i).join() === 'V4,B5,V7,B8' && b8.subs.length === 1 && b8.subs[0].start === 1);
    // HDB3 1 0000 11 0000 0000, last pulse -, even count at start:
    // 1 -> + (count 1, odd) ; 0000 -> 000V with V=+ ; 1 -> - ; 1 -> + (count 2, even) ; 0000 -> B00V, B=-,V=- ;
    // count reset 0 (even) ; 0000 -> B00V, B=+, V=+
    const h3 = M.lineEncode(M.parseBits('100001100000000').bits, 'hdb3', { lastPulse: -1 });
    t.check('HDB3 100001100000000 -> +000+-+-00-+00+', sgn(h3.half.filter((_, i) => i % 2 === 0)) === '+000+-+-00-+00+', sgn(h3.half));
    const vs = h3.marks.filter(m => m.label === 'V').map(m => h3.half[2 * m.i]);
    t.check('HDB3 successive V pulses alternate in polarity', vs.length === 3 && vs[0] === -vs[1] && vs[1] === -vs[2]);
    t.check('HDB3 with odd starting count: 0000 -> 000V', sgn(M.lineEncode([0, 0, 0, 0], 'hdb3', { lastPulse: 1, hdb3Odd: true }).half) === '000000++');
    const st = M.lineStats(M.lineEncode(M.parseBits('01001110').bits, 'pnrz'));
    // polar NRZ 01001110: four 1s, four 0s -> DC 0; level changes at 0|1,1|0,0|1,1|0 = 4
    t.check('Polar NRZ stats: DC 0, 4 transitions, longest flat 3 bits', st.dc === 0 && st.transitions === 4 && st.longestFlat === 3);
    t.check('Bit parser rejects junk and empty', !M.parseBits('0120').ok && !M.parseBits('  ').ok && !M.parseBits('1'.repeat(65)).ok && M.parseBits('01 10').bits.length === 4);

    // Delays: 1000 B = 8000 bit at 200 Mb/s -> 8000/2e8 = 40 us ; 2 km / 2e8 m/s = 10 us
    const nd = M.nodal({ L: 8000, R: 200e6, d: 2000, v: 2e8, tproc: 0, tqueue: 0 });
    t.check('t_trans 40 us, t_prop 10 us', t.rel(nd.tTrans, 40e-6, 1e-12) && t.rel(nd.tProp, 10e-6, 1e-12));
    // Kurose P1.31: 8e6-bit message, 2 switches, 2 Mb/s links, no header/prop. Message switching: 3 x 4 s = 12 s;
    // 800 packets of 1e4 bits: t_pkt = 5 ms, total = (800 + 2) x 5 ms = 4.01 s
    const sw = M.switching({ M: 8e6, P: 1e4, h: 0, S: 2, R: 2e6, tprop: 0, tproc: 0 });
    t.check('Kurose message vs packet switching: 12 s vs 4.01 s', sw.N === 800 && t.rel(sw.totalMsg, 12, 1e-12) && t.rel(sw.total, 4.01, 1e-12), sw.total);
    const sp = { M: 10500, P: 1000, h: 40, S: 3, R: 1e6, tprop: 1e-4, tproc: 2e-5 }, sw2 = M.switching(sp);
    const sizes = []; for (let k = 0; k < sw2.N; k++) sizes.push(Math.min(sp.P, sp.M - k * sp.P) + sp.h);
    t.check('Closed form (short last packet) = event simulation', sw2.N === 11 && t.rel(sw2.total, M.simulatePackets(sizes, 3, 1e6, 1e-4, 2e-5).total, 1e-12));
    // P_opt = sqrt(M h / S): brute-force minimise (M/P + S)(P + h) over integer P for M = 1e5, h = 100, S = 4 -> sqrt(2.5e6) = 1581.1
    let bestP = 1, bestF = Infinity; for (let P = 1; P <= 20000; P++) { const f = (1e5 / P + 4) * (P + 100); if (f < bestF) { bestF = f; bestP = P; } }
    t.check('P_opt = sqrt(Mh/S) matches brute force', Math.abs(M.switching({ M: 1e5, P: 1000, h: 100, S: 4, R: 1, tprop: 0, tproc: 0 }).Popt - bestP) < 1, bestP);
    // Statistical multiplexing vs direct enumeration of all 2^10 activity patterns (N = 10, p = 0.3, M = 3)
    let en = 0; for (let m = 0; m < 1024; m++) { let k = 0; for (let b = 0; b < 10; b++) k += (m >> b) & 1; if (k > 3) en += Math.pow(0.3, k) * Math.pow(0.7, 10 - k); }
    const smx = M.statMux({ C: 3e5, r: 1e5, N: 10, p: 0.3 });
    t.check('Binomial overload = enumeration (N=10, p=0.3, M=3)', smx.M === 3 && t.rel(smx.over, en, 1e-12), smx.over + ' vs ' + en);
    // Kurose: 1 Mb/s link, 100 kb/s users, p = 0.1, 35 users -> P(>10 active) ~ 0.0004
    const k35 = M.statMux({ C: 1e6, r: 1e5, N: 35, p: 0.1 });
    t.check('Kurose 35 users: P(>10) ~ 0.0004', k35.M === 10 && k35.over > 3.5e-4 && k35.over < 4.8e-4, k35.over);
    t.check('pmf sums to 1 (N = 500)', t.near(M.binomPmf(500, 0.37).reduce((a, b) => a + b, 0), 1, 1e-12));
    t.check('Stat-mux rejects N = 0 / p > 1', !M.statMux({ C: 1, r: 1, N: 0, p: 0.1 }).ok && !M.statMux({ C: 1, r: 1, N: 3, p: 1.2 }).ok);
    // Nyquist 3 kHz, 4 levels -> 2*3000*2 = 12 kb/s ; SNR 30 dB -> 1001 -> C = 3000 log2(1001) = 29.9 kb/s ; Mmax = sqrt(1001) = 31.64
    const ns = M.nyquistShannon(3000, 4, 30);
    t.check('Nyquist 12 kb/s, Shannon 29.9 kb/s, Mmax 31.6', ns.Rnyq === 12000 && t.near(ns.C, 29901.7, 0.1) && t.near(ns.Mmax, 31.638, 1e-3), ns.C);

    // ALOHA maxima found by a fine grid search (not from the closed form)
    let pm = 0, pg = 0, sm = 0, sg = 0;
    for (let G = 0; G <= 3; G += 1e-4) { const a = M.alohaPure(G), b = M.alohaSlotted(G); if (a > pm) { pm = a; pg = G; } if (b > sm) { sm = b; sg = G; } }
    t.check('Pure ALOHA max 1/(2e) = 0.1839 at G = 0.5', t.near(pm, 1 / (2 * Math.E), 1e-8) && t.near(pg, 0.5, 1e-3));
    t.check('Slotted ALOHA max 1/e = 0.3679 at G = 1', t.near(sm, 1 / Math.E, 1e-8) && t.near(sg, 1, 1e-3));
    // 3 nodes p = 0.2, 0.3, 0.5: P_A = 0.2*0.7*0.5 = 0.07, P_B = 0.3*0.8*0.5 = 0.12, P_C = 0.5*0.8*0.7 = 0.28; eff 0.47; idle 0.8*0.7*0.5 = 0.28
    const sn = M.slottedNodes([0.2, 0.3, 0.5], 3);
    t.check('Per-node slotted ALOHA (0.2, 0.3, 0.5)', t.near(sn.nodes[0].s, 0.07, 1e-12) && t.near(sn.nodes[1].s, 0.12, 1e-12) && t.near(sn.nodes[2].s, 0.28, 1e-12) && t.near(sn.eff, 0.47, 1e-12) && t.near(sn.idle, 0.28, 1e-12) && t.near(sn.collision, 0.25, 1e-12));
    // first success of A in slot 3: 0.93^2 * 0.07 = 0.060543
    t.check('First success in slot k is geometric', t.near(sn.nodes[0].firstAtK, 0.060543, 1e-9));
    // equal p: N p (1-p)^(N-1) maximised at p = 1/N; N = 4: (3/4)^3 = 0.421875 ; limit 1/e
    let bp = 0, be = 0; for (let p = 0; p <= 1; p += 1e-5) { const e = M.equalEff(4, p); if (e > be) { be = e; bp = p; } }
    t.check('Equal-p optimum p = 1/N, (1-1/N)^(N-1) -> 1/e', t.near(bp, 0.25, 1e-4) && t.near(M.equalMax(4), 0.421875, 1e-12) && t.near(M.equalMax(1e6), 1 / Math.E, 1e-6));
    // CSMA/CD: 10 Mb/s, 2500 m, v = 2e8 -> t_prop 12.5 us, Lmin = 2*1e7*12.5e-6 = 250 bit
    const cd = M.csmacd({ R: 1e7, d: 2500, v: 2e8, L: 12000 });
    // a = 12.5 us / 1.2 ms = 0.0104167 -> 1/(1+0.0520833) = 0.950495
    t.check('CSMA/CD Lmin 250 bit, eff 1/(1+5a)', t.near(cd.Lmin, 250, 1e-9) && t.near(cd.eff, 0.950495, 1e-6));
    const bo = M.backoff(3, 1e7), bo12 = M.backoff(12, 1e7);
    t.check('Backoff k=3 -> K in 0..7 (8 choices), max 7*51.2 us; k=12 capped at 0..1023; k=16 aborts',
      bo.Kmax === 7 && t.rel(bo.maxWait, 7 * 51.2e-6, 1e-12) && bo12.Kmax === 1023 && M.backoff(16, 1e7).abort && !M.backoff(0, 1e7).ok);

    // ARQ (Stallings): W = 7, a = 10 -> 2a+1 = 21 > W ; P = 0.1
    // SW 0.9/21 = 0.0428571 ; SR 7*0.9/21 = 0.3 ; GBN 7*0.9/(21*(0.9+0.7)) = 6.3/33.6 = 0.1875
    const u1 = M.arqU(10, 7, 0.1);
    t.check('ARQ W < 2a+1: SW 0.04286, SR 0.3, GBN 0.1875', t.near(u1.sw, 0.0428571, 1e-6) && t.near(u1.sr, 0.3, 1e-12) && t.near(u1.gbn, 0.1875, 1e-12));
    // W = 127 >= 21: SR 0.9, GBN 0.9/(1+2) = 0.3
    const u2 = M.arqU(10, 127, 0.1);
    t.check('ARQ W >= 2a+1: SR 1-P, GBN (1-P)/(1+2aP)', t.near(u2.sr, 0.9, 1e-12) && t.near(u2.gbn, 0.3, 1e-12));
    const uA = M.arqU(3, 7, 0.05), uB = M.arqU(3, 6.999999, 0.05);
    t.check('GBN/SR formulas continuous at W = 2a+1', t.near(uA.gbn, uB.gbn, 1e-6) && t.near(uA.sr, uB.sr, 1e-6));
    // Simulation: no loss, a = 1. SW: each frame takes 1 + 2a = 3 -> 4 frames done at 12. GBN W=4 >= 2a+1: pipe full, last ACK at 8 + 2 = 10
    const sSW = M.arqSim({ proto: 'sw', n: 4, W: 1, a: 1, lost: -1 }), sG = M.arqSim({ proto: 'gbn', n: 8, W: 4, a: 1, lost: -1 });
    t.check('ARQ sim, no loss: SW 4 frames in 12 t_f, GBN 8 frames in 10 t_f', t.near(sSW.tEnd, 12, 1e-9) && t.near(sG.tEnd, 10, 1e-9) && sG.transmissions === 8);
    const lG = M.arqSim({ proto: 'gbn', n: 8, W: 4, a: 1, lost: 2 }), lS = M.arqSim({ proto: 'sr', n: 8, W: 4, a: 1, lost: 2 }), lW = M.arqSim({ proto: 'sw', n: 5, W: 1, a: 1, lost: 2 });
    const inOrder = s => s.delivered.map(d => d.f).join() === Array.from({ length: s.delivered.length }, (_, i) => i).join();
    t.check('Lost frame: SR resends 1 frame, GBN resends the window, SW resends 1; all deliver in order',
      lS.transmissions === 9 && lG.transmissions > 9 && lW.transmissions === 6 && inOrder(lS) && inOrder(lG) && lS.delivered.length === 8 && lG.delivered.length === 8, lG.transmissions);
    // SW with loss: frame 2 sent at 6, ends 7, timeout at 7 + (2a+1) = 10, resent 10..11, ACK at 13, frames 3, 4 -> 13 + 3 + 3 = 19
    t.check('SW loss timing: finish at 19 t_f', t.near(lW.tEnd, 19, 1e-9), lW.tEnd);
    // TCP (Kurose Fig. 3.53): ssthresh0 = 8, triple dup at round 8 (cwnd 12).
    // rounds 1-8: 1 2 4 8 9 10 11 12 ; Tahoe 9..: 1 2 4 6 7 ; Reno (+3): 9 10 11
    const tl = { 8: '3dup' }, ta = M.tcpTrace({ rounds: 13, ssthresh0: 8, losses: tl, variant: 'tahoe' }).map(r => r.cwnd).join();
    const re = M.tcpTrace({ rounds: 11, ssthresh0: 8, losses: tl, variant: 'reno', plus3: true }).map(r => r.cwnd).join();
    t.check('TCP Tahoe trace (Kurose 3.53)', ta === '1,2,4,8,9,10,11,12,1,2,4,6,7', ta);
    t.check('TCP Reno trace with +3 MSS', re === '1,2,4,8,9,10,11,12,9,10,11', re);
    t.check('Reno timeout -> cwnd 1; loss parser', M.tcpTrace({ rounds: 6, ssthresh0: 64, losses: { 5: 'to' }, variant: 'reno' })[5].cwnd === 1
      && M.parseLosses('8:3dup, 16:to', 30).losses[16] === 'to' && !M.parseLosses('8:x', 30).ok && !M.parseLosses('40:to', 30).ok);

    // RFC 1071 example: 0001 + f203 = f204; + f4f5 = 1e6f9 -> e6fa; + f6f7 = 1ddf1 -> ddf2; checksum ~ddf2 = 220d
    const ck = M.inetChecksum(M.parseHexWords('0001 f203 f4f5 f6f7').words);
    t.check('Internet checksum RFC 1071: sum ddf2, checksum 220d', ck.sum === 0xddf2 && ck.checksum === 0x220d && ck.steps[2].carry === 1);
    t.check('Receiver: data + checksum sums to ffff; a flipped bit fails', M.inetVerify([1, 0xf203, 0xf4f5, 0xf6f7, 0x220d]).ok && !M.inetVerify([1, 0xf203, 0xf4f4, 0xf6f7, 0x220d]).ok);
    t.check('Hex parser: continuous string, odd byte padded, junk rejected', M.parseHexWords('0001f203f4').words.join() === [1, 0xf203, 0xf400].join() && !M.parseHexWords('12g4').ok && !M.parseHexWords('').ok);
    // Tanenbaum: frame 1101011011, G = 10011 (x^4 + x + 1) -> remainder 1110, transmitted 11010110111110
    const g = M.parsePoly('x^4+x+1').bits, cr = M.crcEncode(M.parseBits('1101011011').bits, g);
    t.check('CRC Tanenbaum: x^4+x+1 = 10011, remainder 1110', g.join('') === '10011' && cr.remainder.join('') === '1110' && cr.frame.join('') === '11010110111110');
    // Forouzan: dataword 1001, divisor 1011 -> remainder 110 (by hand: 1001000 ^1011 -> 0010000 -> 0000110)
    t.check('CRC Forouzan: 1001 / 1011 -> 110', M.crcEncode([1, 0, 0, 1], [1, 0, 1, 1]).remainder.join('') === '110');
    // independent check with a shift-register (LFSR) CRC on integers for a random-ish message
    const lfsr = (msg, gen) => { const r = gen.length - 1, gi = parseInt(gen.join(''), 2); let reg = 0; msg.concat(new Array(r).fill(0)).forEach(b => { reg = (reg << 1) | b; if (reg >> r & 1) reg ^= gi; }); return reg; };
    const msgX = [1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1], gX = [1, 0, 0, 0, 0, 0, 1, 1, 1];
    t.check('CRC-8 long division = LFSR', parseInt(M.crcEncode(msgX, gX).remainder.join(''), 2) === lfsr(msgX, gX));
    const fr = cr.frame.slice(); fr[3] ^= 1;
    t.check('CRC check: clean frame -> 0, single error detected', M.crcCheck(cr.frame, g).ok && !M.crcCheck(fr, g).ok);
    t.check('Generator parser rejects degree 0 / junk', !M.parsePoly('1').ok && !M.parsePoly('x^a').ok && M.parsePoly('x^3+1').bits.join('') === '1001');
    // parity: 1011001 has four 1s -> even parity bit 0, odd parity bit 1
    t.check('Single parity even/odd', M.parityBit([1, 0, 1, 1, 0, 0, 1], false) === 0 && M.parityBit([1, 0, 1, 1, 0, 0, 1], true) === 1);
    // 2-D: rows 1011 / 0110 / 1100 -> row parities 1,0,0 ; column parities 0,0,0,1 ; corner 1
    const p2 = M.parity2D([[1, 0, 1, 1], [0, 1, 1, 0], [1, 1, 0, 0]], false);
    const flip = p2.full.map(r => r.slice()); flip[1][2] ^= 1; const c2 = M.check2D(flip, false);
    t.check('2-D parity: parity bits and single-error location', p2.full.map(r => r[4]).join('') === '1001' && p2.full[3].join('') === '00011' && M.check2D(p2.full, false).ok && c2.locatable && c2.badRows[0] === 1 && c2.badCols[0] === 2);
  });

  /* ============================================================ UI */
  const css = n => { try { return getComputedStyle(document.documentElement).getPropertyValue(n).trim() || null; } catch (e) { return null; } };
  function theme() {
    return { bg: css('--panel2') || '#0f141f', text: css('--text') || '#e6edf3', muted: css('--muted') || '#8b98ab', border: css('--border') || '#243047',
      a: css('--accent2') || '#f9a03f', b: css('--ok') || '#2ecc71', c: css('--accent') || '#4cc9f0', pink: css('--pink') || '#f72585', warn: css('--warn') || '#ffb347', bad: css('--bad') || '#ff5c5c', mono: css('--mono') || 'monospace' };
  }
  const E = (x, u, d) => FSP.fmtEng(x, u, d === undefined ? 4 : d);
  const N4 = (x, d) => (fin(x) ? String(+x.toPrecision(d || 4)) : '—');
  const pct = x => (fin(x) ? (100 * x).toFixed(2) + ' %' : '—');
  const tickStr = v => (Math.abs(v) < 1e-12 ? '0' : String(+v.toPrecision(4)));
  function niceTicks(lo, hi, n) {
    const span = hi - lo; if (!(span > 0)) return [lo];
    const raw = span / n, mag = Math.pow(10, Math.floor(Math.log10(raw))), f = raw / mag, step = (f < 1.5 ? 1 : f < 3 ? 2 : f < 7 ? 5 : 10) * mag, out = [];
    for (let v = Math.ceil(lo / step - 1e-9) * step; v <= hi + step * 1e-9; v += step) out.push(Math.abs(v) < step * 1e-9 ? 0 : +v.toPrecision(12));
    return out;
  }
  const logTicks = (lo, hi) => { const out = []; for (let e = Math.ceil(Math.log10(lo) - 1e-9); e <= Math.floor(Math.log10(hi) + 1e-9); e++) out.push(Math.pow(10, e)); return out; };
  const logStr = v => { const e = Math.round(Math.log10(v)); return e >= -2 && e <= 3 ? String(v) : '1e' + e; };

  // Generic axes + series. Returns mapping helpers so callers can draw extras.
  function plot(c, o) {
    const g = c.prep(), ctx = g.ctx, w = g.w, h = g.h, T = theme();
    const L = 54, R = 12, Tp = o.topPad || 12, B = 32, pw = w - L - R, ph = h - Tp - B;
    ctx.clearRect(0, 0, w, h); ctx.fillStyle = T.bg; ctx.fillRect(0, 0, w, h);
    if (pw < 20 || ph < 20) return null;
    const fx = o.xlog ? Math.log10 : v => v, fy = o.ylog ? Math.log10 : v => v;
    const x0 = fx(o.xmin), x1 = fx(o.xmax), y0 = fy(o.ymin), y1 = fy(o.ymax);
    const X = x => L + (fx(x) - x0) / (x1 - x0) * pw, Y = y => Tp + (1 - (fy(Math.max(y, o.ylog ? o.ymin * 1e-3 : -Infinity)) - y0) / (y1 - y0)) * ph;
    ctx.font = '11px ' + T.mono; ctx.lineWidth = 1;
    ctx.textBaseline = 'middle'; ctx.textAlign = 'right';
    (o.ylog ? logTicks(o.ymin, o.ymax) : niceTicks(o.ymin, o.ymax, 5)).forEach(v => { const yy = Y(v); ctx.strokeStyle = T.border; ctx.beginPath(); ctx.moveTo(L, yy); ctx.lineTo(L + pw, yy); ctx.stroke(); ctx.fillStyle = T.muted; ctx.fillText(o.ylog ? logStr(v) : tickStr(v), L - 4, yy); });
    ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    (o.xlog ? logTicks(o.xmin, o.xmax) : niceTicks(o.xmin, o.xmax, Math.max(3, Math.floor(pw / 70)))).forEach(v => { const xx = X(v); ctx.strokeStyle = T.border; ctx.beginPath(); ctx.moveTo(xx, Tp); ctx.lineTo(xx, Tp + ph); ctx.stroke(); ctx.fillStyle = T.muted; ctx.fillText(o.xlog ? logStr(v) : tickStr(v), xx, Tp + ph + 3); });
    ctx.fillStyle = T.muted; ctx.textAlign = 'center'; ctx.textBaseline = 'bottom'; ctx.fillText(o.xlabel || '', L + pw / 2, h - 1);
    ctx.save(); ctx.translate(10, Tp + ph / 2); ctx.rotate(-PI / 2); ctx.textBaseline = 'top'; ctx.fillText(o.ylabel || '', 0, 0); ctx.restore();
    ctx.save(); ctx.beginPath(); ctx.rect(L, Tp, pw, ph); ctx.clip();
    (o.vlines || []).forEach(v => { ctx.strokeStyle = v.color || T.muted; ctx.setLineDash([4, 4]); ctx.beginPath(); ctx.moveTo(X(v.x), Tp); ctx.lineTo(X(v.x), Tp + ph); ctx.stroke(); ctx.setLineDash([]); if (v.label) { ctx.fillStyle = v.color || T.muted; ctx.textAlign = 'left'; ctx.textBaseline = 'top'; ctx.fillText(v.label, X(v.x) + 3, Tp + 2 + (v.dy || 0)); } });
    (o.series || []).forEach(s => {
      ctx.strokeStyle = s.color; ctx.lineWidth = s.width || 2; ctx.setLineDash(s.dash || []); ctx.beginPath(); let pen = false;
      for (let i = 0; i < s.x.length; i++) { const xv = s.x[i], yv = s.y[i]; if (!fin(xv) || !fin(yv)) { pen = false; continue; } const px = X(xv), py = Y(yv); if (pen) ctx.lineTo(px, py); else { ctx.moveTo(px, py); pen = true; } }
      ctx.stroke(); ctx.setLineDash([]);
      if (s.dots) { ctx.fillStyle = s.color; for (let i = 0; i < s.x.length; i++) if (fin(s.y[i])) { ctx.beginPath(); ctx.arc(X(s.x[i]), Y(s.y[i]), 2.5, 0, 2 * PI); ctx.fill(); } }
    });
    if (o.extra) o.extra(ctx, X, Y, T, { L, Tp, pw, ph });
    ctx.restore();
    (o.points || []).forEach(p => {
      if (!fin(p.x) || !fin(p.y)) return; const px = X(p.x), py = Y(p.y);
      ctx.fillStyle = p.color; ctx.strokeStyle = T.bg; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.arc(px, py, 5, 0, 2 * PI); ctx.fill(); ctx.stroke();
      if (p.label) { const right = px < L + pw * 0.6; ctx.textAlign = right ? 'left' : 'right'; ctx.textBaseline = p.below ? 'top' : 'bottom'; ctx.fillStyle = T.text; ctx.fillText(p.label, px + (right ? 8 : -8), py + (p.below ? 6 : -6)); }
    });
    (o.legend || []).forEach((lg, i) => { ctx.font = '11px ' + T.mono; ctx.fillStyle = lg.color; ctx.textAlign = 'right'; ctx.textBaseline = 'top'; ctx.fillText(lg.text, L + pw - 4, Tp + 4 + i * 14); });
    return { ctx, X, Y, T, L, Tp, pw, ph, w, h };
  }
  function blank(c, text) { const g = c.prep(), T = theme(); g.ctx.clearRect(0, 0, g.w, g.h); g.ctx.fillStyle = T.bg; g.ctx.fillRect(0, 0, g.w, g.h); g.ctx.fillStyle = T.muted; g.ctx.font = '12px ' + T.mono; g.ctx.textAlign = 'center'; g.ctx.textBaseline = 'middle'; g.ctx.fillText(text || 'Fix the inputs', g.w / 2, g.h / 2); }

  /* ---------- small DOM helpers ---------- */
  const el = function () { return FSP.ui.el.apply(null, arguments); };
  function layout(root) { const lay = el('div', { class: 'layout' }), ctl = el('div', { class: 'controls' }), stage = el('div', { class: 'stage' }); lay.appendChild(ctl); lay.appendChild(stage); root.appendChild(lay); return { ctl, stage }; }
  function textIn(parent, label, value, onChange, attrs) {
    const inp = el('input', Object.assign({ type: 'text', 'aria-label': label, spellcheck: 'false', size: 14 }, attrs || {})); inp.value = value; inp.style.width = '100%';
    inp.addEventListener('input', () => { onChange(inp.value); FSP.state.touch(); });
    parent.appendChild(el('div', { class: 'ctl' }, el('label', { text: label }), inp, el('span'))); return inp;
  }
  function msgBox(parent, cls) { const m = el('div', { class: 'msg ' + (cls || 'bad'), hidden: '', role: 'alert' }); parent.appendChild(m); return { el: m, show(t) { m.hidden = false; m.textContent = t; }, hide() { m.hidden = true; } }; }
  function sliders(parent, defs, store, onChange) {
    defs.forEach(d => { store[d.k] = FSP.ui.slider(parent, { label: d.l, min: d.min, max: d.max, step: d.step, value: d.v, unit: d.u, log: d.log, digits: d.digits, onInput: () => { onChange(); FSP.state.touch(); } }); });
  }
  function hud(parent) { const h = el('div', { class: 'hud' }); parent.appendChild(h); return h; }
  function hudSet(h, rows) { while (h.firstChild) h.removeChild(h.firstChild); rows.forEach(r => h.appendChild(el('div', null, el('span', { text: r[0] }), el('span', { text: r[1] })))); }
  function canvasIn(parent, label, height) { const wrap = el('div', { class: 'canvas-wrap' }); parent.appendChild(wrap); const c = FSP.ui.canvas(wrap, { height }); c.cv.setAttribute('role', 'img'); c.cv.setAttribute('aria-label', label); return c; }
  function table(parent) { const wrap = el('div', { class: 'np-scroll' }), tb = el('table', { class: 'np-table' }); wrap.appendChild(tb); parent.appendChild(wrap); return tb; }
  function tableSet(tb, head, rows) {
    tb.innerHTML = ''; const th = el('tr'); head.forEach(hh => th.appendChild(el('th', { text: hh }))); tb.appendChild(th);
    rows.forEach(r => { const tr = el('tr'); r.forEach(c => tr.appendChild(el('td', { text: c }))); tb.appendChild(tr); });
  }
  const sv = (st, k) => st[k].get();
  function setSliders(st, o, keys) { keys.forEach(k => { if (o[k] !== undefined && st[k]) { const v = parseFloat(o[k]); if (fin(v)) st[k].set(v, true); } }); }
  const okStr = (v, max) => typeof v === 'string' && v.length <= (max || 200);

  /* ---------- 1. line coding panel ---------- */
  function buildLine(root) {
    const M = FSP.math.net, { ctl, stage } = layout(root), S = { bits: '01001110', man: 'course', prior: '-1', last: '-1', hodd: '0' }, ui = {};
    const f1 = FSP.ui.fieldset(ctl, 'Data');
    ui.bits = textIn(f1, 'Bits', S.bits, v => { S.bits = v; update(); }, { inputmode: 'numeric' });
    ui.msg = msgBox(f1);
    const f2 = FSP.ui.fieldset(ctl, 'Conventions');
    ui.man = FSP.ui.select(f2, 'Manchester', [['course', 'Course: 1 = high→low'], ['ieee', 'IEEE 802.3: 1 = low→high']], S.man, v => { S.man = v; update(); FSP.state.touch(); });
    ui.prior = FSP.ui.select(f2, 'Level before bit 1', [['-1', 'Low (−V)'], ['1', 'High (+V)']], S.prior, v => { S.prior = v; update(); FSP.state.touch(); });
    ui.last = FSP.ui.select(f2, 'Last pulse (AMI)', [['-1', '−V (first 1 → +V)'], ['1', '+V (first 1 → −V)']], S.last, v => { S.last = v; update(); FSP.state.touch(); });
    ui.hodd = FSP.ui.select(f2, 'HDB3 pulses since last sub.', [['0', 'Even (0)'], ['1', 'Odd']], S.hodd, v => { S.hodd = v; update(); FSP.state.touch(); });
    f2.appendChild(el('div', { class: 'note', text: 'Level before bit 1 is used by NRZ-I and differential Manchester (they encode relative to the previous level). "Last pulse" seeds the AMI alternation for AMI/B8ZS/HDB3.' }));
    f2.appendChild(el('div', { class: 'note', text: 'Polar NRZ here: 1 = +V, 0 = −V (Stallings calls the inverted form NRZ-L: 0 = high). Diff. Manchester: always a mid-bit transition; transition at the start of the bit = 0, none = 1.' }));
    ui.cv = canvasIn(stage, 'Line code waveforms', 440); ui.cv.onResize(() => draw());
    ui.tb = table(stage);
    ui.work = FSP.ui.working(stage);
    let last = null;
    function opts() { return { manchester: S.man, prior: +S.prior, lastPulse: +S.last, hdb3Odd: S.hodd === '1' }; }
    function update() {
      const p = M.parseBits(S.bits, 64);
      if (!p.ok) { ui.msg.show(p.error); last = null; blank(ui.cv, 'Enter a valid bit string'); ui.tb.innerHTML = ''; ui.work.set('Fix the input above.'); return; }
      ui.msg.hide();
      const o = opts(), encs = M.SCHEMES.map(s => ({ id: s[0], name: s[1], e: M.lineEncode(p.bits, s[0], o) }));
      last = { bits: p.bits, encs };
      tableSet(ui.tb, ['Scheme', 'DC (×V)', 'Transitions', 'Longest flat (bits)'], encs.map(x => { const s = M.lineStats(x.e); return [x.name, s.dc.toFixed(3), String(s.transitions), String(s.longestFlat)]; }));
      const sym = v => (v > 0 ? '+' : v < 0 ? '−' : '0'), pairs = e => p.bits.map((_, i) => sym(e.half[2 * i]) + sym(e.half[2 * i + 1])).join(' ');
      const rule = {
        unrz: '1 → +V, 0 → 0 for the whole bit',
        pnrz: '1 → +V, 0 → −V for the whole bit',
        nrzi: '1 → invert level at start of bit, 0 → keep level (prior level ' + sym(o.prior) + 'V)',
        ami: '0 → 0 V, each 1 → pulse opposite to the previous 1 (last pulse before data ' + sym(o.lastPulse) + 'V)',
        man: S.man === 'ieee' ? 'IEEE 802.3: 1 → low→high at mid-bit, 0 → high→low' : 'Course: 1 → high→low at mid-bit, 0 → low→high',
        dman: 'mid-bit transition always; 0 → extra transition at bit start, 1 → none (prior level ' + sym(o.prior) + 'V)',
        b8zs: 'AMI, but each run of 8 zeros → 000VB0VB (V = same polarity as preceding pulse, B = AMI-correct pulse)',
        hdb3: 'AMI, but each run of 4 zeros → 000V if an odd number of pulses since the last substitution, else B00V (start parity: ' + (o.hdb3Odd ? 'odd' : 'even') + ')',
      };
      const lines = ['Bits:         ' + p.bits.join('  '), 'Each bit shown as (1st half, 2nd half) level; + = +V, − = −V, 0 = 0 V.', ''];
      encs.forEach(x => {
        lines.push(x.name + ': ' + rule[x.id]); lines.push('   ' + pairs(x.e));
        x.e.subs.forEach(s => lines.push('   substitution ' + s.pattern + ' at bits ' + (s.start + 1) + '–' + (s.start + s.len)));
        const s = M.lineStats(x.e); lines.push('   DC = mean of half-bit levels = ' + s.dc.toFixed(3) + ' V/V; transitions = ' + s.transitions + (x.e.prior !== undefined ? ' (incl. bit-1 start vs prior level)' : ''));
        lines.push('');
      });
      ui.work.set(lines);
      draw();
    }
    function draw() {
      if (!last) return;
      const g = ui.cv.prep(), ctx = g.ctx, w = g.w, h = g.h, T = theme(), n = last.bits.length;
      ctx.clearRect(0, 0, w, h); ctx.fillStyle = T.bg; ctx.fillRect(0, 0, w, h);
      const narrow = w < 520, LW = narrow ? 74 : 118, R = 10, top = 24, rows = last.encs.length, rh = (h - top - 6) / rows, pw = w - LW - R, bw = pw / n;
      const X = t => LW + t * bw;
      ctx.font = (narrow ? '10px ' : '11px ') + T.mono;
      // bit grid + labels
      for (let i = 0; i <= n; i++) { ctx.strokeStyle = T.border; ctx.setLineDash([3, 3]); ctx.beginPath(); ctx.moveTo(X(i), top - 4); ctx.lineTo(X(i), h - 4); ctx.stroke(); }
      ctx.setLineDash([]); ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
      if (bw >= 7) last.bits.forEach((b, i) => { ctx.fillStyle = T.text; ctx.fillText(String(b), X(i + 0.5), 11); });
      const colors = [T.c, T.c, T.a, T.b, T.pink, T.pink, T.warn, T.warn];
      last.encs.forEach((x, r) => {
        const y0 = top + r * rh, mid = y0 + rh / 2, amp = rh * 0.32, uni = x.id === 'unrz';
        const lv = v => (uni ? y0 + rh * 0.8 - v * rh * 0.6 : mid - v * amp);
        ctx.fillStyle = T.muted; ctx.textAlign = 'left'; ctx.textBaseline = 'middle'; ctx.fillText(x.name, 4, mid);
        ctx.strokeStyle = T.border; ctx.beginPath(); ctx.moveTo(LW, lv(0)); ctx.lineTo(LW + pw, lv(0)); ctx.stroke();
        x.e.subs.forEach(s => { ctx.fillStyle = T.warn; ctx.globalAlpha = 0.12; ctx.fillRect(X(s.start), y0 + 2, s.len * bw, rh - 4); ctx.globalAlpha = 1; });
        ctx.strokeStyle = colors[r]; ctx.lineWidth = 2; ctx.beginPath();
        const hv = x.e.half;
        if (x.e.prior !== undefined) { ctx.moveTo(LW - 8, lv(x.e.prior)); ctx.lineTo(LW, lv(x.e.prior)); ctx.lineTo(LW, lv(hv[0])); } else ctx.moveTo(LW, lv(hv[0]));
        for (let k = 0; k < hv.length; k++) { if (k > 0 && hv[k] !== hv[k - 1]) ctx.lineTo(X(k / 2), lv(hv[k])); ctx.lineTo(X((k + 1) / 2), lv(hv[k])); }
        ctx.stroke(); ctx.lineWidth = 1;
        x.e.marks.forEach(m => { const v = hv[2 * m.i]; ctx.fillStyle = m.label === 'V' ? T.bad : T.b; ctx.textAlign = 'center'; ctx.textBaseline = v > 0 ? 'bottom' : 'top'; ctx.font = 'bold 11px ' + T.mono; ctx.fillText(m.label, X(m.i + 0.5), lv(v) + (v > 0 ? -1 : 2)); ctx.font = (narrow ? '10px ' : '11px ') + T.mono; });
      });
    }
    return {
      update, draw,
      get: () => ({ bits: S.bits, man: S.man, prior: S.prior, last: S.last, hodd: S.hodd }),
      set(o) {
        if (okStr(o.bits, 80)) S.bits = o.bits;
        if (o.man === 'course' || o.man === 'ieee') S.man = o.man;
        ['prior', 'last'].forEach(k => { if (o[k] === '1' || o[k] === '-1') S[k] = o[k]; });
        if (o.hodd === '0' || o.hodd === '1') S.hodd = o.hodd;
        ui.bits.value = S.bits; ui.man.value = S.man; ui.prior.value = S.prior; ui.last.value = S.last; ui.hodd.value = S.hodd;
      },
    };
  }

  /* ---------- 2. delays & multiplexing panel ---------- */
  function buildDelay(root) {
    const M = FSP.math.net, { ctl, stage } = layout(root), st = {}, ui = {};
    const f1 = FSP.ui.fieldset(ctl, 'Link (also each hop below)');
    sliders(f1, [
      { k: 'L', l: 'L (bits)', min: 8, max: 1e9, v: 8000, log: true, u: 'bit' },
      { k: 'R', l: 'R', min: 1e3, max: 1e11, v: 200e6, log: true, u: 'b/s' },
      { k: 'd', l: 'd', min: 1, max: 4e7, v: 2000, log: true, u: 'm' },
      { k: 'v', l: 'v', min: 1e7, max: 3e8, v: 2e8, u: 'm/s' },
      { k: 'tproc', l: 't_proc', min: 0, max: 0.01, v: 0, u: 's' },
      { k: 'tq', l: 't_queue', min: 0, max: 1, v: 0, u: 's' }], st, update);
    f1.appendChild(el('div', { class: 'note', text: 'SI suffixes work in the boxes: 200M, 2k, 20u. 1 byte = 8 bits.' }));
    const f2 = FSP.ui.fieldset(ctl, 'Message vs packet switching');
    sliders(f2, [
      { k: 'M', l: 'Message M', min: 100, max: 1e9, v: 8e6, log: true, u: 'bit' },
      { k: 'P', l: 'Payload P', min: 8, max: 1e8, v: 1e4, log: true, u: 'bit' },
      { k: 'h', l: 'Header h', min: 0, max: 2000, v: 0, step: 1, u: 'bit' },
      { k: 'S', l: 'Switches S', min: 0, max: 10, v: 2, step: 1, u: '' }], st, update);
    ui.msgSw = msgBox(f2);
    const f3 = FSP.ui.fieldset(ctl, 'Statistical multiplexing');
    sliders(f3, [
      { k: 'C', l: 'Link C', min: 1e3, max: 1e10, v: 1e6, log: true, u: 'b/s' },
      { k: 'r', l: 'User rate r', min: 1e2, max: 1e9, v: 1e5, log: true, u: 'b/s' },
      { k: 'N', l: 'Users N', min: 1, max: 1000, v: 35, step: 1, u: '' },
      { k: 'p', l: 'P(active)', min: 0, max: 1, v: 0.1, step: 0.001, u: '' }], st, update);
    ui.ylog = FSP.ui.select(f3, 'y scale', [['lin', 'Linear'], ['log', 'Log']], 'lin', () => { update(); FSP.state.touch(); });
    ui.msgMx = msgBox(f3);
    const f4 = FSP.ui.fieldset(ctl, 'Nyquist & Shannon');
    sliders(f4, [
      { k: 'B', l: 'Bandwidth B', min: 100, max: 1e10, v: 3000, log: true, u: 'Hz' },
      { k: 'lev', l: 'Levels M', min: 2, max: 1024, v: 4, step: 1, u: '' },
      { k: 'snr', l: 'SNR', min: -10, max: 60, v: 30, step: 0.1, u: 'dB' }], st, update);
    ui.cvT = canvasIn(stage, 'Space-time diagram of pipelined packets', 320); ui.cvT.onResize(() => update());
    ui.hudD = hud(stage);
    ui.cvB = canvasIn(stage, 'Binomial distribution of active users', 220); ui.cvB.onResize(() => update());
    ui.hudM = hud(stage);
    ui.work = FSP.ui.working(stage);

    function update() {
      const p = { L: sv(st, 'L'), R: sv(st, 'R'), d: sv(st, 'd'), v: sv(st, 'v'), tproc: sv(st, 'tproc'), tqueue: sv(st, 'tq') };
      const nd = M.nodal(p), W = [];
      W.push('— Single link —',
        't_trans = L/R = ' + N4(p.L, 6) + ' bit / ' + E(p.R, 'b/s') + ' = ' + E(nd.tTrans, 's'),
        't_prop  = d/v = ' + E(p.d, 'm') + ' / ' + N4(p.v, 4) + ' m/s = ' + E(nd.tProp, 's'),
        'd_nodal = t_proc + t_queue + t_trans + t_prop = ' + E(p.tproc, 's') + ' + ' + E(p.tqueue, 's') + ' + ' + E(nd.tTrans, 's') + ' + ' + E(nd.tProp, 's') + ' = ' + E(nd.total, 's'),
        'a = t_prop/t_trans = ' + N4(nd.a) + ';  bandwidth-delay product R·t_prop = ' + N4(nd.bdp) + ' bit (bits "in flight" on the link)', '');
      const rows = [['t_trans = L/R', E(nd.tTrans, 's')], ['t_prop = d/v', E(nd.tProp, 's')], ['Nodal delay', E(nd.total, 's')], ['R·t_prop', N4(nd.bdp) + ' bit']];
      // switching
      const sp = { M: sv(st, 'M'), P: sv(st, 'P'), h: sv(st, 'h'), S: Math.round(sv(st, 'S')), R: p.R, tprop: nd.tProp, tproc: p.tproc }, sw = M.switching(sp);
      if (!sw.ok || sw.N > 1e7) { ui.msgSw.show(sw.ok ? 'Too many packets (N = ' + sw.N + '); increase P.' : sw.error); blank(ui.cvT); }
      else {
        ui.msgSw.hide();
        const lk = sp.S + 1, tp = sp.tprop + sp.tproc;
        rows.push(['Packets N = ⌈M/P⌉', String(sw.N) + (sw.equal ? '' : ' (last ' + N4(sw.Plast, 6) + ' bit)')], ['Packet switching', E(sw.total, 's')], ['Message switching', E(sw.totalMsg, 's')], ['Speed-up', N4(sw.totalMsg / sw.total) + '×'], ['P_opt = √(Mh/S)', sp.S > 0 && sp.h > 0 ? N4(sw.Popt, 5) + ' bit' : (sp.h === 0 ? '→ small as possible (h = 0)' : '∞ (S = 0: one packet)')]);
        W.push('— Store-and-forward over S = ' + sp.S + ' switches (' + lk + ' links, each identical to the link above; queueing ignored) —',
          'N = ⌈M/P⌉ = ⌈' + N4(sp.M, 8) + ' / ' + N4(sp.P, 8) + '⌉ = ' + sw.N + (sw.equal ? '' : '  (last packet carries only ' + N4(sw.Plast, 8) + ' bit)'),
          't_pkt = (P + h)/R = (' + N4(sp.P, 8) + ' + ' + sp.h + ') / ' + E(sp.R, 'b/s') + ' = ' + E(sw.tPkt, 's'),
          'Derivation: packet 1 needs S+1 hops, each t_pkt + t_prop + t_proc  → arrives at (S+1)(t_pkt + t_prop + t_proc).',
          '            Links pipeline, so the other N−1 packets arrive one t_pkt apart behind it:',
          '  T_pkt = (S+1)(t_pkt + t_prop + t_proc) + (N−1)·t_pkt = (N+S)·(P+h)/R + (S+1)·t_prop + (S+1)·t_proc',
          '        = (' + sw.N + ' + ' + sp.S + ')·' + E(sw.tPkt, 's') + ' + ' + lk + '·' + E(sp.tprop, 's') + ' + ' + lk + '·' + E(sp.tproc, 's') + ' = ' + E(sw.totalEq, 's'));
        if (!sw.equal) W.push('  Short last packet: T = (N+S−1)·t_pkt + (P_last+h)/R + (S+1)(t_prop+t_proc) = ' + E(sw.total, 's'));
        W.push('Message switching (whole message + one header per hop):',
          '  T_msg = (S+1)·[(M+h)/R + t_prop + t_proc] = ' + lk + '·[' + E((sp.M + sp.h) / sp.R, 's') + ' + ' + E(sp.tprop, 's') + ' + ' + E(sp.tproc, 's') + '] = ' + E(sw.totalMsg, 's'),
          'Optimal payload: with N ≈ M/P, transmission part ∝ (M/P + S)(P + h) = M + Mh/P + SP + Sh.',
          '  d/dP: −Mh/P² + S = 0  →  P_opt = √(M·h/S)' + (sp.S > 0 && sp.h > 0 ? ' = √(' + N4(sp.M, 6) + '·' + sp.h + '/' + sp.S + ') = ' + N4(sw.Popt, 5) + ' bit' : (sp.S === 0 ? '  (S = 0: no pipelining gain, send one packet)' : '  (h = 0: smaller packets always help)')), '');
        drawTiming(sp, sw);
      }
      hudSet(ui.hudD, rows);
      // stat mux
      const mx = { C: sv(st, 'C'), r: sv(st, 'r'), N: Math.round(sv(st, 'N')), p: sv(st, 'p') }, sm = M.statMux(mx);
      if (!sm.ok) { ui.msgMx.show(sm.error); blank(ui.cvB); hudSet(ui.hudM, []); }
      else {
        ui.msgMx.hide();
        hudSet(ui.hudM, [['M = ⌊C/r⌋ (circuit-switched users)', String(sm.M)], ['P(more than M active)', sm.over < 1e-4 ? sm.over.toExponential(3) : N4(sm.over, 5)], ['Mean active Np', N4(sm.mean)], ['Packet-switched users supported', String(mx.N) + ' (vs ' + sm.M + ' circuits)']]);
        W.push('— Statistical multiplexing —', 'M = ⌊C/r⌋ = ⌊' + E(mx.C, 'b/s') + ' / ' + E(mx.r, 'b/s') + '⌋ = ' + sm.M + ' simultaneous users (this is all circuit switching can admit)',
          'n active ~ Binomial(N = ' + mx.N + ', p = ' + mx.p + '):  P(n) = C(N,n) pⁿ (1−p)^(N−n)',
          'P(overload) = Σ_{n=M+1}^{N} C(N,n) pⁿ (1−p)^(N−n) = Σ_{n=' + (sm.M + 1) + '}^{' + mx.N + '} … = ' + (sm.over < 1e-4 ? sm.over.toExponential(4) : N4(sm.over, 6)));
        const terms = []; for (let n = sm.M + 1; n <= Math.min(mx.N, sm.M + 4); n++) terms.push('P(' + n + ') = ' + sm.pmf[n].toExponential(3));
        if (terms.length) W.push('  first terms: ' + terms.join(', ') + (mx.N > sm.M + 4 ? ', …' : ''));
        W.push('Mean Np = ' + N4(sm.mean) + ', σ = √(Np(1−p)) = ' + N4(sm.sd), '');
        drawBinom(mx, sm);
      }
      // Nyquist/Shannon
      const ns = M.nyquistShannon(sv(st, 'B'), Math.round(sv(st, 'lev')), sv(st, 'snr'));
      hudSet(ui.hudM, [].concat(Array.from(ui.hudM.children).map(d => [d.children[0].textContent, d.children[1].textContent]),
        [['Nyquist R = 2B log₂M', E(ns.Rnyq, 'b/s')], ['Shannon C = B log₂(1+SNR)', E(ns.C, 'b/s')], ['Max useful levels √(1+SNR)', N4(ns.Mmax)], ['Achievable rate min(R, C)', E(Math.min(ns.Rnyq, ns.C), 'b/s')]]));
      W.push('— Nyquist & Shannon —', 'Nyquist (noiseless): R = 2B log₂M = 2·' + E(sv(st, 'B'), 'Hz') + '·log₂' + Math.round(sv(st, 'lev')) + ' = ' + E(ns.Rnyq, 'b/s'),
        'SNR = 10^(' + sv(st, 'snr') + '/10) = ' + N4(ns.snr, 6),
        'Shannon: C = B log₂(1 + SNR) = ' + E(sv(st, 'B'), 'Hz') + '·log₂(' + N4(1 + ns.snr, 6) + ') = ' + E(ns.C, 'b/s'),
        'Setting 2B log₂M = B log₂(1+SNR) → M_max = √(1+SNR) = ' + N4(ns.Mmax) + ' levels; more levels than this cannot raise the rate above C.');
      ui.work.set(W);
    }
    function drawTiming(sp, sw) {
      const maxN = 40, n = Math.min(sw.N, maxN), sizes = [];
      for (let k = 0; k < n; k++) sizes.push(Math.min(sp.P, sp.M - k * sp.P) + sp.h);
      const sim = M.simulatePackets(sizes, sp.S, sp.R, sp.tprop, sp.tproc), links = sp.S + 1;
      const g = ui.cvT.prep(), ctx = g.ctx, w = g.w, h = g.h, T = theme();
      ctx.clearRect(0, 0, w, h); ctx.fillStyle = T.bg; ctx.fillRect(0, 0, w, h);
      const L = 56, R = 16, top = 26, bot = h - 10, tEnd = Math.max(n === sw.N ? sw.total : sim.total, sw.totalMsg), xs = j => L + j * (w - L - R) / links, Y = t => top + t / tEnd * (bot - top);
      ctx.font = '11px ' + T.mono; ctx.textAlign = 'center'; ctx.textBaseline = 'bottom';
      for (let j = 0; j <= links; j++) { ctx.strokeStyle = T.muted; ctx.beginPath(); ctx.moveTo(xs(j), top); ctx.lineTo(xs(j), bot); ctx.stroke(); ctx.fillStyle = T.text; ctx.fillText(j === 0 ? 'Src' : j === links ? 'Dst' : 'S' + j, xs(j), top - 6); }
      ctx.textAlign = 'right'; ctx.textBaseline = 'middle'; ctx.fillStyle = T.muted;
      niceTicks(0, tEnd, 6).forEach(t => { ctx.fillText(E(t, 's', 3), L - 4, Y(t)); ctx.strokeStyle = T.border; ctx.beginPath(); ctx.moveTo(L, Y(t)); ctx.lineTo(w - R, Y(t)); ctx.stroke(); });
      const cols = [T.c, T.a, T.b, T.pink];
      sim.hops.forEach((row, i) => row.forEach((hp, j) => {
        ctx.fillStyle = cols[i % 4]; ctx.globalAlpha = 0.55; ctx.beginPath();
        ctx.moveTo(xs(j), Y(hp.start)); ctx.lineTo(xs(j), Y(hp.end)); ctx.lineTo(xs(j + 1), Y(hp.end + sp.tprop)); ctx.lineTo(xs(j + 1), Y(hp.start + sp.tprop)); ctx.closePath(); ctx.fill(); ctx.globalAlpha = 1;
      }));
      const yP = Y(sw.total), yM = Y(sw.totalMsg);
      ctx.setLineDash([5, 4]); ctx.strokeStyle = T.b; ctx.beginPath(); ctx.moveTo(L, yP); ctx.lineTo(w - R, yP); ctx.stroke(); ctx.strokeStyle = T.bad; ctx.beginPath(); ctx.moveTo(L, yM); ctx.lineTo(w - R, yM); ctx.stroke(); ctx.setLineDash([]);
      ctx.textAlign = 'left'; ctx.textBaseline = 'bottom'; ctx.fillStyle = T.b; ctx.fillText('packet sw. ' + E(sw.total, 's', 3), L + 4, yP - 2);
      ctx.fillStyle = T.bad; ctx.fillText('message sw. ' + E(sw.totalMsg, 's', 3), L + 4, Math.max(top + 12, yM - 2));
      if (n < sw.N) { ctx.fillStyle = T.warn; ctx.textAlign = 'right'; ctx.textBaseline = 'top'; ctx.fillText('first ' + maxN + ' of ' + sw.N + ' packets drawn', w - R - 2, top + 2); }
    }
    function drawBinom(mx, sm) {
      const lg = ui.ylog.value === 'log', N = mx.N, ymax = Math.max.apply(null, sm.pmf) * 1.1 || 1;
      const ymin = lg ? Math.max(1e-15, Math.min(1e-3, Math.pow(10, Math.floor(Math.log10(Math.max(1e-15, sm.over || 1e-15)))))) : 0;
      plot(ui.cvB, {
        xmin: -0.5, xmax: N + 0.5, ymin: lg ? ymin : 0, ymax: lg ? 1 : ymax, ylog: lg, xlabel: 'number of active users n', ylabel: 'P(n)',
        vlines: [{ x: sm.M + 0.5, label: 'M = ' + sm.M, color: theme().bad }],
        extra(ctx, X, Y, T) {
          const bw = Math.max(1, (X(1) - X(0)) * 0.8);
          sm.pmf.forEach((v, n) => { if (!(v > 0)) return; ctx.fillStyle = n > sm.M ? T.bad : T.c; const y = Y(Math.max(v, lg ? ymin : 0)); ctx.fillRect(X(n) - bw / 2, y, bw, Y(lg ? ymin : 0) - y); });
        },
        legend: [{ text: 'P(n > M) = ' + (sm.over < 1e-4 ? sm.over.toExponential(2) : N4(sm.over)), color: theme().bad }],
      });
    }
    const keys = ['L', 'R', 'd', 'v', 'tproc', 'tq', 'M', 'P', 'h', 'S', 'C', 'r', 'N', 'p', 'B', 'lev', 'snr'];
    return {
      update,
      get() { const o = {}; keys.forEach(k => { o[k] = sv(st, k); }); o.ylog = ui.ylog.value; return o; },
      set(o) { setSliders(st, o, keys); if (o.ylog === 'lin' || o.ylog === 'log') ui.ylog.value = o.ylog; },
    };
  }

  /* ---------- 3. multiple access panel ---------- */
  function buildMAC(root) {
    const M = FSP.math.net, { ctl, stage } = layout(root), st = {}, ui = {}, S = { ps: '0.2, 0.3, 0.5' };
    const f1 = FSP.ui.fieldset(ctl, 'ALOHA');
    sliders(f1, [{ k: 'G', l: 'Offered load G', min: 0, max: 3, v: 0.5, step: 0.01, u: 'fr/T' }], st, update);
    const f2 = FSP.ui.fieldset(ctl, 'Slotted ALOHA, per node');
    ui.ps = textIn(f2, 'p_i list', S.ps, v => { S.ps = v; update(); }, { inputmode: 'decimal' });
    sliders(f2, [{ k: 'k', l: 'Slot k', min: 1, max: 50, v: 3, step: 1, u: '' }], st, update);
    ui.msgP = msgBox(f2);
    f2.appendChild(el('div', { class: 'note', text: 'Comma-separated transmit probabilities, one per node (A, B, C …).' }));
    const f3 = FSP.ui.fieldset(ctl, 'Slotted ALOHA, equal p');
    sliders(f3, [{ k: 'Ne', l: 'Nodes N', min: 1, max: 100, v: 4, step: 1, u: '' }, { k: 'pe', l: 'p', min: 0, max: 1, v: 0.25, step: 0.001, u: '' }], st, update);
    const f4 = FSP.ui.fieldset(ctl, 'CSMA/CD (Ethernet)');
    sliders(f4, [
      { k: 'cR', l: 'R', min: 1e6, max: 1e11, v: 1e7, log: true, u: 'b/s' },
      { k: 'cd', l: 'd (max)', min: 1, max: 1e5, v: 2500, log: true, u: 'm' },
      { k: 'cv', l: 'v', min: 1e7, max: 3e8, v: 2e8, u: 'm/s' },
      { k: 'cL', l: 'Frame L', min: 64, max: 1e6, v: 12000, log: true, u: 'bit' },
      { k: 'ck', l: 'Collisions k', min: 1, max: 16, v: 3, step: 1, u: '' }], st, update);
    ui.cvA = canvasIn(stage, 'ALOHA throughput versus offered load', 240); ui.cvA.onResize(() => update());
    ui.hudA = hud(stage);
    ui.tb = table(stage);
    ui.cvE = canvasIn(stage, 'Slotted ALOHA efficiency versus p for equal nodes', 200); ui.cvE.onResize(() => update());
    ui.hudC = hud(stage);
    ui.work = FSP.ui.working(stage);
    function update() {
      const G = sv(st, 'G'), sP = M.alohaPure(G), sS = M.alohaSlotted(G), T = theme(), W = [];
      const gx = [], yp = [], ys = []; for (let i = 0; i <= 300; i++) { const g = 3 * i / 300; gx.push(g); yp.push(M.alohaPure(g)); ys.push(M.alohaSlotted(g)); }
      plot(ui.cvA, { xmin: 0, xmax: 3, ymin: 0, ymax: 0.42, xlabel: 'G (attempts per frame time)', ylabel: 'S (throughput)',
        series: [{ x: gx, y: yp, color: T.a }, { x: gx, y: ys, color: T.c }],
        points: [{ x: 0.5, y: 1 / (2 * Math.E), color: T.a, label: '1/2e', below: true }, { x: 1, y: 1 / Math.E, color: T.c, label: '1/e' }, { x: G, y: sP, color: T.text }, { x: G, y: sS, color: T.text }],
        vlines: [{ x: G }], legend: [{ text: 'slotted S = G e^−G', color: T.c }, { text: 'pure S = G e^−2G', color: T.a }] });
      hudSet(ui.hudA, [['Pure S = G e^−2G', N4(sP, 5)], ['Slotted S = G e^−G', N4(sS, 5)], ['Pure max 1/(2e) at G = 0.5', N4(1 / (2 * Math.E), 5)], ['Slotted max 1/e at G = 1', N4(1 / Math.E, 5)]]);
      W.push('— ALOHA (Poisson attempts, rate G per frame time) —',
        'Pure: a frame succeeds if no other starts within its 2-frame vulnerable period → P = e^(−2G), S = G e^(−2G) = ' + N4(G) + '·e^(−' + N4(2 * G) + ') = ' + N4(sP, 5),
        'Slotted: vulnerable period 1 slot → S = G e^(−G) = ' + N4(sS, 5),
        'dS/dG = 0: pure (1 − 2G)e^(−2G) = 0 → G = 1/2, S_max = 1/(2e) = 0.1839; slotted (1 − G)e^(−G) = 0 → G = 1, S_max = 1/e = 0.3679', '');
      // per node
      const toks = S.ps.split(/[,;\s]+/).filter(Boolean), ps = toks.map(Number), k = Math.round(sv(st, 'k'));
      if (!toks.length || toks.length > 26 || ps.some(p => !fin(p) || p < 0 || p > 1)) {
        ui.msgP.show(!toks.length ? 'Enter at least one probability.' : toks.length > 26 ? 'At most 26 nodes.' : 'Each p_i must be a number between 0 and 1.'); ui.tb.innerHTML = '';
      } else {
        ui.msgP.hide();
        const sn = M.slottedNodes(ps, k), nm = i => String.fromCharCode(65 + i);
        tableSet(ui.tb, ['Node', 'p_i', 'P(success in a slot)', 'P(1st success in slot ' + k + ')', 'E[slots]'],
          sn.nodes.map((x, i) => [nm(i), N4(x.p), N4(x.s, 5), N4(x.firstAtK, 5), fin(x.expSlots) ? N4(x.expSlots) : '∞']).concat([['Any', '', N4(sn.eff, 5), N4(sn.anyFirstAtK, 5), sn.eff > 0 ? N4(1 / sn.eff) : '∞']]));
        W.push('— Slotted ALOHA, node probabilities p = [' + ps.join(', ') + '] —');
        sn.nodes.forEach((x, i) => {
          const others = ps.map((q, j) => (j === i ? null : '(1−' + q + ')')).filter(Boolean).join('·');
          W.push('P(' + nm(i) + ' succeeds) = p_' + nm(i) + ' Π_{j≠' + nm(i) + '}(1−p_j) = ' + x.p + (others ? '·' + others : '') + ' = ' + N4(x.s, 6));
          W.push('   first success in slot ' + k + ' = (1 − ' + N4(x.s, 5) + ')^' + (k - 1) + ' · ' + N4(x.s, 5) + ' = ' + N4(x.firstAtK, 6) + '  (geometric; mean 1/P = ' + (fin(x.expSlots) ? N4(x.expSlots) : '∞') + ' slots)');
        });
        W.push('Efficiency = P(exactly one transmits) = Σ_i P(i succeeds) = ' + N4(sn.eff, 6), 'P(idle) = Π(1−p_i) = ' + N4(sn.idle, 6) + ';  P(collision) = 1 − efficiency − idle = ' + N4(sn.collision, 6), '');
      }
      // equal p
      const Ne = Math.round(sv(st, 'Ne')), pe = sv(st, 'pe'), ee = M.equalEff(Ne, pe), em = M.equalMax(Ne);
      const px = [], py = []; for (let i = 0; i <= 400; i++) { px.push(i / 400); py.push(M.equalEff(Ne, i / 400)); }
      plot(ui.cvE, { xmin: 0, xmax: 1, ymin: 0, ymax: 1.05, xlabel: 'p (each of N nodes)', ylabel: 'efficiency', series: [{ x: px, y: py, color: T.b }],
        points: [{ x: 1 / Ne, y: em, color: T.a, label: 'p = 1/N: ' + N4(em, 4) }, { x: pe, y: ee, color: T.text, below: true }],
        vlines: [{ x: 1 / Ne }], extra(ctx, X, Y, TT) { ctx.strokeStyle = TT.pink; ctx.setLineDash([3, 3]); ctx.beginPath(); ctx.moveTo(X(0), Y(1 / Math.E)); ctx.lineTo(X(1), Y(1 / Math.E)); ctx.stroke(); ctx.setLineDash([]); ctx.fillStyle = TT.pink; ctx.textAlign = 'right'; ctx.textBaseline = 'bottom'; ctx.fillText('1/e', X(1) - 4, Y(1 / Math.E) - 2); } });
      W.push('— Slotted ALOHA, N = ' + Ne + ' equal nodes —', 'E(p) = N p (1−p)^(N−1) = ' + Ne + '·' + pe + '·(1−' + pe + ')^' + (Ne - 1) + ' = ' + N4(ee, 6),
        'dE/dp = N(1−p)^(N−2)[(1−p) − (N−1)p] = 0 → p* = 1/N = ' + N4(1 / Ne) + ', E_max = (1 − 1/N)^(N−1) = ' + N4(em, 6) + '  → 1/e = 0.3679 as N → ∞', '');
      // CSMA/CD
      const cp = { R: sv(st, 'cR'), d: sv(st, 'cd'), v: sv(st, 'cv'), L: sv(st, 'cL') }, cd = M.csmacd(cp), bo = M.backoff(Math.round(sv(st, 'ck')), cp.R);
      const rows = [['t_prop = d/v', E(cd.tprop, 's')], ['L_min = 2R·t_prop', N4(cd.Lmin, 6) + ' bit (' + N4(cd.Lmin / 8, 5) + ' B)'], ['Frame ' + N4(cp.L, 6) + ' bit', cp.L >= cd.Lmin ? 'OK (≥ L_min)' : 'TOO SHORT: collision may go undetected'],
        ['a = t_prop/t_frame', N4(cd.a)], ['Efficiency 1/(1+5a)', pct(cd.eff)]];
      if (bo.abort) rows.push(['After 16 attempts', 'frame discarded']); else rows.push(['Backoff after k = ' + bo.k, 'K ∈ {0 … ' + bo.Kmax + '} → wait 0 … ' + E(bo.maxWait, 's')]);
      hudSet(ui.hudC, rows);
      W.push('— CSMA/CD —', 't_prop = d/v = ' + E(cp.d, 'm') + ' / ' + N4(cp.v) + ' m/s = ' + E(cd.tprop, 's'),
        'A sender must still be transmitting when the collision news returns (round trip 2·t_prop):', '  L_min = R · 2t_prop = ' + E(cp.R, 'b/s') + ' · ' + E(2 * cd.tprop, 's') + ' = ' + N4(cd.Lmin, 6) + ' bit  (slot time = 2 t_prop = ' + E(cd.slot, 's') + ')',
        't_frame = L/R = ' + E(cd.ttrans, 's') + ';  a = t_prop/t_frame = ' + N4(cd.a, 5),
        'Efficiency ≈ 1/(1 + 5a) = 1/(1 + ' + N4(5 * cd.a, 5) + ') = ' + N4(cd.eff, 5) + '  (Kurose approximation; the Metcalfe–Boggs form 1/(1 + 2e·a) gives ' + N4(cd.effExact, 5) + ')');
      if (bo.abort) W.push('Binary exponential backoff: after 16 collisions the frame is discarded.');
      else W.push('Binary exponential backoff after collision k = ' + bo.k + ': m = min(k, 10) = ' + bo.m + ', K uniform in {0, …, 2^' + bo.m + ' − 1 = ' + bo.Kmax + '} (' + bo.choices + ' choices)',
        '  wait = K × 512 bit times = K × ' + E(bo.slot, 's') + ' → 0 … ' + E(bo.maxWait, 's') + ' (mean ' + E(bo.meanWait, 's') + ')');
      ui.work.set(W);
    }
    const keys = ['G', 'k', 'Ne', 'pe', 'cR', 'cd', 'cv', 'cL', 'ck'];
    return {
      update,
      get() { const o = {}; keys.forEach(k => { o[k] = sv(st, k); }); o.ps = S.ps; return o; },
      set(o) { setSliders(st, o, keys); if (okStr(o.ps, 300)) { S.ps = o.ps; ui.ps.value = S.ps; } },
    };
  }

  /* ---------- 4. ARQ & TCP panel ---------- */
  function buildARQ(root, host) {
    const M = FSP.math.net, { ctl, stage } = layout(root), st = {}, ui = {}, S = { proto: 'gbn', losses: '8:3dup, 16:to, 24:3dup', reno: 'half' };
    const f1 = FSP.ui.fieldset(ctl, 'Link & window (efficiency)');
    sliders(f1, [
      { k: 'R', l: 'R', min: 1e3, max: 1e10, v: 1e6, log: true, u: 'b/s' },
      { k: 'L', l: 'Frame L', min: 8, max: 1e6, v: 1000, log: true, u: 'bit' },
      { k: 'd', l: 'd', min: 1, max: 4e7, v: 2e6, log: true, u: 'm' },
      { k: 'v', l: 'v', min: 1e7, max: 3e8, v: 2e8, u: 'm/s' },
      { k: 'W', l: 'Window W', min: 1, max: 1024, v: 7, step: 1, u: 'fr' },
      { k: 'P', l: 'Frame error P', min: 0, max: 0.99, v: 0.01, step: 0.0001, u: '' },
      { k: 'kb', l: 'Seq bits k', min: 1, max: 16, v: 3, step: 1, u: 'bit' }], st, update);
    const f2 = FSP.ui.fieldset(ctl, 'Sequence diagram');
    ui.proto = FSP.ui.select(f2, 'Protocol', [['sw', 'Stop-and-wait'], ['gbn', 'Go-Back-N'], ['sr', 'Selective Repeat']], S.proto, v => { S.proto = v; resim(); FSP.state.touch(); });
    sliders(f2, [
      { k: 'dn', l: 'Frames', min: 1, max: 12, v: 8, step: 1, u: '' },
      { k: 'dW', l: 'Window', min: 1, max: 12, v: 4, step: 1, u: 'fr' },
      { k: 'da', l: 'a (diagram)', min: 0.1, max: 4, v: 1, step: 0.05, u: '' },
      { k: 'dl', l: 'Lost frame', min: -1, max: 11, v: 2, step: 1, u: '#' }], st, resim);
    const row = el('div', { class: 'row' }); f2.appendChild(row);
    ui.play = FSP.ui.button(row, '▶ Play', () => togglePlay()); FSP.ui.button(row, 'Step ›', () => step()); FSP.ui.button(row, 'Show all', () => { stop(); cursor = 1; drawSeq(); });
    [ui.play].concat(Array.from(row.children)).forEach(b => b.classList.add('btn'));
    f2.appendChild(el('div', { class: 'note', text: 'Lost frame = −1 for none (lost on its first transmission only). Time unit = t_frame; ACKs have zero length; timeout = 2a + 1 after a frame finishes. SW/GBN ACK n = "expecting n"; SR ACK n acknowledges frame n.' }));
    const f3 = FSP.ui.fieldset(ctl, 'TCP congestion control');
    sliders(f3, [{ k: 'ss0', l: 'ssthresh₀', min: 1, max: 128, v: 16, step: 1, u: 'MSS' }, { k: 'rounds', l: 'RTT rounds', min: 5, max: 80, v: 30, step: 1, u: '' }], st, update);
    ui.losses = textIn(f3, 'Losses', S.losses, v => { S.losses = v; update(); });
    ui.reno = FSP.ui.select(f3, 'Reno on 3 dup', [['half', 'cwnd = ssthresh'], ['plus3', 'cwnd = ssthresh + 3']], S.reno, v => { S.reno = v; update(); FSP.state.touch(); });
    ui.msgT = msgBox(f3);
    f3.appendChild(el('div', { class: 'note', text: 'Losses: round:type, e.g. "8:3dup, 16:to" (3dup = triple duplicate ACK, to = timeout). Tahoe treats both as cwnd → 1.' }));
    ui.cvU = canvasIn(stage, 'ARQ link utilisation versus a', 240); ui.cvU.onResize(() => update());
    ui.hudU = hud(stage);
    ui.cvS = canvasIn(stage, 'ARQ sequence diagram', 460); ui.cvS.onResize(() => drawSeq());
    ui.hudS = hud(stage);
    ui.cvT = canvasIn(stage, 'TCP congestion window versus RTT round', 260); ui.cvT.onResize(() => update());
    ui.work = FSP.ui.working(stage);
    let sim = null, cursor = 1, raf = 0, lastTs = 0, Wtxt = [], Wseq = [];
    function update() {
      const p = { R: sv(st, 'R'), L: sv(st, 'L'), d: sv(st, 'd'), v: sv(st, 'v'), W: Math.round(sv(st, 'W')), P: sv(st, 'P'), kb: Math.round(sv(st, 'kb')) };
      const tf = p.L / p.R, tp = p.d / p.v, a = tp / tf, u = M.arqU(a, p.W, p.P), T = theme();
      const gbnMax = Math.pow(2, p.kb) - 1, srMax = Math.pow(2, p.kb - 1);
      const ax = [], s1 = [], s2 = [], s3 = []; for (let i = 0; i <= 200; i++) { const aa = Math.pow(10, -2 + 5 * i / 200); ax.push(aa); const uu = M.arqU(aa, p.W, p.P); s1.push(uu.sw); s2.push(uu.gbn); s3.push(uu.sr); }
      plot(ui.cvU, { xmin: 0.01, xmax: 1000, xlog: true, ymin: 0, ymax: 1.05, xlabel: 'a = t_prop / t_frame', ylabel: 'utilisation U',
        series: [{ x: ax, y: s1, color: T.a }, { x: ax, y: s2, color: T.c }, { x: ax, y: s3, color: T.b, dash: [6, 3] }], vlines: [{ x: Math.min(1000, Math.max(0.01, a)) }, { x: (p.W - 1) / 2 > 0.01 ? (p.W - 1) / 2 : 0.01, label: 'W = 2a+1', color: T.pink, dy: 14 }],
        points: [{ x: a, y: u.sw, color: T.a }, { x: a, y: u.gbn, color: T.c }, { x: a, y: u.sr, color: T.b }],
        legend: [{ text: 'Selective Repeat', color: T.b }, { text: 'Go-Back-N', color: T.c }, { text: 'Stop-and-wait', color: T.a }] });
      const warn = [];
      if (p.W > gbnMax) warn.push('GBN needs W ≤ 2^k − 1 = ' + gbnMax);
      if (p.W > srMax) warn.push('SR needs W ≤ 2^(k−1) = ' + srMax);
      hudSet(ui.hudU, [['t_frame = L/R', E(tf, 's')], ['t_prop = d/v', E(tp, 's')], ['a', N4(a)], ['2a + 1', N4(2 * a + 1)],
        ['Stop-and-wait U', pct(u.sw) + '  (' + E(u.sw * p.R, 'b/s') + ')'], ['Go-Back-N U', pct(u.gbn) + '  (' + E(u.gbn * p.R, 'b/s') + ')'], ['Selective Repeat U', pct(u.sr) + '  (' + E(u.sr * p.R, 'b/s') + ')'],
        ['Seq. numbers (k = ' + p.kb + ')', warn.length ? '⚠ ' + warn.join('; ') : 'W = ' + p.W + ' OK for GBN (≤' + gbnMax + ') and SR (≤' + srMax + ')']]);
      Wtxt = ['— ARQ efficiency (Stallings, Data & Computer Communications; a = t_prop/t_frame, P = frame error probability) —',
        't_frame = L/R = ' + E(tf, 's') + ',  t_prop = d/v = ' + E(tp, 's') + ',  a = ' + N4(a, 6) + ',  2a+1 = ' + N4(2 * a + 1, 6),
        'Stop-and-wait: U = (1−P)/(1+2a) = ' + N4(1 - p.P) + '/' + N4(1 + 2 * a) + ' = ' + N4(u.sw, 5),
        'W = ' + p.W + (u.big ? ' ≥ 2a+1: the window never closes (sender transmits continuously)' : ' < 2a+1: the sender stalls waiting for ACKs'),
        u.big ? 'Selective Repeat: U = 1 − P = ' + N4(u.sr, 5) : 'Selective Repeat: U = W(1−P)/(2a+1) = ' + p.W + '·' + N4(1 - p.P) + '/' + N4(2 * a + 1) + ' = ' + N4(u.sr, 5),
        u.big ? 'Go-Back-N: U = (1−P)/(1+2aP) = ' + N4(1 - p.P) + '/' + N4(1 + 2 * a * p.P) + ' = ' + N4(u.gbn, 5)
          : 'Go-Back-N: U = W(1−P)/[(2a+1)(1−P+WP)] = ' + N4(p.W * (1 - p.P)) + '/[' + N4(2 * a + 1) + '·' + N4(1 - p.P + p.W * p.P) + '] = ' + N4(u.gbn, 5),
        'Throughput = U·R. Window limits with k-bit sequence numbers: GBN W ≤ 2^k − 1 = ' + gbnMax + ', SR W ≤ 2^(k−1) = ' + srMax + '.', ''];
      // TCP
      const rounds = Math.round(sv(st, 'rounds')), ss0 = Math.round(sv(st, 'ss0')), pl = M.parseLosses(S.losses, rounds), Wt = [];
      if (!pl.ok) { ui.msgT.show(pl.error); blank(ui.cvT, 'Fix the loss list'); }
      else {
        ui.msgT.hide();
        const ta = M.tcpTrace({ rounds, ssthresh0: ss0, losses: pl.losses, variant: 'tahoe' }), re = M.tcpTrace({ rounds, ssthresh0: ss0, losses: pl.losses, variant: 'reno', plus3: S.reno === 'plus3' });
        const xs = ta.map(r => r.round), ymax = Math.max.apply(null, ta.concat(re).map(r => Math.max(r.cwnd, Math.min(r.ssthresh, 4 * ss0)))) * 1.1;
        const stepXY = tr => { const x = [], y = []; tr.forEach(r => { x.push(r.round - 0.5, r.round + 0.5); y.push(r.ssthresh, r.ssthresh); }); return { x, y }; };
        const sR = stepXY(re), sT = stepXY(ta);
        plot(ui.cvT, { xmin: 0.5, xmax: rounds + 0.5, ymin: 0, ymax: Math.max(4, ymax), xlabel: 'transmission round (RTT)', ylabel: 'cwnd (MSS)',
          series: [{ x: sT.x, y: sT.y, color: T.a, width: 1, dash: [2, 3] }, { x: sR.x, y: sR.y, color: T.c, width: 1, dash: [6, 3] }, { x: xs, y: ta.map(r => r.cwnd), color: T.a, dots: true }, { x: xs, y: re.map(r => r.cwnd), color: T.c, dots: true }],
          vlines: Object.keys(pl.losses).map(r => ({ x: +r, label: pl.losses[r] === 'to' ? 'TO' : '3dup', color: T.bad })),
          legend: [{ text: 'Reno cwnd (dashed: ssthresh)', color: T.c }, { text: 'Tahoe cwnd (dotted: ssthresh)', color: T.a }] });
        const f = x => String(+x.toFixed(2));
        Wt.push('— TCP congestion control (per RTT round; loss in round t changes round t+1) —',
          'Slow start (cwnd < ssthresh): cwnd ← min(2·cwnd, ssthresh). Congestion avoidance (cwnd ≥ ssthresh): cwnd ← cwnd + 1 MSS.',
          'Any loss: ssthresh ← cwnd/2. Timeout: cwnd ← 1 (both). Triple dup ACK: Tahoe cwnd ← 1; Reno fast recovery cwnd ← ssthresh' + (S.reno === 'plus3' ? ' + 3 MSS (Kurose Fig. 3.53)' : ' (the +3 MSS inflation is not shown at RTT granularity)') + '.',
          'round | Tahoe cwnd ssthresh phase | Reno cwnd ssthresh phase | event');
        ta.forEach((r, i) => { const q = re[i]; Wt.push(String(r.round).padStart(5) + ' | ' + f(r.cwnd).padStart(10) + ' ' + f(r.ssthresh).padStart(8) + ' ' + r.phase.padStart(5) + ' | ' + f(q.cwnd).padStart(9) + ' ' + f(q.ssthresh).padStart(8) + ' ' + q.phase.padStart(5) + ' | ' + (r.event === 'to' ? 'timeout' : r.event === '3dup' ? 'triple dup ACK' : '')); });
      }
      Wtxt = Wtxt.concat(Wt);
      resim(true);
    }
    function resim(fromUpdate) {
      const n = Math.round(sv(st, 'dn')), lost = Math.round(sv(st, 'dl'));
      sim = M.arqSim({ proto: S.proto, n, W: Math.round(sv(st, 'dW')), a: sv(st, 'da'), lost: lost < n ? lost : -1 });
      sim.a = sv(st, 'da'); sim.n = n;
      const sends = sim.sends.length, retx = sim.sends.filter(s => s.re).length;
      hudSet(ui.hudS, [['Protocol', { sw: 'Stop-and-wait', gbn: 'Go-Back-N (W = ' + sim.W + ')', sr: 'Selective Repeat (W = ' + sim.W + ')' }[S.proto]], ['Frames delivered', sim.delivered.length + ' / ' + n],
        ['Transmissions', sends + ' (' + retx + ' retransmitted)'], ['Done at', fin(sim.tEnd) ? N4(sim.tEnd) + ' t_frame' : '—'], ['Useful fraction n/(t_end)', fin(sim.tEnd) ? pct(n / sim.tEnd) : '—']]);
      Wseq = ['', '— Sequence diagram events (time in t_frame, a = ' + N4(sim.a) + ', timeout ' + N4(sim.TO) + ' after send ends) —'];
      const evs = [];
      sim.sends.forEach(s => evs.push([s.t0, 'send F' + s.f + (s.re ? ' (retransmission)' : '') + (s.lost ? ' — LOST' : '')]));
      sim.acks.forEach(k => evs.push([k.t0, 'receiver gets F' + k.f + ' → ACK' + k.n]));
      sim.timeouts.forEach(o => evs.push([o.t, 'timeout for F' + o.f + (S.proto === 'sr' ? ' → resend F' + o.f : ' → go back, resend from F' + o.f)]));
      evs.sort((x, y) => x[0] - y[0]).forEach(e => Wseq.push('t = ' + N4(e[0]).padStart(6) + ': ' + e[1]));
      ui.work.set(Wtxt.concat(Wseq));
      if (!fromUpdate) { cursor = 1; stop(); FSP.state.touch(); }
      drawSeq();
    }
    function drawSeq() {
      if (!sim) return;
      const g = ui.cvS.prep(), ctx = g.ctx, w = g.w, h = g.h, T = theme(), a = sim.a;
      ctx.clearRect(0, 0, w, h); ctx.fillStyle = T.bg; ctx.fillRect(0, 0, w, h);
      const xs = Math.min(110, w * 0.24), xr = w - Math.min(110, w * 0.24), top = 28, bot = h - 10;
      const tEnd = (fin(sim.tEnd) ? sim.tEnd : Math.max.apply(null, sim.sends.map(s => s.t1 + a).concat([1]))) + 0.3, Y = t => top + t / tEnd * (bot - top), tc = cursor * tEnd;
      ctx.font = '11px ' + T.mono; ctx.fillStyle = T.text; ctx.textAlign = 'center'; ctx.textBaseline = 'bottom';
      ctx.fillText('Sender', xs, top - 8); ctx.fillText('Receiver', xr, top - 8);
      ctx.strokeStyle = T.muted; ctx.lineWidth = 2; ctx.beginPath(); ctx.moveTo(xs, top); ctx.lineTo(xs, bot); ctx.moveTo(xr, top); ctx.lineTo(xr, bot); ctx.stroke(); ctx.lineWidth = 1;
      ctx.textAlign = 'right'; ctx.textBaseline = 'middle'; ctx.fillStyle = T.muted;
      niceTicks(0, tEnd, 8).forEach(t => { ctx.fillText(tickStr(t), xs - 34, Y(t)); ctx.strokeStyle = T.border; ctx.beginPath(); ctx.moveTo(xs - 30, Y(t)); ctx.lineTo(xs, Y(t)); ctx.stroke(); });
      ctx.save(); ctx.beginPath(); ctx.rect(0, top - 1, w, Math.max(0, Y(tc) - top + 1)); ctx.clip();
      sim.sends.forEach(s => {
        const frac = s.lost ? 0.55 : 1, xe = xs + (xr - xs) * frac, col = s.lost ? T.bad : s.re ? T.warn : T.c;
        ctx.fillStyle = col; ctx.globalAlpha = 0.35; ctx.beginPath(); ctx.moveTo(xs, Y(s.t0)); ctx.lineTo(xs, Y(s.t1)); ctx.lineTo(xe, Y(s.t1 + a * frac)); ctx.lineTo(xe, Y(s.t0 + a * frac)); ctx.closePath(); ctx.fill();
        ctx.globalAlpha = 1; ctx.strokeStyle = col; ctx.stroke();
        ctx.fillStyle = col; ctx.textAlign = 'right'; ctx.textBaseline = 'middle'; ctx.fillText('F' + s.f, xs - 4, Y((s.t0 + s.t1) / 2));
        if (s.lost) { const yx = Y(s.t0 + a * frac + 0.5); ctx.strokeStyle = T.bad; ctx.lineWidth = 2; ctx.beginPath(); ctx.moveTo(xe - 6, yx - 6); ctx.lineTo(xe + 6, yx + 6); ctx.moveTo(xe + 6, yx - 6); ctx.lineTo(xe - 6, yx + 6); ctx.stroke(); ctx.lineWidth = 1; }
      });
      sim.acks.forEach(k => {
        ctx.strokeStyle = T.b; ctx.setLineDash([4, 3]); ctx.beginPath(); ctx.moveTo(xr, Y(k.t0)); ctx.lineTo(xs, Y(k.t1)); ctx.stroke(); ctx.setLineDash([]);
        ctx.fillStyle = T.b; ctx.textAlign = 'left'; ctx.textBaseline = 'middle'; ctx.fillText('ACK' + k.n, xr + 4, Y(k.t0));
      });
      sim.timeouts.forEach(o => { ctx.fillStyle = T.bad; ctx.textAlign = 'right'; ctx.textBaseline = 'middle'; ctx.fillText('TO F' + o.f, xs - 4, Y(o.t) - 1); ctx.strokeStyle = T.bad; ctx.beginPath(); ctx.moveTo(xs - 3, Y(o.t)); ctx.lineTo(xs + 8, Y(o.t)); ctx.stroke(); });
      ctx.restore();
      if (cursor < 1) { ctx.strokeStyle = T.pink; ctx.beginPath(); ctx.moveTo(0, Y(tc)); ctx.lineTo(w, Y(tc)); ctx.stroke(); ctx.fillStyle = T.pink; ctx.textAlign = 'left'; ctx.textBaseline = 'bottom'; ctx.fillText('t = ' + tc.toFixed(2), 4, Y(tc) - 2); }
    }
    function stop() { if (raf) { cancelAnimationFrame(raf); raf = 0; } ui.play.textContent = '▶ Play'; }
    function togglePlay() {
      if (raf) { stop(); return; }
      if (cursor >= 1) cursor = 0; ui.play.textContent = '❚❚ Pause'; lastTs = 0;
      const tick = ts => { if (!host.active) { stop(); return; } if (lastTs) cursor = Math.min(1, cursor + (ts - lastTs) / 8000); lastTs = ts; drawSeq(); if (cursor >= 1) { stop(); return; } raf = requestAnimationFrame(tick); };
      raf = requestAnimationFrame(tick);
    }
    function step() {
      stop(); if (!sim) return;
      const tEnd = (fin(sim.tEnd) ? sim.tEnd : 1) + 0.3, now = cursor >= 1 ? -1 : cursor * tEnd;
      const times = [].concat(sim.sends.map(s => s.t1), sim.sends.map(s => s.t1 + sim.a), sim.acks.map(k => k.t1), sim.timeouts.map(o => o.t)).filter(t => t > now + 1e-6).sort((x, y) => x - y);
      cursor = times.length ? Math.min(1, (times[0] + 1e-3) / tEnd) : 1; drawSeq();
    }
    const keys = ['R', 'L', 'd', 'v', 'W', 'P', 'kb', 'dn', 'dW', 'da', 'dl', 'ss0', 'rounds'];
    return {
      update, stop,
      get() { const o = {}; keys.forEach(k => { o[k] = sv(st, k); }); o.proto = S.proto; o.loss = S.losses; o.reno = S.reno; return o; },
      set(o) {
        setSliders(st, o, keys);
        if (o.proto === 'sw' || o.proto === 'gbn' || o.proto === 'sr') { S.proto = o.proto; ui.proto.value = S.proto; }
        if (okStr(o.loss, 300)) { S.losses = o.loss; ui.losses.value = S.losses; }
        if (o.reno === 'half' || o.reno === 'plus3') { S.reno = o.reno; ui.reno.value = S.reno; }
      },
    };
  }

  /* ---------- 5. error detection panel ---------- */
  function buildErr(root) {
    const M = FSP.math.net, { ctl, stage } = layout(root), ui = {};
    const S = { words: '0001 f203 f4f5 f6f7', rx: '0001 f203 f4f5 f6f7 220d', msg: '1101011011', gen: '10011', flip: '', pbits: '1011001', rows: '1011 0110 1100', pflip: '', odd: '0' };
    const f1 = FSP.ui.fieldset(ctl, 'Internet checksum (RFC 1071)');
    ui.words = textIn(f1, 'Data words', S.words, v => { S.words = v; update(); });
    ui.rx = textIn(f1, 'Received', S.rx, v => { S.rx = v; update(); });
    const r1 = el('div', { class: 'row' }); f1.appendChild(r1); FSP.ui.button(r1, 'Received ← data + checksum', () => { if (ui.sentBlock) { S.rx = ui.sentBlock; ui.rx.value = S.rx; update(); FSP.state.touch(); } }).classList.add('btn');
    ui.msgC = msgBox(f1);
    const f2 = FSP.ui.fieldset(ctl, 'CRC');
    ui.msg = textIn(f2, 'Message', S.msg, v => { S.msg = v; update(); }, { inputmode: 'numeric' });
    ui.gen = textIn(f2, 'Generator', S.gen, v => { S.gen = v; update(); });
    ui.flip = textIn(f2, 'Flip bit #', S.flip, v => { S.flip = v; update(); }, { inputmode: 'numeric', placeholder: 'none' });
    f2.appendChild(el('div', { class: 'note', text: 'Generator as bits (10011) or polynomial (x^4+x+1). Flip bit #: 1 = leftmost bit of the transmitted frame; comma-separate several; blank = no error.' }));
    ui.msgR = msgBox(f2);
    const f3 = FSP.ui.fieldset(ctl, 'Parity');
    ui.odd = FSP.ui.select(f3, 'Parity', [['0', 'Even'], ['1', 'Odd']], S.odd, v => { S.odd = v; update(); FSP.state.touch(); });
    ui.pbits = textIn(f3, 'Single', S.pbits, v => { S.pbits = v; update(); }, { inputmode: 'numeric' });
    ui.rows = textIn(f3, '2-D rows', S.rows, v => { S.rows = v; update(); }, { inputmode: 'numeric' });
    ui.pflip = textIn(f3, 'Flip r,c', S.pflip, v => { S.pflip = v; update(); }, { placeholder: 'e.g. 2,3' });
    f3.appendChild(el('div', { class: 'note', text: '2-D rows separated by spaces or "/" (equal lengths). Flip r,c (1-based, may include the parity row/column) to inject a single error.' }));
    ui.msgP = msgBox(f3);
    ui.hudC = hud(stage);
    ui.hudR = hud(stage);
    ui.cv = canvasIn(stage, 'Two-dimensional parity grid', 230); ui.cv.onResize(() => drawGrid());
    ui.hudP = hud(stage);
    ui.work = FSP.ui.working(stage);
    let grid = null;
    const hex = v => ('0000' + v.toString(16)).slice(-4), bin16 = v => ('0000000000000000' + v.toString(2)).slice(-16).replace(/(.{4})(?!$)/g, '$1 ');
    function update() {
      const W = [];
      // checksum
      const pw = M.parseHexWords(S.words), pr = M.parseHexWords(S.rx);
      if (!pw.ok || !pr.ok) { ui.msgC.show(!pw.ok ? 'Data: ' + pw.error : 'Received: ' + pr.error); hudSet(ui.hudC, []); ui.sentBlock = null; }
      else {
        ui.msgC.hide();
        const ck = M.inetChecksum(pw.words), vf = M.inetVerify(pr.words);
        ui.sentBlock = pw.words.map(hex).join(' ') + ' ' + hex(ck.checksum);
        hudSet(ui.hudC, [['One\'s-complement sum', hex(ck.sum)], ['Checksum = ~sum', hex(ck.checksum)], ['Receiver sum', hex(vf.sum) + (vf.ok ? '  = ffff ✓ accept' : '  ≠ ffff ✗ error detected')]]);
        W.push('— Internet checksum: 16-bit one\'s-complement sum, carries out of bit 15 wrap around (end-around carry) —');
        ck.steps.forEach((s, i) => {
          if (i === 0) { W.push('    ' + hex(s.w)); return; }
          W.push('  + ' + hex(s.w) + '  = ' + s.raw.toString(16) + (s.carry ? '  → carry wraps: ' + hex(s.raw & 0xFFFF) + ' + 1 = ' + hex(s.res) : ''));
        });
        W.push('  sum      = ' + hex(ck.sum) + '  = ' + bin16(ck.sum), '  checksum = ~sum = ' + hex(ck.checksum) + '  = ' + bin16(ck.checksum) + '  (flip every bit)',
          'Receiver adds all words including the checksum; result must be ffff (all ones):');
        vf.steps.forEach((s, i) => { W.push(i === 0 ? '    ' + hex(s.w) : '  + ' + hex(s.w) + '  = ' + s.raw.toString(16) + (s.carry ? ' → ' + hex(s.res) : '')); });
        W.push('  result ' + hex(vf.sum) + (vf.ok ? ' = ffff → no error detected' : ' ≠ ffff → error detected (complement = ' + hex((~vf.sum) & 0xFFFF) + ' ≠ 0)'), '');
      }
      // CRC
      const pm = M.parseBits(S.msg, 128), pg = M.parsePoly(S.gen);
      let flips = [], fe = '';
      if (S.flip.trim()) S.flip.split(/[,\s]+/).filter(Boolean).forEach(x => { const v = Number(x); if (!Number.isInteger(v) || v < 1) fe = 'Flip positions must be whole numbers ≥ 1.'; else flips.push(v); });
      if (!pm.ok || !pg.ok || fe) { ui.msgR.show(!pm.ok ? 'Message: ' + pm.error : !pg.ok ? pg.error : fe); hudSet(ui.hudR, []); }
      else {
        const g = pg.bits, cr = M.crcEncode(pm.bits, g), r = cr.r, nF = cr.frame.length;
        if (flips.some(v => v > nF)) { ui.msgR.show('Flip position must be ≤ frame length ' + nF + '.'); hudSet(ui.hudR, []); }
        else {
          ui.msgR.hide();
          const rx = cr.frame.slice(); flips.forEach(v => { rx[v - 1] ^= 1; });
          const chk = M.crcCheck(rx, g), gs = g.join(''), polyStr = g.map((b, i) => (b ? (g.length - 1 - i === 0 ? '1' : g.length - 1 - i === 1 ? 'x' : 'x^' + (g.length - 1 - i)) : null)).filter(Boolean).join(' + ');
          hudSet(ui.hudR, [['G(x)', polyStr + '  (' + gs + ', r = ' + r + ')'], ['Remainder (CRC)', cr.remainder.join('')], ['Transmitted frame', cr.frame.join('')],
            ['Received frame', rx.join('') + (flips.length ? '  (bit ' + flips.join(', ') + ' flipped)' : '')], ['Receiver remainder', chk.remainder.join('') + (chk.ok ? '  = 0 ✓ accept' : '  ≠ 0 ✗ error detected')]]);
          W.push('— CRC: append r = ' + r + ' zeros, divide by G = ' + gs + ' (mod-2, XOR, no borrows) —', 'Dividend = message · x^' + r + ' = ' + pm.bits.join('') + ' ' + '0'.repeat(r), '');
          W.push.apply(W, divLines(cr.dividend, g, cr.steps, cr.quotient));
          W.push('Remainder = ' + cr.remainder.join('') + ' → transmitted frame = message + remainder = ' + pm.bits.join('') + ' ' + cr.remainder.join(''), '');
          W.push('Receiver divides the received frame ' + rx.join('') + ' by G:');
          W.push.apply(W, divLines(rx, g, chk.steps, chk.quotient));
          W.push('Remainder = ' + chk.remainder.join('') + (chk.ok ? ' (all zero) → accept' : ' (non-zero) → error detected') + (chk.ok && flips.length ? '  ⚠ error pattern is a multiple of G(x): undetected!' : ''), '');
        }
      }
      // parity
      const odd = S.odd === '1', ps = M.parseBits(S.pbits, 64), rowsT = S.rows.split(/[\s\/;|]+/).filter(Boolean), rowsP = rowsT.map(x => M.parseBits(x, 16));
      let perr = '';
      if (!ps.ok) perr = 'Single: ' + ps.error;
      else if (!rowsT.length || rowsP.some(x => !x.ok)) perr = '2-D rows: each row must be 1–16 bits of 0/1.';
      else if (rowsP.some(x => x.bits.length !== rowsP[0].bits.length)) perr = '2-D rows must all have the same length.';
      else if (rowsT.length > 12) perr = 'At most 12 rows.';
      let pf = null;
      if (!perr && S.pflip.trim()) { const m = /^\s*(\d+)\s*[, ]\s*(\d+)\s*$/.exec(S.pflip); if (!m) perr = 'Flip r,c must look like 2,3.'; else pf = [+m[1] - 1, +m[2] - 1]; }
      if (perr) { ui.msgP.show(perr); grid = null; blank(ui.cv, 'Fix the parity inputs'); hudSet(ui.hudP, []); }
      else {
        const p2 = M.parity2D(rowsP.map(x => x.bits), odd), rx = p2.full.map(r => r.slice());
        if (pf && (pf[0] < 0 || pf[0] > p2.R || pf[1] < 0 || pf[1] > p2.C)) { ui.msgP.show('Flip position outside the ' + (p2.R + 1) + '×' + (p2.C + 1) + ' block.'); pf = null; } else ui.msgP.hide();
        if (pf) rx[pf[0]][pf[1]] ^= 1;
        const c2 = M.check2D(rx, odd), pb = M.parityBit(ps.bits, odd), ones = ps.bits.reduce((s, b) => s + b, 0);
        grid = { p2, rx, c2, pf };
        hudSet(ui.hudP, [['Single parity bit', String(pb) + '  → send ' + ps.bits.join('') + pb], ['2-D check', c2.ok ? 'all row/column parities OK' : 'bad rows ' + c2.badRows.map(x => x + 1).join(',') + ', bad cols ' + c2.badCols.map(x => x + 1).join(',')],
          ['Diagnosis', c2.ok ? (pf ? 'error not detected?!' : 'no error') : c2.locatable ? 'single error at (' + (c2.badRows[0] + 1) + ',' + (c2.badCols[0] + 1) + ') → correctable' : 'detected, cannot locate']]);
        W.push('— Parity (' + (odd ? 'odd' : 'even') + ') —', 'Single: ' + ps.bits.join('') + ' has ' + ones + ' ones → parity bit = ' + pb + ' so the total number of ones is ' + (odd ? 'odd' : 'even') + '. Detects any odd number of bit errors.',
          '2-D: row parity bits (right column) = ' + p2.full.slice(0, p2.R).map(r => r[p2.C]).join(' ') + ', column parity row = ' + p2.full[p2.R].join(' ') + ' (corner bit covers both).',
          'Receiver recomputes every row and column: a single error flips exactly one row and one column check → its intersection is the bad bit (correctable); 2 errors in one row are detected by columns only; some 4-bit rectangles go undetected.');
      }
      ui.work.set(W);
      drawGrid();
    }
    function divLines(dividend, g, steps, quotient) {
      const out = [], n = dividend.length, pad = '      ';
      out.push(pad + '  quotient ' + quotient.join(''));
      out.push(pad + '  ' + dividend.join(''));
      steps.forEach((s, k) => {
        out.push(pad + '⊕ ' + ' '.repeat(s.pos) + g.join('') + '   step ' + (k + 1) + ': leading 1 at bit ' + (s.pos + 1));
        const after = s.before.map((b, j) => b ^ g[j]), next = k + 1 < steps.length ? steps[k + 1].pos : n - g.length + 1;
        const show = after.join('') + dividend.slice(s.pos + g.length, Math.min(n, next + g.length)).map(() => '·').join('');
        out.push(pad + '  ' + ' '.repeat(s.pos) + '-'.repeat(g.length));
        out.push(pad + '  ' + ' '.repeat(s.pos) + show);
      });
      if (!steps.length) out.push(pad + '  (no leading 1 within the message part: remainder = last r bits)');
      return out;
    }
    function drawGrid() {
      if (!grid) return;
      const g = ui.cv.prep(), ctx = g.ctx, w = g.w, h = g.h, T = theme(), R = grid.p2.R + 1, C = grid.p2.C + 1;
      ctx.clearRect(0, 0, w, h); ctx.fillStyle = T.bg; ctx.fillRect(0, 0, w, h);
      const cs = Math.max(12, Math.min(34, (w - 80) / (C + 1), (h - 40) / (R + 1))), x0 = (w - (C + 1) * cs) / 2 + cs * 0.3, y0 = (h - (R + 1) * cs) / 2 + cs * 0.3;
      ctx.font = Math.round(cs * 0.45) + 'px ' + T.mono; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
      for (let i = 0; i < R; i++) for (let j = 0; j < C; j++) {
        const par = i === R - 1 || j === C - 1, bad = grid.c2.badRows.indexOf(i) >= 0 || grid.c2.badCols.indexOf(j) >= 0, flipped = grid.pf && grid.pf[0] === i && grid.pf[1] === j;
        ctx.fillStyle = flipped ? T.bad : bad ? T.warn : par ? T.c : T.border; ctx.globalAlpha = flipped ? 0.6 : bad ? 0.25 : par ? 0.25 : 0.5;
        ctx.fillRect(x0 + j * cs + 1, y0 + i * cs + 1, cs - 2, cs - 2); ctx.globalAlpha = 1;
        ctx.fillStyle = T.text; ctx.fillText(String(grid.rx[i][j]), x0 + (j + 0.5) * cs, y0 + (i + 0.5) * cs);
      }
      ctx.font = '10px ' + T.mono; ctx.fillStyle = T.muted;
      for (let i = 0; i < R; i++) { const ok = grid.c2.badRows.indexOf(i) < 0; ctx.fillStyle = ok ? T.b : T.bad; ctx.fillText(ok ? '✓' : '✗', x0 + C * cs + cs * 0.4, y0 + (i + 0.5) * cs); }
      for (let j = 0; j < C; j++) { const ok = grid.c2.badCols.indexOf(j) < 0; ctx.fillStyle = ok ? T.b : T.bad; ctx.fillText(ok ? '✓' : '✗', x0 + (j + 0.5) * cs, y0 + R * cs + cs * 0.4); }
      ctx.strokeStyle = T.c; ctx.setLineDash([3, 3]); ctx.beginPath(); ctx.moveTo(x0 + (C - 1) * cs, y0); ctx.lineTo(x0 + (C - 1) * cs, y0 + R * cs); ctx.moveTo(x0, y0 + (R - 1) * cs); ctx.lineTo(x0 + C * cs, y0 + (R - 1) * cs); ctx.stroke(); ctx.setLineDash([]);
    }
    const tk = ['words', 'rx', 'msg', 'gen', 'flip', 'pbits', 'rows', 'pflip'];
    return {
      update,
      get() { const o = {}; tk.forEach(k => { o[k] = S[k]; }); o.odd = S.odd; return o; },
      set(o) { tk.forEach(k => { if (okStr(o[k], 300)) { S[k] = o[k]; ui[k].value = S[k]; } }); if (o.odd === '0' || o.odd === '1') { S.odd = o.odd; ui.odd.value = S.odd; } },
    };
  }

  /* ---------- tab ---------- */
  const host = { active: false, sub: 'line' };
  let subs = null; const panels = {}, btns = {};
  const DEFS = [['line', 'Line coding', 'lc_'], ['delay', 'Delays & multiplexing', 'dl_'], ['mac', 'Multiple access', 'ma_'], ['arq', 'ARQ & TCP', 'aq_'], ['err', 'Error detection', 'ed_']];
  function showSub(id, silent) {
    if (!panels[id]) id = 'line';
    if (host.sub === 'arq' && id !== 'arq' && subs) subs.arq.stop();
    host.sub = id;
    Object.keys(panels).forEach(k => { panels[k].hidden = k !== id; btns[k].setAttribute('aria-pressed', String(k === id)); btns[k].classList.toggle('active', k === id); });
    if (host.active && subs) subs[id].update();
    if (!silent) FSP.state.touch();
  }
  FSP.registerTab({
    id: 'net', title: 'Networks & Protocols',
    init(panel) {
      if (typeof document !== 'undefined' && !document.getElementById('np-style')) {
        const s = document.createElement('style'); s.id = 'np-style';
        s.textContent = '.np-scroll{overflow-x:auto;max-width:100%;background:var(--panel);border:1px solid var(--border);border-radius:10px;padding:6px 8px}' +
          '.np-table{border-collapse:collapse;font:12px var(--mono);width:100%}.np-table th,.np-table td{border-bottom:1px dashed var(--border);padding:2px 8px;text-align:right;white-space:nowrap}' +
          '.np-table th{color:var(--muted);font-weight:400}.np-table th:first-child,.np-table td:first-child{text-align:left}';
        document.head.appendChild(s);
      }
      const bar = FSP.ui.el('div', { class: 'seg row', role: 'group', 'aria-label': 'Networks topic' }); panel.appendChild(bar);
      DEFS.forEach(d => { btns[d[0]] = FSP.ui.el('button', { type: 'button', class: 'seg-btn', 'aria-pressed': 'false', text: d[1], onclick: () => showSub(d[0]) }); bar.appendChild(btns[d[0]]); });
      DEFS.forEach(d => { panels[d[0]] = FSP.ui.el('div', { hidden: '' }); panel.appendChild(panels[d[0]]); });
      subs = { line: buildLine(panels.line), delay: buildDelay(panels.delay), mac: buildMAC(panels.mac), arq: buildARQ(panels.arq, host), err: buildErr(panels.err) };
      FSP.state.bind('net', {
        get() { const o = { sub: host.sub }; DEFS.forEach(d => { const g = subs[d[0]].get(); Object.keys(g).forEach(k => { o[d[2] + k] = g[k]; }); }); return o; },
        set(o) {
          o = o || {};
          DEFS.forEach(d => { const r = {}; Object.keys(o).forEach(k => { if (k.indexOf(d[2]) === 0) r[k.slice(d[2].length)] = o[k]; }); try { subs[d[0]].set(r); } catch (e) { /* ignore bad state */ } });
          showSub(panels[o.sub] ? o.sub : host.sub, true);
        },
      });
      showSub(host.sub, true);
    },
    activate() { host.active = true; if (subs) subs[host.sub].update(); },
    deactivate() { host.active = false; if (subs) subs.arq.stop(); },
  });
})();
