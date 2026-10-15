/* FSP module: EE2005 Electrical Devices & Machines — Transformers & Actuators (tab id 'edm2').
   Sub-topics: multi-limb magnetic circuit (B–H table, forward chain + nonlinear solve), inductance/energy/force/core loss,
   autotransformer, three-phase transformer banks, all-day efficiency.
   Pure math lives in FSP.math.edm2 (no DOM). All DOM work happens inside init/activate. */
(function () {
  'use strict';
  const PI = Math.PI, S3 = Math.sqrt(3), MU0 = 4e-7 * PI;   // exam value of mu0
  const fin = Number.isFinite;
  const E = {}; FSP.math.edm2 = E; E.MU0 = MU0;
  const NAMES = ['left', 'centre', 'right'], SH = ['L', 'C', 'R'];

  /* ---------- number formatting for working text (pure) ---------- */
  const SUPD = { '-': '⁻', '0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴', '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹' };
  function fx(x, d) {
    if (!fin(x)) return '—'; d = d || 4; if (x === 0) return '0';
    const a = Math.abs(x);
    if (a >= 1e5 || a < 1e-3) {
      let e = Math.floor(Math.log10(a)), m = +(x / Math.pow(10, e)).toPrecision(d);
      if (Math.abs(m) >= 10) { m = +(m / 10).toPrecision(d); e++; }
      return m + '×10' + String(e).split('').map(ch => SUPD[ch]).join('');
    }
    return String(+x.toPrecision(d));
  }
  E.fx = fx;
  const arrowOf = phi => (phi > 0 ? '↑' : phi < 0 ? '↓' : '');

  /* =====================================================================
     1. Materials: constant mu_r or B–H table with LINEAR interpolation
     ===================================================================== */
  // The student's past-paper table (origin (0,0) is added implicitly).
  E.DEFAULT_TABLE = { H: [125, 286, 350, 500, 675, 1000, 1200], B: [0.4, 0.75, 0.85, 1.0, 1.08, 1.15, 1.175] };
  E.tableToText = (t, sep) => t.H.map((h, i) => h + ' ' + t.B[i]).join(sep || '\n');

  // Parse rows "H B" (H in A/m, B in T; separators: space, comma, newline, semicolon). Returns {ok, H, B, msg} without the origin.
  E.parseTable = function (text) {
    const rows = []; let junk = null;
    String(text == null ? '' : text).split(/\r?\n|;/).forEach(line => {
      line = line.replace(/#.*$/, '').trim(); if (!line) return;
      const p = line.split(/[\s,]+/).map(Number);
      if (p.length >= 2 && fin(p[0]) && fin(p[1])) rows.push([p[0], p[1]]); else if (junk === null) junk = line;
    });
    if (junk !== null) return { ok: false, msg: 'Cannot read row "' + junk + '" (expected: H B).' };
    const pts = rows.filter(r => !(r[0] === 0 && r[1] === 0)).sort((a, b) => a[0] - b[0]);
    if (pts.length < 2) return { ok: false, msg: 'Need at least two (H, B) points with H > 0.' };
    for (let i = 0; i < pts.length; i++) {
      if (!(pts[i][0] > 0) || !(pts[i][1] > 0)) return { ok: false, msg: 'H and B must be positive (first quadrant).' };
      if (i > 0 && !(pts[i][0] > pts[i - 1][0])) return { ok: false, msg: 'H values must be strictly increasing (duplicate H = ' + pts[i][0] + ').' };
      if (i > 0 && !(pts[i][1] > pts[i - 1][1])) return { ok: false, msg: 'B must strictly increase with H (check H = ' + pts[i][0] + ').' };
    }
    return { ok: true, H: pts.map(r => r[0]), B: pts.map(r => r[1]), msg: '' };
  };

  // All material functions are odd. Table: piecewise-linear through (0,0) and the points; beyond the last point the
  // last segment is extended (flagged by outside()).
  E.makeMaterial = function (o) {
    if (!o || o.mode !== 'table') {
      const mur = o && o.mur > 0 ? o.mur : NaN, mu = MU0 * mur;
      return {
        mode: 'const', mur, ok: mur > 0, Hof: B => B / mu, Bof: H => H * mu, dHdB: () => 1 / mu, outside: () => false,
        wDensity: B => B * B / (2 * mu),
        describeH: B => 'H = B/(μ0·μr) = ' + fx(B) + '/(4π×10⁻⁷ × ' + fx(mur) + ') = ' + fx(B / mu) + ' A/m',
        describeB: H => 'B = μ0·μr·H = 4π×10⁻⁷ × ' + fx(mur) + ' × ' + fx(H) + ' = ' + fx(H * mu) + ' T',
      };
    }
    const Hs = [0].concat(o.H), Bs = [0].concat(o.B), n = Hs.length;
    const segB = b => { for (let i = 0; i < n - 2; i++) if (b <= Bs[i + 1]) return i; return n - 2; };
    const segH = h => { for (let i = 0; i < n - 2; i++) if (h <= Hs[i + 1]) return i; return n - 2; };
    const H1 = b => { const i = segB(b); return Hs[i] + (b - Bs[i]) * (Hs[i + 1] - Hs[i]) / (Bs[i + 1] - Bs[i]); };
    const B1 = h => { const i = segH(h); return Bs[i] + (h - Hs[i]) * (Bs[i + 1] - Bs[i]) / (Hs[i + 1] - Hs[i]); };
    const out = b => Math.abs(b) > Bs[n - 1] * (1 + 1e-12);
    const node = (arr, v) => arr.findIndex(z => Math.abs(z - v) <= 1e-12 * Math.max(1, Math.abs(z)));
    return {
      mode: 'table', ok: true, H: Hs, B: Bs, Bmax: Bs[n - 1],
      Hof: b => (b < 0 ? -H1(-b) : H1(b)),
      Bof: h => (h < 0 ? -B1(-h) : B1(h)),
      dHdB: b => { const i = segB(Math.abs(b)); return (Hs[i + 1] - Hs[i]) / (Bs[i + 1] - Bs[i]); },
      outside: out,
      // energy density integral_0^B H dB (exact for piecewise-linear H(B): trapezoids)
      wDensity(b) {
        const a = Math.abs(b); let w = 0;
        for (let i = 0; i < n - 1; i++) {
          const top = i === n - 2 ? a : Math.min(a, Bs[i + 1]); if (top <= Bs[i]) break;
          w += 0.5 * (Hs[i] + H1(top)) * (top - Bs[i]);
        }
        return w;
      },
      describeH(b) {
        const a = Math.abs(b), k = node(Bs, a); if (k > 0) return 'H = ' + fx(Hs[k]) + ' A/m (table point)';
        const i = segB(a), h = H1(a);
        return 'H = ' + fx(Hs[i]) + ' + (' + fx(a) + ' − ' + fx(Bs[i]) + ')/(' + fx(Bs[i + 1]) + ' − ' + fx(Bs[i]) + ') × (' + fx(Hs[i + 1]) + ' − ' + fx(Hs[i]) + ') = ' + fx(h) + ' A/m' + (out(a) ? '  [B beyond table — extrapolated!]' : '');
      },
      describeB(h) {
        const a = Math.abs(h), k = node(Hs, a); if (k > 0) return 'B = ' + fx(Bs[k]) + ' T (table point)';
        const i = segH(a), bb = B1(a);
        return 'B = ' + fx(Bs[i]) + ' + (' + fx(a) + ' − ' + fx(Hs[i]) + ')/(' + fx(Hs[i + 1]) + ' − ' + fx(Hs[i]) + ') × (' + fx(Bs[i + 1]) + ' − ' + fx(Bs[i]) + ') = ' + fx(bb) + ' T' + (out(bb) ? '  [beyond table — extrapolated!]' : '');
      },
    };
  };

  /* =====================================================================
     2. Three-limb magnetic circuit
     Convention: flux and coil MMF positive UPWARD (↑) in every limb.
     For limb k:  s_k N_k I_k − D_k(Φ_k) = U   (U = MMF rise bottom yoke → top yoke, common to all limbs)
     Flux law at the top junction: Σ Φ_k = 0.
     D_k(Φ) = H(Φ/A_k)·(l_k − l_g) + Φ·l_g/(μ0·A_g).  l_k is the limb's mean path between the two junctions.
     ===================================================================== */
  // s: {on, l (m), w (m), d (m), lg (m), fringe, N, I (A), sense (+1 ↑, −1 ↓)}
  E.makeLimb = function (s) {
    const lg = s.lg > 0 ? s.lg : 0, A = s.w * s.d;
    const Ag = s.fringe && lg > 0 ? (s.w + lg) * (s.d + lg) : A;
    return { on: !!s.on, l: s.l, w: s.w, d: s.d, A, lg, Ag, li: s.l - lg, N: s.N || 0, I: s.I || 0, sense: s.sense < 0 ? -1 : 1, fringe: !!s.fringe && lg > 0 };
  };
  E.limbMMF = L => L.sense * L.N * L.I;
  E.drop = function (L, mat, phi) {
    const B = phi / L.A, H = mat.Hof(B), Fi = H * L.li, Bg = phi / L.Ag, Hg = Bg / MU0, Fg = L.lg > 0 ? Hg * L.lg : 0;
    return { phi, B, H, Fi, Bg, Hg, Fg, F: Fi + Fg, outside: mat.outside(B) };
  };
  E.dropSlope = (L, mat, phi) => mat.dHdB(phi / L.A) / L.A * L.li + (L.lg > 0 ? L.lg / (MU0 * L.Ag) : 0);
  // inverse of the limb drop: Φ such that D(Φ) = F (monotone)
  E.dropInv = function (L, mat, F) {
    if (F === 0) return 0; const s = F < 0 ? -1 : 1, a = Math.abs(F);
    if (!(L.lg > 0)) return s * L.A * mat.Bof(a / L.li);
    if (mat.mode === 'const') return s * a / (L.li / (MU0 * mat.mur * L.A) + L.lg / (MU0 * L.Ag));
    let lo = 0, hi = Math.min(L.A * mat.Bof(a / L.li), a * MU0 * L.Ag / L.lg);
    for (let k = 0; k < 200; k++) { const m = 0.5 * (lo + hi); if (E.drop(L, mat, m).F < a) lo = m; else hi = m; if (hi - lo <= 1e-15 * hi) break; }
    return s * 0.5 * (lo + hi);
  };
  E.mcValidate = function (c) {
    if (!c || !c.mat || !c.mat.ok) return 'Material is invalid (μr must be > 0, or fix the B–H table).';
    const P = [0, 1, 2].filter(k => c.limbs[k].on);
    if (P.length < 2) return 'At least two limbs are needed to form a closed magnetic path.';
    for (const k of P) {
      const L = c.limbs[k], nm = NAMES[k] + ' limb: ';
      if (!(L.l > 0) || !(L.w > 0) || !(L.d > 0)) return nm + 'length, width and depth must be > 0.';
      if (!(L.lg >= 0) || !(L.li > 0)) return nm + 'air gap must be shorter than the limb length.';
      if (!(L.N >= 0) || !fin(L.I)) return nm + 'turns must be ≥ 0 and current finite.';
    }
    return '';
  };
  function mcPack(c, phi, extra) {
    const P = [0, 1, 2].filter(k => c.limbs[k].on), S = c.limbs.map(L => (L.on ? E.limbMMF(L) : 0));
    const drops = c.limbs.map((L, k) => (L.on ? E.drop(L, c.mat, phi[k]) : null));
    const Us = P.map(k => S[k] - drops[k].F), U = Us.reduce((a, v) => a + v, 0) / Us.length;
    const kcl = P.reduce((a, k) => a + phi[k], 0), mmfResid = Math.max.apply(null, Us.map(u => Math.abs(u - U)));
    const outside = P.some(k => drops[k].outside);
    return Object.assign({ ok: true, msg: '', P, S, phi, drops, U, kcl, mmfResid, outside }, extra || {});
  }
  // Independent solver: bisection on U.  Σ_k D_k^{-1}(S_k − U) is decreasing in U and changes sign on [min S, max S].
  E.mcSolveBisect = function (c) {
    const err = E.mcValidate(c); if (err) return { ok: false, msg: err };
    const P = [0, 1, 2].filter(k => c.limbs[k].on), S = c.limbs.map(L => (L.on ? E.limbMMF(L) : 0));
    const f = U => P.reduce((a, k) => a + E.dropInv(c.limbs[k], c.mat, S[k] - U), 0);
    let lo = Math.min.apply(null, P.map(k => S[k])), hi = Math.max.apply(null, P.map(k => S[k])), it = 0;
    if (hi > lo) for (; it < 300; it++) { const m = 0.5 * (lo + hi); if (f(m) > 0) lo = m; else hi = m; if (hi - lo <= 1e-15 * Math.max(Math.abs(hi), 1)) break; }
    const U = 0.5 * (lo + hi), phi = [0, 0, 0]; P.forEach(k => { phi[k] = E.dropInv(c.limbs[k], c.mat, S[k] - U); });
    return mcPack(c, phi, { method: 'bisection on U', iter: it, converged: true });
  };
  // Given all coil currents, find the fluxes: damped Newton on the window (loop) equations, fallback to bisection.
  // Unknowns: the fluxes of the non-reference limbs; reference = centre limb (if present). Φ_ref = −Σ others.
  // Loop j:  [S_j − D_j(Φ_j)] − [S_ref − D_ref(Φ_ref)] = 0   (Kirchhoff's MMF law around window j)
  E.mcSolveCurrents = function (c) {
    const err = E.mcValidate(c); if (err) return { ok: false, msg: err };
    const Ls = c.limbs, mat = c.mat, P = [0, 1, 2].filter(k => Ls[k].on), S = Ls.map(L => (L.on ? E.limbMMF(L) : 0));
    const ref = P.indexOf(1) >= 0 ? 1 : P[P.length - 1], oth = P.filter(k => k !== ref), m = oth.length;
    const phiAll = x => { const p = [0, 0, 0]; let sum = 0; oth.forEach((k, j) => { p[k] = x[j]; sum += x[j]; }); p[ref] = -sum; return p; };
    const resid = x => { const p = phiAll(x), Ur = S[ref] - E.drop(Ls[ref], mat, p[ref]).F; return oth.map(k => S[k] - E.drop(Ls[k], mat, p[k]).F - Ur); };
    const norm = r => Math.max.apply(null, r.map(Math.abs));
    const tol = 1e-10 * Math.max(1, P.reduce((a, k) => a + Math.abs(S[k]), 0));
    let x = oth.map(() => 0), r = resid(x), it = 0, conv = norm(r) <= tol; const hist = [norm(r)];
    while (!conv && it < 100) {
      it++;
      const p = phiAll(x), dr = E.dropSlope(Ls[ref], mat, p[ref]);
      const J = oth.map((k, j) => oth.map((k2, j2) => -(j === j2 ? E.dropSlope(Ls[k], mat, p[k]) : 0) - dr));
      let d;
      if (m === 1) d = [-r[0] / J[0][0]];
      else { const det = J[0][0] * J[1][1] - J[0][1] * J[1][0]; d = [(-r[0] * J[1][1] + r[1] * J[0][1]) / det, (-r[1] * J[0][0] + r[0] * J[1][0]) / det]; }
      if (!d.every(fin)) break;
      let t = 1, xn = x, rn = r; const n0 = norm(r);
      for (let k = 0; k < 50; k++) { xn = x.map((v, j) => v + t * d[j]); rn = resid(xn); if (norm(rn) < n0 * (1 - 1e-4 * t) || norm(rn) <= tol) break; t /= 2; }
      x = xn; r = rn; hist.push(norm(r)); conv = norm(r) <= tol;
      if (!fin(norm(r))) break;
    }
    if (conv) return mcPack(c, phiAll(x), { method: 'Newton (damped) on the window equations', iter: it, converged: true, hist, ref, oth });
    const b = E.mcSolveBisect(c);
    return Object.assign(b, { method: 'bisection on U (Newton did not converge in ' + it + ' iterations)', hist, ref, oth });
  };

  // Text block: MMF drop in one limb for flux phi (magnitudes, tutorial style)
  function dropLines(c, k, phi) {
    const L = c.limbs[k], d = E.drop(L, c.mat, phi), a = Math.abs, out = [];
    out.push(cap(NAMES[k]) + ' limb: |Φ| = ' + fx(a(phi)) + ' Wb ' + arrowOf(phi) + ', A = ' + fx(L.w * 100) + ' × ' + fx(L.d * 100) + ' = ' + fx(L.A * 1e4) + ' cm²  →  B = Φ/A = ' + fx(a(d.B)) + ' T');
    out.push('   iron: ' + c.mat.describeH(a(d.B)) + ';  l_iron = ' + (L.lg > 0 ? 'l − l_g = ' + fx(L.l) + ' − ' + fx(L.lg) + ' = ' : '') + fx(L.li) + ' m;  F_iron = H·l = ' + fx(a(d.Fi)) + ' At');
    if (L.lg > 0) {
      out.push('   gap: A_g = ' + (L.fringe ? '(a + l_g)(b + l_g) = (' + fx(L.w) + ' + ' + fx(L.lg) + ')(' + fx(L.d) + ' + ' + fx(L.lg) + ') = ' : 'A = ') + fx(L.Ag * 1e4) + ' cm²;  B_g = Φ/A_g = ' + fx(a(d.Bg)) + ' T;  H_g = B_g/μ0 = ' + fx(a(d.Hg)) + ' A/m;  F_gap = H_g·l_g = ' + fx(a(d.Fg)) + ' At');
    }
    out.push('   MMF drop D_' + SH[k] + ' = ' + (L.lg > 0 ? fx(a(d.Fi)) + ' + ' + fx(a(d.Fg)) + ' = ' : '') + fx(a(d.F)) + ' At ' + arrowOf(phi));
    return out;
  }
  // Inverse step text: given MMF drop F across limb k, find Φ
  function invLines(c, k, F) {
    const L = c.limbs[k], phi = E.dropInv(L, c.mat, F), out = [];
    if (!(L.lg > 0)) {
      const H = Math.abs(F) / L.li;
      out.push(cap(NAMES[k]) + ' limb: D_' + SH[k] + ' = ' + fx(F) + ' At  →  H = |D|/l = ' + fx(Math.abs(F)) + '/' + fx(L.li) + ' = ' + fx(H) + ' A/m');
      out.push('   ' + c.mat.describeB(H) + ';  |Φ| = B·A = ' + fx(Math.abs(phi) / L.A) + ' × ' + fx(L.A) + ' = ' + fx(Math.abs(phi)) + ' Wb ' + arrowOf(phi));
    } else {
      out.push(cap(NAMES[k]) + ' limb (has a gap): need Φ with H(Φ/A)·l_iron + Φ·l_g/(μ0 A_g) = ' + fx(Math.abs(F)) + ' At;');
      out.push('   solved numerically (bisection): |Φ| = ' + fx(Math.abs(phi)) + ' Wb ' + arrowOf(phi) + ', B = ' + fx(Math.abs(phi) / L.A) + ' T');
    }
    return { phi, lines: out };
  }
  function cap(s) { return s.charAt(0).toUpperCase() + s.slice(1); }
  function sameLimb(a, b) { const r = (u, v) => Math.abs(u - v) <= 1e-12 * Math.max(Math.abs(u), Math.abs(v), 1e-30); return r(a.l, b.l) && r(a.A, b.A) && r(a.lg, b.lg) && r(a.Ag, b.Ag) && r(E.limbMMF(a), E.limbMMF(b)); }

  // Given the flux phiG (signed, ↑ +) in limb g, find the current in the coil on limb x (other coil currents known).
  E.mcSolveForCurrent = function (c, g, phiG, x) {
    const err = E.mcValidate(c); if (err) return { ok: false, msg: err };
    const Ls = c.limbs, mat = c.mat, P = [0, 1, 2].filter(k => Ls[k].on), S = Ls.map(L => (L.on ? E.limbMMF(L) : 0));
    if (!Ls[g] || !Ls[g].on) return { ok: false, msg: 'The limb with the given flux is not present.' };
    if (!Ls[x] || !Ls[x].on) return { ok: false, msg: 'The limb with the unknown coil is not present.' };
    if (!(Ls[x].N > 0)) return { ok: false, msg: 'The unknown coil (' + NAMES[x] + ' limb) needs N > 0 turns.' };
    if (!fin(phiG)) return { ok: false, msg: 'Given flux is not a number.' };
    S[x] = 0; // unknown
    const steps = [], phi = [0, 0, 0]; let U, method = 'forward chain';
    phi[g] = phiG;
    steps.push('Convention: Φ and coil MMF positive ↑ in every limb. Each limb: s·N·I − D(Φ) = U (same U for all limbs, they join the same two yokes); flux law ΣΦ = 0.');
    steps.push('Given: ' + NAMES[g] + ' limb |Φ| = ' + fx(Math.abs(phiG)) + ' Wb ' + arrowOf(phiG) + ' (B = ' + fx(Math.abs(phiG) / Ls[g].A) + ' T).');
    dropLines(c, g, phiG).forEach(s => steps.push(s));
    const others = P.filter(k => k !== g && k !== x);
    if (x !== g) {
      const Dg = E.drop(Ls[g], mat, phiG).F; U = S[g] - Dg;
      steps.push('U = s·N·I(' + SH[g] + ') − D_' + SH[g] + ' = ' + fx(S[g]) + ' − (' + fx(Dg) + ') = ' + fx(U) + ' At');
      let sum = phiG;
      others.forEach(k => {
        const F = S[k] - U; steps.push('MMF law for the ' + NAMES[k] + ' limb: D_' + SH[k] + ' = s·N·I − U = ' + fx(S[k]) + ' − (' + fx(U) + ') = ' + fx(F) + ' At');
        const r = invLines(c, k, F); r.lines.forEach(s => steps.push(s)); phi[k] = r.phi; sum += r.phi;
      });
      phi[x] = -sum;
      steps.push('Flux law: Φ_' + SH[x] + ' = −(' + [g].concat(others).map(k => 'Φ_' + SH[k]).join(' + ') + ') = ' + fx(phi[x]) + ' Wb  (|Φ| = ' + fx(Math.abs(phi[x])) + ' Wb ' + arrowOf(phi[x]) + ')');
      dropLines(c, x, phi[x]).forEach(s => steps.push(s));
    } else {
      const Dg = E.drop(Ls[g], mat, phiG).F;
      if (others.length === 1) {
        const y = others[0]; phi[y] = -phiG;
        steps.push('Only one other limb: Φ_' + SH[y] + ' = −Φ_' + SH[g] + ' = ' + fx(phi[y]) + ' Wb (series path)');
        dropLines(c, y, phi[y]).forEach(s => steps.push(s));
        U = S[y] - E.drop(Ls[y], mat, phi[y]).F; steps.push('U = s·N·I(' + SH[y] + ') − D_' + SH[y] + ' = ' + fx(U) + ' At');
      } else if (sameLimb(Ls[others[0]], Ls[others[1]])) {
        const [y, z] = others; phi[y] = phi[z] = -phiG / 2;
        steps.push('The ' + NAMES[y] + ' and ' + NAMES[z] + ' limbs are identical (same l, A, gap, coil MMF): by symmetry each carries half the flux, Φ = ' + fx(phi[y]) + ' Wb');
        dropLines(c, y, phi[y]).forEach(s => steps.push(s));
        U = S[y] - E.drop(Ls[y], mat, phi[y]).F; steps.push('U = s·N·I(' + SH[y] + ') − D_' + SH[y] + ' = ' + fx(U) + ' At  (the ' + NAMES[z] + ' limb has the same drop)');
        method = 'forward chain (symmetry)';
      } else {
        const [y, z] = others, f = u => E.dropInv(Ls[y], mat, S[y] - u) + E.dropInv(Ls[z], mat, S[z] - u) + phiG;
        let lo = -1, hi = 1, k = 0; while (f(lo) < 0 && k++ < 200) lo *= 2; k = 0; while (f(hi) > 0 && k++ < 200) hi *= 2;
        let it = 0; for (; it < 300; it++) { const mm = 0.5 * (lo + hi); if (f(mm) > 0) lo = mm; else hi = mm; if (hi - lo <= 1e-14 * Math.max(1, Math.abs(hi))) break; }
        U = 0.5 * (lo + hi); phi[y] = E.dropInv(Ls[y], mat, S[y] - U); phi[z] = E.dropInv(Ls[z], mat, S[z] - U);
        steps.push('The ' + NAMES[y] + ' and ' + NAMES[z] + ' limbs are in parallel but not identical: find U with Φ_' + SH[y] + '(U) + Φ_' + SH[z] + '(U) = −Φ_' + SH[g] + ' (bisection, ' + it + ' iterations): U = ' + fx(U) + ' At');
        dropLines(c, y, phi[y]).forEach(s => steps.push(s)); dropLines(c, z, phi[z]).forEach(s => steps.push(s));
        method = 'bisection on U';
      }
    }
    const Dx = E.drop(Ls[x], mat, phi[x]).F, sNI = U + Dx, I = sNI / (Ls[x].sense * Ls[x].N);
    steps.push('Coil on the ' + NAMES[x] + ' limb: s·N·I = U + D_' + SH[x] + ' = ' + fx(U) + ' + (' + fx(Dx) + ') = ' + fx(sNI) + ' At');
    steps.push('I = ' + fx(sNI) + ' / (' + (Ls[x].sense > 0 ? '+' : '−') + fx(Ls[x].N) + ') = ' + fx(I) + ' A' + (I < 0 ? '  (negative: current must flow opposite to the chosen winding sense)' : ''));
    const c2 = { mat, limbs: Ls.map((L, k) => (k === x ? Object.assign({}, L, { I }) : L)) };
    return mcPack(c2, phi, { I, NI: Math.abs(sNI), sNI, steps, method, x, g });
  };

  /* =====================================================================
     3. Inductance, energy, force, torque, sinusoidal excitation, core loss
     ===================================================================== */
  // Single-coil gapped core / relay: {N, lc (iron), A, lg (total gap), Ag, mat}
  E.coreNI = (k, phi) => k.mat.Hof(phi / k.A) * k.lc + (k.lg > 0 ? phi * k.lg / (MU0 * k.Ag) : 0);
  E.coreFlux = function (k, i) {
    const NI = k.N * i; if (NI === 0) return 0; const s = NI < 0 ? -1 : 1, a = Math.abs(NI);
    if (k.mat.mode === 'const') return s * a / (k.lc / (MU0 * k.mat.mur * k.A) + (k.lg > 0 ? k.lg / (MU0 * k.Ag) : 0));
    let lo = 0, hi = k.A * k.mat.Bof(a / k.lc);
    for (let j = 0; j < 200; j++) { const m = 0.5 * (lo + hi); if (E.coreNI(k, m) < a) lo = m; else hi = m; if (hi - lo <= 1e-15 * hi) break; }
    return s * 0.5 * (lo + hi);
  };
  E.coreAnalyze = function (k, i) {
    const phi = E.coreFlux(k, i), lam = k.N * phi, B = phi / k.A, Bg = phi / k.Ag;
    const Wcore = k.mat.wDensity(B) * k.A * k.lc, Wgap = k.lg > 0 ? Bg * Bg / (2 * MU0) * k.Ag * k.lg : 0, W = Wcore + Wgap;
    const Rc = k.mat.mode === 'const' ? k.lc / (MU0 * k.mat.mur * k.A) : (phi !== 0 ? k.mat.Hof(B) * k.lc / phi : NaN);
    const Rg = k.lg > 0 ? k.lg / (MU0 * k.Ag) : 0;
    const h = Math.max(Math.abs(i) * 1e-5, 1e-9), Linc = k.N * (E.coreFlux(k, i + h) - E.coreFlux(k, i - h)) / (2 * h);
    return { phi, lam, B, Bg, H: k.mat.Hof(B), Lsec: i !== 0 ? lam / i : NaN, Linc, Rc, Rg, W, Wcore, Wgap, Wco: lam * i - W, outside: k.mat.outside(B) };
  };
  // co-energy W'(i) = ∫_0^i λ(i') di'  (Simpson, n even)
  E.coEnergy = function (k, i, n) {
    n = n || 200; if (n % 2) n++; if (i === 0) return 0; const h = i / n; let s = 0;
    for (let j = 0; j <= n; j++) s += (j === 0 || j === n ? 1 : j % 2 ? 4 : 2) * k.N * E.coreFlux(k, j * h);
    return s * h / 3;
  };
  // energy at fixed flux linkage: W(λ) = ∫_0^λ i dλ'  (Simpson)
  E.energyAtLambda = function (k, lam, n) {
    n = n || 200; if (n % 2) n++; if (lam === 0) return 0; const h = lam / n; let s = 0;
    for (let j = 0; j <= n; j++) s += (j === 0 || j === n ? 1 : j % 2 ? 4 : 2) * E.coreNI(k, j * h / k.N) / k.N;
    return s * h / 3;
  };
  // Relay / plunger: r = {N, lc, A, ng (number of gaps), mat}; gap length x each, total gap ng·x, pole area A (no fringing)
  E.relayCore = (r, x) => ({ N: r.N, lc: r.lc, A: r.A, lg: r.ng * x, Ag: r.A, mat: r.mat });
  E.relayL = function (r, x) { const R = r.lc / (MU0 * r.mat.mur * r.A) + r.ng * x / (MU0 * r.A); return r.N * r.N / R; };
  E.relaydLdx = function (r, x) { const R = r.lc / (MU0 * r.mat.mur * r.A) + r.ng * x / (MU0 * r.A); return -r.N * r.N * (r.ng / (MU0 * r.A)) / (R * R); };
  // force in +x direction (opening the gap); negative = attraction.  Linear: ½ i² dL/dx.  Nonlinear: ∂W'/∂x at constant i.
  E.relayForce = function (r, i, x) {
    if (r.mat.mode === 'const') return 0.5 * i * i * E.relaydLdx(r, x);
    // Nonlinear: at constant flux linkage the iron energy does not depend on x, only the gap energy does:
    // F = −∂W/∂x|λ = −n_g·B²A/(2μ0) with B = Φ/A the (fixed) gap flux density. Exact, and cheaper and far less noisy than
    // differencing a numerically integrated co-energy across the kinks of a piecewise-linear table.
    const B = E.coreFlux(E.relayCore(r, x), i) / r.A;
    return -r.ng * B * B * r.A / (2 * MU0);
  };
  E.relayForceFDco = function (r, i, x, h, n) { h = h || x * 1e-4; n = n || 400; return (E.coEnergy(E.relayCore(r, x + h), i, n) - E.coEnergy(E.relayCore(r, x - h), i, n)) / (2 * h); };
  E.relayForceFDenergy = function (r, i, x, h) {
    h = h || x * 1e-4; const lam = r.N * E.coreFlux(E.relayCore(r, x), i);
    return -(E.energyAtLambda(E.relayCore(r, x + h), lam, 400) - E.energyAtLambda(E.relayCore(r, x - h), lam, 400)) / (2 * h);
  };
  // pull-in current: |F(i, x0)| = Fs
  E.relayPullIn = function (r, x0, Fs) {
    if (!(Fs > 0) || !(x0 > 0)) return NaN;
    if (r.mat.mode === 'const') return Math.sqrt(2 * Fs / Math.abs(E.relaydLdx(r, x0)));
    let lo = 0, hi = 1e-3, k = 0; while (Math.abs(E.relayForce(r, hi, x0)) < Fs && k++ < 80) hi *= 2;
    if (k >= 80) return NaN;
    for (let j = 0; j < 80; j++) { const m = 0.5 * (lo + hi); if (Math.abs(E.relayForce(r, m, x0)) < Fs) lo = m; else hi = m; }
    return 0.5 * (lo + hi);
  };
  // rotary reluctance actuator: L(θ) = (Ld+Lq)/2 + (Ld−Lq)/2·cos 2θ (θ from the aligned d-axis)
  E.rotL = (Ld, Lq, th) => (Ld + Lq) / 2 + (Ld - Lq) / 2 * Math.cos(2 * th);
  E.rotTorque = (Ld, Lq, i, th) => -0.5 * i * i * (Ld - Lq) * Math.sin(2 * th);

  // sinusoidal excitation: E_rms = (2π/√2) f N Φmax ≈ 4.44 f N Φmax
  E.K444 = 2 * PI / Math.SQRT2;
  E.sineFlux = function (V, f, N, A) { const phim = V / (E.K444 * f * N); return { phim, Bm: phim / A, lamPk: N * phim }; };
  E.coreLoss = function (o) { const Ph = o.kh * o.f * Math.pow(o.Bm, o.n) * o.vol, Pe = o.ke * o.f * o.f * o.Bm * o.Bm * o.vol; return { Ph, Pe, P: Ph + Pe }; };
  // two-frequency separation at the same Bmax:  P/f = a + b·f,  a = k_h Bm^n V,  b = k_e Bm² V
  E.separateLosses = function (o) {
    const ok = [o.f1, o.P1, o.f2, o.P2, o.Bm, o.n, o.vol].every(v => fin(v) && v > 0);
    if (!ok) return { ok: false, msg: 'Frequencies, powers, Bmax, n and volume must all be positive.' };
    if (Math.abs(o.f2 - o.f1) < 1e-9 * Math.max(o.f1, o.f2)) return { ok: false, msg: 'The two test frequencies must differ.' };
    const y1 = o.P1 / o.f1, y2 = o.P2 / o.f2, b = (y2 - y1) / (o.f2 - o.f1), a = y1 - b * o.f1;
    const res = { ok: true, msg: '', y1, y2, a, b, kh: a / (Math.pow(o.Bm, o.n) * o.vol), ke: b / (o.Bm * o.Bm * o.vol) };
    res.Ph1 = a * o.f1; res.Pe1 = b * o.f1 * o.f1; res.Ph2 = a * o.f2; res.Pe2 = b * o.f2 * o.f2;
    if (a < 0 || b < 0) { res.ok = false; res.msg = 'Data inconsistent: ' + (a < 0 ? 'hysteresis term a' : 'eddy term b') + ' comes out negative (P/f must rise linearly with f).'; }
    return res;
  };

  /* =====================================================================
     4. Autotransformer from a two-winding transformer
     ===================================================================== */
  // o: {V1, V2, S (VA), common: 'V1'|'V2' (which winding is the common one), dir: 'down'|'up', x, pf, Pcore, Pcu (W at rated)}
  E.auto = function (o) {
    if (![o.V1, o.V2, o.S].every(v => fin(v) && v > 0)) return { ok: false, msg: 'V1, V2 and S must be positive.' };
    const Vc = o.common === 'V2' ? o.V2 : o.V1, Vs = o.common === 'V2' ? o.V1 : o.V2;
    const VH = Vc + Vs, VL = Vc, Sa = o.S * VH / (VH - VL);
    const IH = Sa / VH, IL = Sa / VL, Icom = IL - IH, Iser = IH;
    const r = { ok: true, msg: '', Vc, Vs, VH, VL, ratio: VH / VL, Sa, Sw: o.S, IH, IL, Icom, Iser, IserRated: o.S / Vs, IcomRated: o.S / Vc,
      Strans: (VH - VL) * IH, Scond: VL * IH, gain: Sa / o.S };
    const x = fin(o.x) && o.x >= 0 ? o.x : 1, pf = fin(o.pf) ? Math.min(1, Math.max(0, o.pf)) : 1, Pc = o.Pcore > 0 ? o.Pcore : 0, Pcu = o.Pcu > 0 ? o.Pcu : 0;
    const loss = Pc + x * x * Pcu, P2w = x * o.S * pf, Pa = x * Sa * pf;
    r.loss = loss; r.P2w = P2w; r.Pauto = Pa;
    r.eta2w = P2w > 0 ? P2w / (P2w + loss) : NaN; r.etaAuto = Pa > 0 ? Pa / (Pa + loss) : NaN;
    return r;
  };

  /* =====================================================================
     5. Three-phase transformer banks (three single-phase units)
     ===================================================================== */
  // o: {conn:'YY'|'YD'|'DY'|'DD'|'VV', S (VA, full 3-unit bank), VHL, VLL (line-line rms), conv:'ANSI'|'IEC11', Req, Xeq (Ω per unit, referred to its HV winding)}
  E.threePhase = function (o) {
    const conns = ['YY', 'YD', 'DY', 'DD', 'VV'];
    if (conns.indexOf(o.conn) < 0) return { ok: false, msg: 'Unknown connection.' };
    if (![o.S, o.VHL, o.VLL].every(v => fin(v) && v > 0)) return { ok: false, msg: 'S and both line voltages must be positive.' };
    const hY = o.conn[0] === 'Y', lY = o.conn === 'VV' ? false : o.conn[1] === 'Y';
    const VHw = hY ? o.VHL / S3 : o.VHL, VLw = lY ? o.VLL / S3 : o.VLL, a = VHw / VLw;
    const Sunit = o.S / 3, Sbank = o.conn === 'VV' ? S3 * Sunit : o.S;
    const IHL = Sbank / (S3 * o.VHL), ILL = Sbank / (S3 * o.VLL);
    const IHw = o.conn === 'VV' ? IHL : (hY ? IHL : IHL / S3), ILw = o.conn === 'VV' ? ILL : (lY ? ILL : ILL / S3);
    const mixed = o.conn === 'YD' || o.conn === 'DY';
    const shift = mixed ? (o.conv === 'IEC11' ? 30 : -30) : 0;           // angle of LV line voltage relative to HV (deg)
    const Zw = { re: o.Req || 0, im: o.Xeq || 0 };
    const kH = (hY || o.conn === 'VV') ? (hY ? 1 : 1 / 3) : 1 / 3;
    const ZYH = { re: Zw.re * kH, im: Zw.im * kH }, kr = (o.VLL / o.VHL) * (o.VLL / o.VHL), ZYL = { re: ZYH.re * kr, im: ZYH.im * kr };
    const Zbase = o.VHL * o.VHL / o.S, Zpu = { re: ZYH.re / Zbase, im: ZYH.im / Zbase };
    return { ok: true, msg: '', hY, lY, VHw, VLw, a, lineRatio: o.VHL / o.VLL, Sunit, Sbank, IHL, ILL, IHw, ILw, shift, ZYH, ZYL, Zbase, Zpu,
      vvFraction: o.conn === 'VV' ? Sbank / o.S : 1, vvUtil: o.conn === 'VV' ? Sbank / (2 * Sunit) : 1 };
  };

  /* =====================================================================
     6. All-day (energy) efficiency over a 24-h profile
     ===================================================================== */
  // rows "hours kW pf"
  E.parseProfile = function (text) {
    const rows = []; let bad = null;
    String(text == null ? '' : text).split(/\r?\n|;/).forEach(line => {
      line = line.replace(/#.*$/, '').trim(); if (!line) return;
      const p = line.split(/[\s,]+/).map(Number);
      if (p.length >= 2 && p.slice(0, Math.min(3, p.length)).every(fin)) rows.push({ h: p[0], kW: p[1], pf: p.length >= 3 ? p[2] : 1 }); else if (bad === null) bad = line;
    });
    if (bad !== null) return { ok: false, msg: 'Cannot read row "' + bad + '" (expected: hours kW pf).' };
    if (!rows.length) return { ok: false, msg: 'Enter at least one row: hours kW pf.' };
    for (const r of rows) {
      if (!(r.h > 0)) return { ok: false, msg: 'Hours must be > 0.' };
      if (!(r.kW >= 0)) return { ok: false, msg: 'Load kW must be ≥ 0.' };
      if (!(r.pf > 0 && r.pf <= 1)) return { ok: false, msg: 'Power factor must be in (0, 1].' };
    }
    const tot = rows.reduce((a, r) => a + r.h, 0);
    if (tot > 24 + 1e-9) return { ok: false, msg: 'Hours add up to ' + fx(tot) + ' h (> 24 h).' };
    return { ok: true, msg: '', rows, hours: tot };
  };
  // o: {S (kVA), Pcore (kW), Pcu (kW at full load), rows}
  E.allDay = function (o) {
    if (![o.S, o.Pcu].every(v => fin(v) && v > 0) || !(o.Pcore >= 0)) return { ok: false, msg: 'S and P_cu,FL must be > 0, P_core ≥ 0.' };
    const rows = o.rows.map(r => { const x = r.kW / r.pf / o.S; return Object.assign({}, r, { x, Eout: r.kW * r.h, Ecu: x * x * o.Pcu * r.h }); });
    const hours = rows.reduce((a, r) => a + r.h, 0);
    const Eout = rows.reduce((a, r) => a + r.Eout, 0), Ecu = rows.reduce((a, r) => a + r.Ecu, 0), Ecore = o.Pcore * 24;
    const eta = Eout > 0 ? Eout / (Eout + Ecu + Ecore) : NaN;
    return { ok: true, msg: '', rows, hours, idle: 24 - hours, Eout, Ecu, Ecore, Eloss: Ecu + Ecore, eta, overload: rows.some(r => r.x > 1 + 1e-9) };
  };
  // max-efficiency load fraction x* = √(P_core/P_cu,FL)
  E.maxEffPoint = function (S, Pcore, Pcu, pf) {
    if (!(S > 0) || !(Pcu > 0) || !(Pcore >= 0)) return { x: NaN, eta: NaN };
    const x = Math.sqrt(Pcore / Pcu), Po = x * S * pf; return { x, Skva: x * S, eta: Po > 0 ? Po / (Po + 2 * Pcore) : NaN };
  };
  E.etaAt = (S, Pcore, Pcu, x, pf) => { const Po = x * S * pf; return Po > 0 ? Po / (Po + Pcore + x * x * Pcu) : NaN; };

  /* =====================================================================
     Tests (expected values derived by hand / closed forms, see comments)
     ===================================================================== */

  FSP.math.edm2._shared = { MU0, NAMES, PI, S3, SH, arrowOf, cap, dropLines, fin, fx };
})();
