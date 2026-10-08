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
  FSP.registerTests('edm2', t => {
    const tbl = E.parseTable(E.tableToText(E.DEFAULT_TABLE)), mat = E.makeMaterial({ mode: 'table', H: tbl.H, B: tbl.B });
    t.check('past-paper B–H table parses (7 points)', tbl.ok && tbl.H.length === 7);
    t.check('parseTable rejects non-monotone B, duplicate H, junk, too few points',
      !E.parseTable('100 0.5\n200 0.4').ok && !E.parseTable('100 0.5;100 0.6;300 0.9').ok && !E.parseTable('100 abc').ok && !E.parseTable('100 0.5').ok);
    // hand: B=0.9 lies between (350,0.85) and (500,1.0): H = 350 + 0.05/0.15*150 = 400
    t.check('linear interpolation: H(0.75)=286 (node), H(0.9)=400, B(600)=1.08−… = 1.0457 (hand)',
      t.near(mat.Hof(0.75), 286, 1e-9) && t.near(mat.Hof(0.9), 400, 1e-9) && t.near(mat.Bof(600), 1.0 + 100 / 175 * 0.08, 1e-12));
    // hand: below the first point, line from origin: H(0.2) = 125*0.2/0.4 = 62.5
    t.check('below first point: straight line from origin, H(0.2)=62.5; odd symmetry', t.near(mat.Hof(0.2), 62.5, 1e-12) && t.near(mat.Hof(-0.9), -400, 1e-9));
    t.check('beyond table flagged as extrapolation', mat.outside(1.2) && !mat.outside(1.175));

    // ---- three-limb worked example (default state). Coil on centre limb N=500 (↑); gap 0.2 mm in the right limb;
    // depth 2 cm; widths L 2 cm, C 4 cm, R 2 cm → A_L=A_R=4 cm², A_C=8 cm²; l_L=l_R=0.30 m, l_C=0.10 m.
    // Given B_R = 0.75 T (↓) → Φ_R = 3.0e-4 Wb.
    // Right limb: H = 286 A/m (table); iron length 0.30−0.0002 = 0.2998 m → 85.7428 At;
    //   gap H_g = 0.75/μ0 = 596 831 A/m; F_g = 596831×2e-4 = 119.366 At → D_R = 205.109 At.
    // Left limb (no coil) in parallel: H_L = 205.109/0.30 = 683.697 A/m → B_L = 1.08 + 8.697/325×0.07 = 1.081873 T
    //   → Φ_L = 1.081873×4e-4 = 4.32749e-4 Wb.
    // Flux law: Φ_C = 4.32749e-4 + 3e-4 = 7.32749e-4 Wb → B_C = 0.915937 T → H_C = 350 + 0.065937/0.15×150 = 415.937 A/m
    //   → F_C = 41.5937 At.  NI = 41.5937 + 205.109 = 246.703 At → I = 0.493405 A ≈ 0.493 A.
    const lim = (on, l, w, lg, N, I, s) => E.makeLimb({ on, l, w, d: 0.02, lg, fringe: false, N, I, sense: s });
    const c0 = { mat, limbs: [lim(true, 0.3, 0.02, 0, 0, 0, 1), lim(true, 0.1, 0.04, 0, 500, 0, 1), lim(true, 0.3, 0.02, 2e-4, 0, 0, 1)] };
    const r0 = E.mcSolveForCurrent(c0, 2, -0.75 * 4e-4, 1);
    t.check('three-limb past-paper chain: I = 0.493 A (hand 0.493405)', r0.ok && t.rel(r0.I, 0.493405328, 1e-6) && r0.I.toPrecision(3) === '0.493', r0.ok ? 'I=' + r0.I : r0.msg);
    t.check('three-limb chain intermediates: B_L=1.08187 T, B_C=0.915937 T, H_C=415.937 A/m, U=205.109 At',
      t.rel(Math.abs(r0.drops[0].B), 1.0818731, 1e-6) && t.rel(r0.drops[1].B, 0.9159366, 1e-6) && t.rel(r0.drops[1].H, 415.93657, 1e-6) && t.rel(r0.U, 205.10901, 1e-6));
    // round trip: put I back and solve for fluxes (Newton) → B_R = 0.75 T ↓
    const c0b = { mat, limbs: [c0.limbs[0], Object.assign({}, c0.limbs[1], { I: r0.I }), c0.limbs[2]] };
    const rb = E.mcSolveCurrents(c0b);
    t.check('given currents (Newton) recovers B_R = −0.75 T, converged', rb.ok && rb.converged && rb.method.indexOf('Newton') === 0 && t.rel(rb.drops[2].B, -0.75, 1e-8), rb.ok ? rb.method + ' iter=' + rb.iter : rb.msg);
    t.check('solution satisfies flux law (ΣΦ≈0) and MMF law (same U on every limb)', Math.abs(rb.kcl) < 1e-15 && rb.mmfResid < 1e-7);
    // symmetric: no gap, given B_C = 1.0 T ↑ (Φ_C = 8e-4) → outer limbs 4e-4 each → B = 1.0 T, H = 500 everywhere.
    // NI = 500×0.1 + 500×0.3 = 200 At → I = 0.4 A
    const cs = { mat, limbs: [lim(true, 0.3, 0.02, 0, 0, 0, 1), lim(true, 0.1, 0.04, 0, 500, 0, 1), lim(true, 0.3, 0.02, 0, 0, 0, 1)] };
    const rs = E.mcSolveForCurrent(cs, 1, 8e-4, 1);
    t.check('symmetric core, B_C = 1.0 T: I = 0.400 A (hand: 200 At / 500)', rs.ok && t.rel(rs.I, 0.4, 1e-12) && rs.method.indexOf('symmetry') > 0);
    // asymmetric (gap in right limb) with given B in centre: x = g → U by bisection; verify by Newton round trip
    const ra = E.mcSolveForCurrent(c0, 1, 8e-4, 1), rra = ra.ok && E.mcSolveCurrents({ mat, limbs: [c0.limbs[0], Object.assign({}, c0.limbs[1], { I: ra.I }), c0.limbs[2]] });
    t.check('given B_C with asymmetric outer limbs: current reproduces Φ_C under the forward solve', ra.ok && rra.ok && t.rel(rra.phi[1], 8e-4, 1e-8), ra.ok ? 'I=' + ra.I : ra.msg);
    // coil on left limb, flux given in right limb: other branch of the chain
    const cl = { mat, limbs: [lim(true, 0.3, 0.02, 0, 400, 0, 1), lim(true, 0.1, 0.04, 0, 0, 0, 1), lim(true, 0.3, 0.02, 2e-4, 0, 0, 1)] };
    const rl = E.mcSolveForCurrent(cl, 2, -2e-4, 0), rlb = rl.ok && E.mcSolveCurrents({ mat, limbs: [Object.assign({}, cl.limbs[0], { I: rl.I }), cl.limbs[1], cl.limbs[2]] });
    t.check('coil on left, given Φ_R: forward solve reproduces Φ_R', rl.ok && rlb.ok && t.rel(rlb.phi[2], -2e-4, 1e-8));
    // two coils, Newton vs bisection-on-U (independent algorithm)
    const c2c = { mat, limbs: [lim(true, 0.3, 0.02, 0, 300, 0.6, 1), lim(true, 0.1, 0.04, 0, 0, 0, 1), lim(true, 0.25, 0.03, 3e-4, 200, 1.5, -1)] };
    const n2 = E.mcSolveCurrents(c2c), b2 = E.mcSolveBisect(c2c);
    t.check('two coils: Newton fluxes == bisection-on-U fluxes (1e-9 rel)', n2.ok && b2.ok && [0, 1, 2].every(k => t.near(n2.phi[k], b2.phi[k], 1e-9 * 1e-3)), n2.method);
    // linear (μr=1000) three-limb vs closed-form reluctance network: Φ_C = NI/(R_C + R_L‖R_R)
    // R_L = 0.3/(μ0·1000·4e-4), R_C = 0.1/(μ0·1000·8e-4), R_R = 0.2998/(μ0·1000·4e-4) + 2e-4/(μ0·4e-4); python: Φ_C = 1.058346e-3 Wb
    const lin = E.makeMaterial({ mode: 'const', mur: 1000 });
    const cL = { mat: lin, limbs: [c0.limbs[0], Object.assign({}, c0.limbs[1], { I: 1 }), c0.limbs[2]] }, rL = E.mcSolveCurrents(cL);
    t.check('constant μr: Φ_C = NI/(R_C + R_L‖R_R) = 1.058346 mWb, Φ_L/Φ_R = R_R/R_L', rL.ok && t.rel(rL.phi[1], 1.058346045e-3, 1e-8) && t.rel(rL.phi[0] / rL.phi[2], 0.6613670334 / 0.3969790117, 1e-8));
    // Chapman Ex. 1-2 style 2-limb series: μr=4000, iron 0.40 m, A=12 cm², gap 0.05 cm, N=400, B=0.5 T, no fringing.
    // R_c = 0.4/(4000·μ0·12e-4) = 66 315 A/Wb; R_g = 5e-4/(μ0·12e-4) = 331 573 A/Wb; I = 0.5·12e-4·397 887/400 = 0.59683 A
    const mu4k = E.makeMaterial({ mode: 'const', mur: 4000 });
    const cS = { mat: mu4k, limbs: [E.makeLimb({ on: true, l: 0.2, w: 0.03, d: 0.04, lg: 0, N: 400, I: 0, sense: 1 }), E.makeLimb({ on: false, l: 0.1, w: 0.03, d: 0.04, lg: 0 }), E.makeLimb({ on: true, l: 0.2005, w: 0.03, d: 0.04, lg: 5e-4, N: 0, I: 0, sense: 1 })] };
    const rS = E.mcSolveForCurrent(cS, 2, -0.5 * 12e-4, 0);
    t.check('2-limb series circuit with gap (textbook): I = 0.5968 A', rS.ok && t.rel(rS.I, 0.5968310366, 1e-8), rS.ok ? 'I=' + rS.I : rS.msg);
    // fringing: a=3 cm, b=4 cm, lg=1 mm → A_g = 3.1×4.1 cm²;  R_g = 1e-3/(μ0·0.031·0.041) = 626 101 A/Wb
    const Lf = E.makeLimb({ on: true, l: 0.2, w: 0.03, d: 0.04, lg: 1e-3, fringe: true, N: 0, I: 0, sense: 1 });
    t.check('fringing: A_g = (a+l_g)(b+l_g) = 12.71 cm², R_g = 626 101 A/Wb', t.rel(Lf.Ag, 0.031 * 0.041, 1e-12) && t.rel(E.drop(Lf, mu4k, 1e-4).Fg / 1e-4, 626101.27, 1e-7));
    t.check('invalid circuits rejected (one limb, gap ≥ length, unknown coil N=0)',
      !E.mcSolveCurrents({ mat, limbs: [c0.limbs[0], Object.assign({}, c0.limbs[1], { on: false }), Object.assign({}, c0.limbs[2], { on: false })] }).ok &&
      !E.mcSolveCurrents({ mat, limbs: [E.makeLimb({ on: true, l: 0.001, w: 0.02, d: 0.02, lg: 0.002 }), c0.limbs[1], c0.limbs[2]] }).ok &&
      !E.mcSolveForCurrent(c0, 2, -3e-4, 0).ok);

    // ---- inductance & energy. N=200, lc=0.2, A=4 cm², μr=2000, lg=1 mm:
    // R_c = 0.2/(2000·μ0·4e-4) = 198 944; R_g = 1e-3/(μ0·4e-4) = 1 989 437; L = 200²/2 188 381 = 0.0182784 H
    const mu2k = E.makeMaterial({ mode: 'const', mur: 2000 }), kc = { N: 200, lc: 0.2, A: 4e-4, lg: 1e-3, Ag: 4e-4, mat: mu2k }, ka = E.coreAnalyze(kc, 3);
    t.check('gapped core L = N²/(R_c+R_g) = 18.278 mH (hand)', t.rel(ka.Lsec, 0.0182783573, 1e-8));
    t.check('linear: W = ½LI², W_gap/W_core = R_g/R_c = 10, co-energy = energy', t.rel(ka.W, 0.5 * ka.Lsec * 9, 1e-10) && t.rel(ka.Wgap / ka.Wcore, 10, 1e-10) && t.rel(ka.Wco, ka.W, 1e-10));
    const kt = { N: 300, lc: 0.3, A: 4e-4, lg: 2e-4, Ag: 4e-4, mat }, kta = E.coreAnalyze(kt, 1.5);
    t.check('B–H table core: field energy (exact trapezoids) == ∫ i dλ (Simpson); W + W\' = λi', kta.B > 0.9 && t.rel(kta.W, E.energyAtLambda(kt, kta.lam, 2000), 1e-5) && t.rel(kta.W + kta.Wco, kta.lam * 1.5, 1e-12) && t.rel(kta.Wco, E.coEnergy(kt, 1.5, 2000), 1e-4), 'B=' + kta.B.toFixed(4));
    // ---- force. Ideal iron (μr→∞), N=1000, i=1 A, 2 gaps of 1 mm, A=4 cm²: B = μ0·1000/2e-3 = 0.6283 T;
    // F = 2·B²A/(2μ0) = 125.66 N (attractive) = ½ i² N² μ0 A/(2 x²)
    const rel1 = { N: 1000, lc: 0.2, A: 4e-4, ng: 2, mat: E.makeMaterial({ mode: 'const', mur: 1e12 }) };
    t.check('relay (ideal iron): F = ½ i² dL/dx = −125.66 N = −2·B²A/(2μ0)', t.rel(E.relayForce(rel1, 1, 1e-3), -0.5 * 1e6 * MU0 * 4e-4 / (2 * 1e-6), 1e-9) && t.rel(E.relayForce(rel1, 1, 1e-3), -2 * Math.pow(MU0 * 1000 / 2e-3, 2) * 4e-4 / (2 * MU0), 1e-9));
    const rel2 = { N: 1000, lc: 0.2, A: 4e-4, ng: 2, mat: mu2k };
    t.check('linear: ½ i² dL/dx == finite difference of co-energy ∂W\'/∂x (1e-5)', t.rel(E.relayForce(rel2, 0.8, 1.5e-3), E.relayForceFDco(rel2, 0.8, 1.5e-3), 1e-5));
    const rel3 = { N: 1000, lc: 0.2, A: 4e-4, ng: 2, mat };
    // Three routes to the same force. Exact: F = −n_g·B²A/(2μ0), B = Φ/A (flux linkage held fixed, so only the gap energy moves).
    // The two finite-difference routes are checked against it. The co-energy route integrates λ(i') (piecewise linear: kinks at the
    // table breakpoints) by Simpson, whose error at a kink is O(Δi²·Δslope); dividing by the small difference 2h amplifies it, so a
    // 0.1 % step and n=400 gave 5e-4 error. With h = 1 % of x and n=2000 the observed error is 3e-5, hence the 1e-4 tolerance.
    const fex = E.relayForce(rel3, 0.6, 0.3e-3), Bx = E.coreFlux(E.relayCore(rel3, 0.3e-3), 0.6) / 4e-4;
    const fco = E.relayForceFDco(rel3, 0.6, 0.3e-3, 3e-6, 2000), fen = E.relayForceFDenergy(rel3, 0.6, 0.3e-3);
    t.check('nonlinear (table): F = −n_g B²A/(2μ0) == ∂W\'/∂x|i == −∂W/∂x|λ (three routes, 1e-4)', fex < 0 && t.rel(fex, -2 * Bx * Bx * 4e-4 / (2 * MU0), 1e-12) && t.rel(fco, fex, 1e-4) && t.rel(fen, fex, 1e-4), fex.toFixed(4) + ' / ' + fco.toFixed(4) + ' / ' + fen.toFixed(4));
    const ip = E.relayPullIn(rel2, 2e-3, 30), ipn = E.relayPullIn(rel3, 2e-3, 30);
    t.check('pull-in current: |F(i_pull, x0)| = F_s (linear closed form and nonlinear bisection)', t.rel(Math.abs(E.relayForce(rel2, ip, 2e-3)), 30, 1e-9) && t.rel(Math.abs(E.relayForce(rel3, ipn, 2e-3)), 30, 1e-6));
    const th = 0.4, h = 1e-6, fdT = 0.5 * 4 * (E.rotL(0.2, 0.05, th + h) - E.rotL(0.2, 0.05, th - h)) / (2 * h);
    t.check('rotary: T = ½ i² dL/dθ = −½i²(Ld−Lq)sin2θ (finite-difference check); zero at aligned/unaligned', t.rel(E.rotTorque(0.2, 0.05, 2, th), fdT, 1e-7) && Math.abs(E.rotTorque(0.2, 0.05, 2, 0)) < 1e-15 && Math.abs(E.rotTorque(0.2, 0.05, 2, PI / 2)) < 1e-15);
    // 4.44 rule: 230 V, 50 Hz, 200 turns: Φmax = 230/(√2·π·50·200) = 5.17682 mWb; 4.44·f·N·Φ within 0.07 %
    const sf = E.sineFlux(230, 50, 200, 25e-4);
    t.check('E = 4.44 f N Φmax: Φmax = 5.1768 mWb, Bmax = 2.0707 T for 25 cm²', t.rel(sf.phim, 5.17681882e-3, 1e-8) && t.rel(4.44 * 50 * 200 * sf.phim, 230, 7e-4) && t.rel(sf.Bm, 5.17681882e-3 / 25e-4, 1e-8));
    // separation, hand: P1=100 W @50 Hz, P2=130 W @60 Hz: P/f = 2, 2.16667 → b = 0.0166667, a = 1.16667;
    // at 50 Hz: P_h = 58.333 W, P_e = 41.667 W
    const sp = E.separateLosses({ f1: 50, P1: 100, f2: 60, P2: 130, Bm: 1.2, n: 1.6, vol: 1e-3 });
    t.check('loss separation (hand): a = 1.16667 J, b = 0.0166667 J·s, P_h(50) = 58.333 W, P_e(50) = 41.667 W', sp.ok && t.rel(sp.a, 7 / 6, 1e-12) && t.rel(sp.b, 1 / 60, 1e-12) && t.rel(sp.Ph1, 175 / 3, 1e-12) && t.rel(sp.Pe1, 125 / 3, 1e-12));
    const kh0 = 120, ke0 = 0.05, gen = f => E.coreLoss({ kh: kh0, ke: ke0, n: 1.7, f, Bm: 1.3, vol: 2e-3 }).P;
    const sp2 = E.separateLosses({ f1: 25, P1: gen(25), f2: 70, P2: gen(70), Bm: 1.3, n: 1.7, vol: 2e-3 });
    t.check('loss separation recovers k_h, k_e from synthetic data', sp2.ok && t.rel(sp2.kh, kh0, 1e-10) && t.rel(sp2.ke, ke0, 1e-10));
    t.check('loss separation: equal frequencies / negative term rejected', !E.separateLosses({ f1: 50, P1: 100, f2: 50, P2: 120, Bm: 1, n: 1.6, vol: 1 }).ok && !E.separateLosses({ f1: 50, P1: 100, f2: 60, P2: 110, Bm: 1, n: 1.6, vol: 1 }).ok);

    // ---- autotransformer: 5 kVA 220/110 V → 330/220 V step-down (common = 220 V winding):
    // S_auto = 5·330/110 = 15 kVA; I_H = 15000/330 = 45.4545 A (= 5000/110, series winding rated);
    // I_L = 15000/220 = 68.1818 A; common winding I_L − I_H = 22.7273 A (= 5000/220); conducted 10 kVA, transformed 5 kVA
    const au = E.auto({ V1: 220, V2: 110, S: 5000, common: 'V1', dir: 'down', x: 1, pf: 0.8, Pcore: 50, Pcu: 100 });
    t.check('auto 5 kVA 220/110 → 330/220: 15 kVA, I_H = 45.45 A, I_L = 68.18 A, I_common = 22.73 A',
      au.ok && au.VH === 330 && au.VL === 220 && t.rel(au.Sa, 15000, 1e-12) && t.rel(au.IH, 5000 / 110, 1e-12) && t.rel(au.IL, 15000 / 220, 1e-12) && t.rel(au.Icom, 5000 / 220, 1e-12));
    t.check('auto power split: transformed 5 kVA + conducted 10 kVA', t.rel(au.Strans, 5000, 1e-12) && t.rel(au.Scond, 10000, 1e-12));
    // common = 110 V winding → 330/110 V: S = 5·330/220 = 7.5 kVA
    const au2 = E.auto({ V1: 220, V2: 110, S: 5000, common: 'V2', dir: 'up' });
    t.check('auto with the 110 V winding common: 330/110 V, 7.5 kVA', au2.ok && au2.VL === 110 && t.rel(au2.Sa, 7500, 1e-12));
    // efficiencies (hand): 2-w: 4000/(4000+150) = 0.963855; auto: 12000/12150 = 0.987654
    t.check('auto efficiency: 96.386 % (2-winding) → 98.765 % (auto), same losses', t.rel(au.eta2w, 4000 / 4150, 1e-12) && t.rel(au.etaAuto, 12000 / 12150, 1e-12));
    t.check('auto rejects zero/negative rating', !E.auto({ V1: 0, V2: 110, S: 5000 }).ok);

    // ---- three-phase: 150 kVA bank Y-Δ, 11 kV / 415 V.
    // HV winding = 11000/√3 = 6350.85 V, LV winding = 415 V → a = 15.3033; line ratio 11000/415 = √3·a.
    // I_HL = 150e3/(√3·11e3) = 7.8730 A (= HV winding); I_LL = 150e3/(√3·415) = 208.68 A; LV winding 120.482 A = 50e3/415
    const yd = E.threePhase({ conn: 'YD', S: 150e3, VHL: 11000, VLL: 415, conv: 'ANSI', Req: 12, Xeq: 30 });
    t.check('Y-Δ: line-voltage ratio = √3·a (a = winding ratio)', yd.ok && t.rel(yd.lineRatio, S3 * yd.a, 1e-12) && t.rel(yd.a, 11000 / S3 / 415, 1e-12));
    t.check('Y-Δ currents: I_HL = 7.873 A, I_LL = 208.68 A, LV winding 120.48 A (= 50 kVA/415 V)', t.rel(yd.IHL, 150e3 / (S3 * 11000), 1e-12) && t.rel(yd.ILw, 50e3 / 415, 1e-12) && t.rel(yd.IHw, yd.IHL, 1e-12));
    t.check('Y-Δ phase shift: −30° (ANSI, LV lags) / +30° (IEC clock 11)', yd.shift === -30 && E.threePhase({ conn: 'YD', S: 150e3, VHL: 11000, VLL: 415, conv: 'IEC11' }).shift === 30);
    const dy = E.threePhase({ conn: 'DY', S: 150e3, VHL: 11000, VLL: 415 }), dd = E.threePhase({ conn: 'DD', S: 150e3, VHL: 11000, VLL: 415 });
    t.check('Δ-Y ratio = a/√3; Δ-Δ and Y-Y: ratio = a, no shift', t.rel(dy.lineRatio, dy.a / S3, 1e-12) && t.rel(dd.lineRatio, dd.a, 1e-12) && dd.shift === 0 && E.threePhase({ conn: 'YY', S: 1, VHL: 11000, VLL: 415 }).shift === 0);
    const vv = E.threePhase({ conn: 'VV', S: 150e3, VHL: 11000, VLL: 415 });
    t.check('open-Δ (V-V): capacity = √3·S_unit = 57.7 % of the Δ-Δ bank, 86.6 % unit utilisation', t.rel(vv.vvFraction, 1 / S3, 1e-12) && t.near(vv.vvFraction * 100, 57.735, 1e-3) && t.rel(vv.vvUtil, S3 / 2, 1e-12));
    // per-phase equivalent: Y on HV → Z_Y(HV) = Z_w = 12 + j30; referred to LV × (415/11000)² = 0.0170805 + j0.0427012
    t.check('per-phase equivalent: Y-HV → Z_Y = Z_w; referred to LV by (V_LL/V_HL)²; Δ-HV → Z_w/3',
      t.rel(yd.ZYH.re, 12, 1e-12) && t.rel(yd.ZYL.im, 30 * (415 / 11000) ** 2, 1e-12) && t.rel(E.threePhase({ conn: 'DY', S: 150e3, VHL: 11000, VLL: 415, Req: 12, Xeq: 30 }).ZYH.re, 4, 1e-12));

    // ---- all-day efficiency (hand): 100 kVA, P_core 1 kW, P_cu,FL 1.5 kW; 6 h 80 kW pf 0.8 (x=1), 10 h 40 kW pf 1 (x=0.4), 8 h idle.
    // E_out = 480+400 = 880 kWh; E_cu = 6·1.5 + 10·1.5·0.16 = 11.4 kWh; E_core = 24 kWh → η = 880/915.4 = 96.133 %
    const pr = E.parseProfile('6 80 0.8\n10 40 1\n8 0 1'), ad = pr.ok && E.allDay({ S: 100, Pcore: 1, Pcu: 1.5, rows: pr.rows });
    t.check('all-day efficiency = 880/915.4 = 96.133 % (hand)', ad && ad.ok && t.rel(ad.eta, 880 / 915.4, 1e-12) && t.rel(ad.Ecu, 11.4, 1e-12) && t.rel(ad.Ecore, 24, 1e-12));
    const me = E.maxEffPoint(100, 1, 1.5, 0.8);
    t.check('max-efficiency load fraction x = √(P_core/P_cu) = 0.8165 and η(x) is the maximum', t.rel(me.x, 0.8164965809, 1e-9) && t.rel(me.eta, E.etaAt(100, 1, 1.5, me.x, 0.8), 1e-12) && E.etaAt(100, 1, 1.5, me.x * 1.05, 0.8) < me.eta && E.etaAt(100, 1, 1.5, me.x * 0.95, 0.8) < me.eta);
    t.check('profile validation: > 24 h, pf > 1, junk rejected', !E.parseProfile('20 10 1\n5 10 1').ok && !E.parseProfile('2 10 1.2').ok && !E.parseProfile('two 10').ok);
  });

  /* =====================================================================
     UI
     ===================================================================== */
  const css = n => { try { return getComputedStyle(document.documentElement).getPropertyValue(n).trim() || null; } catch (e) { return null; } };
  function theme() {
    return { bg: css('--panel2') || '#0f141f', panel: css('--panel') || '#121826', text: css('--text') || '#e6edf3', muted: css('--muted') || '#8b98ab', border: css('--border') || '#243047',
      a: css('--accent2') || '#f9a03f', b: css('--ok') || '#2ecc71', c: css('--accent') || '#4cc9f0', pink: css('--pink') || '#f72585', warn: css('--warn') || '#ffb347', bad: css('--bad') || '#ff5c5c', mono: css('--mono') || 'monospace' };
  }
  const f2 = (x, d) => (fin(x) ? x.toFixed(d === undefined ? 2 : d) : '—');
  const fE = (x, u, d) => (fin(x) ? FSP.fmtEng(x, u, d || 4) : '—');
  const tickStr = v => (Math.abs(v) < 1e-12 ? '0' : Math.abs(v) >= 1e5 || Math.abs(v) < 1e-3 ? v.toExponential(1) : String(+v.toPrecision(4)));
  function niceTicks(lo, hi, n) {
    const span = hi - lo; if (!(span > 0)) return [lo];
    const raw = span / n, mag = Math.pow(10, Math.floor(Math.log10(raw))), f = raw / mag, step = (f < 1.5 ? 1 : f < 3 ? 2 : f < 7 ? 5 : 10) * mag, out = [];
    for (let v = Math.ceil(lo / step - 1e-9) * step; v <= hi + step * 1e-9; v += step) out.push(Math.abs(v) < step * 1e-9 ? 0 : +v.toPrecision(12));
    return out;
  }
  function withAlpha(col, a) {
    const m = /^#([0-9a-f]{6})$/i.exec(col || ''); if (!m) return col;
    const n = parseInt(m[1], 16); return 'rgba(' + (n >> 16) + ',' + ((n >> 8) & 255) + ',' + (n & 255) + ',' + a + ')';
  }
  // generic plot. o: {xmin,xmax,ymin,ymax,y2min,y2max,xlabel,ylabel,y2label,series:[{x,y,color,width,dash,axis}],fills:[{x,y,color}],bars:[{x0,x1,y,color}],
  //   points:[{x,y,label,color,below,axis}],vlines:[{x,label,color}],hlines:[{y,label,color}],legend:[{text,color}]}
  function drawPlot(c, o) {
    const g = c.prep(), ctx = g.ctx, w = g.w, h = g.h, T = theme(), y2 = o.y2min !== undefined;
    const L = 52, R = y2 ? 42 : 10, Tp = 10, B = 32, pw = w - L - R, ph = h - Tp - B;
    ctx.clearRect(0, 0, w, h); ctx.fillStyle = T.bg; ctx.fillRect(0, 0, w, h);
    if (pw < 20 || ph < 20 || !(o.xmax > o.xmin) || !(o.ymax > o.ymin)) return null;
    const X = x => L + (x - o.xmin) / (o.xmax - o.xmin) * pw, Y = y => Tp + (1 - (y - o.ymin) / (o.ymax - o.ymin)) * ph;
    const Y2 = y => Tp + (1 - (y - o.y2min) / (o.y2max - o.y2min)) * ph;
    ctx.font = '11px ' + T.mono; ctx.lineWidth = 1; ctx.textBaseline = 'middle'; ctx.textAlign = 'right';
    niceTicks(o.ymin, o.ymax, 5).forEach(v => { const yy = Y(v); ctx.strokeStyle = T.border; ctx.beginPath(); ctx.moveTo(L, yy); ctx.lineTo(L + pw, yy); ctx.stroke(); ctx.fillStyle = T.muted; ctx.fillText(tickStr(v), L - 4, yy); });
    if (y2) { ctx.textAlign = 'left'; niceTicks(o.y2min, o.y2max, 5).forEach(v => { ctx.fillStyle = T.muted; ctx.fillText(tickStr(v), L + pw + 4, Y2(v)); }); }
    ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    niceTicks(o.xmin, o.xmax, Math.max(3, Math.floor(pw / 70))).forEach(v => { const xx = X(v); ctx.strokeStyle = T.border; ctx.beginPath(); ctx.moveTo(xx, Tp); ctx.lineTo(xx, Tp + ph); ctx.stroke(); ctx.fillStyle = T.muted; ctx.fillText(tickStr(v), xx, Tp + ph + 3); });
    ctx.strokeStyle = T.muted; if (o.ymin < 0 && o.ymax > 0) { ctx.beginPath(); ctx.moveTo(L, Y(0)); ctx.lineTo(L + pw, Y(0)); ctx.stroke(); }
    ctx.fillStyle = T.muted; ctx.textAlign = 'center'; ctx.textBaseline = 'bottom'; ctx.fillText(o.xlabel || '', L + pw / 2, h - 1);
    ctx.save(); ctx.translate(10, Tp + ph / 2); ctx.rotate(-PI / 2); ctx.textBaseline = 'top'; ctx.fillText(o.ylabel || '', 0, 0); ctx.restore();
    if (y2) { ctx.save(); ctx.translate(w - 3, Tp + ph / 2); ctx.rotate(-PI / 2); ctx.textBaseline = 'bottom'; ctx.fillText(o.y2label || '', 0, 0); ctx.restore(); }
    ctx.save(); ctx.beginPath(); ctx.rect(L, Tp, pw, ph); ctx.clip();
    (o.bars || []).forEach(b => { if (!fin(b.y)) return; ctx.fillStyle = b.color; const y0 = Y(Math.max(o.ymin, 0)), y1 = Y(b.y); ctx.fillRect(X(b.x0) + 1, Math.min(y0, y1), Math.max(1, X(b.x1) - X(b.x0) - 2), Math.abs(y1 - y0)); });
    (o.fills || []).forEach(f => { ctx.fillStyle = f.color; ctx.beginPath(); f.x.forEach((xv, i) => { const px = X(xv), py = Y(f.y[i]); if (i) ctx.lineTo(px, py); else ctx.moveTo(px, py); }); ctx.closePath(); ctx.fill(); });
    (o.vlines || []).forEach(v => { if (!fin(v.x)) return; ctx.strokeStyle = v.color || T.muted; ctx.setLineDash([4, 4]); ctx.beginPath(); ctx.moveTo(X(v.x), Tp); ctx.lineTo(X(v.x), Tp + ph); ctx.stroke(); ctx.setLineDash([]); });
    (o.hlines || []).forEach(v => { if (!fin(v.y)) return; ctx.strokeStyle = v.color || T.muted; ctx.setLineDash([4, 4]); ctx.beginPath(); ctx.moveTo(L, Y(v.y)); ctx.lineTo(L + pw, Y(v.y)); ctx.stroke(); ctx.setLineDash([]); });
    (o.series || []).forEach(s => {
      const fy = s.axis === 'y2' ? Y2 : Y; ctx.strokeStyle = s.color; ctx.lineWidth = s.width || 2; ctx.setLineDash(s.dash || []); ctx.beginPath(); let pen = false;
      for (let i = 0; i < s.x.length; i++) { const xv = s.x[i], yv = s.y[i]; if (!fin(xv) || !fin(yv)) { pen = false; continue; } const px = X(xv), py = Math.max(-1e4, Math.min(1e4, fy(yv))); if (pen) ctx.lineTo(px, py); else { ctx.moveTo(px, py); pen = true; } }
      ctx.stroke(); ctx.setLineDash([]);
    });
    ctx.restore();
    ctx.font = '11px ' + T.mono;
    (o.vlines || []).forEach((v, i) => { if (!fin(v.x) || !v.label) return; const xx = X(v.x), right = xx < L + pw * 0.6; ctx.fillStyle = v.color || T.muted; ctx.textAlign = right ? 'left' : 'right'; ctx.textBaseline = 'top'; ctx.fillText(v.label, xx + (right ? 4 : -4), Tp + 3 + (v.row || 0) * 13); });
    (o.hlines || []).forEach(v => { if (!fin(v.y) || !v.label) return; ctx.fillStyle = v.color || T.muted; ctx.textAlign = 'left'; ctx.textBaseline = 'bottom'; ctx.fillText(v.label, L + 4, Y(v.y) - 2); });
    (o.points || []).forEach(p => {
      if (!fin(p.x) || !fin(p.y)) return; const px = X(p.x), py = (p.axis === 'y2' ? Y2 : Y)(p.y);
      if (px < L - 1 || px > L + pw + 1 || py < Tp - 1 || py > Tp + ph + 1) return;
      ctx.fillStyle = p.color; ctx.strokeStyle = T.bg; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.arc(px, py, 5, 0, 2 * PI); ctx.fill(); ctx.stroke();
      if (p.label) { const right = px < L + pw * 0.62; ctx.textAlign = right ? 'left' : 'right'; ctx.textBaseline = p.below ? 'top' : 'bottom'; ctx.fillStyle = T.text; ctx.fillText(p.label, px + (right ? 8 : -8), py + (p.below ? 6 : -6)); }
    });
    (o.legend || []).forEach((lg, i) => { ctx.fillStyle = lg.color; ctx.textAlign = 'right'; ctx.textBaseline = 'top'; ctx.fillText(lg.text, L + pw - 4, Tp + 4 + i * 14); });
    return { X, Y, L, Tp, pw, ph };
  }
  function hudSet(el, rows) {
    while (el.firstChild) el.removeChild(el.firstChild);
    rows.forEach(r => el.appendChild(FSP.ui.el('div', null, FSP.ui.el('span', { text: r[0] }), FSP.ui.el('span', { text: r[1] }))));
  }
  const wrapCanvas = (parent, label, height) => { const wrap = FSP.ui.el('div', { class: 'canvas-wrap' }); parent.appendChild(wrap); const c = FSP.ui.canvas(wrap, { height }); c.cv.setAttribute('role', 'img'); c.cv.setAttribute('aria-label', label); return c; };
  function layout(root) { const lay = FSP.ui.el('div', { class: 'layout' }); root.appendChild(lay); const ctl = FSP.ui.el('div', { class: 'controls' }), stage = FSP.ui.el('div', { class: 'stage' }); lay.appendChild(ctl); lay.appendChild(stage); return { ctl, stage }; }
  // defs: [{k,l,min,max,v,u,log,step}]
  function makeSliders(parent, defs, st, onChange) {
    defs.forEach(d => { st[d.k] = FSP.ui.slider(parent, { label: d.l, min: d.min, max: d.max, step: d.step, value: d.v, unit: d.u, log: d.log, onInput: () => { onChange(); FSP.state.touch(); } }); });
  }
  const sel = (parent, label, opts, v, onChange) => FSP.ui.select(parent, label, opts, v, () => { onChange(); FSP.state.touch(); });
  const rowOf = selEl => selEl.parentNode;
  function stateIO(st, sels, extra) {
    return {
      get() { const o = {}; Object.keys(st).forEach(k => { o[k] = st[k].get(); }); Object.keys(sels).forEach(k => { o[k] = sels[k].value; }); if (extra) Object.assign(o, extra.get()); return o; },
      set(o) {
        o = o || {};
        Object.keys(st).forEach(k => { if (o[k] !== undefined) { const v = parseFloat(o[k]); if (fin(v)) st[k].set(v, true); } });
        Object.keys(sels).forEach(k => { if (o[k] !== undefined && Array.prototype.some.call(sels[k].options, op => op.value === o[k])) sels[k].value = o[k]; });
        if (extra) extra.set(o);
      },
    };
  }
  function msgBox(parent) { const m = FSP.ui.el('div', { class: 'msg bad', hidden: '' }); parent.appendChild(m); return { el: m, show(t) { m.textContent = t; m.hidden = !t; } }; }
  function arrowHead(ctx, x0, y0, x1, y1, s) {
    const a = Math.atan2(y1 - y0, x1 - x0); ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1, y1); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(x1, y1); ctx.lineTo(x1 - s * Math.cos(a - 0.4), y1 - s * Math.sin(a - 0.4)); ctx.lineTo(x1 - s * Math.cos(a + 0.4), y1 - s * Math.sin(a + 0.4)); ctx.closePath(); ctx.fill();
  }

  /* ---------- shared B–H table (used by sub-topics 1 and 2) ---------- */
  const shared = { text: E.tableToText(E.DEFAULT_TABLE), parsed: E.parseTable(E.tableToText(E.DEFAULT_TABLE)), listeners: [] };
  shared.tableMat = () => (shared.parsed.ok ? E.makeMaterial({ mode: 'table', H: shared.parsed.H, B: shared.parsed.B }) : { ok: false });
  shared.setText = (t, silent) => { shared.text = String(t); shared.parsed = E.parseTable(shared.text); if (!silent) shared.listeners.forEach(f => f()); };

  /* ---------- 1. multi-limb magnetic circuit ---------- */
  function buildMC(root) {
    const st = {}, ui = {}, sels = {}, { ctl, stage } = layout(root);
    const fsP = FSP.ui.fieldset(ctl, 'Problem');
    sels.mode = sel(fsP, 'Solve for', [['I', 'coil current (given Φ or B)'], ['phi', 'fluxes (given currents)']], 'I', update);
    sels.g = sel(fsP, 'Given in', [['0', 'left limb'], ['1', 'centre limb'], ['2', 'right limb']], '2', update);
    sels.gq = sel(fsP, 'Given', [['B', 'flux density B'], ['phi', 'flux Φ']], 'B', update);
    makeSliders(fsP, [{ k: 'gB', l: 'B given', min: 0, max: 2.5, v: 0.75, step: 0.0001, u: 'T' }, { k: 'gPhi', l: 'Φ given', min: 0, max: 50, v: 0.3, step: 0.0001, u: 'mWb' }], st, update);
    sels.gd = sel(fsP, 'Direction', [['-1', '↓ down'], ['1', '↑ up']], '-1', update);
    sels.x = sel(fsP, 'Unknown coil', [['0', 'on left limb'], ['1', 'on centre limb'], ['2', 'on right limb']], '1', update);
    const fsM = FSP.ui.fieldset(ctl, 'Core material & geometry');
    sels.mat = sel(fsM, 'Material', [['table', 'B–H table (linear interp.)'], ['const', 'constant μr']], 'table', update);
    makeSliders(fsM, [{ k: 'mur', l: 'μr', min: 1, max: 100000, v: 1000, u: '', log: true }, { k: 'd', l: 'depth', min: 0.1, max: 50, v: 2, u: 'cm', log: true }], st, update);
    sels.fr = sel(fsM, 'Fringing', [['0', 'none (A_g = A)'], ['1', 'A_g = (a+l_g)(b+l_g)']], '0', update);
    ui.tblWrap = FSP.ui.el('div'); fsM.appendChild(ui.tblWrap);
    ui.tblWrap.appendChild(FSP.ui.el('div', { class: 'note', text: 'B–H table, one "H B" pair per row (H in A/m, B in T). (0,0) is implied; straight lines between points.' }));
    ui.tbl = FSP.ui.el('textarea', { rows: '8', spellcheck: 'false', 'aria-label': 'B–H table' }); ui.tbl.value = shared.text;
    ui.tbl.style.cssText = 'width:100%;box-sizing:border-box;background:var(--panel2);color:var(--text);border:1px solid var(--border);border-radius:6px;font:12px var(--mono);padding:4px';
    ui.tbl.addEventListener('change', () => { shared.setText(ui.tbl.value); update(); FSP.state.touch(); }); ui.tblWrap.appendChild(ui.tbl);
    const rr = FSP.ui.el('div', { class: 'row' }); ui.tblWrap.appendChild(rr);
    FSP.ui.button(rr, 'Past-paper table', () => { ui.tbl.value = E.tableToText(E.DEFAULT_TABLE); shared.setText(ui.tbl.value); update(); FSP.state.touch(); }).className = 'btn';
    ui.tblMsg = msgBox(ui.tblWrap);
    const limbDefs = [
      { on: '1', l: 0.3, w: 2, g: 0, N: 0, I: 0, s: '1' }, { on: '1', l: 0.1, w: 4, g: 0, N: 500, I: 0.5, s: '1' }, { on: '1', l: 0.3, w: 2, g: 0.2, N: 0, I: 0, s: '1' }];
    ui.limbRows = [];
    limbDefs.forEach((D, k) => {
      const fs = FSP.ui.fieldset(ctl, cap(NAMES[k]) + ' limb'), p = SH[k];
      sels['on' + p] = sel(fs, 'Present', [['1', 'yes'], ['0', 'no (removed)']], D.on, update);
      makeSliders(fs, [
        { k: 'l' + p, l: 'mean l', min: 0.01, max: 2, v: D.l, u: 'm', log: true }, { k: 'w' + p, l: 'width a', min: 0.1, max: 50, v: D.w, u: 'cm', log: true },
        { k: 'g' + p, l: 'air gap', min: 0, max: 10, v: D.g, step: 0.01, u: 'mm' }, { k: 'N' + p, l: 'turns N', min: 0, max: 5000, v: D.N, step: 1, u: '' },
        { k: 'I' + p, l: 'current I', min: 0, max: 100, v: D.I, step: 0.0001, u: 'A' }], st, update);
      sels['s' + p] = sel(fs, 'Coil drives', [['1', 'flux ↑'], ['-1', 'flux ↓']], D.s, update);
      ui.limbRows.push(fs);
    });
    ctl.appendChild(FSP.ui.el('div', { class: 'note', text: 'Mean length l = path along the limb between the two junctions (include its share of the yokes). Iron length = l − l_g. Area A = width a × depth b.' }));
    ui.msg = msgBox(ctl);
    ui.cv = wrapCanvas(stage, 'Three-limb magnetic circuit with flux arrows', 320); ui.cv.onResize(update);
    ui.hud = FSP.ui.el('div', { class: 'hud' }); stage.appendChild(ui.hud);
    ui.warn = FSP.ui.el('div', { class: 'msg warn', hidden: '' }); stage.appendChild(ui.warn);
    ui.work = FSP.ui.working(stage); ui.work.el.open = true;
    shared.listeners.push(() => { if (ui.tbl.value !== shared.text) ui.tbl.value = shared.text; });

    function circuit() {
      const fringe = sels.fr.value === '1', mat = sels.mat.value === 'table' ? shared.tableMat() : E.makeMaterial({ mode: 'const', mur: st.mur.get() });
      const limbs = SH.map(p => E.makeLimb({ on: sels['on' + p].value === '1', l: st['l' + p].get(), w: st['w' + p].get() / 100, d: st.d.get() / 100, lg: st['g' + p].get() / 1000, fringe, N: Math.round(st['N' + p].get()), I: st['I' + p].get(), sense: +sels['s' + p].value }));
      return { mat, limbs };
    }
    let last = null;
    function update() {
      if (!ui.cv) return;
      const modeI = sels.mode.value === 'I', g = +sels.g.value, x = +sels.x.value;
      [sels.g, sels.gq, sels.gd, sels.x].forEach(s => { rowOf(s).hidden = !modeI; });
      st.gB.el.hidden = !modeI || sels.gq.value !== 'B'; st.gPhi.el.hidden = !modeI || sels.gq.value !== 'phi';
      st.mur.el.hidden = sels.mat.value !== 'const'; ui.tblWrap.hidden = sels.mat.value !== 'table';
      ui.tblMsg.show(shared.parsed.ok ? '' : shared.parsed.msg);
      SH.forEach((p, k) => { const on = sels['on' + p].value === '1'; ['l', 'w', 'g', 'N', 'I'].forEach(q => { st[q + p].el.hidden = !on || (q === 'I' && modeI && k === x); }); rowOf(sels['s' + p]).hidden = !on; });
      const c = circuit(); let r;
      if (modeI) {
        const L = c.limbs[g], dir = +sels.gd.value, phiG = sels.gq.value === 'B' ? dir * st.gB.get() * L.A : dir * st.gPhi.get() * 1e-3;
        r = !(Math.abs(phiG) > 0) ? { ok: false, msg: 'Enter a non-zero given flux / flux density.' } : E.mcSolveForCurrent(c, g, phiG, x);
      } else r = E.mcSolveCurrents(c);
      last = { c, r };
      ui.msg.show(r.ok ? '' : r.msg);
      drawCore(ui.cv, c, r.ok ? r : null, modeI ? x : -1);
      if (!r.ok) { hudSet(ui.hud, []); ui.work.set(r.msg); ui.warn.hidden = true; return; }
      const rows = [];
      if (modeI) rows.push(['coil current I (' + NAMES[x] + ')', fx(r.I, 4) + ' A'], ['coil MMF N·I', fx(Math.abs(r.sNI), 5) + ' At']);
      r.P.forEach(k => { const d = r.drops[k]; rows.push(['Φ ' + NAMES[k], fx(Math.abs(d.phi) * 1e3, 4) + ' mWb ' + arrowOf(d.phi)], ['B / H ' + NAMES[k], fx(Math.abs(d.B), 4) + ' T / ' + fx(Math.abs(d.H), 4) + ' A/m']); });
      rows.push(['U (yoke-to-yoke MMF)', fx(r.U, 5) + ' At'], ['ΣΦ check', fx(r.kcl, 3) + ' Wb']);
      if (!modeI) rows.push(['method', r.method], ['iterations', String(r.iter)]);
      hudSet(ui.hud, rows);
      ui.warn.hidden = !r.outside; ui.warn.textContent = r.outside ? 'A flux density exceeds the last B–H table point: H was extrapolated along the last segment, so the answer is outside the given data.' : '';
      if (modeI) ui.work.set(r.steps);
      else {
        const lines = ['Convention: Φ and coil MMF positive ↑ in every limb.  For each limb: s·N·I − D(Φ) = U (common yoke-to-yoke MMF);  flux law ΣΦ = 0.',
          'Coil MMFs: ' + r.P.map(k => SH[k] + ': ' + fx(r.S[k]) + ' At').join(',  ')];
        if (r.oth && r.oth.length === 2) lines.push('Unknowns Φ_L, Φ_R with Φ_C = −(Φ_L + Φ_R). Window equations (Kirchhoff MMF law):', '  left window:  [S_L − D_L(Φ_L)] − [S_C − D_C(Φ_C)] = 0', '  right window: [S_R − D_R(Φ_R)] − [S_C − D_C(Φ_C)] = 0');
        else lines.push('Two limbs: single loop, Φ_a = −Φ_b, S_a − D_a(Φ_a) = S_b − D_b(Φ_b).');
        lines.push('D(Φ) is nonlinear (B–H table) → solved by ' + r.method + (r.hist ? '; residual history (At): ' + r.hist.map(v => fx(v, 3)).join(' → ') : '') + '.');
        lines.push('Solution check (forward chain with the solved fluxes):');
        r.P.forEach(k => dropLines(c, k, r.phi[k]).forEach(s => lines.push(s)));
        lines.push('U per limb: ' + r.P.map(k => SH[k] + ': ' + fx(r.S[k]) + ' − (' + fx(r.drops[k].F) + ') = ' + fx(r.S[k] - r.drops[k].F)).join(';  ') + ' At');
        lines.push('ΣΦ = ' + fx(r.kcl, 3) + ' Wb (flux law satisfied)');
        ui.work.set(lines);
      }
    }
    // core drawing
    function drawCore(cv, c, r, xUnknown) {
      const gg = cv.prep(), ctx = gg.ctx, w = gg.w, h = gg.h, T = theme();
      ctx.clearRect(0, 0, w, h); ctx.fillStyle = T.bg; ctx.fillRect(0, 0, w, h); if (w < 120) return;
      const x0 = 12, x1 = w - 12, y0 = 46, y1 = h - 46, wmax = Math.max.apply(null, c.limbs.map(L => L.w || 0)) || 1;
      const lw = c.limbs.map(L => Math.max(14, Math.min(Math.min(48, w * 0.1), (L.w / wmax) * Math.min(48, w * 0.1))));
      const yk = Math.max(16, Math.min(30, Math.max.apply(null, lw) * 0.7));
      const cxs = [x0 + lw[0] / 2, (x0 + x1) / 2, x1 - lw[2] / 2];
      ctx.fillStyle = withAlpha(T.muted, 0.35); ctx.strokeStyle = T.muted; ctx.lineWidth = 1;
      ctx.fillRect(x0, y0, x1 - x0, yk); ctx.fillRect(x0, y1 - yk, x1 - x0, yk);
      const gapPx = 7, yGap = y0 + yk + (y1 - y0 - 2 * yk) * 0.3;
      c.limbs.forEach((L, k) => {
        const xl = cxs[k] - lw[k] / 2;
        if (!L.on) { ctx.setLineDash([4, 4]); ctx.strokeRect(xl, y0 + yk, lw[k], y1 - y0 - 2 * yk); ctx.setLineDash([]); return; }
        ctx.fillStyle = withAlpha(T.muted, 0.35); ctx.fillRect(xl, y0 + yk, lw[k], y1 - y0 - 2 * yk);
        if (L.lg > 0) { ctx.fillStyle = T.bg; ctx.fillRect(xl - 1, yGap - gapPx / 2, lw[k] + 2, gapPx); ctx.fillStyle = T.a; ctx.font = '10px ' + T.mono; ctx.textAlign = k === 2 ? 'right' : 'left'; ctx.textBaseline = 'middle'; ctx.fillText('g ' + fx(L.lg * 1e3, 3) + 'mm', k === 2 ? xl - 4 : xl + lw[k] + 4, yGap); }
      });
      ctx.strokeStyle = T.muted; ctx.strokeRect(x0, y0, x1 - x0, y1 - y0);
      const yc0 = y0 + yk + (y1 - y0 - 2 * yk) * 0.52, yc1 = y0 + yk + (y1 - y0 - 2 * yk) * 0.86;
      ctx.font = '10px ' + T.mono;
      c.limbs.forEach((L, k) => {
        if (!L.on) return; const cxk = cxs[k];
        // coil
        if (L.N > 0) {
          const cw = lw[k] + 12; ctx.strokeStyle = k === xUnknown ? T.pink : T.a; ctx.lineWidth = 1.6;
          for (let yy = yc0; yy <= yc1; yy += 5) { ctx.beginPath(); ctx.moveTo(cxk - cw / 2, yy); ctx.lineTo(cxk + cw / 2, yy + 2.5); ctx.stroke(); }
          ctx.fillStyle = k === xUnknown ? T.pink : T.a; const ax = k === 2 ? cxk - cw / 2 - 7 : cxk + cw / 2 + 7, up = L.sense > 0;
          arrowHead(ctx, ax, up ? yc1 : yc0, ax, up ? yc0 : yc1, 6);
          ctx.textAlign = 'center'; ctx.textBaseline = 'top';
          const cur = r && k === xUnknown ? r.I : L.I;
          ctx.fillText('N=' + L.N, cxk, y1 + 4); ctx.fillText((k === xUnknown ? 'I=' + fx(cur, 3) : 'I=' + fx(cur, 3)) + 'A', cxk, y1 + 17);
        }
        // flux arrow and label
        if (r) {
          const d = r.drops[k], f = d.phi; if (Math.abs(f) > 0) {
            ctx.strokeStyle = T.c; ctx.fillStyle = T.c; ctx.lineWidth = 2.4; const ya = y0 + yk + 6, yb = yc0 - 8;
            arrowHead(ctx, cxk, f > 0 ? yb : ya, cxk, f > 0 ? ya : yb, 8);
          }
          ctx.fillStyle = T.text; ctx.textAlign = k === 0 ? 'left' : k === 2 ? 'right' : 'center'; ctx.textBaseline = 'bottom';
          const tx = k === 0 ? x0 : k === 2 ? x1 : cxk;
          ctx.fillText('Φ' + SH[k] + '=' + fx(Math.abs(f) * 1e3, 3) + 'mWb' + arrowOf(f), tx, y0 - 15);
          ctx.fillStyle = d.outside ? T.bad : T.muted; ctx.fillText('B=' + fx(Math.abs(d.B), 3) + 'T', tx, y0 - 3);
        }
      });
    }
    const io = stateIO(st, sels, { get: () => ({ tbl: shared.text.split(/\r?\n/).join(';') }), set(o) { if (typeof o.tbl === 'string' && o.tbl.length < 4000) { shared.setText(o.tbl.split(';').join('\n'), false); ui.tbl.value = shared.text; } } });
    return { update, get: io.get, set: io.set };
  }

  /* ---------- 2. inductance, energy, force, losses ---------- */
  function buildIE(root) {
    const st = {}, sels = {}, ui = {}, { ctl, stage } = layout(root);
    const fs0 = FSP.ui.fieldset(ctl, 'View');
    sels.view = sel(fs0, 'Topic', [['core', 'gapped core: L & energy'], ['relay', 'relay / plunger force'], ['rot', 'rotary actuator torque'], ['sine', 'sinusoidal excitation & core loss']], 'core', update);
    sels.mat = sel(fs0, 'Material', [['const', 'constant μr'], ['table', 'B–H table (sub-topic 1)']], 'const', update);
    makeSliders(fs0, [{ k: 'mur', l: 'μr', min: 10, max: 100000, v: 2000, u: '', log: true }], st, update);
    ui.fsCore = FSP.ui.fieldset(ctl, 'Gapped core');
    makeSliders(ui.fsCore, [{ k: 'N', l: 'turns N', min: 1, max: 5000, v: 500, step: 1, u: '' }, { k: 'lc', l: 'l iron', min: 0.01, max: 2, v: 0.3, u: 'm', log: true },
      { k: 'A', l: 'area A', min: 0.1, max: 500, v: 10, u: 'cm²', log: true }, { k: 'lg', l: 'air gap', min: 0, max: 10, v: 1, step: 0.01, u: 'mm' },
      { k: 'i', l: 'current i', min: 0.001, max: 50, v: 1, u: 'A', log: true }], st, update);
    ui.fsRel = FSP.ui.fieldset(ctl, 'Relay (U-core + armature)');
    makeSliders(ui.fsRel, [{ k: 'rN', l: 'turns N', min: 1, max: 10000, v: 1000, step: 1, u: '' }, { k: 'rlc', l: 'l iron', min: 0.01, max: 2, v: 0.2, u: 'm', log: true },
      { k: 'rA', l: 'pole A', min: 0.1, max: 100, v: 4, u: 'cm²', log: true }, { k: 'rx', l: 'gap x', min: 0.01, max: 10, v: 1, u: 'mm', log: true },
      { k: 'ri', l: 'current i', min: 0.001, max: 20, v: 0.5, u: 'A', log: true }, { k: 'rFs', l: 'spring F', min: 0.01, max: 1000, v: 20, u: 'N', log: true },
      { k: 'rx0', l: 'open gap x0', min: 0.05, max: 10, v: 2, u: 'mm', log: true }], st, update);
    sels.ng = sel(ui.fsRel, 'Gaps', [['2', '2 (armature across both poles)'], ['1', '1 (single gap)']], '2', update);
    ui.fsRot = FSP.ui.fieldset(ctl, 'Rotary reluctance actuator');
    makeSliders(ui.fsRot, [{ k: 'Ld', l: 'L aligned', min: 0.001, max: 5, v: 0.2, u: 'H', log: true }, { k: 'Lq', l: 'L unalign.', min: 0.0005, max: 5, v: 0.05, u: 'H', log: true },
      { k: 'ti', l: 'current i', min: 0, max: 50, v: 2, step: 0.01, u: 'A' }, { k: 'th', l: 'angle θ', min: 0, max: 180, v: 30, step: 0.5, u: '°' }], st, update);
    ui.fsSine = FSP.ui.fieldset(ctl, 'Sinusoidal excitation & losses');
    makeSliders(ui.fsSine, [{ k: 'V', l: 'V rms', min: 1, max: 50000, v: 230, u: 'V', log: true }, { k: 'f', l: 'f', min: 1, max: 2000, v: 50, u: 'Hz', log: true },
      { k: 'sN', l: 'turns N', min: 1, max: 10000, v: 200, step: 1, u: '' }, { k: 'sA', l: 'area A', min: 0.1, max: 1000, v: 25, u: 'cm²', log: true },
      { k: 'vol', l: 'core vol.', min: 1, max: 1e6, v: 3000, u: 'cm³', log: true }, { k: 'n', l: 'Steinmetz n', min: 1, max: 3, v: 1.6, step: 0.01, u: '' },
      { k: 'f1', l: 'test f1', min: 1, max: 2000, v: 50, u: 'Hz', log: true }, { k: 'P1', l: 'test P1', min: 0.001, max: 1e6, v: 100, u: 'W', log: true },
      { k: 'f2', l: 'test f2', min: 1, max: 2000, v: 60, u: 'Hz', log: true }, { k: 'P2', l: 'test P2', min: 0.001, max: 1e6, v: 130, u: 'W', log: true }], st, update);
    ui.fsSine.appendChild(FSP.ui.el('div', { class: 'note', text: 'Both test readings must be at the same Bmax (V/f kept constant). Bmax used = value from V, f, N, A above.' }));
    ui.msg = msgBox(ctl);
    ui.cv = wrapCanvas(stage, 'Inductance, force or loss plot', 300); ui.cv.onResize(update);
    ui.hud = FSP.ui.el('div', { class: 'hud' }); stage.appendChild(ui.hud);
    ui.warn = FSP.ui.el('div', { class: 'msg warn', hidden: '' }); stage.appendChild(ui.warn);
    ui.work = FSP.ui.working(stage);
    shared.listeners.push(() => update());

    function material() { return sels.mat.value === 'table' ? shared.tableMat() : E.makeMaterial({ mode: 'const', mur: st.mur.get() }); }
    function update() {
      if (!ui.cv) return;
      const v = sels.view.value, T = theme();
      ui.fsCore.hidden = v !== 'core'; ui.fsRel.hidden = v !== 'relay'; ui.fsRot.hidden = v !== 'rot'; ui.fsSine.hidden = v !== 'sine';
      const needMat = v === 'core' || v === 'relay'; rowOf(sels.mat).hidden = !needMat; st.mur.el.hidden = !needMat || sels.mat.value !== 'const';
      ui.warn.hidden = true; ui.msg.show('');
      const mat = material();
      if (needMat && !mat.ok) { ui.msg.show('B–H table invalid: ' + shared.parsed.msg + ' (edit it in sub-topic 1).'); drawPlot(ui.cv, {}); hudSet(ui.hud, []); ui.work.set('—'); return; }
      if (v === 'core') return core(mat, T);
      if (v === 'relay') return relay(mat, T);
      if (v === 'rot') return rot(T);
      return sine(T);
    }
    function core(mat, T) {
      const k = { N: Math.round(st.N.get()), lc: st.lc.get(), A: st.A.get() * 1e-4, lg: st.lg.get() * 1e-3, mat }; k.Ag = k.A;
      const i = st.i.get(), a = E.coreAnalyze(k, i);
      const n = 120, is = [], ls = []; for (let j = 0; j <= n; j++) { const ii = 1.5 * i * j / n; is.push(ii); ls.push(k.N * E.coreFlux(k, ii)); }
      const ii0 = is.filter(z => z <= i + 1e-15), ll0 = ls.slice(0, ii0.length);
      const lmax = Math.max.apply(null, ls) * 1.1 || 1;
      drawPlot(ui.cv, {
        xmin: 0, xmax: 1.5 * i, ymin: 0, ymax: lmax, xlabel: 'current i (A)', ylabel: 'flux linkage λ (Wb-t)',
        fills: [{ x: ii0.concat([i]), y: ll0.concat([0]), color: withAlpha(T.b, 0.25) }, { x: [0].concat(ii0), y: [a.lam].concat(ll0), color: withAlpha(T.a, 0.25) }],
        series: [{ x: is, y: ls, color: T.c, width: 2.2 }], points: [{ x: i, y: a.lam, label: 'operating point', color: T.text }],
        legend: [{ text: 'energy W (above curve)', color: T.a }, { text: "co-energy W' (below)", color: T.b }],
      });
      const lin = mat.mode === 'const';
      hudSet(ui.hud, [['flux Φ', fx(a.phi, 4) + ' Wb'], ['B core / B gap', fx(a.B, 4) + ' / ' + fx(a.Bg, 4) + ' T'], ['λ = NΦ', fx(a.lam, 4) + ' Wb-t'],
        ['L = λ/i (secant)', fE(a.Lsec, 'H')], ['L incremental dλ/di', fE(a.Linc, 'H')], ['R core / R gap', fx(a.Rc, 4) + ' / ' + fx(a.Rg, 4) + ' A/Wb'],
        ['energy W', fE(a.W, 'J')], ['  in core / in gap', fE(a.Wcore, 'J') + ' / ' + fE(a.Wgap, 'J')], ["co-energy W'", fE(a.Wco, 'J')], ['gap share of W', f2(100 * a.Wgap / a.W, 1) + ' %']]);
      if (a.outside) { ui.warn.hidden = false; ui.warn.textContent = 'Core B is beyond the B–H table: extrapolated.'; }
      ui.work.set([
        'R_core = ' + (lin ? 'l_c/(μ0 μr A) = ' + fx(k.lc) + '/(4π×10⁻⁷ × ' + fx(mat.mur) + ' × ' + fx(k.A) + ')' : 'H_c·l_c/Φ (secant, nonlinear)') + ' = ' + fx(a.Rc) + ' A/Wb',
        'R_gap = l_g/(μ0 A) = ' + fx(k.lg) + '/(4π×10⁻⁷ × ' + fx(k.A) + ') = ' + fx(a.Rg) + ' A/Wb',
        lin ? 'L = N²/(R_core + R_gap) = ' + k.N + '²/' + fx(a.Rc + a.Rg) + ' = ' + fE(a.Lsec, 'H') : 'Nonlinear: Φ solved from N·i = H(Φ/A)·l_c + Φ·R_gap (bisection);  L = λ/i = ' + fE(a.Lsec, 'H'),
        'Φ = ' + fx(a.phi) + ' Wb,  B = Φ/A = ' + fx(a.B) + ' T,  λ = NΦ = ' + fx(a.lam) + ' Wb-t',
        'W_gap = (B²/2μ0)·A·l_g = ' + fx(a.Bg) + '²/(2 × 4π×10⁻⁷) × ' + fx(k.A) + ' × ' + fx(k.lg) + ' = ' + fE(a.Wgap, 'J'),
        'W_core = ' + (lin ? '(B²/2μ0μr)·A·l_c' : '(∫₀^B H dB)·A·l_c (exact trapezoids on the table)') + ' = ' + fE(a.Wcore, 'J'),
        'W = ∫ i dλ = ' + fE(a.W, 'J') + (lin ? ' = ½ L i² = ½ × ' + fx(a.Lsec) + ' × ' + fx(i) + '² ✓' : ''),
        "Co-energy W' = λ·i − W = " + fE(a.Wco, 'J') + (lin ? ' (= W for a linear core)' : ' (≠ W: saturating core)'),
        'Ratio W_gap/W_core = ' + fx(a.Wgap / a.Wcore) + (lin ? ' = R_gap/R_core: most energy is stored in the gap.' : ''),
      ]);
    }
    function relay(mat, T) {
      const r = { N: Math.round(st.rN.get()), lc: st.rlc.get(), A: st.rA.get() * 1e-4, ng: +sels.ng.value, mat };
      const i = st.ri.get(), x = st.rx.get() * 1e-3, Fs = st.rFs.get(), x0 = st.rx0.get() * 1e-3, lin = mat.mode === 'const';
      const F = E.relayForce(r, i, x), ip = E.relayPullIn(r, x0, Fs);
      const xm = Math.max(x, x0) * 1.5, n = 100, xs = [], fs = [], fp = [];
      for (let j = 1; j <= n; j++) { const xx = xm * j / n; xs.push(xx * 1e3); fs.push(Math.abs(E.relayForce(r, i, xx))); if (fin(ip)) fp.push(Math.abs(E.relayForce(r, ip, xx))); }
      const ymax = Math.min(Math.max(Fs * 3, Math.abs(F) * 1.3), Math.max.apply(null, fs.concat(fp).concat([Fs])) * 1.05) || 1;
      drawPlot(ui.cv, {
        xmin: 0, xmax: xm * 1e3, ymin: 0, ymax, xlabel: 'gap x (mm)', ylabel: 'attractive force |F| (N)',
        series: [{ x: xs, y: fs, color: T.c, width: 2.2 }].concat(fp.length ? [{ x: xs, y: fp, color: T.pink, width: 1.5, dash: [5, 3] }] : []),
        hlines: [{ y: Fs, label: 'spring/load F_s', color: T.warn }], vlines: [{ x: x0 * 1e3, label: 'x0 (open)', color: T.muted }],
        points: [{ x: x * 1e3, y: Math.abs(F), label: f2(Math.abs(F), 2) + ' N', color: T.text }],
        legend: [{ text: 'F at i = ' + fx(i, 3) + ' A', color: T.c }].concat(fp.length ? [{ text: 'F at pull-in ' + fx(ip, 3) + ' A', color: T.pink }] : []),
      });
      const Lx = lin ? E.relayL(r, x) : NaN, dL = lin ? E.relaydLdx(r, x) : NaN, B = E.coreFlux(E.relayCore(r, x), i) / r.A;
      hudSet(ui.hud, [['force F (+ = opening)', fx(F, 4) + ' N'], ['attractive |F|', fx(Math.abs(F), 4) + ' N'], ['gap B', fx(B, 4) + ' T'],
        ['L(x)', lin ? fE(Lx, 'H') : fE(r.N * E.coreFlux(E.relayCore(r, x), i) / i, 'H') + ' (secant)'], ['dL/dx', lin ? fx(dL, 4) + ' H/m' : '— (nonlinear)'],
        ['pull-in current at x0', fin(ip) ? fx(ip, 4) + ' A' : '—'], ['closes at present i?', fin(ip) ? (i >= ip ? 'yes (i ≥ i_pull)' : 'no (i < i_pull)') : '—']]);
      if (mat.outside(B)) { ui.warn.hidden = false; ui.warn.textContent = 'Iron B exceeds the B–H table: extrapolated.'; }
      const Ri = lin ? r.lc / (MU0 * mat.mur * r.A) : NaN, Rg = r.ng * x / (MU0 * r.A);
      ui.work.set(lin ? [
        'Reluctance: R(x) = l_c/(μ0 μr A) + n·x/(μ0 A) = ' + fx(Ri) + ' + ' + r.ng + '×' + fx(x) + '/(4π×10⁻⁷ × ' + fx(r.A) + ') = ' + fx(Ri + Rg) + ' A/Wb',
        'L(x) = N²/R(x) = ' + r.N + '²/' + fx(Ri + Rg) + ' = ' + fE(Lx, 'H'),
        'dL/dx = −N²·(n/(μ0 A))/R(x)² = ' + fx(dL) + ' H/m',
        'F = ½ i² dL/dx = ½ × ' + fx(i) + '² × (' + fx(dL) + ') = ' + fx(F) + ' N  (negative → pulls the armature in, closing the gap)',
        'Check (iron ideal): F ≈ n·B²A/(2μ0) with B = μ0 N i/(n x) — equal when R_iron ≪ R_gap.',
        'Pull-in at x0 = ' + fx(x0 * 1e3) + ' mm: ½ i² |dL/dx(x0)| = F_s → i = √(2F_s/|dL/dx|) = √(2 × ' + fx(Fs) + '/' + fx(Math.abs(E.relaydLdx(r, x0))) + ') = ' + fx(ip) + ' A',
        'Force ∝ i² and ∝ 1/R(x)²: it rises steeply as the gap closes, so once pulled in the relay stays closed (hold current < pull-in current).',
      ] : [
        'Nonlinear core: λ(i, x) from N·i = H(λ/(NA))·l_c + (λ/N)·n x/(μ0 A) (bisection).',
        "Co-energy W'(i, x) = ∫₀^i λ(i', x) di' (Simpson, 120 steps);  F = ∂W'/∂x at constant i (central difference).",
        "F = " + fx(F) + ' N.  Equivalent route: F = −∂W/∂x at constant λ (same value; see tests).',
        'Pull-in current found by bisection on |F(i, x0)| = F_s: i = ' + fx(ip) + ' A',
      ]);
    }
    function rot(T) {
      const Ld = st.Ld.get(), Lq = st.Lq.get(), i = st.ti.get(), th = st.th.get() * PI / 180;
      if (!(Ld > Lq)) { ui.msg.show('Need L_aligned > L_unaligned for a reluctance torque.'); }
      const ths = [], Ts = [], Lsr = []; for (let j = 0; j <= 180; j++) { ths.push(j); Ts.push(E.rotTorque(Ld, Lq, i, j * PI / 180)); Lsr.push(E.rotL(Ld, Lq, j * PI / 180)); }
      const Tm = 0.5 * i * i * Math.abs(Ld - Lq) || 1, Tn = E.rotTorque(Ld, Lq, i, th);
      drawPlot(ui.cv, { xmin: 0, xmax: 180, ymin: -1.2 * Tm, ymax: 1.2 * Tm, y2min: 0, y2max: Math.max(Ld, Lq) * 1.1, xlabel: 'rotor angle θ from aligned axis (deg)', ylabel: 'torque T (N·m)', y2label: 'L(θ) (H)',
        series: [{ x: ths, y: Ts, color: T.c, width: 2.2 }, { x: ths, y: Lsr, color: T.a, width: 1.5, dash: [5, 3], axis: 'y2' }],
        points: [{ x: st.th.get(), y: Tn, label: fx(Tn, 3) + ' N·m', color: T.text }], legend: [{ text: 'T(θ)', color: T.c }, { text: 'L(θ)', color: T.a }] });
      hudSet(ui.hud, [['L(θ)', fE(E.rotL(Ld, Lq, th), 'H')], ['dL/dθ', fx(-(Ld - Lq) * Math.sin(2 * th), 4) + ' H/rad'], ['torque T', fx(Tn, 4) + ' N·m'], ['max |T| (θ = 45°)', fx(0.5 * i * i * (Ld - Lq), 4) + ' N·m']]);
      ui.work.set([
        'L(θ) = (L_d + L_q)/2 + (L_d − L_q)/2 · cos 2θ  (θ measured from the aligned position)',
        'T = ½ i² dL/dθ = −½ i² (L_d − L_q) sin 2θ = −½ × ' + fx(i) + '² × ' + fx(Ld - Lq) + ' × sin(' + fx(2 * st.th.get()) + '°) = ' + fx(Tn) + ' N·m',
        'Negative T for 0 < θ < 90°: the torque pulls the rotor back toward alignment (minimum reluctance). Peak at θ = 45°.',
        'Linear magnetics assumed (L independent of i), so co-energy W\' = ½ L(θ) i² and T = ∂W\'/∂θ at constant i.',
      ]);
    }
    function sine(T) {
      const V = st.V.get(), f = st.f.get(), N = Math.round(st.sN.get()), A = st.sA.get() * 1e-4, vol = st.vol.get() * 1e-6, n = st.n.get();
      const sf = E.sineFlux(V, f, N, A), sp = E.separateLosses({ f1: st.f1.get(), P1: st.P1.get(), f2: st.f2.get(), P2: st.P2.get(), Bm: sf.Bm, n, vol });
      if (!sp.ok) ui.msg.show(sp.msg);
      const fmax = Math.max(st.f1.get(), st.f2.get(), f) * 1.4, fsr = [], yy = [], ph = [], pe = [];
      for (let j = 0; j <= 100; j++) { const ff = fmax * j / 100; fsr.push(ff); yy.push(sp.a + sp.b * ff); }
      const ymax = Math.max(sp.y1, sp.y2, sp.a + sp.b * fmax) * 1.15 || 1;
      drawPlot(ui.cv, { xmin: 0, xmax: fmax, ymin: 0, ymax: fin(ymax) && ymax > 0 ? ymax : 1, xlabel: 'frequency f (Hz)', ylabel: 'P/f (J per cycle)',
        series: sp.ok ? [{ x: fsr, y: yy, color: T.c, width: 2.2 }] : [],
        points: [{ x: st.f1.get(), y: sp.y1, label: 'test 1', color: T.a }, { x: st.f2.get(), y: sp.y2, label: 'test 2', color: T.a, below: true }].concat(sp.ok ? [{ x: 0, y: sp.a, label: 'a = hysteresis/cycle', color: T.b }] : []),
        legend: [{ text: 'P/f = a + b·f', color: T.c }] });
      const atF = sp.ok ? { Ph: sp.a * f, Pe: sp.b * f * f } : { Ph: NaN, Pe: NaN };
      hudSet(ui.hud, [['Φmax', fx(sf.phim * 1e3, 4) + ' mWb'], ['Bmax', fx(sf.Bm, 4) + ' T'], ['a (intercept)', fx(sp.a, 4) + ' J'], ['b (slope)', fx(sp.b, 4) + ' J·s'],
        ['k_h', fx(sp.kh, 4) + ' J/(m³·T^n)'], ['k_e', fx(sp.ke, 4) + ' J·s/(m³·T²)'], ['P_h / P_e at f1', fx(sp.Ph1, 4) + ' / ' + fx(sp.Pe1, 4) + ' W'], ['P_h / P_e at f2', fx(sp.Ph2, 4) + ' / ' + fx(sp.Pe2, 4) + ' W'],
        ['at f = ' + fx(f) + ' Hz (same Bmax)', fx(atF.Ph, 4) + ' + ' + fx(atF.Pe, 4) + ' W']]);
      if (sf.Bm > 2) { ui.warn.hidden = false; ui.warn.textContent = 'Bmax > 2 T: a real silicon-steel core would be deeply saturated — increase N or A.'; }
      ui.work.set([
        'e = N dΦ/dt with Φ = Φmax sin ωt → E_rms = ωNΦmax/√2 = (2π/√2) f N Φmax = 4.443 f N Φmax (≈ 4.44)',
        'Φmax = V/(4.443 f N) = ' + fx(V) + '/(4.443 × ' + fx(f) + ' × ' + N + ') = ' + fx(sf.phim) + ' Wb;  Bmax = Φmax/A = ' + fx(sf.Bm) + ' T',
        'Steinmetz: P_h = k_h f Bmax^n V,   eddy: P_e = k_e f² Bmax² V   (n = ' + fx(n) + ')',
        'At constant Bmax:  P/f = a + b f  with a = k_h Bmax^n V and b = k_e Bmax² V',
        'Test 1: P1/f1 = ' + fx(st.P1.get()) + '/' + fx(st.f1.get()) + ' = ' + fx(sp.y1) + ';   Test 2: P2/f2 = ' + fx(st.P2.get()) + '/' + fx(st.f2.get()) + ' = ' + fx(sp.y2),
        'b = (' + fx(sp.y2) + ' − ' + fx(sp.y1) + ')/(' + fx(st.f2.get()) + ' − ' + fx(st.f1.get()) + ') = ' + fx(sp.b) + ';   a = ' + fx(sp.y1) + ' − b × ' + fx(st.f1.get()) + ' = ' + fx(sp.a),
        'P_h(f1) = a f1 = ' + fx(sp.Ph1) + ' W,  P_e(f1) = b f1² = ' + fx(sp.Pe1) + ' W  (sum = ' + fx(sp.Ph1 + sp.Pe1) + ' W ✓)',
        'k_h = a/(Bmax^n V) = ' + fx(sp.a) + '/(' + fx(sf.Bm) + '^' + fx(n) + ' × ' + fx(vol) + ') = ' + fx(sp.kh) + ';   k_e = b/(Bmax² V) = ' + fx(sp.ke),
      ]);
    }
    const io = stateIO(st, sels);
    return { update, get: io.get, set: io.set };
  }

  /* ---------- 3. autotransformer ---------- */
  function buildAuto(root) {
    const st = {}, sels = {}, ui = {}, { ctl, stage } = layout(root);
    const fs1 = FSP.ui.fieldset(ctl, 'Two-winding transformer');
    makeSliders(fs1, [{ k: 'V1', l: 'V1 winding', min: 1, max: 500000, v: 220, u: 'V', log: true }, { k: 'V2', l: 'V2 winding', min: 1, max: 500000, v: 110, u: 'V', log: true },
      { k: 'S', l: 'S rating', min: 0.01, max: 100000, v: 5, u: 'kVA', log: true }, { k: 'Pc', l: 'P core', min: 0, max: 100000, v: 50, step: 0.1, u: 'W' }, { k: 'Pcu', l: 'P cu (FL)', min: 0, max: 100000, v: 100, step: 0.1, u: 'W' }], st, update);
    const fs2 = FSP.ui.fieldset(ctl, 'Auto connection & load');
    sels.common = sel(fs2, 'Common wdg', [['V1', 'V1 winding'], ['V2', 'V2 winding']], 'V1', update);
    sels.dir = sel(fs2, 'Use as', [['down', 'step-down (H → L)'], ['up', 'step-up (L → H)']], 'down', update);
    makeSliders(fs2, [{ k: 'x', l: 'load', min: 0, max: 1.5, v: 1, step: 0.01, u: '× rated' }, { k: 'pf', l: 'pf', min: 0.05, max: 1, v: 0.8, step: 0.01, u: '' }], st, update);
    fs2.appendChild(FSP.ui.el('div', { class: 'note', text: 'Windings connected series-aiding: V_H = V_common + V_series, V_L = V_common.' }));
    ui.msg = msgBox(ctl);
    ui.cv = wrapCanvas(stage, 'Autotransformer connection with currents', 300); ui.cv.onResize(update);
    ui.hud = FSP.ui.el('div', { class: 'hud' }); stage.appendChild(ui.hud);
    ui.work = FSP.ui.working(stage);
    function update() {
      if (!ui.cv) return;
      const r = E.auto({ V1: st.V1.get(), V2: st.V2.get(), S: st.S.get() * 1e3, common: sels.common.value, dir: sels.dir.value, x: st.x.get(), pf: st.pf.get(), Pcore: st.Pc.get(), Pcu: st.Pcu.get() });
      ui.msg.show(r.ok ? '' : r.msg); if (!r.ok) { hudSet(ui.hud, []); ui.work.set(r.msg); return; }
      const x = st.x.get(), down = sels.dir.value === 'down';
      draw(r, x, down);
      hudSet(ui.hud, [['voltage ratio', fx(r.VH, 5) + ' / ' + fx(r.VL, 5) + ' V  (' + fx(r.ratio, 4) + ':1)'], ['S auto rating', fx(r.Sa / 1e3, 5) + ' kVA (' + fx(r.gain, 4) + '× S)'],
        ['I_H (= series wdg)', fx(r.IH, 5) + ' A'], ['I_L', fx(r.IL, 5) + ' A'], ['common wdg I_L − I_H', fx(r.Icom, 5) + ' A'],
        ['transformed S', fx(r.Strans / 1e3, 5) + ' kVA'], ['conducted S', fx(r.Scond / 1e3, 5) + ' kVA'],
        ['η two-winding @ x', fx(100 * r.eta2w, 5) + ' %'], ['η auto @ x', fx(100 * r.etaAuto, 5) + ' %']]);
      ui.work.set([
        'Common winding: ' + fx(r.Vc) + ' V (rated ' + fx(r.IcomRated) + ' A); series winding: ' + fx(r.Vs) + ' V (rated ' + fx(r.IserRated) + ' A)',
        'V_H = V_common + V_series = ' + fx(r.Vc) + ' + ' + fx(r.Vs) + ' = ' + fx(r.VH) + ' V;   V_L = ' + fx(r.VL) + ' V;   ratio a = ' + fx(r.ratio),
        'The series winding carries the HV-side current: I_H = I_series,rated = S/V_series = ' + fx(r.S || r.Sw) + '/' + fx(r.Vs) + ' = ' + fx(r.IH) + ' A',
        'S_auto = V_H·I_H = V_H·S/(V_H − V_L) = ' + fx(r.VH) + ' × ' + fx(r.Sw) + '/' + fx(r.VH - r.VL) + ' = ' + fx(r.Sa) + ' VA',
        'I_L = S_auto/V_L = ' + fx(r.IL) + ' A;  common winding I_L − I_H = ' + fx(r.Icom) + ' A = S/V_common ✓ (both windings exactly at rating)',
        'Transformed (magnetically) = (V_H − V_L)·I_H = ' + fx(r.Strans) + ' VA = S (two-winding rating);  conducted = V_L·I_H = ' + fx(r.Scond) + ' VA',
        'Losses unchanged (same winding currents and flux): P = P_core + x²P_cu = ' + fx(st.Pc.get()) + ' + ' + fx(x) + '² × ' + fx(st.Pcu.get()) + ' = ' + fx(r.loss) + ' W',
        'η_2w = xS·pf/(xS·pf + P) = ' + fx(r.P2w) + '/' + fx(r.P2w + r.loss) + ' = ' + fx(100 * r.eta2w) + ' %;   η_auto = ' + fx(r.Pauto) + '/' + fx(r.Pauto + r.loss) + ' = ' + fx(100 * r.etaAuto) + ' %',
        'Trade-off: no electrical isolation between HV and LV, and a higher fault level (lower per-unit impedance on the auto base).',
      ]);
    }
    function draw(r, x, down) {
      const g = ui.cv.prep(), ctx = g.ctx, w = g.w, h = g.h, T = theme(); ctx.clearRect(0, 0, w, h); ctx.fillStyle = T.bg; ctx.fillRect(0, 0, w, h); if (w < 160) return;
      const cx = w / 2, yH = 40, yT = h * 0.47, yN = h - 40, xl = Math.max(30, cx - Math.min(170, w * 0.4)), xr = Math.min(w - 30, cx + Math.min(170, w * 0.4));
      ctx.lineWidth = 2; ctx.strokeStyle = T.muted;
      const coil = (ya, yb, col) => { const n = 6, step = (yb - ya) / n; ctx.strokeStyle = col; ctx.lineWidth = 2; ctx.beginPath(); for (let k = 0; k < n; k++) ctx.arc(cx, ya + step * (k + 0.5), step / 2, -PI / 2, PI / 2, false); ctx.stroke(); };
      ctx.strokeStyle = T.muted; ctx.beginPath(); ctx.moveTo(xl, yH); ctx.lineTo(cx, yH); ctx.lineTo(cx, yH + 12); ctx.moveTo(cx, yT - 12); ctx.lineTo(cx, yT + 12); ctx.moveTo(cx, yN - 12); ctx.lineTo(cx, yN); ctx.lineTo(xl, yN); ctx.moveTo(cx, yN); ctx.lineTo(xr, yN); ctx.moveTo(cx, yT); ctx.lineTo(xr, yT); ctx.stroke();
      coil(yH + 12, yT - 12, T.a); coil(yT + 12, yN - 12, T.c);
      // core bar
      ctx.strokeStyle = T.muted; ctx.lineWidth = 3; ctx.beginPath(); ctx.moveTo(cx + 26, yH + 10); ctx.lineTo(cx + 26, yN - 10); ctx.stroke(); ctx.lineWidth = 1.5;
      const T0 = (t, xx, yy, al, col) => { ctx.fillStyle = col || T.text; ctx.font = '11px ' + T.mono; ctx.textAlign = al || 'left'; ctx.textBaseline = 'middle'; ctx.fillText(t, xx, yy); };
      [[xl, yH], [xl, yN], [xr, yT], [xr, yN]].forEach(p => { ctx.fillStyle = T.text; ctx.beginPath(); ctx.arc(p[0], p[1], 3.5, 0, 2 * PI); ctx.fill(); });
      T0('H', xl - 4, yH - 12, 'left'); T0('L', xr + 4, yT - 12, 'right');
      T0('V_H = ' + fx(r.VH, 4) + ' V', xl, (yH + yN) / 2, 'left', T.muted); T0('V_L = ' + fx(r.VL, 4) + ' V', xr, (yT + yN) / 2, 'right', T.muted);
      T0('series ' + fx(r.Vs, 4) + ' V', cx - 14, (yH + yT) / 2, 'right', T.a); T0('common ' + fx(r.Vc, 4) + ' V', cx - 14, (yT + yN) / 2, 'right', T.c);
      T0(down ? 'source' : 'load', xl, yH + 16, 'left', T.muted); T0(down ? 'load' : 'source', xr, yT + 16, 'right', T.muted);
      // currents
      ctx.strokeStyle = T.pink; ctx.fillStyle = T.pink; ctx.lineWidth = 2;
      const mxH = (xl + cx) / 2, mxL = (cx + xr) / 2, IH = x * r.IH, IL = x * r.IL, Ic = x * r.Icom;
      if (down) arrowHead(ctx, mxH - 18, yH, mxH + 18, yH, 7); else arrowHead(ctx, mxH + 18, yH, mxH - 18, yH, 7);
      if (down) arrowHead(ctx, mxL - 18, yT, mxL + 18, yT, 7); else arrowHead(ctx, mxL + 18, yT, mxL - 18, yT, 7);
      const ya = (yT + yN) / 2; if (down) arrowHead(ctx, cx + 12, ya + 16, cx + 12, ya - 16, 7); else arrowHead(ctx, cx + 12, ya - 16, cx + 12, ya + 16, 7);
      const yb = (yH + yT) / 2; if (down) arrowHead(ctx, cx + 12, yb - 16, cx + 12, yb + 16, 7); else arrowHead(ctx, cx + 12, yb + 16, cx + 12, yb - 16, 7);
      ctx.font = '11px ' + T.mono; ctx.textBaseline = 'bottom'; ctx.textAlign = 'center';
      ctx.fillText('I_H=' + fx(IH, 4) + 'A', mxH, yH - 6); ctx.fillText('I_L=' + fx(IL, 4) + 'A', mxL, yT - 6);
      ctx.textAlign = 'left'; ctx.textBaseline = 'middle'; ctx.fillText(fx(Ic, 4) + 'A', cx + 32, ya); ctx.fillText(fx(IH, 4) + 'A', cx + 32, yb);
      T0('S_auto = ' + fx(r.Sa / 1e3, 4) + ' kVA  (load ' + fx(x, 3) + '×)', 8, 12, 'left', T.text);
    }
    const io = stateIO(st, sels);
    return { update, get: io.get, set: io.set };
  }

  /* ---------- 4. three-phase transformers ---------- */
  function buildTP(root) {
    const st = {}, sels = {}, ui = {}, { ctl, stage } = layout(root);
    const fs1 = FSP.ui.fieldset(ctl, 'Bank of three single-phase units');
    sels.conn = sel(fs1, 'HV-LV', [['YY', 'Y – Y'], ['YD', 'Y – Δ'], ['DY', 'Δ – Y'], ['DD', 'Δ – Δ'], ['VV', 'open-Δ (V – V)']], 'YD', update);
    makeSliders(fs1, [{ k: 'S', l: 'S bank (3 units)', min: 0.1, max: 1e6, v: 150, u: 'kVA', log: true }, { k: 'VH', l: 'V line HV', min: 1, max: 1e6, v: 11000, u: 'V', log: true },
      { k: 'VL', l: 'V line LV', min: 1, max: 1e6, v: 415, u: 'V', log: true }], st, update);
    sels.conv = sel(fs1, 'Phase shift', [['ANSI', 'ANSI: LV lags HV by 30° (Yd1/Dy1)'], ['IEC11', 'IEC clock 11: LV leads 30° (Dyn11/YNd11)']], 'ANSI', update);
    const fs2 = FSP.ui.fieldset(ctl, 'Each unit: series impedance');
    makeSliders(fs2, [{ k: 'R', l: 'R_eq (HV wdg)', min: 0, max: 10000, v: 12, step: 0.01, u: 'Ω' }, { k: 'X', l: 'X_eq (HV wdg)', min: 0, max: 10000, v: 30, step: 0.01, u: 'Ω' }], st, update);
    fs2.appendChild(FSP.ui.el('div', { class: 'note', text: 'Impedance of one single-phase unit referred to its own HV winding.' }));
    ui.msg = msgBox(ctl);
    ui.cv = wrapCanvas(stage, 'HV and LV line-voltage phasors', 300); ui.cv.onResize(update);
    ui.hud = FSP.ui.el('div', { class: 'hud' }); stage.appendChild(ui.hud);
    ui.work = FSP.ui.working(stage);
    const nm = { Y: 'Y', D: 'Δ' };
    function update() {
      if (!ui.cv) return;
      const conn = sels.conn.value, r = E.threePhase({ conn, S: st.S.get() * 1e3, VHL: st.VH.get(), VLL: st.VL.get(), conv: sels.conv.value, Req: st.R.get(), Xeq: st.X.get() });
      ui.msg.show(r.ok ? '' : r.msg); if (!r.ok) { hudSet(ui.hud, []); ui.work.set(r.msg); return; }
      rowOf(sels.conv).hidden = !(conn === 'YD' || conn === 'DY');
      draw(r);
      const hn = conn === 'VV' ? 'Δ (open)' : nm[conn[0]], ln = conn === 'VV' ? 'Δ (open)' : nm[conn[1]];
      hudSet(ui.hud, [['HV: line / winding V', fx(st.VH.get(), 5) + ' / ' + fx(r.VHw, 5) + ' V'], ['LV: line / winding V', fx(st.VL.get(), 5) + ' / ' + fx(r.VLw, 5) + ' V'],
        ['winding ratio a', fx(r.a, 5)], ['line-voltage ratio', fx(r.lineRatio, 5)], ['unit rating', fx(r.Sunit / 1e3, 5) + ' kVA'], ['bank capacity', fx(r.Sbank / 1e3, 5) + ' kVA'],
        ['HV line / winding I', fx(r.IHL, 5) + ' / ' + fx(r.IHw, 5) + ' A'], ['LV line / winding I', fx(r.ILL, 5) + ' / ' + fx(r.ILw, 5) + ' A'],
        ['LV vs HV line voltage', (r.shift > 0 ? '+' : '') + r.shift + '°'], ['Z per phase (Y-eq) HV', fx(r.ZYH.re, 4) + ' + j' + fx(r.ZYH.im, 4) + ' Ω'],
        ['Z per phase (Y-eq) LV', fx(r.ZYL.re, 4) + ' + j' + fx(r.ZYL.im, 4) + ' Ω'], ['Z per unit', fx(r.Zpu.re, 4) + ' + j' + fx(r.Zpu.im, 4) + ' pu']]);
      const L = [
        'Connection ' + hn + ' (HV) – ' + ln + ' (LV). Y: V_winding = V_line/√3, I_winding = I_line.  Δ: V_winding = V_line, I_winding = I_line/√3.',
        'HV winding V = ' + fx(r.VHw) + ' V;  LV winding V = ' + fx(r.VLw) + ' V;  turns ratio a = ' + fx(r.a),
        'Line-voltage ratio V_HL/V_LL = ' + fx(r.lineRatio) + (conn === 'YD' ? ' = √3·a' : conn === 'DY' ? ' = a/√3' : ' = a'),
      ];
      if (conn === 'VV') L.push('Open-Δ: two units supply the full 3-phase voltage set; each carries the line current. Capacity = √3·V_L·I_rated = √3·S_unit = ' + fx(r.Sbank / 1e3) + ' kVA = ' + fx(100 * r.vvFraction) + ' % of the Δ-Δ bank (3 S_unit); unit utilisation √3/2 = ' + fx(100 * r.vvUtil) + ' %.');
      else L.push('Unit rating = S/3 = ' + fx(r.Sunit / 1e3) + ' kVA.');
      L.push('I_line = S/(√3 V_line):  HV ' + fx(r.Sbank) + '/(√3 × ' + fx(st.VH.get()) + ') = ' + fx(r.IHL) + ' A;  LV = ' + fx(r.ILL) + ' A;  winding currents ' + fx(r.IHw) + ' A / ' + fx(r.ILw) + ' A');
      if (conn === 'YD' || conn === 'DY') L.push('Phase shift: ' + (sels.conv.value === 'ANSI' ? 'ANSI C57.12 standard — the HV line voltage leads the LV line voltage by 30° (V_ab lags V_AB by 30°); IEC equivalent Yd1 / Dy1.' : 'IEC clock number 11 (e.g. Dyn11): LV line voltage leads HV by 30° (LV phasor at "11 o\'clock" when HV is at 12). Achieved by reversing the Δ connection relative to ANSI.'));
      else L.push('Y-Y / Δ-Δ (and V-V): no phase shift between HV and LV line voltages (clock 0).');
      L.push('Per-phase (Y-equivalent, line-to-neutral) model referred to HV: Z = ' + (r.hY ? 'Z_winding (Y on HV)' : 'Z_winding/3 (Δ→Y)') + ' = ' + fx(r.ZYH.re) + ' + j' + fx(r.ZYH.im) + ' Ω;  referred to LV × (V_LL/V_HL)² = ' + fx(r.ZYL.re) + ' + j' + fx(r.ZYL.im) + ' Ω');
      L.push('Per unit: Z_base(HV) = V_HL²/S = ' + fx(r.Zbase) + ' Ω → Z = ' + fx(r.Zpu.re) + ' + j' + fx(r.Zpu.im) + ' pu (same on either side). Per-phase analysis then uses V_line/√3 and ignores the 30° shift (add it back at the end).');
      ui.work.set(L);
    }
    function draw(r) {
      const g = ui.cv.prep(), ctx = g.ctx, w = g.w, h = g.h, T = theme(); ctx.clearRect(0, 0, w, h); ctx.fillStyle = T.bg; ctx.fillRect(0, 0, w, h); if (w < 120) return;
      const cx = w / 2, cy = h / 2 + 8, R = Math.min(w, h) / 2 - 34;
      ctx.strokeStyle = T.border; ctx.lineWidth = 1; ctx.beginPath(); ctx.arc(cx, cy, R, 0, 2 * PI); ctx.stroke();
      const ph = (ang, len, col, lab, dash) => {
        const a = ang * PI / 180, x1 = cx + len * Math.cos(a), y1 = cy - len * Math.sin(a); ctx.strokeStyle = col; ctx.fillStyle = col; ctx.lineWidth = 2.2; ctx.setLineDash(dash || []);
        arrowHead(ctx, cx, cy, x1, y1, 9); ctx.setLineDash([]); ctx.font = '11px ' + T.mono; ctx.textAlign = 'center'; ctx.textBaseline = 'middle'; ctx.fillText(lab, cx + (len + 16) * Math.cos(a), cy - (len + 14) * Math.sin(a));
      };
      const base = 90; // V_AB drawn at 12 o'clock
      ['AB', 'BC', 'CA'].forEach((n, k) => ph(base - 120 * k, R, T.c, 'V' + n));
      ['ab', 'bc', 'ca'].forEach((n, k) => ph(base - 120 * k + r.shift, R * 0.72, T.a, 'V' + n, [6, 3]));
      ctx.fillStyle = T.text; ctx.font = '11px ' + T.mono; ctx.textAlign = 'left'; ctx.textBaseline = 'top';
      ctx.fillText('HV (solid) at 12 o\'clock; LV (dashed) ' + (r.shift === 0 ? 'in phase' : r.shift < 0 ? 'lags 30°' : 'leads 30°'), 6, 5);
      ctx.fillStyle = T.muted; ctx.fillText('lengths normalised; ABC sequence, CCW rotation', 6, 19);
      if (r.shift !== 0) { ctx.strokeStyle = T.pink; ctx.lineWidth = 1.2; ctx.beginPath(); const a0 = -base * PI / 180, a1 = -(base + r.shift) * PI / 180; ctx.arc(cx, cy, R * 0.45, Math.min(a0, a1), Math.max(a0, a1)); ctx.stroke(); }
    }
    const io = stateIO(st, sels);
    return { update, get: io.get, set: io.set };
  }

  /* ---------- 5. all-day efficiency ---------- */
  function buildAD(root) {
    const st = {}, sels = {}, ui = {}, { ctl, stage } = layout(root);
    ctl.appendChild(FSP.ui.el('div', { class: 'note', text: 'Regulation, OC/SC tests and η(load) at one operating point are in the Magnetics & Transformer tab. This adds the all-day (energy) efficiency over a 24-h load cycle.' }));
    const fs1 = FSP.ui.fieldset(ctl, 'Transformer');
    makeSliders(fs1, [{ k: 'S', l: 'S rated', min: 0.1, max: 1e5, v: 100, u: 'kVA', log: true }, { k: 'Pc', l: 'P core', min: 0, max: 1000, v: 1, step: 0.001, u: 'kW' },
      { k: 'Pcu', l: 'P cu (FL)', min: 0.001, max: 1000, v: 1.5, u: 'kW', log: true }, { k: 'pf', l: 'pf for η(x)', min: 0.05, max: 1, v: 0.8, step: 0.01, u: '' }], st, update);
    const fs2 = FSP.ui.fieldset(ctl, '24-h load profile');
    fs2.appendChild(FSP.ui.el('div', { class: 'note', text: 'One row per period: hours  kW  pf. Remaining hours are no-load (core loss only — the transformer stays energised).' }));
    ui.prof = FSP.ui.el('textarea', { rows: '6', spellcheck: 'false', 'aria-label': 'load profile' }); ui.prof.value = '6 80 0.8\n10 40 1\n8 0 1';
    ui.prof.style.cssText = 'width:100%;box-sizing:border-box;background:var(--panel2);color:var(--text);border:1px solid var(--border);border-radius:6px;font:12px var(--mono);padding:4px';
    ui.prof.addEventListener('change', () => { update(); FSP.state.touch(); }); fs2.appendChild(ui.prof);
    ui.msg = msgBox(ctl);
    ui.cv1 = wrapCanvas(stage, '24-hour load profile', 220); ui.cv1.onResize(update);
    ui.cv2 = wrapCanvas(stage, 'Efficiency versus load fraction', 240); ui.cv2.onResize(update);
    ui.hud = FSP.ui.el('div', { class: 'hud' }); stage.appendChild(ui.hud);
    ui.warn = FSP.ui.el('div', { class: 'msg warn', hidden: '' }); stage.appendChild(ui.warn);
    ui.work = FSP.ui.working(stage);
    function update() {
      if (!ui.cv1) return;
      const T = theme(), S = st.S.get(), Pc = st.Pc.get(), Pcu = st.Pcu.get(), pf = st.pf.get(), pr = E.parseProfile(ui.prof.value);
      const ad = pr.ok ? E.allDay({ S, Pcore: Pc, Pcu, rows: pr.rows }) : pr;
      ui.msg.show(ad.ok ? '' : ad.msg);
      const me = E.maxEffPoint(S, Pc, Pcu, pf);
      // efficiency curve
      const xs = [], es = []; for (let j = 1; j <= 150; j++) { const x = j / 100; xs.push(x); es.push(100 * E.etaAt(S, Pc, Pcu, x, pf)); }
      const emin = Math.max(0, Math.min.apply(null, es.filter(fin)) - 1);
      drawPlot(ui.cv2, { xmin: 0, xmax: 1.5, ymin: emin, ymax: 100, xlabel: 'load fraction x = S_load/S_rated', ylabel: 'efficiency η (%) at pf ' + fx(pf, 3),
        series: [{ x: xs, y: es, color: T.c, width: 2.2 }], vlines: [{ x: me.x, label: 'x* = ' + fx(me.x, 4), color: T.warn }, { x: 1, color: T.muted }],
        points: [{ x: me.x, y: 100 * me.eta, label: 'η max ' + fx(100 * me.eta, 4) + '%', color: T.warn, below: true }].concat(ad.ok ? ad.rows.filter(r => r.kW > 0).map(r => ({ x: r.x, y: 100 * E.etaAt(S, Pc, Pcu, r.x, r.pf), color: T.a })) : []) });
      if (!ad.ok) { drawPlot(ui.cv1, {}); hudSet(ui.hud, [['max-efficiency load x*', fx(me.x, 4)]]); ui.work.set(ad.msg); ui.warn.hidden = true; return; }
      // profile bars (kVA fraction) over 24 h
      const bars = []; let t0 = 0; ad.rows.forEach(r => { bars.push({ x0: t0, x1: t0 + r.h, y: r.x, color: withAlpha(T.c, 0.7) }); t0 += r.h; });
      if (t0 < 24) bars.push({ x0: t0, x1: 24, y: 0, color: T.muted });
      const ymax = Math.max(1.1, Math.max.apply(null, ad.rows.map(r => r.x)) * 1.15);
      drawPlot(ui.cv1, { xmin: 0, xmax: 24, ymin: 0, ymax, xlabel: 'hour of day', ylabel: 'load fraction x', bars,
        hlines: [{ y: me.x, label: 'x* (max η)', color: T.warn }, { y: 1, label: 'rated', color: T.muted }] });
      ui.warn.hidden = !ad.overload; ui.warn.textContent = ad.overload ? 'A period exceeds rated kVA (x > 1).' : '';
      hudSet(ui.hud, [['energy out', fx(ad.Eout, 5) + ' kWh'], ['core loss energy (24 h)', fx(ad.Ecore, 5) + ' kWh'], ['copper loss energy', fx(ad.Ecu, 5) + ' kWh'],
        ['all-day efficiency', fx(100 * ad.eta, 5) + ' %'], ['max-η load fraction x*', fx(me.x, 4) + '  (' + fx(me.Skva, 4) + ' kVA)'], ['η max at pf ' + fx(pf, 3), fx(100 * me.eta, 5) + ' %'],
        ['idle hours', fx(ad.idle, 4) + ' h']]);
      const L = ['All-day (energy) efficiency η_AD = E_out / (E_out + E_core + E_cu) over 24 h.',
        'Core loss runs all 24 h: E_core = P_core × 24 = ' + fx(Pc) + ' × 24 = ' + fx(ad.Ecore) + ' kWh'];
      ad.rows.forEach((r, i) => L.push('  period ' + (i + 1) + ': ' + fx(r.h) + ' h at ' + fx(r.kW) + ' kW, pf ' + fx(r.pf) + ' → S = ' + fx(r.kW / r.pf) + ' kVA, x = ' + fx(r.x) + ';  E_out = ' + fx(r.Eout) + ' kWh;  E_cu = x²·P_cu·h = ' + fx(r.x) + '² × ' + fx(Pcu) + ' × ' + fx(r.h) + ' = ' + fx(r.Ecu) + ' kWh'));
      if (ad.idle > 0) L.push('  remaining ' + fx(ad.idle) + ' h: no load (E_out = 0, E_cu = 0)');
      L.push('E_out = ' + fx(ad.Eout) + ' kWh;  E_cu = ' + fx(ad.Ecu) + ' kWh;  η_AD = ' + fx(ad.Eout) + '/(' + fx(ad.Eout) + ' + ' + fx(ad.Ecore) + ' + ' + fx(ad.Ecu) + ') = ' + fx(100 * ad.eta) + ' %');
      L.push('Maximum efficiency where copper loss = core loss: x*² P_cu = P_core → x* = √(P_core/P_cu,FL) = √(' + fx(Pc) + '/' + fx(Pcu) + ') = ' + fx(me.x) + '  (' + fx(me.Skva) + ' kVA)');
      L.push('η_max = x*·S·pf / (x*·S·pf + 2 P_core) = ' + fx(100 * me.eta) + ' % at pf ' + fx(pf) + '.  Distribution transformers (lightly loaded most of the day) are designed with low core loss, i.e. small x*.');
      ui.work.set(L);
    }
    const io = stateIO(st, sels, { get: () => ({ prof: ui.prof.value.split(/\r?\n/).join(';') }), set(o) { if (typeof o.prof === 'string' && o.prof.length < 4000) ui.prof.value = o.prof.split(';').join('\n'); } });
    return { update, get: io.get, set: io.set };
  }

  /* ---------- tab ---------- */
  const host = { active: false, sub: 'mc' };
  let subs = null; const panels = {}, btns = {};
  function showSub(id, silent) {
    host.sub = id;
    Object.keys(panels).forEach(k => { panels[k].hidden = k !== id; btns[k].setAttribute('aria-pressed', String(k === id)); btns[k].classList.toggle('active', k === id); });
    if (host.active) subs[id].update();
    if (!silent) FSP.state.touch();
  }
  FSP.registerTab({
    id: 'edm2', title: 'Transformers & Actuators',
    init(panel) {
      const bar = FSP.ui.el('div', { class: 'seg row', role: 'group', 'aria-label': 'Transformers and actuators topic' }); panel.appendChild(bar);
      const defs = [['mc', 'Multi-limb magnetic circuit'], ['ie', 'Inductance, force & losses'], ['auto', 'Autotransformer'], ['tp', 'Three-phase banks'], ['ad', 'All-day efficiency']];
      defs.forEach(d => { btns[d[0]] = FSP.ui.el('button', { type: 'button', class: 'seg-btn', 'aria-pressed': 'false', text: d[1], onclick: () => showSub(d[0]) }); bar.appendChild(btns[d[0]]); });
      defs.forEach(d => { panels[d[0]] = FSP.ui.el('div', { hidden: '' }); panel.appendChild(panels[d[0]]); });
      subs = { mc: buildMC(panels.mc), ie: buildIE(panels.ie), auto: buildAuto(panels.auto), tp: buildTP(panels.tp), ad: buildAD(panels.ad) };
      const pre = { mc: 'mc_', ie: 'ie_', auto: 'at_', tp: 'tp_', ad: 'ad_' };
      FSP.state.bind('edm2', {
        get() { const o = { sub: host.sub }; Object.keys(subs).forEach(s => { const g = subs[s].get(); Object.keys(g).forEach(k => { o[pre[s] + k] = g[k]; }); }); return o; },
        set(o) {
          o = o || {};
          Object.keys(subs).forEach(s => { const r = {}; Object.keys(o).forEach(k => { if (k.indexOf(pre[s]) === 0) r[k.slice(pre[s].length)] = o[k]; }); try { subs[s].set(r); } catch (e) { /* ignore bad state */ } });
          showSub(panels[o.sub] ? o.sub : 'mc', true);
        },
      });
      showSub(host.sub, true);
    },
    activate() { host.active = true; if (subs) subs[host.sub].update(); },
    deactivate() { host.active = false; },
  });
})();
