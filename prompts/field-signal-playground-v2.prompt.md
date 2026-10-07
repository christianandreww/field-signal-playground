ROLE
You are a principal computational physicist and senior front-end/DSP engineer extending an existing, working app: "Field & Signal Playground". The v1 app (EM boundary simulator + biquad DSP lab) is in the repo as `index.html` and already passes its 13 built-in self-tests. DO NOT break v1 behavior. Read `index.html` first and reuse its existing math (Complex, makeMedium, solveInterface, rbjCoefs, freqResponse, root utilities, BIQUAD_CORE_SRC/WORKLET_SRC) instead of rewriting it.

STEP 0: REFACTOR TO MULTI-FILE (no behavior change)
Constraint: must open by double-clicking `index.html` (file://). Therefore NO ES modules, NO fetch of local files. Use classic `<script src="...">` tags in dependency order; each file attaches to one global namespace `FSP` (e.g. `FSP.complex`, `FSP.em`, `FSP.dsp`). Math files must also load under Node (`if (typeof module !== 'undefined') module.exports = FSP;`-style guard) with no DOM access at load time.
Layout:
- index.html (shell, tabs, CSS link, script tags)
- css/style.css
- js/core/{complex,util,prng,fft,state}.js   (pure math; util has formatting + SI parsing; state = URL-hash codec)
- js/modules/{em,tline,dsp,filters,audio,comms,machines}.js  (each: pure `math` part + `ui` part registered via `FSP.registerTab({id,title,init,activate,deactivate,tests})`)
- tests/run.js  (Node harness: loads core + modules' math, runs all registered tests, prints PASS/FAIL, exits non-zero on failure)
Tab switching must come from the registry (no hard-coded tab names). `FSP.runSelfTests()` aggregates every module's tests; the in-page "Run self-tests" button and `node tests/run.js` run the same suite. Phase 0 is done only when all 13 v1 tests still pass in both places.

OUTPUT CONTRACT
- Work in phases below. After EACH phase: run `node tests/run.js`, fix failures, then commit. Do not start a phase with failing tests.
- Complete code, no TODOs/stubs/placeholders. Zero dependencies, no network, no CDN, no build step.
- Numerical hygiene: guard divisions, clamp acos/asin arguments, no NaN/Infinity ever reaching the UI (show "—" instead). Canvas must handle devicePixelRatio and resize.
- Every module has a "Show working" collapsible panel that prints the formulas with the live numeric values substituted (intermediate quantities, not just the final answer).
- If git is blocked in your session, say so and skip commits; do not loop on it.
- Final answer: short summary table (phase, files, tests passing) and any known limitations. Do not paste whole files.

PHASE 1: URL STATE + SHARED UI
- `FSP.state`: serialize each tab's parameters to `location.hash` (compact `key=value&...`, per-tab prefix), restore on load, debounce updates (200 ms), "Copy link" button. Unknown/invalid keys are ignored, never throw.
- Shared helpers: slider+numeric input pair with SI-suffix parsing ("4.7k", "10u"), unit labels, tooltip, reset button.
TEST: round-trip encode/decode of a full state object is lossless; garbage hash does not throw.

PHASE 2: EM ADDITIONS (extend existing tab)
- Skin depth δ = sqrt(2/(ωμσ)) and surface impedance Zs = (1+j)/(σδ) readouts for any medium with σ>0 (use the exact complex-k value alongside the good-conductor approximation and show their ratio).
- Power balance display: R+T=1 check shown live with residual.
TESTS: copper (σ=5.8e7, μr=1) at 1 MHz: δ = 66.09 µm (±0.1%). Zs magnitude = sqrt(ωμ/(2σ)) ·√2 consistent with (1+j)/(σδ). Lossless dielectric interface: |R+T−1| < 1e-12 at 20 random angles/polarizations.

PHASE 3: TRANSMISSION LINES + SMITH CHART (tab `tline`)
- Inputs: Z0 (real or complex), load ZL (R+jX), frequency, line length (in λ or meters), velocity factor; optional single shunt stub / quarter-wave transformer matching helper.
- Outputs: Γ_L, VSWR, return loss, mismatch loss, input impedance Zin(l) = Z0 (ZL + jZ0 tan βl)/(Z0 + jZL tan βl) (use the Γ form to avoid tan singularities), V/I standing-wave plot along the line, Smith chart with constant-|Γ| circle, locus of Zin as l varies, click-to-set load, wavelengths-toward-generator scale.
TESTS: ZL=Z0 → Γ=0, VSWR=1. ZL=short → |Γ|=1, Zin at λ/4 is open (|Zin|>1e6·Z0 or reported "∞"). Z0=50, ZL=100: Γ=1/3, VSWR=2, RL=9.542 dB (±0.001). Quarter-wave transformer between 50 Ω and 100 Ω needs Zt=70.71 Ω and gives Γ=0 at design frequency. Zin is periodic with λ/2 (error < 1e-9).

PHASE 4: CONVOLUTION / DFT / WINDOWS LAB (tab `dsp`, new sub-section) + `core/fft.js`
- Implement iterative radix-2 FFT/IFFT (pure functions, in-place on Float64Array pairs) plus a direct O(N²) DFT used only for tests.
- Convolution visualizer: two user-editable sequences (type values or draw), animated flip-and-slide with the running sum, result stem plot; show linear vs circular convolution.
- Window explorer (rect, Hann, Hamming, Blackman) with time-domain shape, magnitude spectrum in dB (zero-padded x64), and a readout of peak sidelobe, ENBW, coherent gain. Spectrum analyzer of a user sinusoid (frequency can be off-bin) to show leakage with/without window.
TESTS (N=1024, periodic windows w[n]=…2πn/N): FFT vs direct DFT max error < 1e-9 for N=64 random input; IFFT(FFT(x))=x within 1e-12; Parseval holds within 1e-9; convolution of [1,2,3] and [0,1,0.5] = [0,1,2.5,4,1.5]; FFT-based linear convolution equals direct within 1e-9. Window figures (tolerance ±0.05 dB / ±0.005): Rect sidelobe −13.26 dB, ENBW 1.000, CG 1.000; Hann −31.47 dB, 1.500, 0.500; Hamming −42.67 dB, 1.363, 0.540; Blackman −58.11 dB, 1.727, 0.420.

PHASE 5: HIGHER-ORDER FILTERS (tab `filters`)
- Butterworth and Chebyshev-I (ripple input) low-pass/high-pass/band-pass/band-stop, order 1–10. Design: analog prototype poles → frequency transform → bilinear transform with prewarping → pair into biquad sections (SOS cascade, ordered for stability/dynamic range). Never form a high-order polynomial transfer function.
- Plots: magnitude (dB) and phase, group delay, step response, s-plane (analog) and z-plane (digital) pole-zero views, per-section response toggles, SOS coefficient table with copy-as-JSON.
- Audio path: feed the SOS cascade through the existing worklet (one biquad per section, same smoothing/ramping) with IIRFilterNode fallback (max 20 feedback coefs → chain several nodes if needed).
TESTS (fs=48 kHz, fc=1 kHz, low-pass): Butterworth order 2/4/6 → −3.010 dB at fc (±0.01) and −24.48 / −48.92 / −73.38 dB at 4 kHz (±0.1 dB). Chebyshev-I order 4, 1 dB ripple: |H| at fc = −1.000 dB (±0.01), −60.58 dB at 4 kHz (±0.3 dB), DC gain −1 dB for even order. All poles inside the unit circle for every order/type/band. Magnitude of cascade equals product of section magnitudes (error < 1e-12).

PHASE 6: AUDIO TOOLS (tab `audio`; reuse the v1 graph)
- Fix known v1 weaknesses first: replace the 2 s looped noise buffer with a continuous noise generator (small AudioWorklet or ScriptProcessor-free approach) so there is no audible loop; replace the compressor-as-limiter with a soft clip (tanh WaveShaper) before the master gain.
- YIN pitch detector (de Cheveigné) on the mic or file input: difference function, cumulative mean normalization, absolute threshold 0.10, parabolic interpolation, fmin 50 Hz / fmax 1000 Hz, return null when unvoiced. Tuner display (note name, cents, needle).
- Scrolling spectrogram (log-frequency axis option, dB range controls) from an AnalyserNode.
- Reference tone generator (sine/square/saw, frequency, level) and audio-file loading (`<input type=file>` → decodeAudioData → looped player) as sources.
- getUserMedia requires a secure context; on file:// show a clear message and keep everything else working.
TESTS (pure YIN function, fs=48 kHz, 4096 samples, tolerance 0.1% unless stated): sine 110/440/880 Hz → same frequency; sawtooth-like harmonic series Σ sin(2π·220k t)/k, k=1..5 → 220 Hz (no octave error); sine 440 Hz + Gaussian noise σ=0.22 (seeded PRNG) → 440 ±2%; pure white noise (seeded) → null. Cents conversion: 440→0¢ on A4, 466.16 Hz → +100 ¢ ±0.1.

PHASE 7: COMMUNICATIONS LAB (tab `comms`) + `core/prng.js`
- Seeded PRNG (mulberry32 or xoshiro128**) and Box–Muller Gaussian; every simulation takes a seed so results are reproducible.
- Constellation/BER: BPSK, QPSK, 16-QAM with Gray mapping, AWGN channel at chosen Eb/N0, live received scatter, Monte Carlo BER with confidence indication (stop at ≥100 errors or N max), overlaid against theoretical curves (BPSK/QPSK: Q(√(2Eb/N0)); 16-QAM: (3/8)·erfc(√(0.4·Eb/N0)) approx).
- Eye diagram: raised-cosine pulse shaping with roll-off β slider, sampling jitter and noise controls, overlaid traces over 2 symbol periods.
- Shannon calculator: C = B·log2(1+SNR), with bidirectional solve (given any two of B, SNR, C find the third).
- CIDR/subnet calculator: IPv4 address + prefix → network, broadcast, mask, wildcard, host range, host count (handle /31 and /32 correctly), plus a split-into-N-subnets helper. Input validation with inline errors.
TESTS: BPSK theoretical BER at Eb/N0 = 0/4/8 dB = 0.07865 / 0.012501 / 1.909e-4 (±0.1%). Monte Carlo BPSK at 4 dB with 2e5 bits and fixed seed lands within 15% of theory and is bit-identical across two runs. Gray mapping: adjacent constellation points differ by exactly 1 bit. Raised-cosine pulse has zeros at nonzero integer symbol times (|p|<1e-9) for β=0.25/0.5/1 (handle the t=±1/(2βT) singularity). Shannon: B=1 MHz, SNR=30 dB → 9.967 Mbit/s (±0.01%). Subnetting: 192.168.1.130/26 → network 192.168.1.128, broadcast 192.168.1.191, 62 hosts; /31 → 2 usable (RFC 3021); /32 → 1; 10.0.0.0/8 → 16,777,214 hosts; invalid inputs (octet 256, prefix 33) rejected.

PHASE 8: MACHINES & MAGNETICS (tab `machines`; largest module, build as 5 sub-tabs, each with its own tests)
8a Magnetic circuit: series core + air gap (optional parallel leg), mean length, area, μr from either a constant or a user-editable B–H table (piecewise-linear in log-H or monotone cubic; ship a default silicon-steel table). Given NI or flux, solve for the other (nonlinear: bisection/Newton on H with bounded iterations, report convergence). Show reluctance network, B, H, flux, inductance L = N²/R_total, stored energy, saturation indicator.
8b Transformer: single-phase equivalent circuit (R1, X1, Rc, Xm, R2', X2'), open-circuit and short-circuit test inputs → derive parameters, then load analysis: voltage regulation, efficiency vs load curve with max-efficiency point (where copper loss = core loss), phasor diagram.
8c Induction motor: per-phase equivalent circuit via Thevenin; torque–speed curve over slip −1…2, markers for starting torque, breakdown torque and slip, rated point; synchronous speed ns = 120f/P; mode labels (motor/generator/braking).
8d DC machine: separately/shunt/series excitation, armature reaction ignored; speed–torque and efficiency curves; back-EMF and power flow.
8e Three-phase & rotating field: three phasors/waveforms, line vs phase quantities (√3 relations, Y/Δ), animated rotating MMF vector from three 120°-spaced windings; direction reversal by swapping two phases.
TESTS: Core N=500, l=0.3 m, A=1e-3 m², μr=2000, gap 1 mm → L = 0.2732 H (±0.2%); with gap=0 L is larger by the correct ratio; B–H solve returns a flux whose forward re-evaluation reproduces NI within 1e-6 relative. Transformer: results from OC/SC test inputs, when fed back into the load analysis at rated load, reproduce the rated secondary voltage within 0.5%; max-efficiency load fraction equals sqrt(Pcore/Pcu,rated) (±1e-6). Induction motor: ns for 4 poles/50 Hz = 1500 rpm; torque at s=0 is 0; torque is positive for 0<s<1 and negative for s<0; T(s) from the Thevenin form equals T computed by solving the full circuit directly (error <1e-9) at 5 slips; breakdown slip = R2'/sqrt(Rth²+(Xth+X2')²). DC: power balance Pin = Pout + losses within 1e-9; separately-excited speed falls linearly with torque. Three-phase: phase sum = 0 at every sample (|ε|<1e-12); resultant MMF magnitude = 1.5×single-phase peak and rotates at the supply frequency; swapping two phases reverses rotation.

PHASE 9: POLISH
- Responsive layout (≥360 px wide), keyboard-focusable controls, dark/light respecting `prefers-color-scheme`, no console errors, no layout shift when switching tabs.
- README update: module list, how to run tests (`node tests/run.js`), how to open (double-click index.html), known limitations (mic needs https/localhost; audio features cannot be exercised headlessly).
- Final: run the full suite, list the count of passing tests per module.

DO NOT
- Do not use ES modules, bundlers, TypeScript, or any third-party library.
- Do not modify the verified v1 Fresnel/biquad math except to move it into the new files.
- Do not weaken a test tolerance to make it pass; if a test fails, debug the code (and say so if you believe the expected value is wrong, with your own derivation).
