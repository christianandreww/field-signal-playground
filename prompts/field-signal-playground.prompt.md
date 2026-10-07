ROLE
You are a principal computational physicist and senior systems software engineer.

DELIVERABLE (OUTPUT CONTRACT)
- Output exactly one fenced code block containing a complete, standalone `index.html`, from <!DOCTYPE html> to </html>. No prose before or after it.
- Zero dependencies: no CDNs, libraries, frameworks, web fonts or external files. HTML, CSS and vanilla JS only, all inline.
- No TODOs, placeholders, stubs or "similar to above" abbreviations. Every equation, complex-number operation, renderer and audio node is fully implemented.
- Write it for correctness first, then polish. Self-check against the ACCEPTANCE TESTS section before finishing, and fix any failure silently.

APP: "Field & Signal Playground"
Dark-mode engineering dashboard, responsive from 360px wide to desktop, two modules in an accessible tab bar (keyboard-navigable, aria-selected). Switching tabs must pause the inactive module's animation loop and stop/suspend its audio. Use requestAnimationFrame and pause when document.hidden. All canvases must handle devicePixelRatio and window resize. Use Pointer Events (not mouse events) so dragging works on touch.

======================================================
MODULE 1: EM WAVE & DIELECTRIC BOUNDARY SIMULATOR
======================================================
CONVENTIONS (state them in a collapsible "Theory" panel)
- Time dependence e^{jωt}; fields written E = Re{Ê e^{jωt}}.
- Interface at z = 0, medium 1 in z < 0, medium 2 in z > 0, plane of incidence = x–z.
- TE (s): E along y. TM (p): H along y.
- r and t are E-field amplitude coefficients. R and T are power ratios.

INPUTS
- Angle of incidence θi: 0–89°, step 0.1°.
- Polarization: TE / TM toggle.
- Frequency f (1 MHz – 1 THz, log slider) and derived free-space wavelength.
- Primary material inputs: εr1, μr1, σ1 and εr2, μr2, σ2. Show derived n = sqrt(μr·εr) read-only (do NOT offer n as a second editable input).
- Display wavelength is normalized for the canvas (about 8–14 wavelengths across the view); show physical λ in the HUD only.

PHYSICS (exact, no small-angle approximations)
- Use complex arithmetic throughout (write a small Complex class or functional helpers: add, sub, mul, div, sqrt with correct branch, conj, abs, arg, exp).
- Complex permittivity ε_c = ε0·εr − jσ/ω. Wave impedance η = sqrt(μ/ε_c). Complex index n_c = sqrt(μr·ε_c/ε0).
- Snell: n_c1 sinθi = n_c2 sinθt, with complex cosθt = sqrt(1 − sin²θt). Choose the branch so the transmitted wave decays (Im(kz2) ≤ 0 in the e^{jωt} convention), which automatically gives the evanescent solution in total internal reflection.
- Impedance-form Fresnel coefficients (valid for μ ≠ μ0):
  TE: r = (η2 cosθi − η1 cosθt)/(η2 cosθi + η1 cosθt); t = 2 η2 cosθi/(η2 cosθi + η1 cosθt)
  TM: r = (η2 cosθt − η1 cosθi)/(η2 cosθt + η1 cosθi); t = 2 η2 cosθi/(η2 cosθt + η1 cosθi)
  (state the sign convention for TM r in the Theory panel and keep it consistent in the phase readout).
- Power: R = |r|². T = Re(Sz,trans)/Re(Sz,inc), the normal Poynting flux ratio computed from the complex fields, so it stays correct in the lossy and evanescent cases. Never use T = |t|².
- Critical angle: θc = asin(n2/n1), only for lossless media with n1 > n2; otherwise display "N/A".
- Brewster angle: find the angle where r_TM = 0 (r_TE = 0 is also possible when μr1 ≠ μr2). Solve analytically from the impedance form; verify numerically; display "N/A" when no real solution exists. For lossy media, display the "pseudo-Brewster" angle (the minimum of |r_TM|) found by a golden-section search and label it as such.
- TIR: for θi > θc, |r| = 1 (lossless), T = 0, with the phase shift of r computed from the complex r. Compute the evanescent decay length δ = 1/|Im(kz2)| and show it in the HUD and as a decay-envelope overlay.
- Numerical stability: clamp sinθ to [−1, 1] before sqrt only in lossless real cases; use a small epsilon to guard denominators at θi → 0 (normal) and the 89° limit; never output NaN or Infinity (display "—" instead); guard the cosθi → 0 case.

FIELDS & POYNTING
- Compute the complex phasor field at every pixel for each medium:
  Medium 1: incident + reflected (their interference gives the standing-wave pattern).
  Medium 2: transmitted (propagating or evanescent).
  Render TE as Ey and TM as Hy (state which one is shown, with a toggle to switch the shown field between E-magnitude and H-magnitude for TM, if both are derived correctly).
- Instantaneous field = Re{Ê(x,z) e^{jωt}}. Animate t with a speed slider and a pause/play button.
- Instantaneous Poynting vector S(x,z,t) = E(x,z,t) × H(x,z,t) from the real fields (so interference cross-terms are included). Also compute time-averaged ⟨S⟩ = ½ Re(E × H*). Provide a toggle between the two for the arrows.
- Render a heat map (diverging colormap, symmetric around 0) using ImageData at a reduced internal resolution upscaled to the display size, with a quality selector. Draw a vector-arrow overlay of S on a coarse grid, normalized and length-capped.

OVERLAYS
- Interface line, normal vector, incident/reflected/transmitted ray paths with angle arcs and labels (θi, θr = θi, θt or "evanescent"), and wavefront lines.
- Medium labels and the ε, μ, σ shown in each half-plane.

HUD (live, monospaced, updated at ≥ 10 Hz)
- r (magnitude ∠ phase) and t (magnitude ∠ phase), Γ = r, τ = t.
- Reflection phase shift (degrees, unwrapped sensibly).
- R, T, and R + T with a pass/fail badge (tolerance 1e-9 lossless; for lossy media, R + T measured at the interface must still be 1 within 1e-9 because T is the flux entering medium 2).
- Brewster angle, critical angle, evanescent decay length.
- Optional small plot: R and T vs θi for both polarizations, with a marker at the current angle.

======================================================
MODULE 2: DSP & BIQUAD FILTER LAB
======================================================
FILTER MATH
- Coefficients from the RBJ Audio EQ Cookbook for Lowpass, Highpass, Bandpass (constant 0 dB peak gain), Notch and Peaking EQ. Normalize so a0 = 1 and display the six values b0, b1, b2, a0 (=1), a1, a2 (and the pre-normalization a0 as a note).
- Inputs: fc (log slider, 10 Hz – 0.49·fs, clamp), fs (default 48000; selectable 8000/16000/22050/44100/48000/96000, and match the AudioContext rate or resample note), Q (0.1–20), gain in dB (−24 to +24, used by the Peaking EQ only).
- Direct Form II Transposed:
  y[n] = b0 x[n] + s1
  s1 = b1 x[n] − a1 y[n] + s2
  s2 = b2 x[n] − a2 y[n]
  Use this one implementation for BOTH the impulse-response display and the audio path.
- Frequency response H(e^{jω}) by evaluating the numerator and denominator polynomials on the unit circle at 512+ log-spaced points (20 Hz – fs/2). Magnitude in dB with a floor at −120 dB to avoid log(0); unwrapped phase in degrees.
- Poles/zeros: solve the quadratics with complex roots. Handle the degenerate cases (b2 = 0, double roots, real roots, b0 = 0).
- Stability: compute |p|; if any |p| ≥ 1 show an "UNSTABLE" badge and do not send those coefficients to audio (keep the last stable set).

POLE–ZERO EDITOR (interactive z-plane)
- Unit circle, axes, grid, poles (×) and zeros (○), with the double-root multiplicity marked.
- Drag poles/zeros with Pointer Events. Enforce conjugate symmetry: dragging one of a conjugate pair moves its mirror; real roots slide along the real axis only; a pair can merge into two real roots, and vice versa, when crossing the axis.
- Clamp |pole| ≤ 0.9999 so dragging cannot create an unstable filter. Zeros may lie anywhere, clamped to |z| ≤ 3.
- Rebuild coefficients from the roots: a = [1, −2Re(p), |p|²] (for a conjugate pair) or from the real roots, b likewise, scaled by a gain K. Keep K so |H| at the reference frequency (DC for LP, fs/2 for HP, fc for BP/Notch/Peak) equals the same level as before the drag, and offer a "Normalize gain" toggle.
- While dragging, show a "Custom" preset indicator and decouple the sliders; moving a slider or preset button re-derives the roots from the cookbook coefficients and snaps the editor back.
- Keyboard-nudge support for the selected root (arrow keys).

PLOTS
- Magnitude (dB vs log frequency, −60 to +30 dB default with auto-range option), phase (degrees vs log frequency), and impulse response h[n] as a stem plot (first 128 samples by default, with an N slider up to 1024; auto-trim when the response decays below 1e-6). Draw a marker at fc on the Bode plots. Hover/touch cursor readout of frequency, dB and phase.

AUDIO (Web Audio API)
- Create the AudioContext only after a user gesture (a "Start audio" button), and handle the suspended state.
- Filter node: an AudioWorkletProcessor created from a Blob URL (inline code string) implementing the DF-II-T biquad above, with coefficients sent through port.postMessage (use per-sample or per-block coefficient smoothing, such as a short linear ramp, to avoid zipper noise). State a fallback: if AudioWorklet is unavailable, use IIRFilterNode recreated on change with a short crossfade, and show a notice.
- Sources (selector): white noise, pink noise (Paul Kellet's refined filter or Voss-McCartney), logarithmic frequency sweep (looping, 20 Hz → 20 kHz over an adjustable 2–20 s), and microphone via getUserMedia (echoCancellation/noiseSuppression/autoGainControl off). Handle permission-denied, insecure-context and no-device errors with visible messages. Show a feedback warning when the mic is selected and the output goes to the speakers.
- Safety: master gain default 0.2, and a DynamicsCompressor or limiter before the destination. Include a Bypass toggle (A/B).
- Wiring: source → [input AnalyserNode] → worklet filter → [output AnalyserNode] → master gain → limiter → destination. Disconnect/stop everything on tab switch or Stop, release the mic tracks, and avoid leaked nodes.
- Spectrum analyzer: fftSize 4096, smoothingTimeConstant ≈ 0.8, log-frequency x-axis, dB y-axis (−100 to 0), input (dim) and filtered output (bright) overlaid, with the theoretical |H| curve optionally overlaid. Animate at 60 FPS and compute the bins once per frame.
- Audio sample-rate handling: set coefficient math using the context's actual sampleRate when audio is running (show a note if it differs from the selected fs).

======================================================
ACCEPTANCE TESTS (verify internally; include a hidden self-test function runSelfTests() callable from the console, which logs PASS/FAIL)
EM:
1. n1 = 1, n2 = 1.5 (μr = 1): Brewster ≈ 56.3099°; at θi = 0 r_TE = −0.2, t_TE = 0.8, R = 0.04, T = 0.96.
2. n1 = 1.5, n2 = 1: critical angle ≈ 41.8103°; for θi > θc, |r| = 1, T = 0 (≤ 1e-9), and the decay length is finite.
3. R + T = 1 within 1e-9 for lossless media at 500 random (θ, ε, μ, polarization) combinations, with and without TIR.
4. θi = 0 and θi = 89° produce no NaN/Infinity for any polarization.
5. With σ2 > 0, R + T (flux through the interface) still equals 1 within 1e-9, and the pseudo-Brewster angle is reported and labelled.
DSP:
6. Lowpass at fc with Q = 0.7071: |H(fc)| = −3.01 dB ± 0.02; DC gain 0 dB; Nyquist < −30 dB at fc = fs/10.
7. Peaking EQ with +12 dB: |H(fc)| = +12.00 dB ± 0.01. Notch: |H(fc)| < −80 dB.
8. Impulse response via DF-II-T matches the inverse DFT of H(e^{jω}) to within 1e-6 for 256 samples.
9. Roots → coefficients → roots round trip within 1e-9, including a double-root case.
10. Dragging a pole to |p| ≥ 1 never yields |p| > 0.9999 in the output coefficients.
