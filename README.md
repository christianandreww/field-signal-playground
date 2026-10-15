# Field & Signal Playground

A zero-dependency, build-free engineering sandbox that runs by double-clicking `index.html` (works from `file://`). Thirteen tools grouped by subject (Electromagnetics, Signals & DSP, Communications, Machines), all with a "Show working" panel and shareable URL state ("Copy link").

| Tab | What it does |
|---|---|
| EM Wave & Boundary | Oblique incidence on a planar interface: Fresnel (impedance form, lossy/magnetic media), TIR, Brewster, live field map, Poynting flux, skin depth and surface impedance |
| DSP & Biquad Lab | RBJ biquads, draggable pole-zero editor, Bode/impulse, live audio through an AudioWorklet |
| Convolution & Spectra | Flip-and-slide convolution, FFT, window explorer (sidelobe/ENBW/coherent gain), spectral leakage |
| Transmission Lines | Smith chart, VSWR/return loss, Zin along the line, stub and quarter-wave matching |
| Butterworth / Chebyshev | Order 1–10 LP/HP/BP/BS as second-order-section cascades (bilinear with prewarping), s/z-plane views, audio path |
| Audio Tools | YIN pitch tuner, spectrogram, tone generator, audio-file player |
| Comms & Networks | BPSK/QPSK/16-QAM BER Monte Carlo vs theory, raised-cosine eye diagram, Shannon calculator, CIDR/subnet and VLSM calculator |
| Magnetics & Transformer | Magnetic circuit with air gap and B–H saturation; transformer equivalent circuit, regulation, efficiency |
| Rotating Machines | Induction motor torque–speed and speed control (voltage, rotor resistance, V/f), DC machine incl. compound, three-phase rotating field |
| Networks & Protocols | Line coding (NRZ, NRZ-I, AMI, Manchester, diff. Manchester, B8ZS, HDB3), delays and store-and-forward, statistical multiplexing, ALOHA and CSMA/CD, ARQ efficiency and TCP congestion control, Internet checksum, CRC and parity (IE3017) |
| z-Transform & FIR Design | Sampling and aliasing, difference equations, pole-zero, inverse z by partial fractions with ROC, window and frequency-sampling FIR design, filter structures with block diagrams, DFT properties and quantisation (IE3014) |
| Fields & Polarization | Charge distributions (closed form vs numeric), Biot–Savart and Ampère, Faraday and displacement current, plane waves and polarization, normal incidence and standing waves (EE3101) |
| Transformers & Actuators | Three-limb magnetic circuits with a B–H table, inductance, energy and force, core losses, autotransformers, three-phase transformer connections, all-day efficiency (EE2005) |

## Run

- Open `index.html` in a browser. No server, no install.
- Microphone input needs `https` or `localhost` (browser rule); everything else works from `file://`.

## Your workspace

Choose a subject in the desktop sidebar or a tool in the mobile picker. The four subject identities are EE3101 (teal), IE3014 (violet), IE3017 (emerald), and EE2005 (copper). Each tool includes a **What to try** prompt.

The appearance selector offers Light, Dark, and System independently of the subject. This preference stays on your device when storage is available; shared links continue to describe the experiment, rather than the appearance. Appearance changes redraw plots without starting audio or resetting a filter.

The wave tool leads with its visualization and three key results. Start with an example, adjust the angle and polarization, and expand **Material properties**, **Display options**, or **Detailed results & angle sweep** when needed. Field arrows are off initially and can be enabled in Display options. On mobile, the graph comes first, controls follow, and detailed theory stays below them. Reduced-motion preferences start the wave animation paused.

## Test

```
node tests/run.js        # all math/self-tests, no browser needed (Node 18+)
node tests/browser.js    # browser regressions; requires the setup below
```

In the page, append `?selftest` to the URL (or run `runSelfTests()` in the console) to run the same suite.

For browser tests (Node 20+; development only; opening the app still needs no install):

```
npm ci
npx playwright install chromium
npm run test:all
```

Browser tests fail if Playwright or Chromium is missing. They exercise visible desktop navigation, keyboard focus across subjects, the mobile picker, custom-filter saved links, hostile URL state, and delayed microphone/worklet cancellation. Permission timing is mocked; no real microphone is accessed. `PW_MODULE` and `CHROMIUM` can point to an existing Playwright module or Chromium executable. GitHub Actions runs numerical checks on Node 18 and 22 and browser checks on Node 22 for pushes and pull requests.

## Layout

```
index.html            shell + script tags (order matters; no ES modules, so file:// works)
css/style.css         base controls and module styles
css/studio.css        responsive Lab Studio shell, appearance and subject palettes
js/core/              registries, safe URL state, UI helpers, model notes, PRNG, FFT, original EM/biquad core
js/math/              extracted numerical code for fields, DSP design, transformers/actuators, rotating machines
js/selftests/         extracted numerical checks and URL-state regressions (Node and browser)
js/modules/           tab interfaces; older modules still include their numerical code and self-tests
js/v1-ui.js           EM and DSP tab UI, tab switching, state restore
js/studio-ui.js       course identity, tool introductions, guided prompts, wave layout
tests/                run.js (Node), browser.js and studio.js (Playwright regressions)
prompts/              the prompts used to specify and generate this project
```

## Notes and limits

- Every tool includes a collapsible **Model assumptions & limits** note. These are teaching models; check the stated geometry, material data, sample rate, and operating assumptions before comparing against a real system.
- Browser tests exercise real Web Audio contexts and graphs with controlled asynchronous permission/worklet responses. Actual microphone hardware, listening quality, and browser permission dialogs still need manual checks.
- v1 noise comes from a never-repeating AudioWorklet generator (a 10 s crossfaded loop is used only when AudioWorklet is unavailable); the master output is a tanh soft clip (0.9 ceiling) rather than a compressor.
- Machines models are teaching-level: induction motor has no separate core-loss branch, DC machine ignores armature reaction.

See [Architecture and maintenance](docs/ARCHITECTURE.md) for script ordering, saved-state validation, audio cleanup, and adding a tool.
