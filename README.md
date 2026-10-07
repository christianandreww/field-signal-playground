# Field & Signal Playground

A zero-dependency, build-free engineering sandbox that runs by double-clicking `index.html` (works from `file://`). Nine tabs, all with a "Show working" panel and shareable URL state ("Copy link").

| Tab | What it does |
|---|---|
| EM Wave & Boundary | Oblique incidence on a planar interface: Fresnel (impedance form, lossy/magnetic media), TIR, Brewster, live field map, Poynting flux, skin depth and surface impedance |
| DSP & Biquad Lab | RBJ biquads, draggable pole-zero editor, Bode/impulse, live audio through an AudioWorklet |
| Convolution & Spectra | Flip-and-slide convolution, FFT, window explorer (sidelobe/ENBW/coherent gain), spectral leakage |
| Transmission Lines | Smith chart, VSWR/return loss, Zin along the line, stub and quarter-wave matching |
| Butterworth / Chebyshev | Order 1–10 LP/HP/BP/BS as second-order-section cascades (bilinear with prewarping), s/z-plane views, audio path |
| Audio Tools | YIN pitch tuner, spectrogram, tone generator, audio-file player |
| Comms & Networks | BPSK/QPSK/16-QAM BER Monte Carlo vs theory, raised-cosine eye diagram, Shannon calculator, CIDR/subnet calculator |
| Magnetics & Transformer | Magnetic circuit with air gap and B–H saturation; transformer equivalent circuit, regulation, efficiency |
| Rotating Machines | Induction motor torque–speed, DC machine, three-phase rotating field |

## Run

- Open `index.html` in a browser. No server, no install.
- Microphone input needs `https` or `localhost` (browser rule); everything else works from `file://`.

## Test

```
node tests/run.js        # all math/self-tests, no browser needed (Node 18+)
node tests/browser.js    # optional: headless smoke test of every tab (needs Playwright + Chromium)
```

In the page, append `?selftest` to the URL (or run `runSelfTests()` in the console) to run the same suite.

## Layout

```
index.html            shell + script tags (order matters; no ES modules, so file:// works)
css/style.css
js/core/              fsp.js (registry, URL state, UI helpers), prng.js, fft.js, v1-math.js (EM + biquad core, v1 tests)
js/modules/           one file per tab: pure math in FSP.math.<name>, UI in init(), tests via FSP.registerTests
js/v1-ui.js           EM and DSP tab UI, tab switching, state restore
tests/                run.js (Node), browser.js (Playwright smoke)
prompts/              the prompts used to specify and generate this project
```

## Notes and limits

- Audio (worklet/IIR chains, mic, tuner) cannot be exercised headlessly; those paths were checked for load/render and the pure math is unit tested.
- v1 noise comes from a never-repeating AudioWorklet generator (a 10 s crossfaded loop is used only when AudioWorklet is unavailable); the master output is a tanh soft clip (0.9 ceiling) rather than a compressor.
- Machines models are teaching-level: induction motor has no separate core-loss branch, DC machine ignores armature reaction.
