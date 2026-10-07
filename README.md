# Field & Signal Playground

A single-file, zero-dependency interactive engineering simulation web app with two modules:

1. **EM Wave & Dielectric Boundary Simulator**: Fresnel coefficients, Brewster and critical angles, total internal reflection with evanescent decay, and Poynting vector visualization for plane waves at a planar interface.
2. **DSP & Biquad Filter Lab**: RBJ cookbook biquads (LP/HP/BP/Notch/Peaking), an interactive z-plane pole-zero editor, Bode and impulse-response plots, and a Web Audio harness with a live FFT analyzer.

## Status

`index.html` is generated and committed. Open it directly in a browser; run `runSelfTests()` in the console (or load `index.html?selftest=1`) to print the acceptance checks. All 13 checks pass in Chromium and in Node (core math only).

## Repo layout

```
.
├── README.md
├── .gitignore
├── prompts/
│   └── field-signal-playground.prompt.md   # the full spec to give the model
└── index.html                              # generated output goes here
```

## Workflow

1. Paste `prompts/field-signal-playground.prompt.md` into Claude Fable.
2. Save the single code block it returns as `index.html` in the repo root.
3. Open `index.html` in a browser (no build step, no server needed for Module 1; use `python -m http.server` or `localhost` if the microphone is blocked on `file://`).
4. Open the browser console and run `runSelfTests()` to check the acceptance tests.
5. Commit the result.

## Notes

- Web Audio, AudioWorklet and `getUserMedia` need a secure context (`https://` or `localhost`) in most browsers.
- The audio chain defaults to a low master gain and a limiter. Keep volume low when testing with white noise and high-Q peaking filters.
- If the model's output truncates, ask it to continue from the exact last line rather than regenerating.
