# Architecture and maintenance

The app stays build-free and works from `file://`. Classic script order in `index.html` is the dependency list for both the page and the numerical test runner. Runtime code has no package dependencies; Playwright is a development dependency only.

## Layers

`js/core/fsp.js` owns tab/test registration, URL state, and small interface helpers. `model-notes.js` supplies the assumptions shown inside each tool. FFT, PRNG, and the original EM/biquad core are separate scripts.

The four largest subject modules now have three files each:

| Family | Math | Self-tests | Interface |
|---|---|---|---|
| Fields | `js/math/fields.js` | `js/selftests/fields.js` | `js/modules/fields.js` |
| DSP design | `js/math/dsp2.js` | `js/selftests/dsp2.js` | `js/modules/dsp2.js` |
| Transformers/actuators | `js/math/edm2.js` | `js/selftests/edm2.js` | `js/modules/edm2.js` |
| Rotating machines | `js/math/machines2.js` | `js/selftests/machines2.js` | `js/modules/machines2.js` |

Math remains available under the existing `FSP.math` names. A module's `_shared` property contains private closure helpers needed by its own interface and tests; it is not a supported API. Load math before its self-tests and interface. Interface scripts are marked `data-dom` so the Node runner skips them. Other subject modules retain their original single-file structure; they can be extracted using the same pattern when changed.

`js/v1-ui.js` retains the original EM/DSP interfaces and global navigation. Subject groups must be revealed before focus moves to a tab in another group.

## Saved links

Hash values are strings and are untrusted. The decoder uses dictionaries without prototypes, rejects `__proto__`, `constructor`, and `prototype` as tab or parameter names, and restores only registered own bindings. Each tool still validates its own values.

Custom biquad links save the normalized five coefficients `[b0,b1,b2,a1,a2]`, custom mode, normalization option, impulse length, and magnitude auto-range. Restore requires finite bounded coefficients, stable poles within the editor limit, and zeros within the editor limit. Invalid custom data falls back to the selected cookbook filter. Restore uses the saved coefficients directly so it does not normalize the filter a second time. Older cookbook links remain valid. Links never grant permission or start audio.

## Audio lifecycle

Audio Tools uses a start-generation counter. Stop, deactivation, and source changes invalidate pending work. A stream returned to an obsolete request is stopped immediately, before it can be attached. A pending start exposes a Cancel button. File decodes use a separate generation so selecting a second file cannot be overwritten by an earlier decode.

The original DSP lab tracks both context and source generations. Butterworth/Chebyshev tracks context generations. Worklet completions must match both the original context and current generation before changing a graph or status. This also protects stop/restart sequences where an old promise resolves after a new context exists.

## Verification

Run `node tests/run.js` for numerical and saved-state checks. `--with js/modules/fields.js` (or `dsp2.js`, `edm2.js`, `netproto.js`) retains the focused-test interface, including extracted math and self-tests.

Run `npm ci`, `npx playwright install chromium`, and `npm run test:browser` for actual visible navigation, all top-level subtopics, keyboard and mobile selection, filter-response round trips, malformed links, and controlled audio races. The test does not use a physical microphone or verify perceptual audio quality. Missing tooling fails rather than producing a successful skip.

When adding a tool, register its numerical API and self-tests without accessing the DOM at load time. Register the interface with `FSP.registerTab`, bind validated state, add model assumptions to `FSP.modelNotes`, and place its ID in the appropriate subject group. Update the browser tool count and verify the new interface through visible controls.
