/* Teaching-model scope, kept beside each tool via FSP.ui.modelNotes(). */
FSP.modelNotes = {
  em: [
    'Uniform, isotropic media separated by an infinite planar boundary; time convention exp(jωt). Material values are constant at the chosen frequency.',
    'The field map uses normalized wavelength and display time. Read physical wavelength, flux ratios, and decay length from the numerical readouts.',
  ],
  dsp: [
    'A linear second-order digital filter with real coefficients. Pole edits are limited to radius 0.9999; zeros to radius 3.',
    'Plots describe the filter before output soft clipping. Clipping is nonlinear, so high-level audio can differ from the plotted response.',
    'Saved links retain custom coefficients and editor options. They do not start audio or grant microphone access.',
  ],
  dsplab: [
    'Finite sampled sequences are used. Circular convolution wraps at the selected length; linear convolution includes the full output length.',
    'FFT resolution and leakage depend on record length and window choice. Zero padding adds plot samples, not new resolving power.',
  ],
  tline: [
    'A uniform transmission line with specified characteristic impedance and attenuation. Length is measured in wavelengths.',
    'Stub and quarter-wave matching are single-frequency teaching designs; practical bandwidth and component parasitics are not modeled.',
  ],
  filters: [
    'Butterworth and Chebyshev designs are transformed to digital second-order sections using a bilinear transform with prewarping.',
    'Audio is redesigned at the actual AudioContext sample rate. A design whose edges cannot fit that rate cannot be played.',
  ],
  audio: [
    'The YIN tuner searches 50 Hz to 1 kHz for a dominant periodic pitch. Silence, noise, and multiple simultaneous notes may not yield a reliable note.',
    'The spectrogram uses the browser analyser and actual audio sample rate. Microphone monitoring is muted to avoid speaker feedback.',
  ],
  comms: [
    'BER simulation uses ideal symbol decisions with additive Gaussian noise. Finite runs have uncertainty; zero observed errors does not prove zero BER.',
    'Shannon capacity is an ideal channel bound. CIDR /31 links use both addresses; VLSM uses the classic network/broadcast reservation.',
  ],
  magnetics: [
    'Lumped magnetic and transformer equivalent circuits are teaching approximations. Geometry, material data, and equivalent-circuit parameters determine their accuracy.',
    'Saturation is represented by the selected material model, rather than a measured dynamic hysteresis loop.',
  ],
  rotating: [
    'Induction-machine results use a steady-state equivalent circuit without a separate core-loss branch. DC-machine models ignore armature reaction.',
    'Speed-control and load intersections are steady-state operating points; mechanical and electrical startup transients are not simulated.',
  ],
  net: [
    'Access, delay, ARQ, and congestion-control examples use the displayed textbook assumptions; they are not packet captures or full network-stack simulations.',
    'Checksums, CRC, and parity detect specific error patterns. An accepted frame is not proof that every possible error was detected.',
  ],
  dspdesign: [
    'Transfer functions use normalized a0 = 1. An inverse z-transform depends on the selected region of convergence, not just its poles.',
    'FIR specifications are checked numerically on a frequency grid. Quantization and finite record length can change the observed response or SNR.',
  ],
  fields: [
    'Charge, filament, sheet, and wave models use the displayed ideal geometry and uniform material properties. Fields on ideal sources can be singular.',
    'Numerical quadrature has a finite tolerance and evaluation budget. Compare against closed form where available and check convergence near sources.',
  ],
  edm2: [
    'Magnetic circuits are lumped. The B–H table uses linear interpolation with extrapolation flagged outside the supplied data; this is not a dynamic hysteresis model.',
    'Energy, force, transformer-bank, and all-day efficiency calculations use the chosen equivalent circuits and load profile. Unmodeled construction losses require measured data.',
  ],
};
