# Changelog

## mizerEvolution 0.1.0

- Upgraded package to be compatible with `mizer` 3.4.
- Switched to S3 `MizerParams` and `MizerSim` classes and modernized
  slot access.
- Registered extension via
  [`recordExtension()`](https://sizespectrum.org/mizer/reference/recordExtension.html).
- Renamed `min_w_inf` and `max_w_inf` in
  [`evoParams()`](https://sizespectrum.org/mizerEvolution/reference/evoParams.md)
  to `min_w_max` and `max_w_max`, and `R_factor` to `reproduction_level`
  (with lifecycle deprecation warnings for backwards compatibility).
- Updated rate function calls
  ([`getMort()`](https://sizespectrum.org/mizer/reference/getMort.html),
  [`getERepro()`](https://sizespectrum.org/mizer/reference/getERepro.html)).
- Fixed dimensions and data structures in
  [`evoProject()`](https://sizespectrum.org/mizerEvolution/reference/evoProject.md),
  `finalTouch()`, and all plotting functions.
- Added comprehensive unit test suite covering constructors, simulation
  stepping, and plotting functions.
