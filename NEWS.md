# mizerEvolution 0.1.0

* Upgraded package to be compatible with `mizer` 3.4.
* Switched to S3 `MizerParams` and `MizerSim` classes and modernized slot access.
* Registered extension via `recordExtension()`.
* Renamed `min_w_inf` and `max_w_inf` in `evoParams()` to `min_w_max` and `max_w_max`, and `R_factor` to `reproduction_level` (with lifecycle deprecation warnings for backwards compatibility).
* Updated rate function calls (`getMort()`, `getERepro()`).
* Fixed dimensions and data structures in `evoProject()`, `finalTouch()`, and all plotting functions.
* Added comprehensive unit test suite covering constructors, simulation stepping, and plotting functions.
