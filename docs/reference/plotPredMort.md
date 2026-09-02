# Plot predation mortality rate of each species against size

After running a projection, plot the predation mortality rate of each
species by size.

## Usage

``` r
plotPredMort(object, species = NULL, time_range, highlight = NULL)
```

## Arguments

- object:

  An object of class
  [mizer::MizerSim](https://sizespectrum.org/mizer/reference/MizerSim.html).

- species:

  Plot only a specific species and its phenotypes. Default is NULL.

- time_range:

  The time range (either a vector of values, a vector of min and max
  time, or a single value) to average the abundances over. Default is
  the entire simulation time range.

- highlight:

  Vector of species names to highlight.

## Value

A ggplot object.

## See also

[mizer::plotting_functions](https://sizespectrum.org/mizer/reference/plotting_functions.html),
[`mizer::getPredMort()`](https://sizespectrum.org/mizer/reference/getPredMort.html)

Other plotting functions:
[`plotDynamics()`](https://sizespectrum.org/mizerEvolution/reference/plotDynamics.md),
[`plotPredRate()`](https://sizespectrum.org/mizerEvolution/reference/plotPredRate.md),
[`plotSS()`](https://sizespectrum.org/mizerEvolution/reference/plotSS.md),
[`plotSpawn()`](https://sizespectrum.org/mizerEvolution/reference/plotSpawn.md),
[`plotStarvation()`](https://sizespectrum.org/mizerEvolution/reference/plotStarvation.md),
[`plotevoFeeding()`](https://sizespectrum.org/mizerEvolution/reference/plotevoFeeding.md),
[`plotevoGrowth()`](https://sizespectrum.org/mizerEvolution/reference/plotevoGrowth.md),
[`plotevoMortality()`](https://sizespectrum.org/mizerEvolution/reference/plotevoMortality.md),
[`plotevoTrait()`](https://sizespectrum.org/mizerEvolution/reference/plotevoTrait.md)

## Examples

``` r
# \donttest{
params <- evoParams(no_sp = 5)
#> ℹ No `a` column so using a = 0.01 in w = a l^b, with w in g and l in cm.
#> ℹ No `b` column so using the isometric default b = 3 in w = a l^b.
#> ℹ Using f0, h, lambda, kappa and the predation kernel to calculate gamma.
#> Warning: Some column names in your species parameter data frame are very close to standard parameter names: `zeta` (did you mean `beta`?). Please check for mis-spellings.
sim <- evoProject(params, t_max = 20)
plotPredMort(sim)

# }
```
