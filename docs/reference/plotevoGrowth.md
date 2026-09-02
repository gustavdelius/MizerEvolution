# Plot growth

Plots the instantaneous growth rates of species against body size.

## Usage

``` r
plotevoGrowth(
  object,
  time_range = max(as.numeric(dimnames(object@n)$time)),
  species = TRUE,
  print_it = TRUE,
  returnData = FALSE,
  save_it = FALSE,
  ylim = c(NA, NA),
  nameSave = "Growth.png",
  ...
)
```

## Arguments

- object:

  An object of class
  [mizer::MizerSim](https://sizespectrum.org/mizer/reference/MizerSim.html).

- time_range:

  The time range (either a vector of values, a vector of min and max
  time, or a single value) to average the abundances over. Default is
  the entire simulation time range.

- species:

  Plot only a specific species and its phenotypes. Default is NULL.

- print_it:

  A boolean value that displays the plot. Default is TRUE.

- returnData:

  A boolean value that controls the output. Either a dataframe ready for
  ggplot (TRUE) or the plot itself (FALSE). Default is FALSE.

- save_it:

  A boolean value to save the plot. Default is FALSE.

- ylim:

  Limits for the y axis.

- nameSave:

  File name for saving the plot. Default is `"Biomass.png"`.

- ...:

  Additional arguments passed to methods.

## Value

A ggplot object, or data frame if `returnData = TRUE`.

## See also

[mizer::plotting_functions](https://sizespectrum.org/mizer/reference/plotting_functions.html)

Other plotting functions:
[`plotDynamics()`](https://sizespectrum.org/mizerEvolution/reference/plotDynamics.md),
[`plotPredMort()`](https://sizespectrum.org/mizerEvolution/reference/plotPredMort.md),
[`plotPredRate()`](https://sizespectrum.org/mizerEvolution/reference/plotPredRate.md),
[`plotSS()`](https://sizespectrum.org/mizerEvolution/reference/plotSS.md),
[`plotSpawn()`](https://sizespectrum.org/mizerEvolution/reference/plotSpawn.md),
[`plotStarvation()`](https://sizespectrum.org/mizerEvolution/reference/plotStarvation.md),
[`plotevoFeeding()`](https://sizespectrum.org/mizerEvolution/reference/plotevoFeeding.md),
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
plotevoGrowth(sim)
#> Warning: log-10 transformation introduced infinite values.

# }
```
