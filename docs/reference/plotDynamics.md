# Plot biomass variations through time of species and phenotypes

Plots the total biomass of each species and optionally its individual
phenotypes over the course of the simulation.

## Usage

``` r
plotDynamics(
  object,
  phenotype = TRUE,
  species = NULL,
  trait = NULL,
  SpIdx = NULL,
  time_range = c(min(as.numeric(dimnames(object@n)$time)),
    max(as.numeric(dimnames(object@n)$time))),
  print_it = TRUE,
  returnData = FALSE,
  save_it = FALSE,
  nameSave = "Biomass.png",
  ylimit = c(NA, NA)
)
```

## Arguments

- object:

  An object of class
  [mizer::MizerSim](https://sizespectrum.org/mizer/reference/MizerSim.html).

- phenotype:

  A boolean value to display the phenotype's abundance on the plot.
  Default is TRUE.

- species:

  Plot only a specific species and its phenotypes. Default is NULL.

- trait:

  Display the trait value as color gradient. Works only if `species` is
  not NULL.

- SpIdx:

  Display a subset of species.

- time_range:

  The time range (either a vector of values, a vector of min and max
  time, or a single value) to average the abundances over. Default is
  the entire simulation time range.

- print_it:

  A boolean value that displays the plot. Default is TRUE.

- returnData:

  A boolean value that controls the output. Either a dataframe ready for
  ggplot (TRUE) or the plot itself (FALSE). Default is FALSE.

- save_it:

  A boolean value to save the plot. Default is FALSE.

- nameSave:

  File name for saving the plot. Default is `"Biomass.png"`.

- ylimit:

  To edit the y axis scale on the plot.

## Value

A ggplot object, or a list of data frames if `returnData = TRUE`.

## See also

[mizer::plotting_functions](https://sizespectrum.org/mizer/reference/plotting_functions.html)

Other plotting functions:
[`plotPredMort()`](https://sizespectrum.org/mizerEvolution/reference/plotPredMort.md),
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
plotDynamics(sim)

# }
```
