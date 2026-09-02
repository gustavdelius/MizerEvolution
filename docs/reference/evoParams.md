# Set up parameters for an evolutionary model

Sets up a `MizerParams` object for an evolutionary or invasion model.

## Usage

``` r
evoParams(
  no_sp = 11,
  min_w_max = 10,
  max_w_max = 10^4,
  min_w = 10^(-3),
  max_w = max_w_max,
  eta = 10^(-0.6),
  min_w_mat = min_w_max * eta,
  no_w = log10(max_w_max/min_w) * 20 + 1,
  min_w_pp = 1e-10,
  w_pp_cutoff = min_w_mat,
  n = 2/3,
  p = n,
  lambda = 2.05,
  r_pp = 0.1,
  kappa = 0.005,
  alpha = 0.4,
  h = 40,
  beta = 100,
  sigma = 1.3,
  f0 = 0.6,
  fc = 0.25,
  ks = NA,
  gamma = NA,
  zeta = 0.2,
  ext_mort_prop = 0,
  reproduction_level = 1/4,
  R_factor = deprecated(),
  gear_names = "knife_edge_gear",
  knife_edge_size = 1000,
  RDD = "extinctionRDD",
  egg_size_scaling = FALSE,
  resource_scaling = FALSE,
  perfect_scaling = FALSE,
  updateParams = NULL,
  ea_int = 0,
  ca_int = 0,
  ed_int = 0,
  cd_int = 0,
  t_ref = 0,
  t_d = 0,
  min_w_inf = deprecated(),
  max_w_inf = deprecated(),
  info_level = default_info_level()
)
```

## Arguments

- no_sp:

  The number of species in the model.

- min_w_max:

  The maximum size of the smallest species in the community. This will
  be rounded to lie on a grid point.

- max_w_max:

  The maximum size of the largest species in the community. This will be
  rounded to lie on a grid point.

- min_w:

  The size of the the egg of the smallest species. This also defines the
  start of the community size spectrum.

- max_w:

  The largest size in the model. By default this is set to the largest
  maximum size `max_w_max`. Setting it to something larger only makes
  sense if you plan to add larger species to the model later.

- eta:

  Ratio between maturity size and maximum size of a species. Ignored if
  `min_w_mat` is supplied. Default is 10^(-0.6), approximately 1/4.

- min_w_mat:

  The maturity size of the smallest species. Default value is
  `eta * min_w_max`. This will be rounded to lie on a grid point.

- no_w:

  The number of size bins in the community spectrum. These bins will be
  equally spaced on a logarithmic scale. Default value is such that
  there are 50 bins for each factor of 10 in weight.

- min_w_pp:

  The smallest size of the resource spectrum. By default this is set to
  the smallest value at which any of the consumers can feed.

- w_pp_cutoff:

  The largest size of the resource spectrum. Default value is max_w_max
  unless `perfect_scaling = TRUE` when it is Inf.

- n:

  Scaling exponent of the maximum intake rate.

- p:

  Scaling exponent of the standard metabolic rate. By default this is
  equal to the exponent `n`.

- lambda:

  Exponent of the abundance power law.

- r_pp:

  Growth rate parameter for the resource spectrum.

- kappa:

  Coefficient in abundance power law.

- alpha:

  The assimilation efficiency of the community.

- h:

  Maximum food intake rate.

- beta:

  Preferred predator prey mass ratio.

- sigma:

  Width of prey size preference.

- f0:

  Expected average feeding level. Used to set `gamma`, the coefficient
  in the search rate. Ignored if `gamma` is given explicitly.

- fc:

  Critical feeding level. Used to determine `ks` if it is not given
  explicitly.

- ks:

  Standard metabolism coefficient. If not provided, default will be
  calculated from critical feeding level argument `fc`.

- gamma:

  Volumetric search rate. If not provided, default is determined by
  [`mizer::get_gamma_default()`](https://sizespectrum.org/mizer/reference/get_gamma_default.html)
  using the value of `f0`.

- zeta:

  Amplitude of lognormal distribution around trait when generating new
  phenotypes.

- ext_mort_prop:

  The proportion of the total mortality that comes from external
  mortality, i.e., from sources not explicitly modelled. A number in the
  interval \[0, 1).

- reproduction_level:

  Level of density dependence in reproduction.

- R_factor:

  **\[deprecated\]** Use `reproduction_level = 1 / R_factor` instead.

- gear_names:

  The names of the fishing gears for each species. A character vector,
  the same length as the number of species.

- knife_edge_size:

  The minimum size at which the gear or gears select fish. A single
  value for each gear or a vector with one value for each gear.

- RDD:

  Stock recruitment relationship function name. Default is
  `"extinctionRDD"`.

- egg_size_scaling:

  If TRUE, the egg size is a constant fraction of the maximum size of
  each species. This fraction is `min_w / min_w_max`. If FALSE, all
  species have the egg size `w_min`.

- resource_scaling:

  If TRUE, the carrying capacity for larger resource is reduced to
  compensate for the fact that fish eggs and larvae are present in the
  same size range.

- perfect_scaling:

  If TRUE then parameters are set so that the community abundance,
  growth before reproduction and death are perfect power laws. In
  particular all other scaling corrections are turned on.

- updateParams:

  Optional `MizerParams` object to convert into evolutionary parameters.

- ea_int:

  Temperature parameter Ea for intake.

- ca_int:

  Temperature parameter ca for intake.

- ed_int:

  Temperature parameter Ed for metabolism.

- cd_int:

  Temperature parameter cd for metabolism.

- t_ref:

  Reference temperature.

- t_d:

  Deactivation temperature.

- min_w_inf:

  **\[deprecated\]** Use `min_w_max` instead.

- max_w_inf:

  **\[deprecated\]** Use `max_w_max` instead.

- info_level:

  Controls the amount of information messages shown.

## Value

An object of type `MizerParams`

## Examples

``` r
if (FALSE) { # \dontrun{
params <- evoParams(no_sp = 5)
sim <- evoProject(params, t_max = 5)
plotDynamics(sim)
} # }
```
