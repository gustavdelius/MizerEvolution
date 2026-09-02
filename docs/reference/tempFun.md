# Calculate temperature scalar

Calculates a physiological temperature scalar based on the Padfield
(2016) temperature performance equation.

## Usage

``` r
tempFun(w, temperature, t_ref, Ea, c_a)
```

## Arguments

- w:

  Vector of body weights.

- temperature:

  Ambient temperature in Celsius.

- t_ref:

  Reference temperature in Celsius.

- Ea:

  Activation energy.

- c_a:

  Temperature scaling parameter.

## Value

A numeric matrix of temperature scalars.
