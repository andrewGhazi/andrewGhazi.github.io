# Draw a sample from a multivariate normal distribution

Draw a sample from a multivariate normal distribution with specified
mean and covariance.

## Usage

``` r
draw_mvnorm(mean = NULL, Sigma = NULL, L = NULL)
```

## Arguments

- mean:

  mean vector, defaults to 0 vector

- Sigma:

  covariance matrix

- L:

  optional cholesky factor of `Sigma`

## Value

A vector giving one draw from the specified multivariate normal
distribution.

## Details

One of `Sigma` or `L` must be provided. Pre-computing L with
`L = t(chol(Sigma))` will make calling this function repeatedly much
faster.
