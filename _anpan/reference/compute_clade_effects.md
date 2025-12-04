# Compute a clade effect

This function computes the average phylogenetic effect of members of a
user-defined clade against that of all non-members.

## Usage

``` r
compute_clade_effects(
  clade_members,
  anpan_pglmm_result,
  plot_difference = TRUE
)
```

## Arguments

- clade_members:

  a character vector listing members of the clade of interest

- anpan_pglmm_result:

  a result list from
  [`anpan_pglmm()`](https://andrewghazi.github.io/anpan/reference/anpan_pglmm.md)

- plot_difference:

  logical indicating whether to plot 50 three computed variables

## Value

a list of two tibbles
