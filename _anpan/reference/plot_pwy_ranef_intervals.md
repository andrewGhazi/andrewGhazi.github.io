# Plot a pathway random effects result

Plot a pathway random effects result

## Usage

``` r
plot_pwy_ranef_intervals(pwy_ranef_res, group_ind = "crc", ncol = 1)
```

## Arguments

- pwy_ranef_res:

  a result from
  [`anpan_pwy_ranef()`](https://andrewghazi.github.io/anpan/reference/anpan_pwy_ranef.md)
  or
  [`anpan_pwy_ranef_batch()`](https://andrewghazi.github.io/anpan/reference/anpan_pwy_ranef_batch.md)

- group_ind:

  a character giving the name of the column for the 0/1 group indicator
  variable in `bug_pwy_dat`

- ncol:

  The number of columns to use if faceting multiple bugs.
