# Plot a pathway random effects result

Plot a pathway random effects result

## Usage

``` r
plot_pwy_ranef(
  bug_pwy_dat,
  pwy_ranef_res,
  group_ind = "crc",
  group_labels = c("ctrl", "case"),
  bug_name = NULL,
  max_pwy = 20,
  post_draws = 30,
  verbose = TRUE
)
```

## Arguments

- bug_pwy_dat:

  a data frame with a row for each observation and columns "pwy",
  "log10_species_abd", "log10_pwy_abd", and a group indicator column
  named according to the `group_ind` argument

- group_ind:

  a character giving the name of the column for the 0/1 group indicator
  variable in `bug_pwy_dat`

- group_labels:

  labels for the 0/1 indicator to use on the plots

- bug_name:

  name of the bug (if using a result from
  [`anpan_pwy_ranef()`](https://andrewghazi.github.io/anpan/reference/anpan_pwy_ranef.md))

- max_pwy:

  the maximum number of bug:pwy facets to include

- post_draws:

  number of post draws to draw in each facet

- verbose:

  logical for verbosity

## Details

This function plots bug:pwy data alongside posterior draws. If there's a
strong group effect in a particular bug:pwy combination, you will see
wide separation of the red and blue posterior lines.

If `bug_name` is specified, `bug_pwy_dat` is first filtered to just data
from that bug.

If specified, `bug_name` must exactly match the corresponding entries in
`bug_pwy_dat`.
