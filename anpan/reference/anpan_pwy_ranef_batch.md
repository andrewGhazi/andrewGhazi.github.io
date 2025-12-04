# Fit the pathway random effects model for multiple bugs

Fit the pathway random effects model for multiple bugs

## Usage

``` r
anpan_pwy_ranef_batch(
  bug_pwy_dat,
  group_ind = "crc",
  out_dir = NULL,
  group_exp_rate = 3,
  ...
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

- out_dir:

  output directory

- group_exp_rate:

  rate parameter of the exponential distribution prior on group effects

- ...:

  other arguments to pass to cmdstanr::sample()

## Value

a tibble of row-binded anpan_pwy_ranef results

## Details

In addition to the column requirements of
[`anpan_pwy_ranef()`](https://andrewghazi.github.io/anpan/reference/anpan_pwy_ranef.md),
the input data frame `bug_pwy_dat` here must also contain a variable
called "bug" which gives a unique identifier for each bug.

## See also

[`plot_pwy_ranef()`](https://andrewghazi.github.io/anpan/reference/plot_pwy_ranef.md)
,
[`plot_pwy_ranef_intervals()`](https://andrewghazi.github.io/anpan/reference/plot_pwy_ranef_intervals.md),
[`cmdstanr::CmdStanMCMC`](https://mc-stan.org/cmdstanr/reference/CmdStanMCMC.html)

## Examples

``` r
if (FALSE) { # \dontrun{
library(tidyverse)
library(anpan)

set.seed(123)
input_dat = tibble(bug = rep(paste0("bug", 1:5), each = 200),
                   pwy = rep(paste0('pwy', 1:5), times = 200),
                   log10_species_abd = rnorm(1000),
                   log10_pwy_abd = rnorm(1000, mean = .8*log10_species_abd),
                   group = sample(c(0,1), size = 1000, replace = TRUE))
 # ^ the pathway and and group are NOT related in any pathway.

res = anpan_pwy_ranef_batch(input_dat, group_ind = "group")

# Examine the summary
pwy_group_res = res |>
  dplyr::select(bug, summary_df) |>             # select the two main columns
  tidyr::unnest(c(summary_df)) |>               # unnest
  dplyr::filter(grepl("^pwy_eff", variable)) |> # get just the pwy:group terms
  dplyr::arrange(-abs(mean))                    # sort by decreasing effect size

print(pwy_group_res)

pwy_group_res |> filter(hit)
# ^ Here, there are no hits because we simulated with no dependence
} # }
```
