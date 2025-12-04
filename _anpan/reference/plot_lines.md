# Make a filter-labelled line plot

Make a filter-labelled line plot

## Usage

``` r
plot_lines(
  bug_file = NULL,
  meta_file = NULL,
  covariates,
  outcome,
  genomes_stats,
  omit_na = FALSE,
  fgf = NULL,
  bug_name = NULL,
  plot_ext = "pdf",
  plot_dir = NULL,
  subset_line = 200
)
```

## Arguments

- bug_file:

  a path to a gene family file

- meta_file:

  a path to the corresponding metadata file

- covariates:

  covariates to account for (as a vector of strings)

- outcome:

  the name of the outcome variable

- omit_na:

  logical indicating whether to omit incomplete cases of the metadata

- fgf:

  a filtered gene family data frame

- bug_name:

  name of the bug

- plot_ext:

  extension to use for plots

- subset_line:

  integer gives the number of points to take along the lines

## Details

The required input is either

- the gene family file and the metadata file

- OR a pre-filtered gene family file

`subset_line` is used to make saving plots faster. Plotting thousands of
lines each with tens of thousands of points along them is too much
visual detail and makes saving the plot very slow. Set `subset_line` to
0 to turn off subsetting.
