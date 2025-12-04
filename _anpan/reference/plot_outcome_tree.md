# Plot a tree file showing the outcome variable

Plot a tree file, and show the outcome variable as a colored dot on the
end of each tip.

## Usage

``` r
plot_outcome_tree(
  tree_file,
  meta_file,
  covariates = c("age", "gender"),
  offset = NULL,
  outcome = "crc",
  omit_na = FALSE,
  ladderize = TRUE,
  color_bars = FALSE,
  verbose = TRUE,
  trim_pattern = NULL,
  return_tree_df = FALSE
)
```

## Arguments

- tree_file:

  either a path to a tree file readable by ape::read.tree() or an object
  of class "phylo" that is already read into R.

- meta_file:

  either a data frame of metadata or a path to file containing the
  metadata

- covariates:

  covariates to account for (as a vector of strings)

- offset:

  a variable to include as an offset

- outcome:

  the name of the outcome variable

- omit_na:

  logical indicating whether to omit incomplete cases of the metadata

- ladderize:

  logical indicating whether to run
  [`ape::ladderize()`](https://rdrr.io/pkg/ape/man/ladderize.html) on
  the tree before running the model

- color_bars:

  if true, show color bars below the plot showing the covariates and
  outcome variables.

- trim_pattern:

  optional pattern to trim from tip labels of the tree

- return_tree_df:

  if true, return a list containing 1) the plot, 2) the segment data
  frame, and 3) the labelled terminal segment data frame. Otherwise,
  just return the plot.

## Details

Showing the covariates as color bar annotations isn't supported yet.
