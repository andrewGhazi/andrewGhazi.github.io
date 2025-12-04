# Plot a tree and the PGLMM posterior predictive

Plot a tree and the PGLMM posterior predictive

## Usage

``` r
plot_tree_with_post_pred(
  tree_file,
  meta_file,
  covariates = c("age", "gender"),
  offset = NULL,
  outcome = "crc",
  omit_na = FALSE,
  ladderize = TRUE,
  color_bars = FALSE,
  verbose = TRUE,
  fit,
  labels,
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

- fit:

  a pglmm fit from
  [`anpan_pglmm()`](https://andrewghazi.github.io/anpan/reference/anpan_pglmm.md)

- labels:

  the ordered tip labels from the tree

- trim_pattern:

  optional pattern to trim from tip labels of the tree

- return_tree_df:

  if true, return a list containing 1) the plot, 2) the segment data
  frame, and 3) the labelled terminal segment data frame. Otherwise,
  just return the plot.

## Value

either the plot or (if return_tree_df = TRUE) a list containing the
plot, the segment df, the terminal segment df, and the yrep df.

## Details

The whiskers of each box plot are the 90 box is the 50 case of binary
outcomes, the dot for each leaf represents the mean of the posterior
predictions (which is a proportion).
