# Plot a tree and the PGLMM posterior on phylogenetic effects

Plot a tree and the PGLMM posterior on phylogenetic effects

## Usage

``` r
plot_tree_with_post(
  tree_file,
  meta_file,
  fit,
  labels,
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

- fit:

  a pglmm fit from
  [`anpan_pglmm()`](https://andrewghazi.github.io/anpan/reference/anpan_pglmm.md)

- labels:

  the ordered tip labels from the tree

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

## Value

either the plot or (if return_tree_df = TRUE) a list containing the
plot, the segment df, the terminal segment df, and the yrep df.

## Details

The whiskers of each box plot are the 90% posterior intervals, the box
is the 50% interval, and the middle line is the posterior mean.

The `labels` needs to contain the leaf labels *in the order produced by
the Cholesky factorization of the correlation matrix* (which is how the
data are passed to the sampler). This is not necessarily the order of
the leaves on the x-axis of the tree, nor the order of the samples in
the metadata. The simplest way to get the sample_ids in the proper order
is to take the `sample_id` column from the model_input result from
anpan_pglmm().

## Examples

``` r
if (FALSE) { # \dontrun{
# Using the result simulated in the vignette:
plot_tree_with_post(tr, metadata,
                    fit        = result$pglmm_fit,
                    labels     = result$model_input$sample_id,
                    outcome    = 'outcome',
                    covariate  = 'covariate',
                    color_bars = TRUE)
} # }
```
