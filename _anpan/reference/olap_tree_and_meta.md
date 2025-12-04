# Overlap a tree and metadata

Overlap a tree's leaves with the observations in metadata, returning an
intersected tree and metadata data.table.

## Usage

``` r
olap_tree_and_meta(
  tree_file,
  meta_file,
  covariates,
  offset,
  outcome,
  omit_na = FALSE,
  ladderize = TRUE,
  trim_pattern = NULL,
  verbose = TRUE
)
```

## Arguments

- tree_file:

  either a path to a tree file readable by
  [`ape::read.tree()`](https://rdrr.io/pkg/ape/man/read.tree.html) or an
  object of class "phylo" that is already read into R. Ignored if
  `cor_mat` is supplied.

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

- trim_pattern:

  optional pattern to trim from tip labels of the tree

## Value

A list with two elements: the tree and metadata, both cut down to
samples that are present in the other.
