# Filter a batch of files

This function applies anpan::read_and_filter() to a set of files.

## Usage

``` r
filter_batch(
  bug_dir,
  meta_file,
  filter_stats_dir,
  pivot_wide = TRUE,
  minmax_thresh = NULL,
  covariates = NULL,
  outcome = NULL,
  filtering_method = "kmeans",
  discretize_inputs = TRUE,
  discard_poorly_covered_samples = TRUE,
  omit_na = FALSE,
  plot_ext = "pdf",
  verbose = TRUE
)
```

## Arguments

- bug_dir:

  a directory of gene family files

- meta_file:

  path to a metadata tsv

- filter_stats_dir:

  directory to save filtering statistics to

- pivot_wide:

  logical indicating whether to return data in wide format

- minmax_thresh:

  genes must have at least this many (or N - this many) non-zero
  observations or else be discarded. NULL defaults to
  `floor(0.005*nrow(metadata))`.

- covariates:

  covariates to account for (as a vector of strings)

- outcome:

  the name of the outcome variable

- filtering_method:

  either "kmeans" or "none"

- discretize_inputs:

  logical indicating whether to discretize the input abundance
  measurements (0/nonzero –\> FALSE/TRUE) before passing them to the
  modelling function

- discard_poorly_covered_samples:

  logical indicating whether to discard samples where the genes of a bug
  are poorly covered

- omit_na:

  logical indicating whether to omit incomplete cases of the metadata

- plot_ext:

  extension to use for plots

## Value

A list of filtered data frames
