# Read and filter a gene family file

Read and filter a gene family file

## Usage

``` r
read_and_filter(
  bug_file,
  metadata,
  pivot_wide = TRUE,
  minmax_thresh,
  covariates = NULL,
  outcome = NULL,
  genomes_file = NULL,
  filtering_method = "kmeans",
  discretize_inputs = TRUE,
  discard_poorly_covered_samples = TRUE,
  save_filter_stats = TRUE,
  filter_stats_dir = NULL,
  plot_ext = "pdf",
  verbose = TRUE
)
```

## Arguments

- bug_file:

  path to a gene family file (usually from HUMAnN)

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

- genomes_file:

  optional file giving gene presence/absence of representative isolate
  genomes

- filtering_method:

  either "kmeans" or "none"

- discretize_inputs:

  logical indicating whether to discretize the input abundance
  measurements (0/nonzero –\> FALSE/TRUE) before passing them to the
  modelling function

- discard_poorly_covered_samples:

  logical indicating whether to discard samples where the genes of a bug
  are poorly covered

- save_filter_stats:

  logical indicating whether to save filter statistics

- filter_stats_dir:

  directory to save filtering statistics to

- plot_ext:

  extension to use for plots
