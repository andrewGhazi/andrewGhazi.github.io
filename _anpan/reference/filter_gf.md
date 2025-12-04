# Filter a gene family file

This is the function that anpan uses to filter an input genefamily
data.table. This can be useful if you want to play around with different
filtering options without repeatedly re-reading or checking the file.

## Usage

``` r
filter_gf(
  gf,
  samp_stats,
  filtering_method = "kmeans",
  covariates = NULL,
  outcome = NULL,
  genomes_file = NULL,
  discard_poorly_covered_samples = TRUE,
  save_filter_stats = FALSE,
  filter_stats_dir = NULL,
  plot_ext = "pdf",
  bug_name = NULL
)
```

## Arguments

- gf:

  a gene family data.table

- samp_stats:

  a data.table of sample statistics

- filtering_method:

  method to use for filtering samples. Either "kmeans" or "none"

- covariates:

  covariates to account for (as a vector of strings)

- outcome:

  the name of the outcome variable

- genomes_file:

  optional file giving gene presence/absence of representative isolate
  genomes

- discard_poorly_covered_samples:

  logical indicating whether to discard samples where the genes of a bug
  are poorly covered

- save_filter_stats:

  logical indicating whether to save filter statistics

- filter_stats_dir:

  directory to save filtering statistics to

- plot_ext:

  extension to use for plots
