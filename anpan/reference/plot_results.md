# Plot the data for top results

This funciton makes a tile plot of the top results of a fit alongside
another tile plot showing the covariates included. Optional annotations
can be included.

## Usage

``` r
plot_results(
  res,
  covariates,
  outcome,
  model_input,
  bug_name = NULL,
  discretize_inputs = TRUE,
  plot_dir = NULL,
  annotation_file = NULL,
  cluster = "none",
  show_trees = FALSE,
  n_top = 50,
  q_threshold = 0.1,
  beta_threshold = 1,
  show_intervals = TRUE,
  stars = TRUE,
  max_anno_width = 70,
  width = NULL,
  height = NULL,
  plot_ext = "pdf"
)
```

## Arguments

- res:

  a data frame of model results (from `anpan`) for the genes of a single
  bug (i.e. the output written to \*gene_terms.tsv.gz)

- covariates:

  character string of the covariates to show

- outcome:

  character string of the outcome variable

- model_input:

  data frame of the model input

- bug_name:

  character string giving the name to use in the title/output file

- plot_dir:

  directory to write output to

- annotation_file:

  optional path file giving annotations

- cluster:

  axis to cluster. either "none", "samples", "genes", or "both"

- show_trees:

  logical to show the trees for the samples (if clustered)

- n_top:

  number of top elements to show from the results

- q_threshold:

  FDR threshold to use for inclusion in the plot.

- beta_threshold:

  Regression coefficient threshold to use for inclusion in the plot. Set
  to 0 to include everything.

- show_intervals:

  logical indicating whether to show the interval plot of estimates on
  the left

- stars:

  logical indicating whether to show significance stars on the

- width:

  width of saved plot in inches

- height:

  height of saved plot in inches

- plot_ext:

  extension to use for plots

## Details

If included, `annotation_file` must be a tsv with two columns: "gene"
and "annotation".

`n_top` is ignored if `q_threshold` is specified.

When `cluster = "none"`, the samples are ordered by metadata and the
genes are ordered by statistical significance.

When significance stars are shown, they encode the following (fairly
standard) significance thresholds: p.value \< .001 ~ \*\*\*, p.value \<
.01 ~ \*\*, p.value \< .05 ~ \*, p.value \< .1 ~ ., p.value \< 1 ~ " "

If applicable, the Q-value used to color the dot on the interval panel
is q_global if present in the input and q_bug_wise otherwise. That means
that you'll get different results if you compare the output of
[`anpan_batch()`](https://andrewghazi.github.io/anpan/reference/anpan_batch.md)
and a manual call to `plot_results()` using the bug-wise results from
the `model_stats/` output directory. If you'd like to replicate the
[`anpan_batch()`](https://andrewghazi.github.io/anpan/reference/anpan_batch.md)
plots exactly, read in the `all_bug_gene_terms.tsv.gz` result from the
top level output directory, then filter it to the bug of interest.

Note that the beta_threshold uses the value of the estimate column
directly, so it is interpreted according to the units of your outcome
variable with a continuous outcome, and on the log-odds scale with a
binary outcome. So the default value of 1 is pretty big for a binary
outcome, but if the spread of your continuous outcome variable is ~5 the
default value of 1 won't exclude very much.
