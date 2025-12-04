# Use repeated measures to refine the gene model

Use repeated measures to refine the gene model

## Usage

``` r
anpan_repeated_measures(
  subject_sample_map,
  bug_dir,
  meta_file,
  out_dir,
  genomes_dir = NULL,
  model_type = "fastglm",
  covariates = c("age", "gender"),
  outcome = "crc",
  omit_na = FALSE,
  filtering_method = "kmeans",
  minmax_thresh = NULL,
  discard_poorly_covered_samples = TRUE,
  skip_large = TRUE,
  save_fit = TRUE,
  annotation_file = NULL,
  save_filter_stats = TRUE,
  verbose = TRUE,
  plot_result = TRUE,
  plot_ext = "pdf",
  q_threshold = 0.1,
  n_top = 50,
  width = 10,
  height = 8,
  ...
)
```

## Arguments

- subject_sample_map:

  a data frame between the mapping between subject_id and sample_id

- bug_dir:

  a directory of gene family files

- meta_file:

  path to a metadata tsv

- out_dir:

  path to the desired output directory

- genomes_dir:

  an optional directory of genome files

- model_type:

  either "horseshoe" or "fastglm"

- covariates:

  character vector of covariates to include in the model

- outcome:

  character string of the outcome variable

- omit_na:

  logical indicating whether to omit incomplete cases of the metadata

- filtering_method:

  method to use for filtering samples. Either "kmeans" or "none"

- minmax_thresh:

  genes must have at least this many (or N - this many) non-zero
  observations or else be discarded. NULL defaults to
  `floor(0.005*nrow(metadata))`.

- discard_poorly_covered_samples:

  logical indicating whether to discard samples where the genes of a bug
  are poorly covered

- skip_large:

  logical indicating whether to skip bugs with over 5k genes. Only used
  when model_type = "horseshoe".

- save_fit:

  logical indicating whether to save horseshoe fit objects. Only used
  when model_type = "horseshoe".

- annotation_file:

  a path to a file giving annotations for each gene

- save_filter_stats:

  logical indicating whether to save filter statistics

- plot_result:

  logical indicating whether or not to plot the results

- plot_ext:

  extension to use for plots

- q_threshold:

  FDR threshold to use for inclusion in the plot.

- n_top:

  number of top elements to show from the results

- width:

  width of saved plot in inches

- height:

  height of saved plot in inches

- ...:

  arguments to pass to \[cmdstanr::sample()\] if applicable

## Details

This function performs the standard anpan filtering on all samples, then
uses the subject-sample map to compute the proportion of samples with
the bug. This gives a gene \_proportion\_ matrix (instead of a
presence/absence matrix) which is then passed to
`anpan_batch(filtering_method = "none", discretize_inputs = FALSE)`.
Subjects that do not have the bug present in at least half their samples
are dropped.

In cases where subject metadata varies by sample, the mean is taken if
the variable is numeric, otherwise it is tabulated and the most frequent
category is selected as the subject-level metadata value. This
tabulation will respect factor ordering if you'd like to alter the value
selected in the event of ties.
