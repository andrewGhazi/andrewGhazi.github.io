# Apply anpan to a many bugs

This function calls anpan() on each gene family file in the `bug_dir`
directory and makes a composite data + results plot for each.

## Usage

``` r
anpan_batch(
  bug_dir,
  meta_file,
  out_dir,
  genomes_dir = NULL,
  prefiltered_dir = NULL,
  model_type = "fastglm",
  covariates = c("age", "gender"),
  outcome = "crc",
  omit_na = FALSE,
  filtering_method = "kmeans",
  discretize_inputs = TRUE,
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
  beta_threshold = 1,
  n_top = 50,
  width = 10,
  height = 8,
  ...
)
```

## Arguments

- bug_dir:

  a directory of gene family files

- meta_file:

  path to a metadata tsv

- out_dir:

  path to the desired output directory

- genomes_dir:

  an optional directory of genome files

- prefiltered_dir:

  an optional directory to pre-filtered data from an earlier run to skip
  the filtering step

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

- discretize_inputs:

  logical indicating whether to discretize the input abundance
  measurements (0/nonzero –\> FALSE/TRUE) before passing them to the
  modelling function

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

- beta_threshold:

  Regression coefficient threshold to use for inclusion in the plot. Set
  to 0 to include everything.

- n_top:

  number of top elements to show from the results

- width:

  width of saved plot in inches

- height:

  height of saved plot in inches

- ...:

  arguments to pass to \[cmdstanr::sample()\] if applicable

## Value

a data.table of model statistics for each bug:gene combination

## Details

`bug_dir` should be a directory of gene (or SNV or pathway) abundance
files, one for each bug.

`annotation` file must have two columns named "gene" and "annotation"

See `?anpan()` for the format / usage if providing genome files. If
provided, genomes_dir must contain ONLY the genome files themselves.

## See also

\[anpan()\]
