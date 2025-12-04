# Run anpan

Run the anpan gene model on a single bug

## Usage

``` r
anpan(
  bug_file,
  meta_file,
  out_dir,
  genomes_file = NULL,
  prefiltered_dir = NULL,
  model_type = "fastglm",
  covariates = c("age", "gender"),
  outcome = "crc",
  omit_na = FALSE,
  filtering_method = "kmeans",
  discretize_inputs = TRUE,
  minmax_thresh = NULL,
  skip_large = TRUE,
  save_fit = TRUE,
  discard_poorly_covered_samples = TRUE,
  plot_ext = "pdf",
  save_filter_stats = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- bug_file:

  path to a gene family file (usually from HUMAnN)

- meta_file:

  path to a metadata tsv

- out_dir:

  path to the desired output directory

- genomes_file:

  optional file giving gene presence/absence of representative isolate
  genomes

- prefiltered_dir:

  an optional directory to pre-filtered data from an earlier run to skip
  the filtering step

- model_type:

  either "horseshoe" or "fastglm"

- covariates:

  covariates to account for (as a vector of strings)

- outcome:

  the name of the outcome variable

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

- skip_large:

  logical indicating whether to skip bugs with over 5k genes. Only used
  when model_type = "horseshoe".

- save_fit:

  logical indicating whether to save horseshoe fit objects. Only used
  when model_type = "horseshoe".

- discard_poorly_covered_samples:

  logical indicating whether to discard samples where the genes of a bug
  are poorly covered

- plot_ext:

  extension to use for plots

- save_filter_stats:

  logical indicating whether to save filter statistics

- ...:

  arguments to pass to \[cmdstanr::sample()\] if applicable

## Value

a data.table of model statistics for each gene

## Details

The specified metadata file must contain columns matching "sample_id"
and the specified covariates and outcome variables.

If provided, `genomes_file` is used to refine the filtering process. The
format must be genes as rows, with the first column giving the gene id
(usually a UniRef90 identifier), and subsequent columns representing
isolate genomes. The entries of the isolate genome columns should give
0/1 indicators of whether or not the gene is present in the isolate. The
gene counts present in these isolates are used to establish the typical
number of genes present in a strain of the species and a lower threshold
on the number of acceptable gene observations. If \>=5 isolate genomes
are available, the lower threshold is 2 standard deviations below the
mean, otherwise it is 2/3 of the mean.

## See also

\[anpan_batch()\]
