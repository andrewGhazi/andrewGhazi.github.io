# Fit a subject-wise PGLMM

If you have a dataset where multiple samples come from the same
individual, running a PGLMM on all the samples without accounting for
this will return a strong spurious signal because samples from the same
individual will (usually) be right next to each other on the tree and
all show the same outcome value. This function uses a subject-sample map
input to aggregate samples by subject, derive a subject-wise correlation
matrix, then run a PGLMM on that.

## Usage

``` r
anpan_subjectwise_pglmm(
  tree_file,
  meta_file,
  subject_sample_map,
  outcome,
  offset = NULL,
  covariates = NULL,
  out_dir = NULL,
  trim_pattern = NULL,
  omit_na = FALSE,
  ladderize = TRUE,
  family = "gaussian",
  show_plot_cor_mat = TRUE,
  show_plot_tree = TRUE,
  show_post = TRUE,
  show_yrep = FALSE,
  save_object = FALSE,
  verbose = TRUE,
  loo_comparison = TRUE,
  reg_noise = TRUE,
  reg_gamma_params = c(1, 2),
  plot_ext = "pdf",
  beta_sd = NULL,
  sigma_phylo_scale = 0.333,
  ...
)
```

## Arguments

- tree_file:

  either a path to a tree file readable by
  [`ape::read.tree()`](https://rdrr.io/pkg/ape/man/read.tree.html) or an
  object of class "phylo" that is already read into R. Ignored if
  `cor_mat` is supplied.

- meta_file:

  either a data frame or a file. Can provide covariate and outcome
  variables for either subjects or samples.

- subject_sample_map:

  a data frame giving a mapping between sample_id (which must match the
  leaves of the tree up to the trim_pattern) and subject_id

- outcome:

  the name of the outcome variable

- offset:

  a variable to include as an offset

- covariates:

  covariates to account for (as a vector of strings)

- out_dir:

  if saving, directory where to save

- trim_pattern:

  optional pattern to trim from tip labels of the tree

- omit_na:

  logical indicating whether to omit incomplete cases of the metadata

- ladderize:

  logical indicating whether to run
  [`ape::ladderize()`](https://rdrr.io/pkg/ape/man/ladderize.html) on
  the tree before running the model

- family:

  string giving the name of the distribution of the outcome variable
  (usually "gaussian" or "binomial")

- show_plot_cor_mat:

  show a plot of the correlation matrix derived from the tree

- show_plot_tree:

  show a plot of the tree overlaid with the outcome.

- show_post:

  show a plot of the tree overlaid with the outcome and posterior
  distribution on phylogenetic effects.

- show_yrep:

  show a plot of the tree overlaid with the outcome and the posterior
  predictive distribution for each observation if plotting the tree

- save_object:

  logical indicating whether to save the model fit object

- loo_comparison:

  logical indicating whether to compare the phylogenetic model against a
  base model (without the phylogenetic term) using
  [`loo::loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html)

- reg_noise:

  logical indicating whether to regularize the ratio of sigma_phylo to
  sigma_resid with a Gamma prior

- reg_gamma_params:

  the shape and rate parameters of the Gamma prior on the noise term
  ratio. Default: c(1,2)

- plot_ext:

  extension to use when saving plots

- beta_sd:

  prior standard deviation parameters on the normal distribution for
  each covariate in the GLM component

- sigma_phylo_scale:

  standard deviation of half-normal prior on `sigma_phylo` for logistic
  PGLMMs when `family = 'binomial'`. Increasing this value can easily
  lead to overfitting.

- ...:

  other arguments to pass to
  [`cmdstanr::sample()`](https://mc-stan.org/cmdstanr/reference/model-method-sample.html)

## Details

The `meta_file` must contain at least one column named "sample_id" or
"subject_id". If the metadata is inferred to be provided by sample,
representative covariate and outcome variable values are selected for
each subject in the manner described in `?anpan_repeated_measures()`

## See also

[`anpan_pglmm()`](https://andrewghazi.github.io/anpan/reference/anpan_pglmm.md)
[`anpan_repeated_measures()`](https://andrewghazi.github.io/anpan/reference/anpan_repeated_measures.md)
