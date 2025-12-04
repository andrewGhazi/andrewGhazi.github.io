# Run PGLMMs on a batch of tree files

This function fits phylogenetic generalized linear mixed models on a
batch of tree files, using the same outcome and covariate arguments.

## Usage

``` r
anpan_pglmm_batch(
  meta_file,
  tree_dir,
  outcome,
  covariates = NULL,
  offset = NULL,
  out_dir = NULL,
  trim_pattern = NULL,
  omit_na = FALSE,
  ladderize = TRUE,
  family = "gaussian",
  show_plot_cor_mat = TRUE,
  show_plot_tree = TRUE,
  save_object = FALSE,
  verbose = TRUE,
  loo_comparison = TRUE,
  run_diagnostics = FALSE,
  reg_noise = TRUE,
  plot_ext = "pdf",
  show_yrep = FALSE,
  show_post = TRUE,
  reg_gamma_params = c(1, 2),
  beta_sd = NULL,
  sigma_phylo_scale = 0.333,
  seed = 123,
  ...
)
```

## Arguments

- meta_file:

  either a data frame of metadata or a path to file containing the
  metadata

- tree_dir:

  string giving the path to a directory of tree files

- outcome:

  the name of the outcome variable

- covariates:

  covariates to account for (as a vector of strings)

- offset:

  a variable to include as an offset

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

- save_object:

  logical indicating whether to save the model fit object

- loo_comparison:

  logical indicating whether to compare the phylogenetic model against a
  base model (without the phylogenetic term) using
  [`loo::loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html)

- run_diagnostics:

  logical indicating whether to run
  [`cmdstanr::cmdstan_diagnose()`](https://mc-stan.org/cmdstanr/reference/fit-method-cmdstan_summary.html)
  and
  [`loo::pareto_k_table()`](https://mc-stan.org/loo/reference/pareto-k-diagnostic.html)
  to check the MCMC and loo diagnostics respectively.

- reg_noise:

  logical indicating whether to regularize the ratio of sigma_phylo to
  sigma_resid with a Gamma prior

- plot_ext:

  extension to use when saving plots

- show_yrep:

  show a plot of the tree overlaid with the outcome and the posterior
  predictive distribution for each observation if plotting the tree

- show_post:

  show a plot of the tree overlaid with the outcome and posterior
  distribution on phylogenetic effects.

- reg_gamma_params:

  the shape and rate parameters of the Gamma prior on the noise term
  ratio. Default: c(1,2)

- beta_sd:

  prior standard deviation parameters on the normal distribution for
  each covariate in the GLM component

- sigma_phylo_scale:

  standard deviation of half-normal prior on `sigma_phylo` for logistic
  PGLMMs when `family = 'binomial'`. Increasing this value can easily
  lead to overfitting.

- seed:

  random seed to pass to furrr_options()

- ...:

  other arguments to pass to
  [`cmdstanr::sample()`](https://mc-stan.org/cmdstanr/reference/model-method-sample.html)

## Value

a tibble listing results for each tree file in input directory that fit
successfully. Columns give the number of leaves on the tree, diagnostic
values, loo comparison values, formatted input data, correlation
matrices, PGLMM and "base" model fits, and loo objects (in list columns
where appropriate).

## Details

See
[`anpan_pglmm()`](https://andrewghazi.github.io/anpan/reference/anpan_pglmm.md)
for details on most of the arguments.

`tree_dir` must contain ONLY tree files readable by ape::read.tree()

If any trees cause an error while fitting, these are saved into a data
frame in a file `pglmm_errors.RData` in the output directory.

The Stan model fitting can't be parallelized via futures, so the most
effective way to parallelize the model fitting AND the importance weight
calculations is a nested future topology (e.g.
`plan(list(sequential, tweak(multisession, workers = 4)))` ) and set
parallel_chains = 4 . This will run sequentially over the trees, running
the model fits with 4 parallel chains for each tree, then compute the
importance weights in the future multisession for each tree.

The tibble result from this function contains a lot of large objects in
list columns, so it can be pretty big (several GBs) when saved to disk
in an RData file (and pretty ugly when not printed as a tibble). So be
careful if you try to save the whole thing.

## See also

[`ape::read.tree`](https://rdrr.io/pkg/ape/man/read.tree.html),
[`ape::write.tree`](https://rdrr.io/pkg/ape/man/write.tree.html),
[`anpan_pglmm()`](https://andrewghazi.github.io/anpan/reference/anpan_pglmm.md)
