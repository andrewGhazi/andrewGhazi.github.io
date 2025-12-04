# Run a phylogenetic generalized linear mixed model

Run a phylogenetic generalized linear mixed model

## Usage

``` r
anpan_pglmm(
  meta_file,
  tree_file = NULL,
  cor_mat = NULL,
  outcome,
  covariates = NULL,
  offset = NULL,
  out_dir = NULL,
  trim_pattern = NULL,
  bug_name = NULL,
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
  run_diagnostics = TRUE,
  reg_noise = TRUE,
  reg_gamma_params = c(1, 2),
  plot_ext = "pdf",
  beta_sd = NULL,
  int_prior_scale = 1,
  sigma_phylo_scale = 0.333,
  ...
)
```

## Arguments

- meta_file:

  either a data frame of metadata or a path to file containing the
  metadata

- tree_file:

  either a path to a tree file readable by
  [`ape::read.tree()`](https://rdrr.io/pkg/ape/man/read.tree.html) or an
  object of class "phylo" that is already read into R. Ignored if
  `cor_mat` is supplied.

- cor_mat:

  a correlation matrix provided as an alternative to a tree.

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

- run_diagnostics:

  logical indicating whether to run
  [`cmdstanr::cmdstan_diagnose()`](https://mc-stan.org/cmdstanr/reference/fit-method-cmdstan_summary.html)
  and
  [`loo::pareto_k_table()`](https://mc-stan.org/loo/reference/pareto-k-diagnostic.html)
  to check the MCMC and loo diagnostics respectively.

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

- int_prior_scale:

  standard deviation of the 0-centered normal prior for the intercept of
  the model (with centered covariates)

- sigma_phylo_scale:

  standard deviation of half-normal prior on `sigma_phylo` for logistic
  PGLMMs when `family = 'binomial'`. Increasing this value can easily
  lead to overfitting.

- ...:

  other arguments to pass to
  [`cmdstanr::sample()`](https://mc-stan.org/cmdstanr/reference/model-method-sample.html)

## Value

A list containing the model input (in the order passed to the model),
estimated correlation matrix, the pglmm fit object, and (if
`loo_comparison` is on) the base fit object and the associated loo
objects.

## Details

the tip labels of the tree must be the sample ids from the metadata. You
can use the `trim_pattern` argument to automatically trim off any
consistent pattern from the tip labels if necessary.

The default error distribution for the outcome is "gaussian". You could
change this to a phylogenetic logistic regression by changing `family`
to "binomial" for example.

The prior for the intercept is a normal distribution centered on the
mean of the outcome variable with a standard deviation of 3\*sd(outcome
variable).

If a prior scale on covariate effects isn't specified (i.e. `beta_sd` is
left at the default NULL), the prior scale is set to 1 / (1 standard
deviation) for each covariate. The values will be shown in a message.
It's preferable to set `beta_sd` with domain knowledge. Consider the
units/scale of your outcome variable if using `family = "gaussian"` or
the range of plausible coefficient values on the log-odds scale if using
`family = "binomial"`.

If supplying `beta_sd` with a categorical predictor with \>2 levels,
only specify a single corresponding element in beta_sd. This appropriate
element will get repeated as necessary.

`sigma_phylo_scale` is used to set the width of the half-normal prior on
sigma_phylo when `family = "binomial"` (no effect when
`family = "gaussian"`). Increasing this value can easily lead to
overfitting. If sigma_phylo is high, then leafwise phylogenetic effects
can get very far from zero (disregarding the correlation structure and
covariate effects), which lets the model perfectly fit each observation.
The default value keeps most of the prior mass on leafwise phylogenetic
effects between +/- 1, which is still pretty liberal on the log-odds
scale.

The model indexes the leaves according to the Cholesky factor of the
correlation matrix. This likely won't be the same order as the input
tree. The `model_input` that is returned as part of the output will have
the samples in a factor with the model index order if needed.

It's normal to see some warnings during warmup, particularly about
"Scale vector is inf".

This function tries to get the `bug_name` argument from the tree file,
but if it's not provided you may need to set it yourself.

Common cmdstanr options that one might want to pass in via ... include
`refresh = 500` (show fewer MCMC progress updates), `adapt_delta = .98`
(avoid divergences at the cost of possibly needing more iterations to
get convergence), `show_messages = FALSE` , and `parallel_chains = 4`
(run the MCMC chains in parallel).

If you want to use the PGLMM log-likelihood data frame with functions
from the `loo` package, you'll need to convert it to a matrix with
[`as.matrix()`](https://rdrr.io/r/base/matrix.html). It's converted to a
tibble internally so that it prints nicely.

If `int_prior_scale` isn't specified, it defaults to 1 for binary
outcomes and 1 standard deviation of the outcome for gaussian outcomes.

If `offset` is categorical, the offset value used is the empirical mean
(continuous outcomes) or proportion (binary outcome) for each category.
Using a categorical offset variable with a continuous outcome is strange
and you probably shouldn't try it unless you know what you're doing.

## See also

[`anpan_pglmm_batch()`](https://andrewghazi.github.io/anpan/reference/anpan_pglmm_batch.md),
[`loo::loo()`](https://mc-stan.org/loo/reference/loo.html),
[`cmdstanr::sample()`](https://mc-stan.org/cmdstanr/reference/model-method-sample.html)

## Examples

``` r
meta = data.frame(x = rnorm(100), sample_id = paste0("t", 1:100))
tr = ape::rtree(100)
anpan_pglmm(meta, tr,
outcome = "x",
iter_sampling = 10, # Use more in practice!
iter_warmup = 10,
show_plot_cor_mat = FALSE,
show_plot_tree = FALSE)
#> (1/4) Checking inputs.
#> (2/4) Fitting model(s).
#> Init values were only set for a subset of parameters. 
#> Missing init values for the following parameters:
#>  - chain 1: beta, centered_cov_intercept, std_phylo_effects
#>  - chain 2: beta, centered_cov_intercept, std_phylo_effects
#>  - chain 3: beta, centered_cov_intercept, std_phylo_effects
#>  - chain 4: beta, centered_cov_intercept, std_phylo_effects
#> 
#> To disable this message use options(cmdstanr_warn_inits = FALSE).
#> Running MCMC with 4 sequential chains...
#> 
#> Chain 1 WARNING: No variance estimation is 
#> Chain 1          performed for num_warmup < 20 
#> Chain 1 Iteration:  1 / 20 [  5%]  (Warmup) 
#> Chain 1 Iteration: 11 / 20 [ 55%]  (Sampling) 
#> Chain 1 Iteration: 20 / 20 [100%]  (Sampling) 
#> Chain 1 finished in 0.0 seconds.
#> Chain 2 WARNING: No variance estimation is 
#> Chain 2          performed for num_warmup < 20 
#> Chain 2 Iteration:  1 / 20 [  5%]  (Warmup) 
#> Chain 2 Iteration: 11 / 20 [ 55%]  (Sampling) 
#> Chain 2 Iteration: 20 / 20 [100%]  (Sampling) 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: normal_id_glm_lpdf: Scale vector is inf, but must be positive finite! (in '/tmp/Rtmp0Zrdfj/model-9ab81ee95ea.stan', line 38, column 2 to column 62)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 finished in 0.0 seconds.
#> Chain 3 WARNING: No variance estimation is 
#> Chain 3          performed for num_warmup < 20 
#> Chain 3 Iteration:  1 / 20 [  5%]  (Warmup) 
#> Chain 3 Iteration: 11 / 20 [ 55%]  (Sampling) 
#> Chain 3 Iteration: 20 / 20 [100%]  (Sampling) 
#> Chain 3 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 3 Exception: normal_id_glm_lpdf: Scale vector is inf, but must be positive finite! (in '/tmp/Rtmp0Zrdfj/model-9ab81ee95ea.stan', line 38, column 2 to column 62)
#> Chain 3 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 3 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 3 
#> Chain 3 finished in 0.0 seconds.
#> Chain 4 WARNING: No variance estimation is 
#> Chain 4          performed for num_warmup < 20 
#> Chain 4 Iteration:  1 / 20 [  5%]  (Warmup) 
#> Chain 4 Iteration: 11 / 20 [ 55%]  (Sampling) 
#> Chain 4 Iteration: 20 / 20 [100%]  (Sampling) 
#> Chain 4 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 4 Exception: normal_id_glm_lpdf: Scale vector is inf, but must be positive finite! (in '/tmp/Rtmp0Zrdfj/model-9ab81ee95ea.stan', line 38, column 2 to column 62)
#> Chain 4 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 4 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 4 
#> Chain 4 finished in 0.0 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 0.0 seconds.
#> Total execution time: 0.6 seconds.
#> 
#> Init values were only set for a subset of parameters. 
#> Missing init values for the following parameters:
#>  - chain 1: beta, centered_cov_intercept
#>  - chain 2: beta, centered_cov_intercept
#>  - chain 3: beta, centered_cov_intercept
#>  - chain 4: beta, centered_cov_intercept
#> 
#> To disable this message use options(cmdstanr_warn_inits = FALSE).
#> Running MCMC with 4 sequential chains...
#> 
#> Chain 1 WARNING: No variance estimation is 
#> Chain 1          performed for num_warmup < 20 
#> Chain 1 Iteration:  1 / 20 [  5%]  (Warmup) 
#> Chain 1 Iteration: 11 / 20 [ 55%]  (Sampling) 
#> Chain 1 Iteration: 20 / 20 [100%]  (Sampling) 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: normal_id_glm_lpdf: Scale vector is inf, but must be positive finite! (in '/tmp/Rtmp0Zrdfj/model-9ab866d133c.stan', line 28, column 2 to column 95)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 finished in 0.0 seconds.
#> Chain 2 WARNING: No variance estimation is 
#> Chain 2          performed for num_warmup < 20 
#> Chain 2 Iteration:  1 / 20 [  5%]  (Warmup) 
#> Chain 2 Iteration: 11 / 20 [ 55%]  (Sampling) 
#> Chain 2 Iteration: 20 / 20 [100%]  (Sampling) 
#> Chain 2 finished in 0.0 seconds.
#> Chain 3 WARNING: No variance estimation is 
#> Chain 3          performed for num_warmup < 20 
#> Chain 3 Iteration:  1 / 20 [  5%]  (Warmup) 
#> Chain 3 Iteration: 11 / 20 [ 55%]  (Sampling) 
#> Chain 3 Iteration: 20 / 20 [100%]  (Sampling) 
#> Chain 3 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 3 Exception: normal_id_glm_lpdf: Matrix of independent variables is inf, but must be finite! (in '/tmp/Rtmp0Zrdfj/model-9ab866d133c.stan', line 28, column 2 to column 95)
#> Chain 3 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 3 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 3 
#> Chain 3 finished in 0.0 seconds.
#> Chain 4 WARNING: No variance estimation is 
#> Chain 4          performed for num_warmup < 20 
#> Chain 4 Iteration:  1 / 20 [  5%]  (Warmup) 
#> Chain 4 Iteration: 11 / 20 [ 55%]  (Sampling) 
#> Chain 4 Iteration: 20 / 20 [100%]  (Sampling) 
#> Chain 4 finished in 0.0 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 0.0 seconds.
#> Total execution time: 0.5 seconds.
#> 
#> Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
#> ℹ Please use tidy evaluation idioms with `aes()`.
#> ℹ See also `vignette("ggplot2-in-packages")` for more information.
#> ℹ The deprecated feature was likely used in the anpan package.
#>   Please report the issue to the authors.

#> (3/4) Evaluating loo comparison.
#> - 0/2 preparing loo inputs
#> - 1/2 precomputing conditional covariance arrays
#> - 2/2 computing integrated importance weights for loo CV
#> Warning: Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.
#> Warning: Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.
#> loo comparison: 
#>           elpd_diff se_diff
#> base_fit   0.0       0.0   
#> pglmm_fit -0.2       0.5   
#> The phylogenetic model seems to fit worse, but the difference doesn't seem clear (less than 2 standard errors difference in ELPD).
#> (4/4) Running diagnostics:
#> Checking sampler transitions treedepth.
#> Treedepth satisfactory for all transitions.
#> 
#> Checking sampler transitions for divergences.
#> No divergent transitions found.
#> 
#> Checking E-BFMI - sampler transitions HMC potential energy.
#> E-BFMI satisfactory.
#> 
#> Rank-normalized split effective sample size satisfactory for all parameters.
#> 
#> The following parameters had rank-normalized split R-hat greater than 1.01:
#>   centered_cov_intercept, sigma_resid, sigma_phylo, std_phylo_effects[1], std_phylo_effects[2], std_phylo_effects[3], std_phylo_effects[4], std_phylo_effects[6], std_phylo_effects[8], std_phylo_effects[9], std_phylo_effects[10], std_phylo_effects[11], std_phylo_effects[12], std_phylo_effects[14], std_phylo_effects[15], std_phylo_effects[18], std_phylo_effects[19], std_phylo_effects[20], std_phylo_effects[22], std_phylo_effects[23], std_phylo_effects[25], std_phylo_effects[26], std_phylo_effects[27], std_phylo_effects[28], std_phylo_effects[29], std_phylo_effects[31], std_phylo_effects[32], std_phylo_effects[33], std_phylo_effects[34], std_phylo_effects[35], std_phylo_effects[36], std_phylo_effects[37], std_phylo_effects[38], std_phylo_effects[39], std_phylo_effects[40], std_phylo_effects[41], std_phylo_effects[42], std_phylo_effects[43], std_phylo_effects[44], std_phylo_effects[46], std_phylo_effects[47], std_phylo_effects[48], std_phylo_effects[49], std_phylo_effects[50], std_phylo_effects[52], std_phylo_effects[53], std_phylo_effects[54], std_phylo_effects[55], std_phylo_effects[56], std_phylo_effects[57], std_phylo_effects[58], std_phylo_effects[60], std_phylo_effects[61], std_phylo_effects[62], std_phylo_effects[63], std_phylo_effects[64], std_phylo_effects[65], std_phylo_effects[66], std_phylo_effects[67], std_phylo_effects[68], std_phylo_effects[69], std_phylo_effects[70], std_phylo_effects[71], std_phylo_effects[72], std_phylo_effects[74], std_phylo_effects[75], std_phylo_effects[76], std_phylo_effects[77], std_phylo_effects[78], std_phylo_effects[79], std_phylo_effects[80], std_phylo_effects[81], std_phylo_effects[82], std_phylo_effects[83], std_phylo_effects[84], std_phylo_effects[86], std_phylo_effects[88], std_phylo_effects[89], std_phylo_effects[91], std_phylo_effects[92], std_phylo_effects[94], std_phylo_effects[96], std_phylo_effects[97], std_phylo_effects[98], std_phylo_effects[99], std_phylo_effects[100], phylo_effect[1], phylo_effect[2], phylo_effect[3], phylo_effect[4], phylo_effect[5], phylo_effect[6], phylo_effect[7], phylo_effect[8], phylo_effect[9], phylo_effect[10], phylo_effect[11], phylo_effect[12], phylo_effect[15], phylo_effect[16], phylo_effect[17], phylo_effect[18], phylo_effect[19], phylo_effect[20], phylo_effect[21], phylo_effect[22], phylo_effect[24], phylo_effect[25], phylo_effect[26], phylo_effect[28], phylo_effect[29], phylo_effect[30], phylo_effect[31], phylo_effect[32], phylo_effect[33], phylo_effect[34], phylo_effect[35], phylo_effect[36], phylo_effect[37], phylo_effect[38], phylo_effect[39], phylo_effect[40], phylo_effect[41], phylo_effect[44], phylo_effect[46], phylo_effect[47], phylo_effect[48], phylo_effect[49], phylo_effect[50], phylo_effect[51], phylo_effect[53], phylo_effect[54], phylo_effect[55], phylo_effect[56], phylo_effect[57], phylo_effect[58], phylo_effect[59], phylo_effect[60], phylo_effect[61], phylo_effect[62], phylo_effect[63], phylo_effect[65], phylo_effect[66], phylo_effect[67], phylo_effect[68], phylo_effect[70], phylo_effect[71], phylo_effect[72], phylo_effect[74], phylo_effect[75], phylo_effect[76], phylo_effect[77], phylo_effect[78], phylo_effect[79], phylo_effect[80], phylo_effect[81], phylo_effect[82], phylo_effect[83], phylo_effect[84], phylo_effect[85], phylo_effect[86], phylo_effect[87], phylo_effect[88], phylo_effect[89], phylo_effect[91], phylo_effect[92], phylo_effect[93], phylo_effect[94], phylo_effect[95], phylo_effect[96], phylo_effect[97], phylo_effect[98], phylo_effect[99], phylo_effect[100], intercept, yrep[1], yrep[3], yrep[4], yrep[5], yrep[6], yrep[7], yrep[10], yrep[11], yrep[14], yrep[15], yrep[17], yrep[19], yrep[22], yrep[25], yrep[26], yrep[27], yrep[30], yrep[31], yrep[33], yrep[34], yrep[35], yrep[38], yrep[39], yrep[40], yrep[43], yrep[46], yrep[47], yrep[51], yrep[54], yrep[55], yrep[56], yrep[57], yrep[58], yrep[59], yrep[60], yrep[61], yrep[63], yrep[64], yrep[65], yrep[66], yrep[67], yrep[69], yrep[70], yrep[72], yrep[74], yrep[76], yrep[77], yrep[78], yrep[79], yrep[80], yrep[81], yrep[82], yrep[83], yrep[84], yrep[86], yrep[89], yrep[91], yrep[94], yrep[95], yrep[96], yrep[97], yrep[100], lin_pred[1], lin_pred[2], lin_pred[3], lin_pred[4], lin_pred[5], lin_pred[6], lin_pred[7], lin_pred[8], lin_pred[9], lin_pred[10], lin_pred[11], lin_pred[12], lin_pred[13], lin_pred[14], lin_pred[15], lin_pred[16], lin_pred[17], lin_pred[18], lin_pred[19], lin_pred[20], lin_pred[21], lin_pred[22], lin_pred[23], lin_pred[24], lin_pred[25], lin_pred[26], lin_pred[27], lin_pred[28], lin_pred[29], lin_pred[30], lin_pred[31], lin_pred[32], lin_pred[33], lin_pred[34], lin_pred[35], lin_pred[36], lin_pred[37], lin_pred[38], lin_pred[39], lin_pred[40], lin_pred[41], lin_pred[42], lin_pred[43], lin_pred[44], lin_pred[45], lin_pred[46], lin_pred[47], lin_pred[48], lin_pred[49], lin_pred[50], lin_pred[51], lin_pred[52], lin_pred[53], lin_pred[54], lin_pred[55], lin_pred[56], lin_pred[57], lin_pred[58], lin_pred[59], lin_pred[60], lin_pred[61], lin_pred[62], lin_pred[63], lin_pred[64], lin_pred[65], lin_pred[66], lin_pred[67], lin_pred[68], lin_pred[69], lin_pred[70], lin_pred[71], lin_pred[72], lin_pred[73], lin_pred[76], lin_pred[77], lin_pred[78], lin_pred[79], lin_pred[80], lin_pred[81], lin_pred[82], lin_pred[83], lin_pred[87], lin_pred[91], lin_pred[93], lin_pred[94], lin_pred[97], lin_pred[98], lin_pred[99]
#> Such high values indicate incomplete mixing and biased estimation.
#> You should consider regularizating your model with additional prior information or a more effective parameterization.
#> 
#> Processing complete.
#> Pareto k diagnostic values:
#>                           Count Pct.    Min. ESS
#> (-Inf, 0.38]   (good)     67    67.0%   26      
#>    (0.38, 1]   (bad)      33    33.0%   <NA>    
#>     (1, Inf)   (very bad)  0     0.0%   <NA>    
#> $model_input
#>     sample_id            x offset_val
#> 1         t34  0.112038083          0
#> 2         t42 -0.049964899          0
#> 3         t50  0.862086482          0
#> 4         t37 -0.279237242          0
#> 5         t19 -0.914074827          0
#> 6         t84 -0.548257264          0
#> 7         t94  0.486148920          0
#> 8         t78  0.172181715          0
#> 9          t6 -1.821817661          0
#> 10        t51 -0.243236740          0
#> 11        t88  0.433889790          0
#> 12        t40  0.070034850          0
#> 13        t86 -2.612334333          0
#> 14        t72  1.337320413          0
#> 15        t14  0.512426950          0
#> 16        t67  0.131670635          0
#> 17        t18  0.542996343          0
#> 18        t15 -1.863011492          0
#> 19        t32  0.243685465          0
#> 20        t24  1.888504929          0
#> 21        t39  1.067307879          0
#> 22        t53  0.019177592          0
#> 23        t87 -0.155693776          0
#> 24        t26 -0.935847354          0
#> 25        t97  0.946347886          0
#> 26        t16 -0.522012515          0
#> 27        t91  1.063101996          0
#> 28        t33  1.623548883          0
#> 29        t29 -1.512399651          0
#> 30        t25 -0.097445104          0
#> 31       t100 -0.387213575          0
#> 32        t74  1.318293384          0
#> 33        t31  0.176488611          0
#> 34        t77 -0.109935672          0
#> 35         t3 -0.005571287          0
#> 36        t17 -0.052601910          0
#> 37        t63 -0.245896412          0
#> 38        t73  0.236696283          0
#> 39        t41 -0.639123324          0
#> 40        t48  0.118194874          0
#> 41        t99 -0.296640025          0
#> 42        t60  1.074345882          0
#> 43        t43 -0.251483443          0
#> 44        t21  0.362951256          0
#> 45        t95  1.672882611          0
#> 46        t59  0.213355750          0
#> 47        t27 -0.015950311          0
#> 48        t81  1.298392759          0
#> 49        t70 -1.470736306          0
#> 50        t20  0.468154420          0
#> 51         t5  1.148411606          0
#> 52        t12  2.065024895          0
#> 53        t44  0.444797116          0
#> 54        t23  0.737776321          0
#> 55        t45  2.755417575          0
#> 56        t69 -1.699450568          0
#> 57         t4  0.621552721          0
#> 58        t96 -0.354361164          0
#> 59        t57  2.682557184          0
#> 60        t22 -1.304543545          0
#> 61        t55  0.549827542          0
#> 62        t58 -0.361221255          0
#> 63        t90  0.424187575          0
#> 64        t28 -0.826788954          0
#> 65        t46  0.046531380          0
#> 66        t71  0.284150344          0
#> 67        t89 -0.381951112          0
#> 68        t92  1.048712620          0
#> 69        t68  0.488628809          0
#> 70         t7 -0.247325302          0
#> 71        t82  0.748791268          0
#> 72        t85  1.110534893          0
#> 73        t64 -1.177563309          0
#> 74        t98  1.316826356          0
#> 75        t49 -1.911720491          0
#> 76         t8 -0.244199607          0
#> 77        t11  0.628982042          0
#> 78        t75  0.523909788          0
#> 79        t79 -0.090327287          0
#> 80         t2 -2.437263611          0
#> 81        t10 -0.553699384          0
#> 82        t62  1.113952419          0
#> 83        t35 -0.133997013          0
#> 84        t65 -0.975850616          0
#> 85        t54  0.029560754          0
#> 86        t83  0.556224329          0
#> 87        t47  0.577709069          0
#> 88        t93 -0.038102895          0
#> 89        t13 -1.630989402          0
#> 90        t76  0.606748047          0
#> 91         t9 -0.282705449          0
#> 92         t1  0.255317055          0
#> 93        t36 -1.910087468          0
#> 94        t56 -2.274114857          0
#> 95        t61 -0.665088249          0
#> 96        t52 -0.206087195          0
#> 97        t80  1.924343341          0
#> 98        t66  1.065057320          0
#> 99        t38 -0.313445978          0
#> 100       t30  0.935363190          0
#> 
#> $cor_mat
#>             t34        t42        t50        t37        t19        t84
#> t34  1.00000000 0.90418569 0.77039910 0.76131233 0.84584644 0.60113262
#> t42  0.90418569 1.00000000 0.73632806 0.72764314 0.80843872 0.57454742
#> t50  0.77039910 0.73632806 1.00000000 0.87214806 0.78008975 0.55440016
#> t37  0.76131233 0.72764314 0.87214806 1.00000000 0.77088867 0.54786106
#> t19  0.84584644 0.80843872 0.78008975 0.77088867 1.00000000 0.69611524
#> t84  0.60113262 0.57454742 0.55440016 0.54786106 0.69611524 1.00000000
#> t94  0.60798619 0.58109788 0.56072092 0.55410728 0.70405171 0.82569687
#> t78  0.46433270 0.44379750 0.42823515 0.42318416 0.53770010 0.47881914
#> t6   0.41329136 0.39501347 0.38116181 0.37666604 0.47859392 0.42618540
#> t51  0.37454446 0.35798016 0.34542711 0.34135284 0.43372477 0.38622965
#> t88  0.36967694 0.35332791 0.34093800 0.33691668 0.42808816 0.38121028
#> t40  0.35203137 0.33646272 0.32466421 0.32083483 0.40765448 0.36301419
#> t86  0.34306477 0.32789266 0.31639468 0.31266283 0.39727110 0.35376785
#> t72  0.34141449 0.32631537 0.31487269 0.31115880 0.39536006 0.35206608
#> t14  0.48131206 0.46002594 0.44389452 0.43865883 0.55736230 0.49632823
#> t67  0.48859375 0.46698560 0.45061013 0.44529523 0.56579454 0.50383710
#> t18  0.48437380 0.46295227 0.44671824 0.44144925 0.56090781 0.49948549
#> t15  0.41159908 0.39339604 0.37960109 0.37512373 0.47663425 0.42444032
#> t32  0.35751669 0.34170544 0.32972309 0.32583404 0.41400651 0.36867064
#> t24  0.34671393 0.33138044 0.31976015 0.31598861 0.40149685 0.35753086
#> t39  0.33194507 0.31726474 0.30613943 0.30252855 0.38439442 0.34230123
#> t53  0.31982534 0.30568100 0.29496189 0.29148285 0.37035969 0.32980338
#> t87  0.33993709 0.32490331 0.31351015 0.30981233 0.39364923 0.35054259
#> t26  0.34599087 0.33068936 0.31909330 0.31532963 0.40065954 0.35678524
#> t97  0.34677049 0.33143450 0.31981232 0.31604016 0.40156235 0.35758918
#> t16  0.38127035 0.36440859 0.35163013 0.34748268 0.44151339 0.39316538
#> t91  0.32507113 0.31069480 0.29979988 0.29626377 0.37643435 0.33521283
#> t33  0.34621695 0.33090544 0.31930180 0.31553567 0.40092134 0.35701837
#> t29  0.36832043 0.35203139 0.33968695 0.33568038 0.42651731 0.37981145
#> t25  0.33821437 0.32325678 0.31192136 0.30824228 0.39165431 0.34876613
#> t100 0.43466697 0.41544374 0.40087566 0.39614737 0.50334700 0.44822789
#> t74  0.40195433 0.38417783 0.37070613 0.36633369 0.46546557 0.41449467
#> t31  0.36203728 0.34602612 0.33389225 0.32995403 0.41924138 0.37333227
#> t77  0.35833347 0.34248611 0.33047638 0.32657845 0.41495235 0.36951291
#> t3   0.36130964 0.34533065 0.33322118 0.32929087 0.41839876 0.37258192
#> t17  0.36958286 0.35323799 0.34085123 0.33683093 0.42797920 0.38111326
#> t63  0.42835555 0.40941144 0.39505489 0.39039526 0.49603834 0.44171956
#> t73  0.36796106 0.35168791 0.33935552 0.33535285 0.42610115 0.37944086
#> t41  0.40797292 0.38993024 0.37625683 0.37181892 0.47243513 0.42070103
#> t48  0.40332539 0.38548825 0.37197059 0.36758324 0.46705326 0.41590850
#> t99  0.35204059 0.33647153 0.32467272 0.32084323 0.40766516 0.36302370
#> t60  0.37839640 0.36166175 0.34897960 0.34486343 0.43818534 0.39020177
#> t43  0.42216934 0.40349882 0.38934960 0.38475727 0.48887467 0.43534035
#> t21  0.40953801 0.39142611 0.37770024 0.37324530 0.47424751 0.42231494
#> t95  0.42435422 0.40558708 0.39136463 0.38674853 0.49140478 0.43759340
#> t59  0.40317596 0.38534543 0.37183278 0.36744706 0.46688022 0.41575441
#> t27  0.39328180 0.37588883 0.36270780 0.35842970 0.45542272 0.40555156
#> t81  0.37005034 0.35368480 0.34128237 0.33725698 0.42852055 0.38159533
#> t70  0.37705931 0.36038379 0.34774646 0.34364482 0.43663698 0.38882296
#> t20  0.38731981 0.37019052 0.35720930 0.35299605 0.44851870 0.39940357
#> t5   0.38970273 0.37246805 0.35940697 0.35516780 0.45127814 0.40186084
#> t12  0.39217615 0.37483208 0.36168810 0.35742203 0.45414237 0.40441142
#> t44  0.46321383 0.44272811 0.42720327 0.42216445 0.53640444 0.47766537
#> t23  0.38293860 0.36600307 0.35316869 0.34900310 0.44344524 0.39488568
#> t45  0.33435422 0.31956734 0.30836129 0.30472420 0.38718423 0.34478554
#> t69  0.35749315 0.34168295 0.32970139 0.32581259 0.41397925 0.36864637
#> t4   0.21653808 0.20696164 0.19970426 0.19734877 0.25075241 0.22329372
#> t96  0.12399654 0.11851277 0.11435696 0.11300813 0.14358874 0.12786503
#> t57  0.12456515 0.11905623 0.11488137 0.11352635 0.14424720 0.12845138
#> t22  0.12086087 0.11551577 0.11146506 0.11015034 0.13995762 0.12463154
#> t55  0.10030913 0.09587294 0.09251103 0.09141987 0.11615858 0.10343861
#> t58  0.09801896 0.09368405 0.09039889 0.08933265 0.11350655 0.10107699
#> t90  0.09002105 0.08603985 0.08302275 0.08204350 0.10424492 0.09282956
#> t28  0.08025197 0.07670281 0.07401313 0.07314015 0.09293227 0.08275571
#> t46  0.07908423 0.07558671 0.07293616 0.07207589 0.09158002 0.08155153
#> t71  0.08437713 0.08064553 0.07781759 0.07689974 0.09770923 0.08700956
#> t89  0.09006528 0.08608213 0.08306354 0.08208382 0.10429614 0.09287517
#> t92  0.09920298 0.09481571 0.09149087 0.09041174 0.11487765 0.10229795
#> t68  0.10494563 0.10030439 0.09678708 0.09564549 0.12152768 0.10821977
#> t7   0.10790871 0.10313643 0.09951981 0.09834599 0.12495894 0.11127529
#> t82  0.10057496 0.09612702 0.09275620 0.09166215 0.11646642 0.10371274
#> t85  0.10914165 0.10431484 0.10065690 0.09946967 0.12638670 0.11254670
#> t64  0.09797249 0.09363964 0.09035604 0.08929030 0.11345275 0.10102908
#> t98  0.10248363 0.09795127 0.09451648 0.09340167 0.11867667 0.10568096
#> t49  0.09690479 0.09261916 0.08937135 0.08831722 0.11221634 0.09992807
#> t8   0.09440787 0.09023266 0.08706853 0.08604157 0.10932488 0.09735324
#> t11  0.09779061 0.09346580 0.09018830 0.08912454 0.11324212 0.10084152
#> t75  0.10124533 0.09676774 0.09337445 0.09227311 0.11724271 0.10440403
#> t79  0.17063874 0.16309221 0.15737316 0.15551696 0.19760070 0.17596240
#> t2   0.12754555 0.12190483 0.11763007 0.11624264 0.14769853 0.13152477
#> t10  0.13741481 0.13133761 0.12673209 0.12523729 0.15912719 0.14170193
#> t62  0.13316802 0.12727864 0.12281545 0.12136685 0.15420938 0.13732265
#> t35  0.13772981 0.13163868 0.12702260 0.12552438 0.15949196 0.14202676
#> t65  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t54  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t83  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t47  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t93  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t13  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t76  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t9   0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t1   0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t36  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t56  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t61  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t52  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t80  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t66  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t38  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t30  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#>             t94        t78         t6        t51        t88        t40
#> t34  0.60798619 0.46433270 0.41329136 0.37454446 0.36967694 0.35203137
#> t42  0.58109788 0.44379750 0.39501347 0.35798016 0.35332791 0.33646272
#> t50  0.56072092 0.42823515 0.38116181 0.34542711 0.34093800 0.32466421
#> t37  0.55410728 0.42318416 0.37666604 0.34135284 0.33691668 0.32083483
#> t19  0.70405171 0.53770010 0.47859392 0.43372477 0.42808816 0.40765448
#> t84  0.82569687 0.47881914 0.42618540 0.38622965 0.38121028 0.36301419
#> t94  1.00000000 0.48427820 0.43104437 0.39063309 0.38555649 0.36715295
#> t78  0.48427820 1.00000000 0.86354495 0.78258587 0.77241552 0.73554626
#> t6   0.43104437 0.86354495 1.00000000 0.80860204 0.79809359 0.75999865
#> t51  0.39063309 0.78258587 0.80860204 1.00000000 0.95611851 0.79036298
#> t88  0.38555649 0.77241552 0.79809359 0.95611851 1.00000000 0.78009156
#> t40  0.36715295 0.73554626 0.75999865 0.79036298 0.78009156 1.00000000
#> t86  0.35780118 0.71681114 0.74064070 0.77023162 0.76022182 0.87491235
#> t72  0.35608002 0.71336299 0.73707793 0.76652650 0.75656485 0.87070367
#> t14  0.50198691 0.85837711 0.76402081 0.69239231 0.68339410 0.65077405
#> t67  0.50958139 0.70480233 0.62732760 0.56851436 0.56112605 0.53434215
#> t18  0.50518017 0.69871500 0.62190941 0.56360414 0.55627964 0.52972707
#> t15  0.42927940 0.49526799 0.44082612 0.39949778 0.39430598 0.37548480
#> t32  0.37287389 0.43019185 0.38290341 0.34700544 0.34249582 0.32614767
#> t24  0.36160710 0.41719314 0.37133357 0.33652030 0.33214694 0.31629276
#> t39  0.34620384 0.39942210 0.35551600 0.32218565 0.31799858 0.30281974
#> t53  0.33356350 0.38483869 0.34253566 0.31042225 0.30638805 0.29176341
#> t87  0.35453916 0.40903872 0.36407552 0.32994270 0.32565482 0.31011053
#> t26  0.36085298 0.41632309 0.37055917 0.33581850 0.33145426 0.31563315
#> t97  0.36166609 0.41726120 0.37139415 0.33657520 0.33220112 0.31634436
#> t16  0.39764789 0.45877410 0.40834378 0.37006073 0.36525148 0.34781715
#> t91  0.33903462 0.39115083 0.34815395 0.31551381 0.31141345 0.29654893
#> t33  0.36108877 0.41659513 0.37080130 0.33603793 0.33167084 0.31583939
#> t29  0.38414171 0.44319175 0.39447432 0.35749155 0.35284565 0.33600348
#> t25  0.35274244 0.40696581 0.36223047 0.32827063 0.32400448 0.30853897
#> t100 0.45333817 0.52302506 0.46553202 0.42188745 0.41640468 0.39652868
#> t74  0.41922036 0.48366268 0.43049651 0.39013659 0.38506644 0.36668629
#> t31  0.37758867 0.43563138 0.38774501 0.35139313 0.34682649 0.33027162
#> t77  0.37372576 0.43117467 0.38377820 0.34779822 0.34327829 0.32689279
#> t3   0.37682976 0.43475582 0.38696569 0.35068688 0.34612941 0.32960782
#> t17  0.38545836 0.44471080 0.39582638 0.35871686 0.35405504 0.33715514
#> t63  0.44675564 0.51543067 0.45877244 0.41576159 0.41035843 0.39077103
#> t73  0.38376690 0.44275933 0.39408943 0.35714275 0.35250138 0.33567564
#> t41  0.42549748 0.49090471 0.43694247 0.39597823 0.39083217 0.37217680
#> t48  0.42065031 0.48531244 0.43196492 0.39146734 0.38637989 0.36793705
#> t99  0.36716257 0.42360259 0.37703847 0.34169035 0.33724980 0.32115206
#> t60  0.39465049 0.45531594 0.40526577 0.36727128 0.36249829 0.34519537
#> t43  0.44030371 0.50798694 0.45214695 0.40975726 0.40443213 0.38512761
#> t21  0.42712979 0.49278794 0.43861869 0.39749730 0.39233149 0.37360456
#> t95  0.44258244 0.51061596 0.45448698 0.41187791 0.40652521 0.38712078
#> t59  0.42049447 0.48513264 0.43180489 0.39132230 0.38623675 0.36780073
#> t27  0.41017529 0.47322721 0.42120815 0.38171903 0.37675828 0.35877470
#> t81  0.38594593 0.44527332 0.39632706 0.35917060 0.35450288 0.33758160
#> t70  0.39325596 0.45370705 0.40383372 0.36597350 0.36121737 0.34397559
#> t20  0.40395721 0.46605328 0.41482281 0.37593233 0.37104678 0.35333582
#> t5   0.40644249 0.46892060 0.41737494 0.37824519 0.37332959 0.35550966
#> t12  0.40902215 0.47189680 0.42002399 0.38064589 0.37569908 0.35776606
#> t44  0.48311127 0.55737486 0.49610595 0.44959502 0.44375216 0.42257080
#> t23  0.39938780 0.46078147 0.41013050 0.37167994 0.36684965 0.34933903
#> t45  0.34871647 0.40232097 0.35809621 0.32452397 0.32030651 0.30501751
#> t69  0.37284934 0.43016353 0.38287821 0.34698260 0.34247327 0.32612620
#> t4   0.22583951 0.26055543 0.23191412 0.21017170 0.20744035 0.19753872
#> t96  0.12932283 0.14920226 0.13280134 0.12035095 0.11878689 0.11311690
#> t57  0.12991587 0.14988646 0.13341033 0.12090284 0.11933161 0.11363563
#> t22  0.12605247 0.14542918 0.12944302 0.11730747 0.11578297 0.11025636
#> t55  0.10461793 0.12069974 0.10743193 0.09735997 0.09609470 0.09150786
#> t58  0.10222938 0.11794402 0.10497914 0.09513713 0.09390074 0.08941863
#> t90  0.09388792 0.10832032 0.09641331 0.08737436 0.08623886 0.08212247
#> t28  0.08369921 0.09656541 0.08595055 0.07789251 0.07688023 0.07321055
#> t46  0.08248131 0.09516028 0.08469988 0.07675909 0.07576155 0.07214527
#> t71  0.08800156 0.10152911 0.09036862 0.08189638 0.08083207 0.07697376
#> t89  0.09393405 0.10837354 0.09646068 0.08741730 0.08628124 0.08216282
#> t92  0.10346426 0.11936873 0.10624723 0.09628634 0.09503502 0.09049877
#> t68  0.10945359 0.12627873 0.11239766 0.10186015 0.10053640 0.09573755
#> t7   0.11254395 0.12984414 0.11557115 0.10473612 0.10337498 0.09844065
#> t82  0.10489518 0.12101961 0.10771664 0.09761799 0.09634936 0.09175037
#> t85  0.11382985 0.13132771 0.11689164 0.10593281 0.10455613 0.09956541
#> t64  0.10218092 0.11788812 0.10492937 0.09509203 0.09385623 0.08937625
#> t98  0.10688583 0.12331626 0.10976084 0.09947054 0.09817784 0.09349157
#> t49  0.10106736 0.11660338 0.10378586 0.09405572 0.09283339 0.08840223
#> t8   0.09846317 0.11359888 0.10111163 0.09163221 0.09044137 0.08612439
#> t11  0.10199122 0.11766926 0.10473457 0.09491549 0.09368199 0.08921032
#> t75  0.10559434 0.12182625 0.10843461 0.09826865 0.09699156 0.09236192
#> t79  0.17796856 0.20532578 0.18275554 0.16562184 0.16346945 0.15566665
#> t2   0.13302430 0.15347271 0.13660237 0.12379562 0.12218680 0.11635452
#> t10  0.14331749 0.16534817 0.14717243 0.13337472 0.13164140 0.12535784
#> t62  0.13888828 0.16023811 0.14262409 0.12925279 0.12757304 0.12148367
#> t35  0.14364602 0.16572720 0.14750979 0.13368045 0.13194316 0.12564520
#> t65  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t54  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t83  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t47  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t93  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t13  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t76  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t9   0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t1   0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t36  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t56  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t61  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t52  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t80  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t66  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t38  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t30  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#>             t86        t72        t14       t67       t18        t15        t32
#> t34  0.34306477 0.34141449 0.48131206 0.4885937 0.4843738 0.41159908 0.35751669
#> t42  0.32789266 0.32631537 0.46002594 0.4669856 0.4629523 0.39339604 0.34170544
#> t50  0.31639468 0.31487269 0.44389452 0.4506101 0.4467182 0.37960109 0.32972309
#> t37  0.31266283 0.31115880 0.43865883 0.4452952 0.4414492 0.37512373 0.32583404
#> t19  0.39727110 0.39536006 0.55736230 0.5657945 0.5609078 0.47663425 0.41400651
#> t84  0.35376785 0.35206608 0.49632823 0.5038371 0.4994855 0.42444032 0.36867064
#> t94  0.35780118 0.35608002 0.50198691 0.5095814 0.5051802 0.42927940 0.37287389
#> t78  0.71681114 0.71336299 0.85837711 0.7048023 0.6987150 0.49526799 0.43019185
#> t6   0.74064070 0.73707793 0.76402081 0.6273276 0.6219094 0.44082612 0.38290341
#> t51  0.77023162 0.76652650 0.69239231 0.5685144 0.5636041 0.39949778 0.34700544
#> t88  0.76022182 0.75656485 0.68339410 0.5611260 0.5562796 0.39430598 0.34249582
#> t40  0.87491235 0.87070367 0.65077405 0.5343421 0.5297271 0.37548480 0.32614767
#> t86  1.00000000 0.95255618 0.63419817 0.5207319 0.5162344 0.36592081 0.31784035
#> t72  0.95255618 1.00000000 0.63114742 0.5182270 0.5137511 0.36416059 0.31631141
#> t14  0.63419817 0.63114742 1.00000000 0.7305750 0.7242651 0.51337857 0.44592277
#> t67  0.52073190 0.51822698 0.73057500 1.0000000 0.9877540 0.52114539 0.45266907
#> t18  0.51623438 0.51375108 0.72426507 0.9877540 1.0000000 0.51664429 0.44875939
#> t15  0.36592081 0.36416059 0.51337857 0.5211454 0.5166443 1.00000000 0.69259180
#> t32  0.31784035 0.31631141 0.44592277 0.4526691 0.4487594 0.69259180 1.00000000
#> t24  0.30823646 0.30675372 0.43244873 0.4389912 0.4351996 0.67166438 0.82366227
#> t39  0.29510661 0.29368703 0.41402785 0.4202916 0.4166616 0.64305371 0.78857699
#> t53  0.28433190 0.28296415 0.39891117 0.4049462 0.4014487 0.61957501 0.75978505
#> t87  0.30221170 0.30075794 0.42399613 0.4304107 0.4266933 0.65853609 0.80756304
#> t26  0.30759364 0.30611400 0.43154688 0.4380757 0.4342920 0.67026365 0.81046129
#> t97  0.30828674 0.30680376 0.43251928 0.4390628 0.4352706 0.67177396 0.81228750
#> t16  0.33895789 0.33732736 0.47555019 0.4827447 0.4785753 0.73860808 0.80700434
#> t91  0.28899552 0.28760534 0.40545413 0.4115882 0.4080333 0.62973731 0.68805197
#> t33  0.30779463 0.30631402 0.43182885 0.4383619 0.4345758 0.67070161 0.73280963
#> t29  0.32744512 0.32586998 0.45939805 0.4663482 0.4623204 0.71352113 0.77959430
#> t25  0.30068016 0.29923377 0.42184742 0.4282295 0.4245309 0.65519879 0.71587122
#> t100 0.38642868 0.38456981 0.54215064 0.5503527 0.5455994 0.81545036 0.70830360
#> t74  0.35734642 0.35562744 0.50134888 0.5089337 0.5045381 0.75408032 0.65499732
#> t31  0.32185926 0.32031099 0.45156121 0.4583928 0.4544337 0.60952476 0.52943576
#> t77  0.31856649 0.31703406 0.44694153 0.4537032 0.4497846 0.60328904 0.52401939
#> t3   0.32121237 0.31966721 0.45065363 0.4574715 0.4535203 0.60829970 0.52837167
#> t17  0.32856745 0.32698691 0.46097264 0.4679466 0.4639050 0.62222847 0.54047025
#> t63  0.38081769 0.37898580 0.53427854 0.5423616 0.5376772 0.55666963 0.48352558
#> t73  0.32712563 0.32555203 0.45894981 0.4658932 0.4618693 0.47818395 0.41535259
#> t41  0.36269708 0.36095236 0.50885573 0.5165541 0.5120927 0.53018139 0.46051777
#> t48  0.35856531 0.35684047 0.50305897 0.5106697 0.5062590 0.52414168 0.45527166
#> t99  0.31297198 0.31146646 0.43909256 0.4457355 0.4418857 0.45749450 0.39738164
#> t60  0.33640289 0.33478465 0.47196558 0.4791059 0.4749679 0.49174520 0.42713194
#> t43  0.37531801 0.37351258 0.52656262 0.5345289 0.5299122 0.54863035 0.47654262
#> t21  0.36408847 0.36233706 0.51080783 0.5185358 0.5140572 0.53221529 0.46228443
#> t95  0.37726042 0.37544565 0.52928777 0.5372953 0.5326547 0.55146971 0.47900890
#> t59  0.35843247 0.35670827 0.50287259 0.5104805 0.5060715 0.52394749 0.45510299
#> t27  0.34963633 0.34795444 0.49053181 0.4979530 0.4936522 0.51108952 0.44393450
#> t81  0.32898305 0.32740051 0.46155573 0.4685385 0.4644918 0.48089908 0.41771096
#> t70  0.33521418 0.33360166 0.47029785 0.4774129 0.4732895 0.49000758 0.42562264
#> t20  0.34433599 0.34267960 0.48309555 0.4904042 0.4861686 0.50334162 0.43720464
#> t5   0.34645447 0.34478788 0.48606772 0.4934214 0.4891597 0.50643835 0.43989447
#> t12  0.34865339 0.34697623 0.48915276 0.4965531 0.4922644 0.50965268 0.44268645
#> t44  0.41180748 0.40982653 0.57775652 0.5864973 0.5814318 0.52293069 0.45421979
#> t23  0.34044100 0.33880335 0.47763097 0.4848570 0.4806693 0.43230649 0.37550323
#> t45  0.29724840 0.29581852 0.41703273 0.4233420 0.4196856 0.37745868 0.32786219
#> t69  0.31781943 0.31629059 0.44589341 0.4526393 0.4487298 0.40358065 0.35055184
#> t4   0.19250720 0.19158117 0.27008322 0.2741693 0.2718013 0.24445385 0.21233364
#> t96  0.11023570 0.10970542 0.15465818 0.1569980 0.1556420 0.13998199 0.12158894
#> t57  0.11074121 0.11020850 0.15536739 0.1577179 0.1563557 0.14062391 0.12214651
#> t22  0.10744802 0.10693115 0.15074713 0.1530278 0.1517061 0.13644208 0.11851416
#> t55  0.08917706 0.08874809 0.12511339 0.1270062 0.1259093 0.11324084 0.09836147
#> t58  0.08714105 0.08672187 0.12225691 0.1241065 0.1230346 0.11065542 0.09611576
#> t90  0.08003073 0.07964575 0.11228129 0.1139800 0.1129955 0.10162643 0.08827314
#> t28  0.07134580 0.07100260 0.10009654 0.1016109 0.1007333 0.09059794 0.07869375
#> t46  0.07030765 0.06996944 0.09864003 0.1001323 0.0992675 0.08927965 0.07754868
#> t71  0.07501316 0.07465231 0.10524175 0.1068339 0.1059112 0.09525490 0.08273881
#> t89  0.08007005 0.07968488 0.11233646 0.1140360 0.1130511 0.10167636 0.08831652
#> t92  0.08819367 0.08776943 0.12373371 0.1256057 0.1245208 0.11199208 0.09727679
#> t68  0.09329902 0.09285021 0.13089640 0.1328767 0.1317291 0.11847507 0.10290794
#> t7   0.09593326 0.09547179 0.13459218 0.1366284 0.1354484 0.12182014 0.10581349
#> t82  0.08941340 0.08898328 0.12544496 0.1273428 0.1262429 0.11354094 0.09862214
#> t85  0.09702938 0.09656263 0.13613000 0.1381895 0.1369960 0.12321204 0.10702249
#> t64  0.08709974 0.08668076 0.12219896 0.1240477 0.1229763 0.11060297 0.09607020
#> t98  0.09111024 0.09067197 0.12782560 0.1297595 0.1286387 0.11569567 0.10049374
#> t49  0.08615053 0.08573612 0.12086724 0.1226958 0.1216361 0.10939762 0.09502323
#> t8   0.08393071 0.08352697 0.11775287 0.1195343 0.1185019 0.10657879 0.09257478
#> t11  0.08693804 0.08651984 0.12197209 0.1238174 0.1227480 0.11039763 0.09589185
#> t75  0.09000937 0.08957639 0.12628110 0.1281916 0.1270844 0.11429773 0.09927949
#> t79  0.15170166 0.15097192 0.21283398 0.2160539 0.2141879 0.19263724 0.16732551
#> t2   0.11339085 0.11284540 0.15908479 0.1614916 0.1600968 0.14398854 0.12506905
#> t10  0.12216484 0.12157718 0.17139449 0.1739875 0.1724848 0.15513013 0.13474668
#> t62  0.11838936 0.11781986 0.16609757 0.1686104 0.1671542 0.15033585 0.13058235
#> t35  0.12244489 0.12185588 0.17178739 0.1743863 0.1728802 0.15548573 0.13505556
#> t65  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.00000000 0.00000000
#> t54  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.00000000 0.00000000
#> t83  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.00000000 0.00000000
#> t47  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.00000000 0.00000000
#> t93  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.00000000 0.00000000
#> t13  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.00000000 0.00000000
#> t76  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.00000000 0.00000000
#> t9   0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.00000000 0.00000000
#> t1   0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.00000000 0.00000000
#> t36  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.00000000 0.00000000
#> t56  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.00000000 0.00000000
#> t61  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.00000000 0.00000000
#> t52  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.00000000 0.00000000
#> t80  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.00000000 0.00000000
#> t66  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.00000000 0.00000000
#> t38  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.00000000 0.00000000
#> t30  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.00000000 0.00000000
#>             t24        t39        t53        t87        t26        t97
#> t34  0.34671393 0.33194507 0.31982534 0.33993709 0.34599087 0.34677049
#> t42  0.33138044 0.31726474 0.30568100 0.32490331 0.33068936 0.33143450
#> t50  0.31976015 0.30613943 0.29496189 0.31351015 0.31909330 0.31981232
#> t37  0.31598861 0.30252855 0.29148285 0.30981233 0.31532963 0.31604016
#> t19  0.40149685 0.38439442 0.37035969 0.39364923 0.40065954 0.40156235
#> t84  0.35753086 0.34230123 0.32980338 0.35054259 0.35678524 0.35758918
#> t94  0.36160710 0.34620384 0.33356350 0.35453916 0.36085298 0.36166609
#> t78  0.41719314 0.39942210 0.38483869 0.40903872 0.41632309 0.41726120
#> t6   0.37133357 0.35551600 0.34253566 0.36407552 0.37055917 0.37139415
#> t51  0.33652030 0.32218565 0.31042225 0.32994270 0.33581850 0.33657520
#> t88  0.33214694 0.31799858 0.30638805 0.32565482 0.33145426 0.33220112
#> t40  0.31629276 0.30281974 0.29176341 0.31011053 0.31563315 0.31634436
#> t86  0.30823646 0.29510661 0.28433190 0.30221170 0.30759364 0.30828674
#> t72  0.30675372 0.29368703 0.28296415 0.30075794 0.30611400 0.30680376
#> t14  0.43244873 0.41402785 0.39891117 0.42399613 0.43154688 0.43251928
#> t67  0.43899118 0.42029161 0.40494623 0.43041070 0.43807568 0.43906279
#> t18  0.43519964 0.41666158 0.40144874 0.42669327 0.43429205 0.43527064
#> t15  0.67166438 0.64305371 0.61957501 0.65853609 0.67026365 0.67177396
#> t32  0.82366227 0.78857699 0.75978505 0.80756304 0.81046129 0.81228750
#> t24  1.00000000 0.91016274 0.87693154 0.88100042 0.78597232 0.78774335
#> t39  0.91016274 1.00000000 0.93029348 0.84347273 0.75249251 0.75418810
#> t53  0.87693154 0.93029348 1.00000000 0.81267647 0.72501806 0.72665174
#> t87  0.88100042 0.84347273 0.81267647 1.00000000 0.77060977 0.77234619
#> t26  0.78597232 0.75249251 0.72501806 0.77060977 1.00000000 0.83718428
#> t97  0.78774335 0.75418810 0.72665174 0.77234619 0.83718428 1.00000000
#> t16  0.78261983 0.74928282 0.72192556 0.76732281 0.78098770 0.78274750
#> t91  0.66726173 0.63883860 0.61551379 0.65421949 0.66587018 0.66737059
#> t33  0.71066699 0.68039494 0.65555286 0.69677636 0.70918492 0.71078293
#> t29  0.75603801 0.72383330 0.69740523 0.74126056 0.75446132 0.75616135
#> t25  0.69424040 0.66466806 0.64040019 0.68067083 0.69279258 0.69435365
#> t100 0.68690144 0.65764172 0.63363039 0.67347532 0.68546893 0.68701349
#> t74  0.63520587 0.60814821 0.58594395 0.62279019 0.63388117 0.63530949
#> t31  0.51343828 0.49156752 0.47361977 0.50340266 0.51236753 0.51352204
#> t77  0.50818557 0.48653856 0.46877442 0.49825262 0.50712577 0.50826847
#> t3   0.51240634 0.49057954 0.47266785 0.50239089 0.51133774 0.51248993
#> t17  0.52413936 0.50181276 0.48349094 0.51389457 0.52304628 0.52422486
#> t63  0.46891533 0.44894109 0.43254968 0.45974995 0.46793742 0.46899182
#> t73  0.40280226 0.38564422 0.37156386 0.39492912 0.40196223 0.40286797
#> t41  0.44660273 0.42757894 0.41196749 0.43787347 0.44567135 0.44667558
#> t48  0.44151513 0.42270805 0.40727445 0.43288532 0.44059437 0.44158716
#> t99  0.38537433 0.36895866 0.35548751 0.37784183 0.38457064 0.38543719
#> t60  0.41422569 0.39658105 0.38210137 0.40612927 0.41336184 0.41429326
#> t43  0.46214336 0.44245759 0.42630290 0.45311035 0.46117958 0.46221875
#> t21  0.44831600 0.42921923 0.41354789 0.43955326 0.44738106 0.44838914
#> t95  0.46453512 0.44474747 0.42850918 0.45545536 0.46356635 0.46461091
#> t59  0.44135156 0.42255145 0.40712356 0.43272494 0.44043113 0.44142356
#> t27  0.43052054 0.41218179 0.39713251 0.42210562 0.42962270 0.43059077
#> t81  0.40508937 0.38783390 0.37367360 0.39717152 0.40424457 0.40515545
#> t70  0.41276199 0.39517970 0.38075119 0.40469418 0.41190119 0.41282932
#> t20  0.42399403 0.40593329 0.39111215 0.41570668 0.42310980 0.42406320
#> t5   0.42660258 0.40843073 0.39351841 0.41826425 0.42571292 0.42667218
#> t12  0.42931020 0.41102301 0.39601604 0.42091894 0.42841489 0.42938024
#> t44  0.44049504 0.42173142 0.40633347 0.43188517 0.43957641 0.44056690
#> t23  0.36415699 0.34864511 0.33591564 0.35703921 0.36339756 0.36421640
#> t45  0.31795548 0.30441164 0.29329718 0.31174075 0.31729240 0.31800735
#> t69  0.33995954 0.32547839 0.31359476 0.33331472 0.33925057 0.34001500
#> t4   0.20591775 0.19714634 0.18994827 0.20189290 0.20548831 0.20595134
#> t96  0.11791500 0.11289221 0.10877038 0.11561024 0.11766909 0.11793424
#> t57  0.11845572 0.11340990 0.10926917 0.11614040 0.11820869 0.11847505
#> t22  0.11493313 0.11003736 0.10601975 0.11268665 0.11469344 0.11495187
#> t55  0.09538937 0.09132610 0.08799167 0.09352489 0.09519044 0.09540493
#> t58  0.09321152 0.08924102 0.08598271 0.09138961 0.09301713 0.09322672
#> t90  0.08560587 0.08195935 0.07896691 0.08393263 0.08542735 0.08561984
#> t28  0.07631593 0.07306513 0.07039743 0.07482427 0.07615678 0.07632838
#> t46  0.07520546 0.07200196 0.06937308 0.07373550 0.07504862 0.07521773
#> t71  0.08023877 0.07682086 0.07401603 0.07867043 0.08007143 0.08025185
#> t89  0.08564794 0.08199962 0.07900571 0.08397387 0.08546932 0.08566191
#> t92  0.09433747 0.09031901 0.08702134 0.09249356 0.09414073 0.09435286
#> t68  0.09979846 0.09554738 0.09205883 0.09784781 0.09959034 0.09981475
#> t7   0.10261622 0.09824511 0.09465806 0.10061049 0.10240222 0.10263296
#> t82  0.09564216 0.09156813 0.08822486 0.09377275 0.09544271 0.09565777
#> t85  0.10378869 0.09936764 0.09573960 0.10176005 0.10357224 0.10380562
#> t64  0.09316733 0.08919872 0.08594196 0.09134629 0.09297304 0.09318253
#> t98  0.09745722 0.09330586 0.08989915 0.09555233 0.09725397 0.09747312
#> t49  0.09215200 0.08822663 0.08500537 0.09035081 0.09195982 0.09216703
#> t8   0.08977754 0.08595331 0.08281505 0.08802275 0.08959031 0.08979218
#> t11  0.09299437 0.08903312 0.08578241 0.09117671 0.09280043 0.09300954
#> t75  0.09627965 0.09217846 0.08881291 0.09439778 0.09607886 0.09629536
#> t79  0.16226959 0.15535744 0.14968515 0.15909788 0.16193118 0.16229606
#> t2   0.12128995 0.11612340 0.11188359 0.11891923 0.12103701 0.12130974
#> t10  0.13067516 0.12510883 0.12054095 0.12812100 0.13040264 0.13069648
#> t62  0.12663666 0.12124236 0.11681565 0.12416143 0.12637257 0.12665732
#> t35  0.13097471 0.12539562 0.12081727 0.12841469 0.13070157 0.13099608
#> t65  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t54  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t83  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t47  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t93  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t13  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t76  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t9   0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t1   0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t36  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t56  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t61  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t52  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t80  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t66  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t38  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t30  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#>             t16        t91        t33        t29        t25       t100
#> t34  0.38127035 0.32507113 0.34621695 0.36832043 0.33821437 0.43466697
#> t42  0.36440859 0.31069480 0.33090544 0.35203139 0.32325678 0.41544374
#> t50  0.35163013 0.29979988 0.31930180 0.33968695 0.31192136 0.40087566
#> t37  0.34748268 0.29626377 0.31553567 0.33568038 0.30824228 0.39614737
#> t19  0.44151339 0.37643435 0.40092134 0.42651731 0.39165431 0.50334700
#> t84  0.39316538 0.33521283 0.35701837 0.37981145 0.34876613 0.44822789
#> t94  0.39764789 0.33903462 0.36108877 0.38414171 0.35274244 0.45333817
#> t78  0.45877410 0.39115083 0.41659513 0.44319175 0.40696581 0.52302506
#> t6   0.40834378 0.34815395 0.37080130 0.39447432 0.36223047 0.46553202
#> t51  0.37006073 0.31551381 0.33603793 0.35749155 0.32827063 0.42188745
#> t88  0.36525148 0.31141345 0.33167084 0.35284565 0.32400448 0.41640468
#> t40  0.34781715 0.29654893 0.31583939 0.33600348 0.30853897 0.39652868
#> t86  0.33895789 0.28899552 0.30779463 0.32744512 0.30068016 0.38642868
#> t72  0.33732736 0.28760534 0.30631402 0.32586998 0.29923377 0.38456981
#> t14  0.47555019 0.40545413 0.43182885 0.45939805 0.42184742 0.54215064
#> t67  0.48274471 0.41158819 0.43836192 0.46634821 0.42822948 0.55035275
#> t18  0.47857528 0.40803333 0.43457582 0.46232039 0.42453089 0.54559939
#> t15  0.73860808 0.62973731 0.67070161 0.71352113 0.65519879 0.81545036
#> t32  0.80700434 0.68805197 0.73280963 0.77959430 0.71587122 0.70830360
#> t24  0.78261983 0.66726173 0.71066699 0.75603801 0.69424040 0.68690144
#> t39  0.74928282 0.63883860 0.68039494 0.72383330 0.66466806 0.65764172
#> t53  0.72192556 0.61551379 0.65555286 0.69740523 0.64040019 0.63363039
#> t87  0.76732281 0.65421949 0.69677636 0.74126056 0.68067083 0.67347532
#> t26  0.78098770 0.66587018 0.70918492 0.75446132 0.69279258 0.68546893
#> t97  0.78274750 0.66737059 0.71078293 0.75616135 0.69435365 0.68701349
#> t16  1.00000000 0.81690203 0.87004136 0.84150767 0.77272387 0.75536378
#> t91  0.81690203 1.00000000 0.90480403 0.71746951 0.65882444 0.64402323
#> t33  0.87004136 0.90480403 1.00000000 0.76414078 0.70168085 0.68591682
#> t29  0.84150767 0.71746951 0.76414078 1.00000000 0.91566337 0.72970772
#> t25  0.77272387 0.65882444 0.70168085 0.91566337 1.00000000 0.67006231
#> t100 0.75536378 0.64402323 0.68591682 0.72970772 0.67006231 1.00000000
#> t74  0.69851580 0.59555463 0.63429536 0.67479059 0.61963404 0.90382893
#> t31  0.56461184 0.48138810 0.51270232 0.54543470 0.50085154 0.64368530
#> t77  0.55883560 0.47646328 0.50745713 0.53985465 0.49572760 0.63710010
#> t3   0.56347705 0.48042058 0.51167185 0.54433845 0.49984490 0.64239158
#> t17  0.57637947 0.49142119 0.52338805 0.55680264 0.51129028 0.65710098
#> t63  0.51565135 0.43964439 0.46824318 0.49813716 0.45742005 0.58786793
#> t73  0.44294890 0.37765827 0.40222488 0.42790407 0.39292771 0.50498356
#> t41  0.49111489 0.41872460 0.44596256 0.47443409 0.43565443 0.55989516
#> t48  0.48552022 0.41395459 0.44088226 0.46902945 0.43069156 0.55351696
#> t99  0.42378396 0.36131824 0.38482192 0.40939006 0.37592702 0.48313458
#> t60  0.45551089 0.38836862 0.41363193 0.44003938 0.40407111 0.51930484
#> t43  0.50820444 0.43329515 0.46148092 0.49094319 0.45081410 0.57937809
#> t21  0.49299893 0.42033093 0.44767338 0.47625414 0.43732571 0.56204306
#> t95  0.51083458 0.43553761 0.46386925 0.49348400 0.45314722 0.58237658
#> t59  0.48534035 0.41380122 0.44071892 0.46885568 0.43053199 0.55331189
#> t27  0.47342982 0.40364630 0.42990343 0.45734969 0.41996649 0.53973331
#> t81  0.44546396 0.37980261 0.40450871 0.43033370 0.39515876 0.50785085
#> t70  0.45390130 0.38699629 0.41217033 0.43848447 0.40264329 0.51746984
#> t20  0.46625282 0.39752720 0.42338627 0.45041646 0.41359998 0.53155118
#> t5   0.46912137 0.39997292 0.42599109 0.45318758 0.41614458 0.53482146
#> t12  0.47209885 0.40251152 0.42869482 0.45606393 0.41878583 0.53821593
#> t44  0.48439847 0.41299818 0.43986363 0.46794579 0.42969648 0.55223811
#> t23  0.40045193 0.34142535 0.36363501 0.38685051 0.35522982 0.45653491
#> t45  0.34964559 0.29810786 0.31749972 0.33776981 0.31016093 0.39861319
#> t69  0.37384275 0.31873836 0.33947224 0.36114512 0.33162556 0.42619915
#> t4   0.22644123 0.19306382 0.20562258 0.21875012 0.20086975 0.25815416
#> t96  0.12966740 0.11055443 0.11774598 0.12526322 0.11502436 0.14782722
#> t57  0.13026201 0.11106140 0.11828593 0.12583764 0.11555183 0.14850511
#> t22  0.12638832 0.10775869 0.11476838 0.12209552 0.11211558 0.14408891
#> t55  0.10489667 0.08943491 0.09525263 0.10133384 0.09305093 0.11958737
#> t58  0.10250176 0.08739301 0.09307791 0.09902027 0.09092647 0.11685705
#> t90  0.09413807 0.08026213 0.08548317 0.09094066 0.08350728 0.10732204
#> t28  0.08392222 0.07155209 0.07620654 0.08107179 0.07444508 0.09567546
#> t46  0.08270107 0.07051094 0.07509766 0.07989211 0.07336183 0.09428329
#> t71  0.08823603 0.07523005 0.08012375 0.08523908 0.07827174 0.10059342
#> t89  0.09418433 0.08030157 0.08552517 0.09098534 0.08354831 0.10737477
#> t92  0.10373993 0.08844867 0.09420224 0.10021639 0.09202482 0.11826863
#> t68  0.10974521 0.09356878 0.09965541 0.10601770 0.09735194 0.12511495
#> t7   0.11284381 0.09621064 0.10246913 0.10901105 0.10010062 0.12864750
#> t82  0.10517466 0.08967193 0.09550507 0.10160239 0.09329753 0.11990429
#> t85  0.11413314 0.09730992 0.10363992 0.11025659 0.10124435 0.13011740
#> t64  0.10245317 0.08735158 0.09303379 0.09897333 0.09088337 0.11680166
#> t98  0.10717062 0.09137368 0.09731752 0.10353055 0.09506809 0.12217978
#> t49  0.10133664 0.08639963 0.09201991 0.09789472 0.08989293 0.11552876
#> t8   0.09872552 0.08417339 0.08964885 0.09537229 0.08757668 0.11255195
#> t11  0.10226297 0.08718941 0.09286107 0.09878959 0.09071465 0.11658482
#> t75  0.10587569 0.09026962 0.09614164 0.10227960 0.09391939 0.12070350
#> t79  0.17844273 0.15214029 0.16203699 0.17238190 0.15829161 0.20343350
#> t2   0.13337872 0.11371871 0.12111609 0.12884849 0.11831657 0.15205831
#> t10  0.14369934 0.12251807 0.13048785 0.13881857 0.12747171 0.16382432
#> t62  0.13925833 0.11873167 0.12645514 0.13452840 0.12353221 0.15876136
#> t35  0.14402875 0.12279892 0.13078697 0.13913678 0.12776391 0.16419986
#> t65  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t54  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t83  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t47  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t93  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t13  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t76  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t9   0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t1   0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t36  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t56  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t61  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t52  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t80  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t66  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t38  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t30  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#>             t74        t31        t77         t3        t17        t63
#> t34  0.40195433 0.36203728 0.35833347 0.36130964 0.36958286 0.42835555
#> t42  0.38417783 0.34602612 0.34248611 0.34533065 0.35323799 0.40941144
#> t50  0.37070613 0.33389225 0.33047638 0.33322118 0.34085123 0.39505489
#> t37  0.36633369 0.32995403 0.32657845 0.32929087 0.33683093 0.39039526
#> t19  0.46546557 0.41924138 0.41495235 0.41839876 0.42797920 0.49603834
#> t84  0.41449467 0.37333227 0.36951291 0.37258192 0.38111326 0.44171956
#> t94  0.41922036 0.37758867 0.37372576 0.37682976 0.38545836 0.44675564
#> t78  0.48366268 0.43563138 0.43117467 0.43475582 0.44471080 0.51543067
#> t6   0.43049651 0.38774501 0.38377820 0.38696569 0.39582638 0.45877244
#> t51  0.39013659 0.35139313 0.34779822 0.35068688 0.35871686 0.41576159
#> t88  0.38506644 0.34682649 0.34327829 0.34612941 0.35405504 0.41035843
#> t40  0.36668629 0.33027162 0.32689279 0.32960782 0.33715514 0.39077103
#> t86  0.35734642 0.32185926 0.31856649 0.32121237 0.32856745 0.38081769
#> t72  0.35562744 0.32031099 0.31703406 0.31966721 0.32698691 0.37898580
#> t14  0.50134888 0.45156121 0.44694153 0.45065363 0.46097264 0.53427854
#> t67  0.50893371 0.45839281 0.45370324 0.45747150 0.46794662 0.54236155
#> t18  0.50453808 0.45443370 0.44978463 0.45352035 0.46390500 0.53767721
#> t15  0.75408032 0.60952476 0.60328904 0.60829970 0.62222847 0.55666963
#> t32  0.65499732 0.52943576 0.52401939 0.52837167 0.54047025 0.48352558
#> t24  0.63520587 0.51343828 0.50818557 0.51240634 0.52413936 0.46891533
#> t39  0.60814821 0.49156752 0.48653856 0.49057954 0.50181276 0.44894109
#> t53  0.58594395 0.47361977 0.46877442 0.47266785 0.48349094 0.43254968
#> t87  0.62279019 0.50340266 0.49825262 0.50239089 0.51389457 0.45974995
#> t26  0.63388117 0.51236753 0.50712577 0.51133774 0.52304628 0.46793742
#> t97  0.63530949 0.51352204 0.50826847 0.51248993 0.52422486 0.46899182
#> t16  0.69851580 0.56461184 0.55883560 0.56347705 0.57637947 0.51565135
#> t91  0.59555463 0.48138810 0.47646328 0.48042058 0.49142119 0.43964439
#> t33  0.63429536 0.51270232 0.50745713 0.51167185 0.52338805 0.46824318
#> t29  0.67479059 0.54543470 0.53985465 0.54433845 0.55680264 0.49813716
#> t25  0.61963404 0.50085154 0.49572760 0.49984490 0.51129028 0.45742005
#> t100 0.90382893 0.64368530 0.63710010 0.64239158 0.65710098 0.58786793
#> t74  1.00000000 0.59524214 0.58915254 0.59404578 0.60764817 0.54362554
#> t31  0.59524214 1.00000000 0.88759688 0.72345241 0.74001793 0.48963948
#> t77  0.58915254 0.88759688 1.00000000 0.71605115 0.73244720 0.48463024
#> t3   0.59404578 0.72345241 0.71605115 1.00000000 0.83692455 0.48865537
#> t17  0.60764817 0.74001793 0.73244720 0.83692455 1.00000000 0.49984454
#> t63  0.54362554 0.48963948 0.48463024 0.48865537 0.49984454 1.00000000
#> t73  0.46697896 0.42060448 0.41630150 0.41975912 0.42937071 0.84208330
#> t41  0.51775797 0.46634076 0.46156987 0.46540348 0.47606022 0.93365092
#> t48  0.51185979 0.46102831 0.45631178 0.46010170 0.47063705 0.86785239
#> t99  0.44677432 0.40240631 0.39828950 0.40159752 0.41079325 0.70004438
#> t60  0.48022244 0.43253278 0.42810777 0.43166345 0.44154763 0.75245377
#> t43  0.53577463 0.48256823 0.47763133 0.48159833 0.49262590 0.83949771
#> t21  0.51974422 0.46812975 0.46334057 0.46718888 0.47788651 0.81437989
#> t95  0.53854746 0.48506569 0.48010325 0.48409078 0.49517542 0.68638768
#> t59  0.51167015 0.46085751 0.45614272 0.45993124 0.47046269 0.65213211
#> t27  0.49911348 0.44954780 0.44494872 0.44864427 0.45891727 0.63612842
#> t81  0.46963047 0.42299267 0.41866526 0.42214251 0.43180868 0.59855183
#> t70  0.47852553 0.43100439 0.42659502 0.43013813 0.43998738 0.60988875
#> t20  0.49154713 0.44273284 0.43820348 0.44184301 0.45196028 0.62648498
#> t5   0.49457129 0.44545669 0.44089946 0.44456138 0.45474089 0.63033933
#> t12  0.49771030 0.44828397 0.44369782 0.44738298 0.45762710 0.63434005
#> t44  0.51067718 0.45996314 0.45525750 0.45903867 0.46954968 0.54421954
#> t23  0.42217652 0.38025125 0.37636111 0.37948700 0.38817644 0.44990597
#> t45  0.36861393 0.33200782 0.32861123 0.33134053 0.33892753 0.39282527
#> t69  0.39412379 0.35498437 0.35135272 0.35427090 0.36238295 0.42001068
#> t4   0.23872571 0.21501848 0.21281874 0.21458632 0.21949989 0.25440573
#> t96  0.13670188 0.12312637 0.12186673 0.12287890 0.12569256 0.14568075
#> t57  0.13732875 0.12369099 0.12242557 0.12344239 0.12626895 0.14634880
#> t22  0.13324491 0.12001270 0.11878492 0.11977149 0.12251401 0.14199672
#> t55  0.11058733 0.09960519 0.09858618 0.09940500 0.10168116 0.11785094
#> t58  0.10806250 0.09733109 0.09633535 0.09713547 0.09935966 0.11516027
#> t90  0.09924508 0.08938931 0.08847481 0.08920965 0.09125236 0.10576371
#> t28  0.08847501 0.07968879 0.07887354 0.07952862 0.08134966 0.09428624
#> t46  0.08718762 0.07852924 0.07772585 0.07837140 0.08016594 0.09291428
#> t71  0.09302285 0.08378499 0.08292783 0.08361660 0.08553124 0.09913279
#> t89  0.09929384 0.08943323 0.08851828 0.08925348 0.09129719 0.10581568
#> t92  0.10936784 0.09850680 0.09749903 0.09830882 0.10055988 0.11655135
#> t68  0.11569891 0.10420915 0.10314304 0.10399971 0.10638108 0.12329826
#> t7   0.11896561 0.10715144 0.10605523 0.10693608 0.10938469 0.12677952
#> t82  0.11088041 0.09986916 0.09884745 0.09966844 0.10195063 0.11816327
#> t85  0.12032488 0.10837573 0.10726700 0.10815791 0.11063449 0.12822808
#> t64  0.10801127 0.09728495 0.09628968 0.09708942 0.09931256 0.11510568
#> t98  0.11298465 0.10176443 0.10072333 0.10155990 0.10388540 0.12040572
#> t49  0.10683417 0.09622475 0.09524032 0.09603135 0.09823026 0.11385127
#> t8   0.10408140 0.09374534 0.09278628 0.09355692 0.09569918 0.11091768
#> t11  0.10781075 0.09710434 0.09611092 0.09690918 0.09912819 0.11489199
#> t75  0.11161946 0.10053482 0.09950630 0.10033276 0.10263017 0.11895087
#> t79  0.18812329 0.16944125 0.16770778 0.16910069 0.17297274 0.20047962
#> t2   0.14061454 0.12665048 0.12535478 0.12639592 0.12929012 0.14985041
#> t10  0.15149506 0.13645047 0.13505452 0.13617623 0.13929437 0.16144557
#> t62  0.14681312 0.13223349 0.13088068 0.13196772 0.13498950 0.15645612
#> t35  0.15184233 0.13676326 0.13536411 0.13648838 0.13961368 0.16181566
#> t65  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t54  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t83  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t47  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t93  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t13  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t76  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t9   0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t1   0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t36  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t56  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t61  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t52  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t80  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t66  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t38  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t30  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#>             t73        t41        t48        t99        t60        t43
#> t34  0.36796106 0.40797292 0.40332539 0.35204059 0.37839640 0.42216934
#> t42  0.35168791 0.38993024 0.38548825 0.33647153 0.36166175 0.40349882
#> t50  0.33935552 0.37625683 0.37197059 0.32467272 0.34897960 0.38934960
#> t37  0.33535285 0.37181892 0.36758324 0.32084323 0.34486343 0.38475727
#> t19  0.42610115 0.47243513 0.46705326 0.40766516 0.43818534 0.48887467
#> t84  0.37944086 0.42070103 0.41590850 0.36302370 0.39020177 0.43534035
#> t94  0.38376690 0.42549748 0.42065031 0.36716257 0.39465049 0.44030371
#> t78  0.44275933 0.49090471 0.48531244 0.42360259 0.45531594 0.50798694
#> t6   0.39408943 0.43694247 0.43196492 0.37703847 0.40526577 0.45214695
#> t51  0.35714275 0.39597823 0.39146734 0.34169035 0.36727128 0.40975726
#> t88  0.35250138 0.39083217 0.38637989 0.33724980 0.36249829 0.40443213
#> t40  0.33567564 0.37217680 0.36793705 0.32115206 0.34519537 0.38512761
#> t86  0.32712563 0.36269708 0.35856531 0.31297198 0.33640289 0.37531801
#> t72  0.32555203 0.36095236 0.35684047 0.31146646 0.33478465 0.37351258
#> t14  0.45894981 0.50885573 0.50305897 0.43909256 0.47196558 0.52656262
#> t67  0.46589319 0.51655413 0.51066966 0.44573552 0.47910588 0.53452890
#> t18  0.46186930 0.51209269 0.50625904 0.44188573 0.47496787 0.52991221
#> t15  0.47818395 0.53018139 0.52414168 0.45749450 0.49174520 0.54863035
#> t32  0.41535259 0.46051777 0.45527166 0.39738164 0.42713194 0.47654262
#> t24  0.40280226 0.44660273 0.44151513 0.38537433 0.41422569 0.46214336
#> t39  0.38564422 0.42757894 0.42270805 0.36895866 0.39658105 0.44245759
#> t53  0.37156386 0.41196749 0.40727445 0.35548751 0.38210137 0.42630290
#> t87  0.39492912 0.43787347 0.43288532 0.37784183 0.40612927 0.45311035
#> t26  0.40196223 0.44567135 0.44059437 0.38457064 0.41336184 0.46117958
#> t97  0.40286797 0.44667558 0.44158716 0.38543719 0.41429326 0.46221875
#> t16  0.44294890 0.49111489 0.48552022 0.42378396 0.45551089 0.50820444
#> t91  0.37765827 0.41872460 0.41395459 0.36131824 0.38836862 0.43329515
#> t33  0.40222488 0.44596256 0.44088226 0.38482192 0.41363193 0.46148092
#> t29  0.42790407 0.47443409 0.46902945 0.40939006 0.44003938 0.49094319
#> t25  0.39292771 0.43565443 0.43069156 0.37592702 0.40407111 0.45081410
#> t100 0.50498356 0.55989516 0.55351696 0.48313458 0.51930484 0.57937809
#> t74  0.46697896 0.51775797 0.51185979 0.44677432 0.48022244 0.53577463
#> t31  0.42060448 0.46634076 0.46102831 0.40240631 0.43253278 0.48256823
#> t77  0.41630150 0.46156987 0.45631178 0.39828950 0.42810777 0.47763133
#> t3   0.41975912 0.46540348 0.46010170 0.40159752 0.43166345 0.48159833
#> t17  0.42937071 0.47606022 0.47063705 0.41079325 0.44154763 0.49262590
#> t63  0.84208330 0.93365092 0.86785239 0.70004438 0.75245377 0.83949771
#> t73  1.00000000 0.87449119 0.74549259 0.60134408 0.64636420 0.72113568
#> t41  0.87449119 1.00000000 0.82655700 0.66673387 0.71664945 0.79955153
#> t48  0.74549259 0.82655700 1.00000000 0.65913859 0.70848554 0.79044322
#> t99  0.60134408 0.66673387 0.65913859 1.00000000 0.88772048 0.77387449
#> t60  0.64636420 0.71664945 0.70848554 0.88772048 1.00000000 0.83181124
#> t43  0.72113568 0.79955153 0.79044322 0.77387449 0.83181124 1.00000000
#> t21  0.69955926 0.77562890 0.76679311 0.72051243 0.77445419 0.86404313
#> t95  0.58961286 0.65372701 0.64627990 0.56410225 0.60633423 0.67647504
#> t59  0.56018704 0.62110143 0.61402599 0.53594958 0.57607389 0.64271418
#> t27  0.54643973 0.60585925 0.59895744 0.52279708 0.56193671 0.62694161
#> t81  0.51416112 0.57007069 0.56357657 0.49191506 0.52874268 0.58990769
#> t70  0.52389962 0.58086815 0.57425103 0.50123221 0.53875737 0.60108088
#> t20  0.53815593 0.59667468 0.58987750 0.51487169 0.55341798 0.61743744
#> t5   0.54146685 0.60034562 0.59350662 0.51803936 0.55682280 0.62123612
#> t12  0.54490350 0.60415598 0.59727357 0.52132732 0.56035691 0.62517906
#> t44  0.46748921 0.51832371 0.51241908 0.44726249 0.48074716 0.53636005
#> t23  0.38647306 0.42849790 0.42361655 0.36975164 0.39743339 0.44340853
#> t45  0.33744025 0.37413330 0.36987126 0.32284032 0.34701002 0.38715218
#> t69  0.36079275 0.40002513 0.39546814 0.34518243 0.37102480 0.41394499
#> t4   0.21853668 0.24230023 0.23954000 0.20908132 0.22473436 0.25073166
#> t96  0.12514100 0.13874876 0.13716816 0.11972656 0.12868999 0.14357686
#> t57  0.12571486 0.13938502 0.13779718 0.12027559 0.12928012 0.14423526
#> t22  0.12197639 0.13524003 0.13369941 0.11669887 0.12543563 0.13994604
#> t55  0.10123496 0.11224319 0.11096454 0.09685486 0.10410598 0.11614897
#> t58  0.09892365 0.10968055 0.10843110 0.09464355 0.10172912 0.11349715
#> t90  0.09085193 0.10073111 0.09958360 0.08692106 0.09342848 0.10423629
#> t28  0.08099268 0.08979978 0.08877680 0.07748840 0.08328963 0.09292458
#> t46  0.07981416 0.08849310 0.08748501 0.07636086 0.08207768 0.09157244
#> t71  0.08515591 0.09441571 0.09334015 0.08147150 0.08757093 0.09770114
#> t89  0.09089657 0.10078060 0.09963253 0.08696377 0.09347438 0.10428751
#> t92  0.10011860 0.11100544 0.10974089 0.09578680 0.10295796 0.11486815
#> t68  0.10591426 0.11743131 0.11609356 0.10133170 0.10891798 0.12151762
#> t7   0.10890469 0.12074692 0.11937140 0.10419274 0.11199321 0.12494860
#> t82  0.10150325 0.11254065 0.11125862 0.09711154 0.10438188 0.11645678
#> t85  0.11014901 0.12212655 0.12073531 0.10538322 0.11327283 0.12637624
#> t64  0.09887676 0.10962856 0.10837970 0.09459869 0.10168090 0.11344336
#> t98  0.10342954 0.11467640 0.11337003 0.09895448 0.10636279 0.11866685
#> t49  0.09779921 0.10843383 0.10719858 0.09356776 0.10057278 0.11220706
#> t8   0.09527923 0.10563984 0.10443641 0.09115681 0.09798134 0.10931584
#> t11  0.09869320 0.10942504 0.10817849 0.09442307 0.10149213 0.11323275
#> t75  0.10217981 0.11329078 0.11200020 0.09775882 0.10507762 0.11723301
#> t79  0.17221370 0.19094012 0.18876497 0.16476258 0.17709767 0.19758434
#> t2   0.12872277 0.14272001 0.14109418 0.12315336 0.13237334 0.14768630
#> t10  0.13868312 0.15376344 0.15201180 0.13268276 0.14261616 0.15911402
#> t62  0.13439714 0.14901140 0.14731390 0.12858222 0.13820863 0.15419662
#> t35  0.13900103 0.15411592 0.15236026 0.13298691 0.14294308 0.15947876
#> t65  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t54  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t83  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t47  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t93  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t13  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t76  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t9   0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t1   0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t36  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t56  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t61  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t52  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t80  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t66  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t38  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t30  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#>             t21        t95        t59        t27        t81        t70
#> t34  0.40953801 0.42435422 0.40317596 0.39328180 0.37005034 0.37705931
#> t42  0.39142611 0.40558708 0.38534543 0.37588883 0.35368480 0.36038379
#> t50  0.37770024 0.39136463 0.37183278 0.36270780 0.34128237 0.34774646
#> t37  0.37324530 0.38674853 0.36744706 0.35842970 0.33725698 0.34364482
#> t19  0.47424751 0.49140478 0.46688022 0.45542272 0.42852055 0.43663698
#> t84  0.42231494 0.43759340 0.41575441 0.40555156 0.38159533 0.38882296
#> t94  0.42712979 0.44258244 0.42049447 0.41017529 0.38594593 0.39325596
#> t78  0.49278794 0.51061596 0.48513264 0.47322721 0.44527332 0.45370705
#> t6   0.43861869 0.45448698 0.43180489 0.42120815 0.39632706 0.40383372
#> t51  0.39749730 0.41187791 0.39132230 0.38171903 0.35917060 0.36597350
#> t88  0.39233149 0.40652521 0.38623675 0.37675828 0.35450288 0.36121737
#> t40  0.37360456 0.38712078 0.36780073 0.35877470 0.33758160 0.34397559
#> t86  0.36408847 0.37726042 0.35843247 0.34963633 0.32898305 0.33521418
#> t72  0.36233706 0.37544565 0.35670827 0.34795444 0.32740051 0.33360166
#> t14  0.51080783 0.52928777 0.50287259 0.49053181 0.46155573 0.47029785
#> t67  0.51853576 0.53729528 0.51048047 0.49795299 0.46853853 0.47741291
#> t18  0.51405720 0.53265470 0.50607148 0.49365220 0.46449179 0.47328953
#> t15  0.53221529 0.55146971 0.52394749 0.51108952 0.48089908 0.49000758
#> t32  0.46228443 0.47900890 0.45510299 0.44393450 0.41771096 0.42562264
#> t24  0.44831600 0.46453512 0.44135156 0.43052054 0.40508937 0.41276199
#> t39  0.42921923 0.44474747 0.42255145 0.41218179 0.38783390 0.39517970
#> t53  0.41354789 0.42850918 0.40712356 0.39713251 0.37367360 0.38075119
#> t87  0.43955326 0.45545536 0.43272494 0.42210562 0.39717152 0.40469418
#> t26  0.44738106 0.46356635 0.44043113 0.42962270 0.40424457 0.41190119
#> t97  0.44838914 0.46461091 0.44142356 0.43059077 0.40515545 0.41282932
#> t16  0.49299893 0.51083458 0.48534035 0.47342982 0.44546396 0.45390130
#> t91  0.42033093 0.43553761 0.41380122 0.40364630 0.37980261 0.38699629
#> t33  0.44767338 0.46386925 0.44071892 0.42990343 0.40450871 0.41217033
#> t29  0.47625414 0.49348400 0.46885568 0.45734969 0.43033370 0.43848447
#> t25  0.43732571 0.45314722 0.43053199 0.41996649 0.39515876 0.40264329
#> t100 0.56204306 0.58237658 0.55331189 0.53973331 0.50785085 0.51746984
#> t74  0.51974422 0.53854746 0.51167015 0.49911348 0.46963047 0.47852553
#> t31  0.46812975 0.48506569 0.46085751 0.44954780 0.42299267 0.43100439
#> t77  0.46334057 0.48010325 0.45614272 0.44494872 0.41866526 0.42659502
#> t3   0.46718888 0.48409078 0.45993124 0.44864427 0.42214251 0.43013813
#> t17  0.47788651 0.49517542 0.47046269 0.45891727 0.43180868 0.43998738
#> t63  0.81437989 0.68638768 0.65213211 0.63612842 0.59855183 0.60988875
#> t73  0.69955926 0.58961286 0.56018704 0.54643973 0.51416112 0.52389962
#> t41  0.77562890 0.65372701 0.62110143 0.60585925 0.57007069 0.58086815
#> t48  0.76679311 0.64627990 0.61402599 0.59895744 0.56357657 0.57425103
#> t99  0.72051243 0.56410225 0.53594958 0.52279708 0.49191506 0.50123221
#> t60  0.77445419 0.60633423 0.57607389 0.56193671 0.52874268 0.53875737
#> t43  0.86404313 0.67647504 0.64271418 0.62694161 0.58990769 0.60108088
#> t21  1.00000000 0.65623486 0.62348413 0.60818347 0.57225762 0.58309650
#> t95  0.65623486 1.00000000 0.87679959 0.82871540 0.77976256 0.71518844
#> t59  0.62348413 0.87679959 1.00000000 0.78735668 0.74084692 0.67949551
#> t27  0.60818347 0.82871540 0.78735668 1.00000000 0.93484102 0.66282031
#> t81  0.57225762 0.77976256 0.74084692 0.93484102 1.00000000 0.62366701
#> t70  0.58309650 0.71518844 0.67949551 0.66282031 0.62366701 1.00000000
#> t20  0.59896367 0.73465008 0.69798588 0.68085691 0.64063818 0.91342126
#> t5   0.60264870 0.73916990 0.70228013 0.68504577 0.64457960 0.87151928
#> t12  0.60647367 0.74386136 0.70673745 0.68939371 0.64867070 0.85941972
#> t44  0.52031212 0.53913591 0.51222924 0.49965884 0.47014362 0.47904840
#> t23  0.43014172 0.44570334 0.42345961 0.41306767 0.38866745 0.39602903
#> t45  0.37556856 0.38915584 0.36973422 0.36066074 0.33935623 0.34578383
#> t69  0.40155973 0.41608730 0.39532162 0.38562020 0.36284133 0.36971375
#> t4   0.24322975 0.25202929 0.23945125 0.23357498 0.21977753 0.22394024
#> t96  0.13928103 0.14431992 0.13711734 0.13375241 0.12585155 0.12823525
#> t57  0.13991973 0.14498173 0.13774612 0.13436576 0.12642867 0.12882330
#> t22  0.13575884 0.14067031 0.13364987 0.13037003 0.12266897 0.12499239
#> t55  0.11267379 0.11675008 0.11092343 0.10820131 0.10180978 0.10373811
#> t58  0.11010131 0.11408454 0.10839092 0.10573095 0.09948534 0.10136965
#> t90  0.10111754 0.10477575 0.09954671 0.09710378 0.09136778 0.09309834
#> t28  0.09014427 0.09340550 0.08874391 0.08656608 0.08145256 0.08299532
#> t46  0.08883258 0.09204636 0.08745260 0.08530646 0.08026734 0.08178765
#> t71  0.09477792 0.09820678 0.09330557 0.09101580 0.08563943 0.08726149
#> t89  0.10116722 0.10482724 0.09959562 0.09715149 0.09141268 0.09314408
#> t92  0.11143128 0.11546263 0.10970023 0.10700813 0.10068708 0.10259415
#> t68  0.11788180 0.12214652 0.11605055 0.11320260 0.10651564 0.10853311
#> t7   0.12121013 0.12559526 0.11932717 0.11639881 0.10952305 0.11159748
#> t82  0.11297239 0.11705949 0.11121740 0.10848806 0.10207959 0.10401303
#> t85  0.12259505 0.12703028 0.12069058 0.11772876 0.11077444 0.11287257
#> t64  0.11004912 0.11403047 0.10833954 0.10568083 0.09943818 0.10132160
#> t98  0.11511633 0.11928099 0.11332803 0.11054690 0.10401681 0.10598695
#> t49  0.10884981 0.11278777 0.10715887 0.10452913 0.09835451 0.10021740
#> t8   0.10604510 0.10988159 0.10439772 0.10183574 0.09582023 0.09763512
#> t11  0.10984482 0.11381877 0.10813841 0.10548463 0.09925358 0.10113350
#> t75  0.11372539 0.11783973 0.11195870 0.10921117 0.10275998 0.10470632
#> t79  0.19167261 0.19860692 0.18869504 0.18406435 0.17319153 0.17647188
#> t2   0.14326752 0.14845063 0.14104190 0.13758065 0.12945366 0.13190559
#> t10  0.15435331 0.15993749 0.15195548 0.14822641 0.13947056 0.14211222
#> t62  0.14958305 0.15499464 0.14725932 0.14364549 0.13516025 0.13772026
#> t35  0.15470714 0.16030412 0.15230381 0.14856619 0.13979028 0.14243798
#> t65  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t54  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t83  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t47  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t93  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t13  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t76  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t9   0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t1   0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t36  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t56  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t61  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t52  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t80  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t66  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t38  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t30  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#>             t20         t5        t12       t44       t23       t45       t69
#> t34  0.38731981 0.38970273 0.39217615 0.4632138 0.3829386 0.3343542 0.3574932
#> t42  0.37019052 0.37246805 0.37483208 0.4427281 0.3660031 0.3195673 0.3416829
#> t50  0.35720930 0.35940697 0.36168810 0.4272033 0.3531687 0.3083613 0.3297014
#> t37  0.35299605 0.35516780 0.35742203 0.4221645 0.3490031 0.3047242 0.3258126
#> t19  0.44851870 0.45127814 0.45414237 0.5364044 0.4434452 0.3871842 0.4139792
#> t84  0.39940357 0.40186084 0.40441142 0.4776654 0.3948857 0.3447855 0.3686464
#> t94  0.40395721 0.40644249 0.40902215 0.4831113 0.3993878 0.3487165 0.3728493
#> t78  0.46605328 0.46892060 0.47189680 0.5573749 0.4607815 0.4023210 0.4301635
#> t6   0.41482281 0.41737494 0.42002399 0.4961060 0.4101305 0.3580962 0.3828782
#> t51  0.37593233 0.37824519 0.38064589 0.4495950 0.3716799 0.3245240 0.3469826
#> t88  0.37104678 0.37332959 0.37569908 0.4437522 0.3668496 0.3203065 0.3424733
#> t40  0.35333582 0.35550966 0.35776606 0.4225708 0.3493390 0.3050175 0.3261262
#> t86  0.34433599 0.34645447 0.34865339 0.4118075 0.3404410 0.2972484 0.3178194
#> t72  0.34267960 0.34478788 0.34697623 0.4098265 0.3388033 0.2958185 0.3162906
#> t14  0.48309555 0.48606772 0.48915276 0.5777565 0.4776310 0.4170327 0.4458934
#> t67  0.49040423 0.49342136 0.49655307 0.5864973 0.4848570 0.4233420 0.4526393
#> t18  0.48616864 0.48915971 0.49226438 0.5814318 0.4806693 0.4196856 0.4487298
#> t15  0.50334162 0.50643835 0.50965268 0.5229307 0.4323065 0.3774587 0.4035807
#> t32  0.43720464 0.43989447 0.44268645 0.4542198 0.3755032 0.3278622 0.3505518
#> t24  0.42399403 0.42660258 0.42931020 0.4404950 0.3641570 0.3179555 0.3399595
#> t39  0.40593329 0.40843073 0.41102301 0.4217314 0.3486451 0.3044116 0.3254784
#> t53  0.39111215 0.39351841 0.39601604 0.4063335 0.3359156 0.2932972 0.3135948
#> t87  0.41570668 0.41826425 0.42091894 0.4318852 0.3570392 0.3117408 0.3333147
#> t26  0.42310980 0.42571292 0.42841489 0.4395764 0.3633976 0.3172924 0.3392506
#> t97  0.42406320 0.42667218 0.42938024 0.4405669 0.3642164 0.3180074 0.3400150
#> t16  0.46625282 0.46912137 0.47209885 0.4843985 0.4004519 0.3496456 0.3738428
#> t91  0.39752720 0.39997292 0.40251152 0.4129982 0.3414253 0.2981079 0.3187384
#> t33  0.42338627 0.42599109 0.42869482 0.4398636 0.3636350 0.3174997 0.3394722
#> t29  0.45041646 0.45318758 0.45606393 0.4679458 0.3868505 0.3377698 0.3611451
#> t25  0.41359998 0.41614458 0.41878583 0.4296965 0.3552298 0.3101609 0.3316256
#> t100 0.53155118 0.53482146 0.53821593 0.5522381 0.4565349 0.3986132 0.4261991
#> t74  0.49154713 0.49457129 0.49771030 0.5106772 0.4221765 0.3686139 0.3941238
#> t31  0.44273284 0.44545669 0.44828397 0.4599631 0.3802513 0.3320078 0.3549844
#> t77  0.43820348 0.44089946 0.44369782 0.4552575 0.3763611 0.3286112 0.3513527
#> t3   0.44184301 0.44456138 0.44738298 0.4590387 0.3794870 0.3313405 0.3542709
#> t17  0.45196028 0.45474089 0.45762710 0.4695497 0.3881764 0.3389275 0.3623830
#> t63  0.62648498 0.63033933 0.63434005 0.5442195 0.4499060 0.3928253 0.4200107
#> t73  0.53815593 0.54146685 0.54490350 0.4674892 0.3864731 0.3374403 0.3607927
#> t41  0.59667468 0.60034562 0.60415598 0.5183237 0.4284979 0.3741333 0.4000251
#> t48  0.58987750 0.59350662 0.59727357 0.5124191 0.4236166 0.3698713 0.3954681
#> t99  0.51487169 0.51803936 0.52132732 0.4472625 0.3697516 0.3228403 0.3451824
#> t60  0.55341798 0.55682280 0.56035691 0.4807472 0.3974334 0.3470100 0.3710248
#> t43  0.61743744 0.62123612 0.62517906 0.5363600 0.4434085 0.3871522 0.4139450
#> t21  0.59896367 0.60264870 0.60647367 0.5203121 0.4301417 0.3755686 0.4015597
#> t95  0.73465008 0.73916990 0.74386136 0.5391359 0.4457033 0.3891558 0.4160873
#> t59  0.69798588 0.70228013 0.70673745 0.5122292 0.4234596 0.3697342 0.3953216
#> t27  0.68085691 0.68504577 0.68939371 0.4996588 0.4130677 0.3606607 0.3856202
#> t81  0.64063818 0.64457960 0.64867070 0.4701436 0.3886674 0.3393562 0.3628413
#> t70  0.91342126 0.87151928 0.85941972 0.4790484 0.3960290 0.3457838 0.3697137
#> t20  1.00000000 0.89523498 0.88280617 0.4920842 0.4068057 0.3551933 0.3797744
#> t5   0.89523498 1.00000000 0.88823749 0.4951117 0.4093085 0.3573785 0.3821109
#> t12  0.88280617 0.88823749 1.00000000 0.4982541 0.4119064 0.3596468 0.3845361
#> t44  0.49208422 0.49511169 0.49825413 1.0000000 0.8171209 0.6765864 0.7234094
#> t23  0.40680574 0.40930854 0.41190639 0.8171209 1.0000000 0.5593336 0.5980422
#> t45  0.35519327 0.35737854 0.35964680 0.6765864 0.5593336 1.0000000 0.8748011
#> t69  0.37977436 0.38211086 0.38453610 0.7234094 0.5980422 0.8748011 1.0000000
#> t4   0.23003409 0.23144933 0.23291833 0.3825893 0.3162863 0.2761583 0.2952698
#> t96  0.13172478 0.13253519 0.13337638 0.2190827 0.1811155 0.1581370 0.1690808
#> t57  0.13232883 0.13314296 0.13398801 0.2200873 0.1819461 0.1588621 0.1698562
#> t22  0.12839367 0.12918360 0.13000352 0.2135425 0.1765354 0.1541379 0.1648050
#> t55  0.10656102 0.10721662 0.10789712 0.1772307 0.1465165 0.1279276 0.1367808
#> t58  0.10412811 0.10476874 0.10543370 0.1731843 0.1431714 0.1250069 0.1336580
#> t90  0.09563172 0.09622008 0.09683078 0.1590533 0.1314892 0.1148069 0.1227521
#> t28  0.08525378 0.08577829 0.08632272 0.1417928 0.1172200 0.1023480 0.1094310
#> t46  0.08401325 0.08453013 0.08506664 0.1397296 0.1155144 0.1008588 0.1078387
#> t71  0.08963604 0.09018751 0.09075992 0.1490813 0.1232454 0.1076090 0.1150561
#> t89  0.09567871 0.09626736 0.09687836 0.1591314 0.1315538 0.1148633 0.1228124
#> t92  0.10538593 0.10603430 0.10670729 0.1752763 0.1449008 0.1265169 0.1352725
#> t68  0.11148650 0.11217240 0.11288435 0.1854227 0.1532888 0.1338407 0.1431031
#> t7   0.11463426 0.11533953 0.11607158 0.1906580 0.1576169 0.1376196 0.1471436
#> t82  0.10684343 0.10750076 0.10818306 0.1777004 0.1469048 0.1282666 0.1371433
#> t85  0.11594404 0.11665737 0.11739779 0.1928364 0.1594178 0.1391920 0.1488248
#> t64  0.10407875 0.10471908 0.10538373 0.1731022 0.1431035 0.1249476 0.1335946
#> t98  0.10887105 0.10954086 0.11023611 0.1810727 0.1496927 0.1307008 0.1397460
#> t49  0.10294451 0.10357786 0.10423526 0.1712158 0.1415440 0.1235859 0.1321387
#> t8   0.10029196 0.10090899 0.10154945 0.1668041 0.1378968 0.1204015 0.1287339
#> t11  0.10388553 0.10452467 0.10518808 0.1727809 0.1428379 0.1247157 0.1333466
#> t75  0.10755558 0.10821729 0.10890414 0.1788848 0.1478840 0.1291216 0.1380574
#> t79  0.18127401 0.18238927 0.18354688 0.3014923 0.2492435 0.2176213 0.2326818
#> t2   0.13549499 0.13632860 0.13719387 0.2253533 0.1862994 0.1626631 0.1739202
#> t10  0.14597936 0.14687747 0.14780970 0.2427907 0.2007149 0.1752497 0.1873779
#> t62  0.14146789 0.14233825 0.14324166 0.2352873 0.1945119 0.1698337 0.1815870
#> t35  0.14631399 0.14721416 0.14814852 0.2433473 0.2011750 0.1756515 0.1878074
#> t65  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t54  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t83  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t47  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t93  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t13  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t76  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t9   0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t1   0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t36  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t56  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t61  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t52  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t80  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t66  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t38  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t30  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000
#>             t4       t96       t57       t22        t55        t58        t90
#> t34  0.2165381 0.1239965 0.1245651 0.1208609 0.10030913 0.09801896 0.09002105
#> t42  0.2069616 0.1185128 0.1190562 0.1155158 0.09587294 0.09368405 0.08603985
#> t50  0.1997043 0.1143570 0.1148814 0.1114651 0.09251103 0.09039889 0.08302275
#> t37  0.1973488 0.1130081 0.1135264 0.1101503 0.09141987 0.08933265 0.08204350
#> t19  0.2507524 0.1435887 0.1442472 0.1399576 0.11615858 0.11350655 0.10424492
#> t84  0.2232937 0.1278650 0.1284514 0.1246315 0.10343861 0.10107699 0.09282956
#> t94  0.2258395 0.1293228 0.1299159 0.1260525 0.10461793 0.10222938 0.09388792
#> t78  0.2605554 0.1492023 0.1498865 0.1454292 0.12069974 0.11794402 0.10832032
#> t6   0.2319141 0.1328013 0.1334103 0.1294430 0.10743193 0.10497914 0.09641331
#> t51  0.2101717 0.1203509 0.1209028 0.1173075 0.09735997 0.09513713 0.08737436
#> t88  0.2074403 0.1187869 0.1193316 0.1157830 0.09609470 0.09390074 0.08623886
#> t40  0.1975387 0.1131169 0.1136356 0.1102564 0.09150786 0.08941863 0.08212247
#> t86  0.1925072 0.1102357 0.1107412 0.1074480 0.08917706 0.08714105 0.08003073
#> t72  0.1915812 0.1097054 0.1102085 0.1069312 0.08874809 0.08672187 0.07964575
#> t14  0.2700832 0.1546582 0.1553674 0.1507471 0.12511339 0.12225691 0.11228129
#> t67  0.2741693 0.1569980 0.1577179 0.1530278 0.12700621 0.12410651 0.11397998
#> t18  0.2718013 0.1556420 0.1563557 0.1517061 0.12590927 0.12303461 0.11299554
#> t15  0.2444539 0.1399820 0.1406239 0.1364421 0.11324084 0.11065542 0.10162643
#> t32  0.2123336 0.1215889 0.1221465 0.1185142 0.09836147 0.09611576 0.08827314
#> t24  0.2059177 0.1179150 0.1184557 0.1149331 0.09538937 0.09321152 0.08560587
#> t39  0.1971463 0.1128922 0.1134099 0.1100374 0.09132610 0.08924102 0.08195935
#> t53  0.1899483 0.1087704 0.1092692 0.1060198 0.08799167 0.08598271 0.07896691
#> t87  0.2018929 0.1156102 0.1161404 0.1126867 0.09352489 0.09138961 0.08393263
#> t26  0.2054883 0.1176691 0.1182087 0.1146934 0.09519044 0.09301713 0.08542735
#> t97  0.2059513 0.1179342 0.1184750 0.1149519 0.09540493 0.09322672 0.08561984
#> t16  0.2264412 0.1296674 0.1302620 0.1263883 0.10489667 0.10250176 0.09413807
#> t91  0.1930638 0.1105544 0.1110614 0.1077587 0.08943491 0.08739301 0.08026213
#> t33  0.2056226 0.1177460 0.1182859 0.1147684 0.09525263 0.09307791 0.08548317
#> t29  0.2187501 0.1252632 0.1258376 0.1220955 0.10133384 0.09902027 0.09094066
#> t25  0.2008698 0.1150244 0.1155518 0.1121156 0.09305093 0.09092647 0.08350728
#> t100 0.2581542 0.1478272 0.1485051 0.1440889 0.11958737 0.11685705 0.10732204
#> t74  0.2387257 0.1367019 0.1373288 0.1332449 0.11058733 0.10806250 0.09924508
#> t31  0.2150185 0.1231264 0.1236910 0.1200127 0.09960519 0.09733109 0.08938931
#> t77  0.2128187 0.1218667 0.1224256 0.1187849 0.09858618 0.09633535 0.08847481
#> t3   0.2145863 0.1228789 0.1234424 0.1197715 0.09940500 0.09713547 0.08920965
#> t17  0.2194999 0.1256926 0.1262690 0.1225140 0.10168116 0.09935966 0.09125236
#> t63  0.2544057 0.1456807 0.1463488 0.1419967 0.11785094 0.11516027 0.10576371
#> t73  0.2185367 0.1251410 0.1257149 0.1219764 0.10123496 0.09892365 0.09085193
#> t41  0.2423002 0.1387488 0.1393850 0.1352400 0.11224319 0.10968055 0.10073111
#> t48  0.2395400 0.1371682 0.1377972 0.1336994 0.11096454 0.10843110 0.09958360
#> t99  0.2090813 0.1197266 0.1202756 0.1166989 0.09685486 0.09464355 0.08692106
#> t60  0.2247344 0.1286900 0.1292801 0.1254356 0.10410598 0.10172912 0.09342848
#> t43  0.2507317 0.1435769 0.1442353 0.1399460 0.11614897 0.11349715 0.10423629
#> t21  0.2432298 0.1392810 0.1399197 0.1357588 0.11267379 0.11010131 0.10111754
#> t95  0.2520293 0.1443199 0.1449817 0.1406703 0.11675008 0.11408454 0.10477575
#> t59  0.2394513 0.1371173 0.1377461 0.1336499 0.11092343 0.10839092 0.09954671
#> t27  0.2335750 0.1337524 0.1343658 0.1303700 0.10820131 0.10573095 0.09710378
#> t81  0.2197775 0.1258516 0.1264287 0.1226690 0.10180978 0.09948534 0.09136778
#> t70  0.2239402 0.1282353 0.1288233 0.1249924 0.10373811 0.10136965 0.09309834
#> t20  0.2300341 0.1317248 0.1323288 0.1283937 0.10656102 0.10412811 0.09563172
#> t5   0.2314493 0.1325352 0.1331430 0.1291836 0.10721662 0.10476874 0.09622008
#> t12  0.2329183 0.1333764 0.1339880 0.1300035 0.10789712 0.10543370 0.09683078
#> t44  0.3825893 0.2190827 0.2200873 0.2135425 0.17723072 0.17318433 0.15905327
#> t23  0.3162863 0.1811155 0.1819461 0.1765354 0.14651653 0.14317139 0.13148924
#> t45  0.2761583 0.1581370 0.1588621 0.1541379 0.12792761 0.12500687 0.11480687
#> t69  0.2952698 0.1690808 0.1698562 0.1648050 0.13678082 0.13365795 0.12275206
#> t4   1.0000000 0.3066162 0.3080223 0.2988624 0.24804248 0.24237938 0.22260231
#> t96  0.3066162 1.0000000 0.8298474 0.5648920 0.46883518 0.45813113 0.42074968
#> t57  0.3080223 0.8298474 1.0000000 0.5674825 0.47098512 0.46023199 0.42267912
#> t22  0.2988624 0.5648920 0.5674825 1.0000000 0.71142533 0.69518266 0.63845886
#> t55  0.2480425 0.4688352 0.4709851 0.7114253 1.00000000 0.90078125 0.82728153
#> t58  0.2423794 0.4581311 0.4602320 0.6951827 0.90078125 1.00000000 0.90419403
#> t90  0.2226023 0.4207497 0.4226791 0.6384589 0.82728153 0.90419403 1.00000000
#> t28  0.1984455 0.3750900 0.3768100 0.5691734 0.64288648 0.62820863 0.57694962
#> t46  0.1955579 0.3696320 0.3713271 0.5608913 0.63353184 0.61906758 0.56855443
#> t71  0.2086461 0.3943705 0.3961790 0.5984303 0.67593246 0.66050014 0.60660628
#> t89  0.2227117 0.4209564 0.4228868 0.6387726 0.72149940 0.70502673 0.64749970
#> t92  0.2453072 0.4636651 0.4657914 0.7035801 0.79470013 0.77655621 0.71319269
#> t68  0.2595075 0.4905057 0.4927551 0.4821301 0.40014650 0.39101070 0.35910597
#> t7   0.2668346 0.5043549 0.5066677 0.4957428 0.41144442 0.40205066 0.36924513
#> t82  0.2486998 0.4700777 0.4722333 0.4620509 0.38348163 0.37472630 0.34415030
#> t85  0.2698834 0.5101175 0.5124568 0.5014071 0.41614550 0.40664441 0.37346404
#> t64  0.2422645 0.4579140 0.4600138 0.4500949 0.37355868 0.36502991 0.33524510
#> t98  0.2534195 0.4789986 0.4811951 0.4708195 0.39075917 0.38183768 0.35068143
#> t49  0.2396243 0.4529236 0.4550006 0.4451898 0.36948766 0.36105183 0.33159161
#> t8   0.2334499 0.4412532 0.4432767 0.4337187 0.35996714 0.35174867 0.32304755
#> t11  0.2418147 0.4570639 0.4591598 0.4492593 0.37286517 0.36435223 0.33462272
#> t75  0.2503575 0.4732109 0.4753809 0.4651306 0.38603767 0.37722398 0.34644419
#> t79  0.4219522 0.5259734 0.5283853 0.5126724 0.42549519 0.41578065 0.38185480
#> t2   0.3153922 0.3931438 0.3949466 0.3832018 0.31804044 0.31077921 0.28542101
#> t10  0.3397967 0.4235646 0.4255069 0.4128533 0.34264986 0.33482677 0.30750640
#> t62  0.3292953 0.4104744 0.4123567 0.4000942 0.33206032 0.32447900 0.29800296
#> t35  0.3405756 0.4245355 0.4264823 0.4137997 0.34343533 0.33559430 0.30821130
#> t65  0.0000000 0.0000000 0.0000000 0.0000000 0.00000000 0.00000000 0.00000000
#> t54  0.0000000 0.0000000 0.0000000 0.0000000 0.00000000 0.00000000 0.00000000
#> t83  0.0000000 0.0000000 0.0000000 0.0000000 0.00000000 0.00000000 0.00000000
#> t47  0.0000000 0.0000000 0.0000000 0.0000000 0.00000000 0.00000000 0.00000000
#> t93  0.0000000 0.0000000 0.0000000 0.0000000 0.00000000 0.00000000 0.00000000
#> t13  0.0000000 0.0000000 0.0000000 0.0000000 0.00000000 0.00000000 0.00000000
#> t76  0.0000000 0.0000000 0.0000000 0.0000000 0.00000000 0.00000000 0.00000000
#> t9   0.0000000 0.0000000 0.0000000 0.0000000 0.00000000 0.00000000 0.00000000
#> t1   0.0000000 0.0000000 0.0000000 0.0000000 0.00000000 0.00000000 0.00000000
#> t36  0.0000000 0.0000000 0.0000000 0.0000000 0.00000000 0.00000000 0.00000000
#> t56  0.0000000 0.0000000 0.0000000 0.0000000 0.00000000 0.00000000 0.00000000
#> t61  0.0000000 0.0000000 0.0000000 0.0000000 0.00000000 0.00000000 0.00000000
#> t52  0.0000000 0.0000000 0.0000000 0.0000000 0.00000000 0.00000000 0.00000000
#> t80  0.0000000 0.0000000 0.0000000 0.0000000 0.00000000 0.00000000 0.00000000
#> t66  0.0000000 0.0000000 0.0000000 0.0000000 0.00000000 0.00000000 0.00000000
#> t38  0.0000000 0.0000000 0.0000000 0.0000000 0.00000000 0.00000000 0.00000000
#> t30  0.0000000 0.0000000 0.0000000 0.0000000 0.00000000 0.00000000 0.00000000
#>             t28        t46        t71        t89        t92        t68
#> t34  0.08025197 0.07908423 0.08437713 0.09006528 0.09920298 0.10494563
#> t42  0.07670281 0.07558671 0.08064553 0.08608213 0.09481571 0.10030439
#> t50  0.07401313 0.07293616 0.07781759 0.08306354 0.09149087 0.09678708
#> t37  0.07314015 0.07207589 0.07689974 0.08208382 0.09041174 0.09564549
#> t19  0.09293227 0.09158002 0.09770923 0.10429614 0.11487765 0.12152768
#> t84  0.08275571 0.08155153 0.08700956 0.09287517 0.10229795 0.10821977
#> t94  0.08369921 0.08248131 0.08800156 0.09393405 0.10346426 0.10945359
#> t78  0.09656541 0.09516028 0.10152911 0.10837354 0.11936873 0.12627873
#> t6   0.08595055 0.08469988 0.09036862 0.09646068 0.10624723 0.11239766
#> t51  0.07789251 0.07675909 0.08189638 0.08741730 0.09628634 0.10186015
#> t88  0.07688023 0.07576155 0.08083207 0.08628124 0.09503502 0.10053640
#> t40  0.07321055 0.07214527 0.07697376 0.08216282 0.09049877 0.09573755
#> t86  0.07134580 0.07030765 0.07501316 0.08007005 0.08819367 0.09329902
#> t72  0.07100260 0.06996944 0.07465231 0.07968488 0.08776943 0.09285021
#> t14  0.10009654 0.09864003 0.10524175 0.11233646 0.12373371 0.13089640
#> t67  0.10161088 0.10013234 0.10683394 0.11403598 0.12560566 0.13287671
#> t18  0.10073327 0.09926750 0.10591122 0.11305106 0.12452081 0.13172906
#> t15  0.09059794 0.08927965 0.09525490 0.10167636 0.11199208 0.11847507
#> t32  0.07869375 0.07754868 0.08273881 0.08831652 0.09727679 0.10290794
#> t24  0.07631593 0.07520546 0.08023877 0.08564794 0.09433747 0.09979846
#> t39  0.07306513 0.07200196 0.07682086 0.08199962 0.09031901 0.09554738
#> t53  0.07039743 0.06937308 0.07401603 0.07900571 0.08702134 0.09205883
#> t87  0.07482427 0.07373550 0.07867043 0.08397387 0.09249356 0.09784781
#> t26  0.07615678 0.07504862 0.08007143 0.08546932 0.09414073 0.09959034
#> t97  0.07632838 0.07521773 0.08025185 0.08566191 0.09435286 0.09981475
#> t16  0.08392222 0.08270107 0.08823603 0.09418433 0.10373993 0.10974521
#> t91  0.07155209 0.07051094 0.07523005 0.08030157 0.08844867 0.09356878
#> t33  0.07620654 0.07509766 0.08012375 0.08552517 0.09420224 0.09965541
#> t29  0.08107179 0.07989211 0.08523908 0.09098534 0.10021639 0.10601770
#> t25  0.07444508 0.07336183 0.07827174 0.08354831 0.09202482 0.09735194
#> t100 0.09567546 0.09428329 0.10059342 0.10737477 0.11826863 0.12511495
#> t74  0.08847501 0.08718762 0.09302285 0.09929384 0.10936784 0.11569891
#> t31  0.07968879 0.07852924 0.08378499 0.08943323 0.09850680 0.10420915
#> t77  0.07887354 0.07772585 0.08292783 0.08851828 0.09749903 0.10314304
#> t3   0.07952862 0.07837140 0.08361660 0.08925348 0.09830882 0.10399971
#> t17  0.08134966 0.08016594 0.08553124 0.09129719 0.10055988 0.10638108
#> t63  0.09428624 0.09291428 0.09913279 0.10581568 0.11655135 0.12329826
#> t73  0.08099268 0.07981416 0.08515591 0.09089657 0.10011860 0.10591426
#> t41  0.08979978 0.08849310 0.09441571 0.10078060 0.11100544 0.11743131
#> t48  0.08877680 0.08748501 0.09334015 0.09963253 0.10974089 0.11609356
#> t99  0.07748840 0.07636086 0.08147150 0.08696377 0.09578680 0.10133170
#> t60  0.08328963 0.08207768 0.08757093 0.09347438 0.10295796 0.10891798
#> t43  0.09292458 0.09157244 0.09770114 0.10428751 0.11486815 0.12151762
#> t21  0.09014427 0.08883258 0.09477792 0.10116722 0.11143128 0.11788180
#> t95  0.09340550 0.09204636 0.09820678 0.10482724 0.11546263 0.12214652
#> t59  0.08874391 0.08745260 0.09330557 0.09959562 0.10970023 0.11605055
#> t27  0.08656608 0.08530646 0.09101580 0.09715149 0.10700813 0.11320260
#> t81  0.08145256 0.08026734 0.08563943 0.09141268 0.10068708 0.10651564
#> t70  0.08299532 0.08178765 0.08726149 0.09314408 0.10259415 0.10853311
#> t20  0.08525378 0.08401325 0.08963604 0.09567871 0.10538593 0.11148650
#> t5   0.08577829 0.08453013 0.09018751 0.09626736 0.10603430 0.11217240
#> t12  0.08632272 0.08506664 0.09075992 0.09687836 0.10670729 0.11288435
#> t44  0.14179282 0.13972960 0.14908133 0.15913142 0.17527632 0.18542269
#> t23  0.11722004 0.11551437 0.12324545 0.13155385 0.14490083 0.15328883
#> t45  0.10234804 0.10085878 0.10760899 0.11486328 0.12651690 0.13384069
#> t69  0.10943102 0.10783869 0.11505605 0.12281238 0.13527248 0.14310312
#> t4   0.19844553 0.19555795 0.20864613 0.22271168 0.24530721 0.25950752
#> t96  0.37508997 0.36963204 0.39437054 0.42095642 0.46366514 0.49050573
#> t57  0.37681002 0.37132706 0.39617901 0.42288680 0.46579137 0.49275505
#> t22  0.56917338 0.56089134 0.59843032 0.63877257 0.70358014 0.48213013
#> t55  0.64288648 0.63353184 0.67593246 0.72149940 0.79470013 0.40014650
#> t58  0.62820863 0.61906758 0.66050014 0.70502673 0.77655621 0.39101070
#> t90  0.57694962 0.56855443 0.60660628 0.64749970 0.71319269 0.35910597
#> t28  1.00000000 0.90815600 0.85773321 0.83193705 0.74720343 0.32013583
#> t46  0.90815600 1.00000000 0.84525234 0.81983154 0.73633088 0.31547754
#> t71  0.85773321 0.84525234 1.00000000 0.87470072 0.78561157 0.33659162
#> t89  0.83193705 0.81983154 0.87470072 1.00000000 0.83857236 0.35928242
#> t92  0.74720343 0.73633088 0.78561157 0.83857236 1.00000000 0.39573392
#> t68  0.32013583 0.31547754 0.33659162 0.35928242 0.39573392 1.00000000
#> t7   0.32917468 0.32438487 0.34609509 0.36942656 0.40690724 0.79212647
#> t82  0.30680315 0.30233886 0.32257361 0.34431940 0.37925281 0.73829158
#> t85  0.33293577 0.32809122 0.35004950 0.37364755 0.41155648 0.75692751
#> t64  0.29886434 0.29451556 0.31422672 0.33540982 0.36943930 0.57813004
#> t98  0.31262552 0.30807651 0.32869526 0.35085374 0.38645010 0.60475000
#> t49  0.29560733 0.29130595 0.31080229 0.33175454 0.36541316 0.57182961
#> t8   0.28799047 0.28379992 0.30279390 0.32320628 0.35599763 0.55709538
#> t11  0.29830949 0.29396880 0.31364335 0.33478713 0.36875343 0.57705674
#> t75  0.30884810 0.30435406 0.32472367 0.34661442 0.38178067 0.59744287
#> t79  0.34041596 0.33546257 0.35791420 0.38204243 0.42080308 0.44516248
#> t2   0.25444715 0.25074470 0.26752638 0.28556126 0.31453327 0.33274094
#> t10  0.27413584 0.27014689 0.28822711 0.30765750 0.33887132 0.35848786
#> t62  0.26566371 0.26179804 0.27931950 0.29814939 0.32839855 0.34740885
#> t35  0.27476425 0.27076616 0.28888782 0.30836275 0.33964812 0.35930963
#> t65  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t54  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t83  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t47  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t93  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t13  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t76  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t9   0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t1   0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t36  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t56  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t61  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t52  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t80  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t66  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t38  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t30  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#>              t7        t82        t85        t64        t98        t49
#> t34  0.10790871 0.10057496 0.10914165 0.09797249 0.10248363 0.09690479
#> t42  0.10313643 0.09612702 0.10431484 0.09363964 0.09795127 0.09261916
#> t50  0.09951981 0.09275620 0.10065690 0.09035604 0.09451648 0.08937135
#> t37  0.09834599 0.09166215 0.09946967 0.08929030 0.09340167 0.08831722
#> t19  0.12495894 0.11646642 0.12638670 0.11345275 0.11867667 0.11221634
#> t84  0.11127529 0.10371274 0.11254670 0.10102908 0.10568096 0.09992807
#> t94  0.11254395 0.10489518 0.11382985 0.10218092 0.10688583 0.10106736
#> t78  0.12984414 0.12101961 0.13132771 0.11788812 0.12331626 0.11660338
#> t6   0.11557115 0.10771664 0.11689164 0.10492937 0.10976084 0.10378586
#> t51  0.10473612 0.09761799 0.10593281 0.09509203 0.09947054 0.09405572
#> t88  0.10337498 0.09634936 0.10455613 0.09385623 0.09817784 0.09283339
#> t40  0.09844065 0.09175037 0.09956541 0.08937625 0.09349157 0.08840223
#> t86  0.09593326 0.08941340 0.09702938 0.08709974 0.09111024 0.08615053
#> t72  0.09547179 0.08898328 0.09656263 0.08668076 0.09067197 0.08573612
#> t14  0.13459218 0.12544496 0.13613000 0.12219896 0.12782560 0.12086724
#> t67  0.13662840 0.12734280 0.13818949 0.12404768 0.12975945 0.12269582
#> t18  0.13544835 0.12624295 0.13699596 0.12297629 0.12863873 0.12163610
#> t15  0.12182014 0.11354094 0.12321204 0.11060297 0.11569567 0.10939762
#> t32  0.10581349 0.09862214 0.10702249 0.09607020 0.10049374 0.09502323
#> t24  0.10261622 0.09564216 0.10378869 0.09316733 0.09745722 0.09215200
#> t39  0.09824511 0.09156813 0.09936764 0.08919872 0.09330586 0.08822663
#> t53  0.09465806 0.08822486 0.09573960 0.08594196 0.08989915 0.08500537
#> t87  0.10061049 0.09377275 0.10176005 0.09134629 0.09555233 0.09035081
#> t26  0.10240222 0.09544271 0.10357224 0.09297304 0.09725397 0.09195982
#> t97  0.10263296 0.09565777 0.10380562 0.09318253 0.09747312 0.09216703
#> t16  0.11284381 0.10517466 0.11413314 0.10245317 0.10717062 0.10133664
#> t91  0.09621064 0.08967193 0.09730992 0.08735158 0.09137368 0.08639963
#> t33  0.10246913 0.09550507 0.10363992 0.09303379 0.09731752 0.09201991
#> t29  0.10901105 0.10160239 0.11025659 0.09897333 0.10353055 0.09789472
#> t25  0.10010062 0.09329753 0.10124435 0.09088337 0.09506809 0.08989293
#> t100 0.12864750 0.11990429 0.13011740 0.11680166 0.12217978 0.11552876
#> t74  0.11896561 0.11088041 0.12032488 0.10801127 0.11298465 0.10683417
#> t31  0.10715144 0.09986916 0.10837573 0.09728495 0.10176443 0.09622475
#> t77  0.10605523 0.09884745 0.10726700 0.09628968 0.10072333 0.09524032
#> t3   0.10693608 0.09966844 0.10815791 0.09708942 0.10155990 0.09603135
#> t17  0.10938469 0.10195063 0.11063449 0.09931256 0.10388540 0.09823026
#> t63  0.12677952 0.11816327 0.12822808 0.11510568 0.12040572 0.11385127
#> t73  0.10890469 0.10150325 0.11014901 0.09887676 0.10342954 0.09779921
#> t41  0.12074692 0.11254065 0.12212655 0.10962856 0.11467640 0.10843383
#> t48  0.11937140 0.11125862 0.12073531 0.10837970 0.11337003 0.10719858
#> t99  0.10419274 0.09711154 0.10538322 0.09459869 0.09895448 0.09356776
#> t60  0.11199321 0.10438188 0.11327283 0.10168090 0.10636279 0.10057278
#> t43  0.12494860 0.11645678 0.12637624 0.11344336 0.11866685 0.11220706
#> t21  0.12121013 0.11297239 0.12259505 0.11004912 0.11511633 0.10884981
#> t95  0.12559526 0.11705949 0.12703028 0.11403047 0.11928099 0.11278777
#> t59  0.11932717 0.11121740 0.12069058 0.10833954 0.11332803 0.10715887
#> t27  0.11639881 0.10848806 0.11772876 0.10568083 0.11054690 0.10452913
#> t81  0.10952305 0.10207959 0.11077444 0.09943818 0.10401681 0.09835451
#> t70  0.11159748 0.10401303 0.11287257 0.10132160 0.10598695 0.10021740
#> t20  0.11463426 0.10684343 0.11594404 0.10407875 0.10887105 0.10294451
#> t5   0.11533953 0.10750076 0.11665737 0.10471908 0.10954086 0.10357786
#> t12  0.11607158 0.10818306 0.11739779 0.10538373 0.11023611 0.10423526
#> t44  0.19065800 0.17770041 0.19283642 0.17310224 0.18107272 0.17121578
#> t23  0.15761686 0.14690482 0.15941775 0.14310352 0.14969271 0.14154399
#> t45  0.13761961 0.12826664 0.13919202 0.12494762 0.13070082 0.12358595
#> t69  0.14714355 0.13714331 0.14882478 0.13359460 0.13974595 0.13213869
#> t4   0.26683457 0.24869983 0.26988337 0.24226449 0.25341954 0.23962430
#> t96  0.50435489 0.47007767 0.51011754 0.45791398 0.47899859 0.45292365
#> t57  0.50666771 0.47223331 0.51245679 0.46001384 0.48119514 0.45500063
#> t22  0.49574281 0.46205089 0.50140706 0.45009489 0.47081948 0.44518978
#> t55  0.41144442 0.38348163 0.41614550 0.37355868 0.39075917 0.36948766
#> t58  0.40205066 0.37472630 0.40664441 0.36502991 0.38183768 0.36105183
#> t90  0.36924513 0.34415030 0.37346404 0.33524510 0.35068143 0.33159161
#> t28  0.32917468 0.30680315 0.33293577 0.29886434 0.31262552 0.29560733
#> t46  0.32438487 0.30233886 0.32809122 0.29451556 0.30807651 0.29130595
#> t71  0.34609509 0.32257361 0.35004950 0.31422672 0.32869526 0.31080229
#> t89  0.36942656 0.34431940 0.37364755 0.33540982 0.35085374 0.33175454
#> t92  0.40690724 0.37925281 0.41155648 0.36943930 0.38645010 0.36541316
#> t68  0.79212647 0.73829158 0.75692751 0.57813004 0.60475000 0.57182961
#> t7   1.00000000 0.92620949 0.77829893 0.59445322 0.62182477 0.58797489
#> t82  0.92620949 1.00000000 0.72540379 0.55405269 0.57956401 0.54801465
#> t85  0.77829893 0.72540379 1.00000000 0.60124532 0.62892961 0.59469297
#> t64  0.59445322 0.55405269 0.60124532 1.00000000 0.85739381 0.81072041
#> t98  0.62182477 0.57956401 0.62892961 0.85739381 1.00000000 0.90844188
#> t49  0.58797489 0.54801465 0.59469297 0.81072041 0.90844188 1.00000000
#> t8   0.57282465 0.53389406 0.57936963 0.78983072 0.88503422 0.85850628
#> t11  0.59334962 0.55302409 0.60012911 0.67166077 0.70258734 0.66434104
#> t75  0.61431133 0.57256120 0.62133033 0.69538905 0.72740818 0.68781073
#> t79  0.45773139 0.42662282 0.46296134 0.41558356 0.43471908 0.41105455
#> t2   0.34213569 0.31888329 0.34604486 0.31063190 0.32493492 0.30724664
#> t10  0.36860957 0.34355794 0.37282122 0.33466806 0.35007783 0.33102086
#> t62  0.35721774 0.33294033 0.36129924 0.32432519 0.33925872 0.32079071
#> t35  0.36945454 0.34434548 0.37367585 0.33543523 0.35088032 0.33177967
#> t65  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t54  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t83  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t47  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t93  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t13  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t76  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t9   0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t1   0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t36  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t56  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t61  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t52  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t80  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t66  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t38  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> t30  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#>              t8        t11        t75       t79        t2       t10       t62
#> t34  0.09440787 0.09779061 0.10124533 0.1706387 0.1275456 0.1374148 0.1331680
#> t42  0.09023266 0.09346580 0.09676774 0.1630922 0.1219048 0.1313376 0.1272786
#> t50  0.08706853 0.09018830 0.09337445 0.1573732 0.1176301 0.1267321 0.1228154
#> t37  0.08604157 0.08912454 0.09227311 0.1555170 0.1162426 0.1252373 0.1213669
#> t19  0.10932488 0.11324212 0.11724271 0.1976007 0.1476985 0.1591272 0.1542094
#> t84  0.09735324 0.10084152 0.10440403 0.1759624 0.1315248 0.1417019 0.1373227
#> t94  0.09846317 0.10199122 0.10559434 0.1779686 0.1330243 0.1433175 0.1388883
#> t78  0.11359888 0.11766926 0.12182625 0.2053258 0.1534727 0.1653482 0.1602381
#> t6   0.10111163 0.10473457 0.10843461 0.1827555 0.1366024 0.1471724 0.1426241
#> t51  0.09163221 0.09491549 0.09826865 0.1656218 0.1237956 0.1333747 0.1292528
#> t88  0.09044137 0.09368199 0.09699156 0.1634694 0.1221868 0.1316414 0.1275730
#> t40  0.08612439 0.08921032 0.09236192 0.1556667 0.1163545 0.1253578 0.1214837
#> t86  0.08393071 0.08693804 0.09000937 0.1517017 0.1133909 0.1221648 0.1183894
#> t72  0.08352697 0.08651984 0.08957639 0.1509719 0.1128454 0.1215772 0.1178199
#> t14  0.11775287 0.12197209 0.12628110 0.2128340 0.1590848 0.1713945 0.1660976
#> t67  0.11953434 0.12381739 0.12819158 0.2160539 0.1614916 0.1739875 0.1686104
#> t18  0.11850193 0.12274799 0.12708440 0.2141879 0.1600968 0.1724848 0.1671542
#> t15  0.10657879 0.11039763 0.11429773 0.1926372 0.1439885 0.1551301 0.1503359
#> t32  0.09257478 0.09589185 0.09927949 0.1673255 0.1250690 0.1347467 0.1305824
#> t24  0.08977754 0.09299437 0.09627965 0.1622696 0.1212900 0.1306752 0.1266367
#> t39  0.08595331 0.08903312 0.09217846 0.1553574 0.1161234 0.1251088 0.1212424
#> t53  0.08281505 0.08578241 0.08881291 0.1496851 0.1118836 0.1205410 0.1168157
#> t87  0.08802275 0.09117671 0.09439778 0.1590979 0.1189192 0.1281210 0.1241614
#> t26  0.08959031 0.09280043 0.09607886 0.1619312 0.1210370 0.1304026 0.1263726
#> t97  0.08979218 0.09300954 0.09629536 0.1622961 0.1213097 0.1306965 0.1266573
#> t16  0.09872552 0.10226297 0.10587569 0.1784427 0.1333787 0.1436993 0.1392583
#> t91  0.08417339 0.08718941 0.09026962 0.1521403 0.1137187 0.1225181 0.1187317
#> t33  0.08964885 0.09286107 0.09614164 0.1620370 0.1211161 0.1304878 0.1264551
#> t29  0.09537229 0.09878959 0.10227960 0.1723819 0.1288485 0.1388186 0.1345284
#> t25  0.08757668 0.09071465 0.09391939 0.1582916 0.1183166 0.1274717 0.1235322
#> t100 0.11255195 0.11658482 0.12070350 0.2034335 0.1520583 0.1638243 0.1587614
#> t74  0.10408140 0.10781075 0.11161946 0.1881233 0.1406145 0.1514951 0.1468131
#> t31  0.09374534 0.09710434 0.10053482 0.1694412 0.1266505 0.1364505 0.1322335
#> t77  0.09278628 0.09611092 0.09950630 0.1677078 0.1253548 0.1350545 0.1308807
#> t3   0.09355692 0.09690918 0.10033276 0.1691007 0.1263959 0.1361762 0.1319677
#> t17  0.09569918 0.09912819 0.10263017 0.1729727 0.1292901 0.1392944 0.1349895
#> t63  0.11091768 0.11489199 0.11895087 0.2004796 0.1498504 0.1614456 0.1564561
#> t73  0.09527923 0.09869320 0.10217981 0.1722137 0.1287228 0.1386831 0.1343971
#> t41  0.10563984 0.10942504 0.11329078 0.1909401 0.1427200 0.1537634 0.1490114
#> t48  0.10443641 0.10817849 0.11200020 0.1887650 0.1410942 0.1520118 0.1473139
#> t99  0.09115681 0.09442307 0.09775882 0.1647626 0.1231534 0.1326828 0.1285822
#> t60  0.09798134 0.10149213 0.10507762 0.1770977 0.1323733 0.1426162 0.1382086
#> t43  0.10931584 0.11323275 0.11723301 0.1975843 0.1476863 0.1591140 0.1541966
#> t21  0.10604510 0.10984482 0.11372539 0.1916726 0.1432675 0.1543533 0.1495830
#> t95  0.10988159 0.11381877 0.11783973 0.1986069 0.1484506 0.1599375 0.1549946
#> t59  0.10439772 0.10813841 0.11195870 0.1886950 0.1410419 0.1519555 0.1472593
#> t27  0.10183574 0.10548463 0.10921117 0.1840644 0.1375807 0.1482264 0.1436455
#> t81  0.09582023 0.09925358 0.10275998 0.1731915 0.1294537 0.1394706 0.1351602
#> t70  0.09763512 0.10113350 0.10470632 0.1764719 0.1319056 0.1421122 0.1377203
#> t20  0.10029196 0.10388553 0.10755558 0.1812740 0.1354950 0.1459794 0.1414679
#> t5   0.10090899 0.10452467 0.10821729 0.1823893 0.1363286 0.1468775 0.1423382
#> t12  0.10154945 0.10518808 0.10890414 0.1835469 0.1371939 0.1478097 0.1432417
#> t44  0.16680410 0.17278088 0.17888484 0.3014923 0.2253533 0.2427907 0.2352873
#> t23  0.13789685 0.14283785 0.14788399 0.2492435 0.1862994 0.2007149 0.1945119
#> t45  0.12040153 0.12471565 0.12912158 0.2176213 0.1626631 0.1752497 0.1698337
#> t69  0.12873390 0.13334658 0.13805742 0.2326818 0.1739202 0.1873779 0.1815870
#> t4   0.23344995 0.24181473 0.25035750 0.4219522 0.3153922 0.3397967 0.3292953
#> t96  0.44125325 0.45706386 0.47321090 0.5259734 0.3931438 0.4235646 0.4104744
#> t57  0.44327671 0.45915982 0.47538091 0.5283853 0.3949466 0.4255069 0.4123567
#> t22  0.43371866 0.44925929 0.46513062 0.5126724 0.3832018 0.4128533 0.4000942
#> t55  0.35996714 0.37286517 0.38603767 0.4254952 0.3180404 0.3426499 0.3320603
#> t58  0.35174867 0.36435223 0.37722398 0.4157806 0.3107792 0.3348268 0.3244790
#> t90  0.32304755 0.33462272 0.34644419 0.3818548 0.2854210 0.3075064 0.2980030
#> t28  0.28799047 0.29830949 0.30884810 0.3404160 0.2544472 0.2741358 0.2656637
#> t46  0.28379992 0.29396880 0.30435406 0.3354626 0.2507447 0.2701469 0.2617980
#> t71  0.30279390 0.31364335 0.32472367 0.3579142 0.2675264 0.2882271 0.2793195
#> t89  0.32320628 0.33478713 0.34661442 0.3820424 0.2855613 0.3076575 0.2981494
#> t92  0.35599763 0.36875343 0.38178067 0.4208031 0.3145333 0.3388713 0.3283986
#> t68  0.55709538 0.57705674 0.59744287 0.4451625 0.3327409 0.3584879 0.3474088
#> t7   0.57282465 0.59334962 0.61431133 0.4577314 0.3421357 0.3686096 0.3572177
#> t82  0.53389406 0.55302409 0.57256120 0.4266228 0.3188833 0.3435579 0.3329403
#> t85  0.57936963 0.60012911 0.62133033 0.4629613 0.3460449 0.3728212 0.3612992
#> t64  0.78983072 0.67166077 0.69538905 0.4155836 0.3106319 0.3346681 0.3243252
#> t98  0.88503422 0.70258734 0.72740818 0.4347191 0.3249349 0.3500778 0.3392587
#> t49  0.85850628 0.66434104 0.68781073 0.4110546 0.3072466 0.3310209 0.3207907
#> t8   1.00000000 0.64722309 0.67008804 0.4004630 0.2993299 0.3224915 0.3125250
#> t11  0.64722309 1.00000000 0.89097663 0.4148120 0.3100552 0.3340467 0.3237231
#> t75  0.67008804 0.89097663 1.00000000 0.4294664 0.3210088 0.3458479 0.3351595
#> t79  0.40046298 0.41481203 0.42946641 1.0000000 0.5450939 0.5872723 0.5691227
#> t2   0.29932987 0.31005521 0.32100876 0.5450939 1.0000000 0.8690441 0.6910910
#> t10  0.32249151 0.33404675 0.34584787 0.5872723 0.8690441 1.0000000 0.7445664
#> t62  0.31252496 0.32372308 0.33515950 0.5691227 0.6910910 0.7445664 1.0000000
#> t35  0.32323076 0.33481249 0.34664067 0.5886185 0.7147649 0.7700722 0.7632920
#> t65  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t54  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t83  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t47  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t93  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t13  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t76  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t9   0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t1   0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t36  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t56  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t61  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t52  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t80  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t66  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t38  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t30  0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000
#>            t35       t65       t54       t83       t47       t93       t13
#> t34  0.1377298 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t42  0.1316387 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t50  0.1270226 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t37  0.1255244 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t19  0.1594920 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t84  0.1420268 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t94  0.1436460 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t78  0.1657272 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t6   0.1475098 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t51  0.1336805 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t88  0.1319432 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t40  0.1256452 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t86  0.1224449 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t72  0.1218559 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t14  0.1717874 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t67  0.1743863 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t18  0.1728802 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t15  0.1554857 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t32  0.1350556 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t24  0.1309747 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t39  0.1253956 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t53  0.1208173 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t87  0.1284147 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t26  0.1307016 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t97  0.1309961 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t16  0.1440287 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t91  0.1227989 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t33  0.1307870 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t29  0.1391368 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t25  0.1277639 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t100 0.1641999 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t74  0.1518423 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t31  0.1367633 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t77  0.1353641 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t3   0.1364884 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t17  0.1396137 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t63  0.1618157 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t73  0.1390010 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t41  0.1541159 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t48  0.1523603 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t99  0.1329869 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t60  0.1429431 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t43  0.1594788 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t21  0.1547071 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t95  0.1603041 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t59  0.1523038 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t27  0.1485662 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t81  0.1397903 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t70  0.1424380 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t20  0.1463140 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t5   0.1472142 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t12  0.1481485 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t44  0.2433473 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t23  0.2011750 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t45  0.1756515 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t69  0.1878074 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t4   0.3405756 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t96  0.4245355 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t57  0.4264823 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t22  0.4137997 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t55  0.3434353 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t58  0.3355943 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t90  0.3082113 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t28  0.2747642 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t46  0.2707662 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t71  0.2888878 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t89  0.3083627 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t92  0.3396481 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t68  0.3593096 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t7   0.3694545 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t82  0.3443455 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t85  0.3736758 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t64  0.3354352 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t98  0.3508803 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t49  0.3317797 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t8   0.3232308 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t11  0.3348125 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t75  0.3466407 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t79  0.5886185 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t2   0.7147649 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t10  0.7700722 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t62  0.7632920 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t35  1.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t65  0.0000000 1.0000000 0.9611232 0.6290814 0.6127259 0.5524875 0.5328777
#> t54  0.0000000 0.9611232 1.0000000 0.6227429 0.6065522 0.5469207 0.5275085
#> t83  0.0000000 0.6290814 0.6227429 1.0000000 0.8049402 0.5051510 0.4872214
#> t47  0.0000000 0.6127259 0.6065522 0.8049402 1.0000000 0.4920176 0.4745540
#> t93  0.0000000 0.5524875 0.5469207 0.5051510 0.4920176 1.0000000 0.8750496
#> t13  0.0000000 0.5328777 0.5275085 0.4872214 0.4745540 0.8750496 1.0000000
#> t76  0.0000000 0.5486492 0.5431211 0.5016416 0.4885994 0.7965983 0.7683241
#> t9   0.0000000 0.5366451 0.5312379 0.4906660 0.4779091 0.7791692 0.7515136
#> t1   0.0000000 0.5548061 0.5492160 0.5072710 0.4940824 0.5519206 0.5323309
#> t36  0.0000000 0.5595380 0.5539002 0.5115975 0.4982964 0.5566279 0.5368711
#> t56  0.0000000 0.5373807 0.5319662 0.4913386 0.4785643 0.5345858 0.5156114
#> t61  0.0000000 0.3454229 0.3419425 0.3158275 0.3076163 0.2874348 0.2772327
#> t52  0.0000000 0.3381286 0.3347217 0.3091581 0.3011203 0.2813650 0.2713783
#> t80  0.0000000 0.2881715 0.2852679 0.2634813 0.2566310 0.2397945 0.2312833
#> t66  0.0000000 0.3122174 0.3090715 0.2854670 0.2780451 0.2598037 0.2505823
#> t38  0.0000000 0.3411552 0.3377178 0.3119255 0.3038157 0.2838835 0.2738075
#> t30  0.0000000 0.3600546 0.3564268 0.3292056 0.3206465 0.2996102 0.2889759
#>            t76        t9        t1       t36       t56       t61       t52
#> t34  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t42  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t50  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t37  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t19  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t84  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t94  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t78  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t6   0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t51  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t88  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t40  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t86  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t72  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t14  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t67  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t18  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t15  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t32  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t24  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t39  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t53  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t87  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t26  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t97  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t16  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t91  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t33  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t29  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t25  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t100 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t74  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t31  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t77  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t3   0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t17  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t63  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t73  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t41  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t48  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t99  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t60  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t43  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t21  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t95  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t59  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t27  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t81  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t70  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t20  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t5   0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t12  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t44  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t23  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t45  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t69  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t4   0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t96  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t57  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t22  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t55  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t58  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t90  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t28  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t46  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t71  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t89  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t92  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t68  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t7   0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t82  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t85  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t64  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t98  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t49  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t8   0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t11  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t75  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t79  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t2   0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t10  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t62  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t35  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> t65  0.5486492 0.5366451 0.5548061 0.5595380 0.5373807 0.3454229 0.3381286
#> t54  0.5431211 0.5312379 0.5492160 0.5539002 0.5319662 0.3419425 0.3347217
#> t83  0.5016416 0.4906660 0.5072710 0.5115975 0.4913386 0.3158275 0.3091581
#> t47  0.4885994 0.4779091 0.4940824 0.4982964 0.4785643 0.3076163 0.3011203
#> t93  0.7965983 0.7791692 0.5519206 0.5566279 0.5345858 0.2874348 0.2813650
#> t13  0.7683241 0.7515136 0.5323309 0.5368711 0.5156114 0.2772327 0.2713783
#> t76  1.0000000 0.8093455 0.5480863 0.5527609 0.5308720 0.2854379 0.2794103
#> t9   0.8093455 1.0000000 0.5360944 0.5406668 0.5192568 0.2791927 0.2732970
#> t1   0.5480863 0.5360944 1.0000000 0.7781681 0.7473532 0.2886411 0.2825458
#> t36  0.5527609 0.5406668 0.7781681 1.0000000 0.8139063 0.2911029 0.2849556
#> t56  0.5308720 0.5192568 0.7473532 0.8139063 1.0000000 0.2795754 0.2736716
#> t61  0.2854379 0.2791927 0.2886411 0.2911029 0.2795754 1.0000000 0.7216418
#> t52  0.2794103 0.2732970 0.2825458 0.2849556 0.2736716 0.7216418 1.0000000
#> t80  0.2381286 0.2329185 0.2408008 0.2428546 0.2332378 0.3868650 0.3786955
#> t66  0.2579987 0.2523539 0.2608940 0.2631191 0.2526998 0.4191462 0.4102950
#> t38  0.2819113 0.2757433 0.2850749 0.2875063 0.2761213 0.4579947 0.4483232
#> t30  0.2975288 0.2910190 0.3008676 0.3034337 0.2914179 0.4833669 0.4731595
#>            t80       t66       t38       t30
#> t34  0.0000000 0.0000000 0.0000000 0.0000000
#> t42  0.0000000 0.0000000 0.0000000 0.0000000
#> t50  0.0000000 0.0000000 0.0000000 0.0000000
#> t37  0.0000000 0.0000000 0.0000000 0.0000000
#> t19  0.0000000 0.0000000 0.0000000 0.0000000
#> t84  0.0000000 0.0000000 0.0000000 0.0000000
#> t94  0.0000000 0.0000000 0.0000000 0.0000000
#> t78  0.0000000 0.0000000 0.0000000 0.0000000
#> t6   0.0000000 0.0000000 0.0000000 0.0000000
#> t51  0.0000000 0.0000000 0.0000000 0.0000000
#> t88  0.0000000 0.0000000 0.0000000 0.0000000
#> t40  0.0000000 0.0000000 0.0000000 0.0000000
#> t86  0.0000000 0.0000000 0.0000000 0.0000000
#> t72  0.0000000 0.0000000 0.0000000 0.0000000
#> t14  0.0000000 0.0000000 0.0000000 0.0000000
#> t67  0.0000000 0.0000000 0.0000000 0.0000000
#> t18  0.0000000 0.0000000 0.0000000 0.0000000
#> t15  0.0000000 0.0000000 0.0000000 0.0000000
#> t32  0.0000000 0.0000000 0.0000000 0.0000000
#> t24  0.0000000 0.0000000 0.0000000 0.0000000
#> t39  0.0000000 0.0000000 0.0000000 0.0000000
#> t53  0.0000000 0.0000000 0.0000000 0.0000000
#> t87  0.0000000 0.0000000 0.0000000 0.0000000
#> t26  0.0000000 0.0000000 0.0000000 0.0000000
#> t97  0.0000000 0.0000000 0.0000000 0.0000000
#> t16  0.0000000 0.0000000 0.0000000 0.0000000
#> t91  0.0000000 0.0000000 0.0000000 0.0000000
#> t33  0.0000000 0.0000000 0.0000000 0.0000000
#> t29  0.0000000 0.0000000 0.0000000 0.0000000
#> t25  0.0000000 0.0000000 0.0000000 0.0000000
#> t100 0.0000000 0.0000000 0.0000000 0.0000000
#> t74  0.0000000 0.0000000 0.0000000 0.0000000
#> t31  0.0000000 0.0000000 0.0000000 0.0000000
#> t77  0.0000000 0.0000000 0.0000000 0.0000000
#> t3   0.0000000 0.0000000 0.0000000 0.0000000
#> t17  0.0000000 0.0000000 0.0000000 0.0000000
#> t63  0.0000000 0.0000000 0.0000000 0.0000000
#> t73  0.0000000 0.0000000 0.0000000 0.0000000
#> t41  0.0000000 0.0000000 0.0000000 0.0000000
#> t48  0.0000000 0.0000000 0.0000000 0.0000000
#> t99  0.0000000 0.0000000 0.0000000 0.0000000
#> t60  0.0000000 0.0000000 0.0000000 0.0000000
#> t43  0.0000000 0.0000000 0.0000000 0.0000000
#> t21  0.0000000 0.0000000 0.0000000 0.0000000
#> t95  0.0000000 0.0000000 0.0000000 0.0000000
#> t59  0.0000000 0.0000000 0.0000000 0.0000000
#> t27  0.0000000 0.0000000 0.0000000 0.0000000
#> t81  0.0000000 0.0000000 0.0000000 0.0000000
#> t70  0.0000000 0.0000000 0.0000000 0.0000000
#> t20  0.0000000 0.0000000 0.0000000 0.0000000
#> t5   0.0000000 0.0000000 0.0000000 0.0000000
#> t12  0.0000000 0.0000000 0.0000000 0.0000000
#> t44  0.0000000 0.0000000 0.0000000 0.0000000
#> t23  0.0000000 0.0000000 0.0000000 0.0000000
#> t45  0.0000000 0.0000000 0.0000000 0.0000000
#> t69  0.0000000 0.0000000 0.0000000 0.0000000
#> t4   0.0000000 0.0000000 0.0000000 0.0000000
#> t96  0.0000000 0.0000000 0.0000000 0.0000000
#> t57  0.0000000 0.0000000 0.0000000 0.0000000
#> t22  0.0000000 0.0000000 0.0000000 0.0000000
#> t55  0.0000000 0.0000000 0.0000000 0.0000000
#> t58  0.0000000 0.0000000 0.0000000 0.0000000
#> t90  0.0000000 0.0000000 0.0000000 0.0000000
#> t28  0.0000000 0.0000000 0.0000000 0.0000000
#> t46  0.0000000 0.0000000 0.0000000 0.0000000
#> t71  0.0000000 0.0000000 0.0000000 0.0000000
#> t89  0.0000000 0.0000000 0.0000000 0.0000000
#> t92  0.0000000 0.0000000 0.0000000 0.0000000
#> t68  0.0000000 0.0000000 0.0000000 0.0000000
#> t7   0.0000000 0.0000000 0.0000000 0.0000000
#> t82  0.0000000 0.0000000 0.0000000 0.0000000
#> t85  0.0000000 0.0000000 0.0000000 0.0000000
#> t64  0.0000000 0.0000000 0.0000000 0.0000000
#> t98  0.0000000 0.0000000 0.0000000 0.0000000
#> t49  0.0000000 0.0000000 0.0000000 0.0000000
#> t8   0.0000000 0.0000000 0.0000000 0.0000000
#> t11  0.0000000 0.0000000 0.0000000 0.0000000
#> t75  0.0000000 0.0000000 0.0000000 0.0000000
#> t79  0.0000000 0.0000000 0.0000000 0.0000000
#> t2   0.0000000 0.0000000 0.0000000 0.0000000
#> t10  0.0000000 0.0000000 0.0000000 0.0000000
#> t62  0.0000000 0.0000000 0.0000000 0.0000000
#> t35  0.0000000 0.0000000 0.0000000 0.0000000
#> t65  0.2881715 0.3122174 0.3411552 0.3600546
#> t54  0.2852679 0.3090715 0.3377178 0.3564268
#> t83  0.2634813 0.2854670 0.3119255 0.3292056
#> t47  0.2566310 0.2780451 0.3038157 0.3206465
#> t93  0.2397945 0.2598037 0.2838835 0.2996102
#> t13  0.2312833 0.2505823 0.2738075 0.2889759
#> t76  0.2381286 0.2579987 0.2819113 0.2975288
#> t9   0.2329185 0.2523539 0.2757433 0.2910190
#> t1   0.2408008 0.2608940 0.2850749 0.3008676
#> t36  0.2428546 0.2631191 0.2875063 0.3034337
#> t56  0.2332378 0.2526998 0.2761213 0.2914179
#> t61  0.3868650 0.4191462 0.4579947 0.4833669
#> t52  0.3786955 0.4102950 0.4483232 0.4731595
#> t80  1.0000000 0.9077854 0.7045673 0.7082639
#> t66  0.9077854 1.0000000 0.7633585 0.7673635
#> t38  0.7045673 0.7633585 1.0000000 0.8384865
#> t30  0.7082639 0.7673635 0.8384865 1.0000000
#> 
#> $pglmm_fit
#>                variable    mean  median   sd  mad      q5     q95 rhat ess_bulk
#>  lp__                   -290.28 -290.60 5.59 7.08 -299.87 -282.39 1.00       20
#>  centered_cov_intercept    0.07    0.04 0.16 0.14   -0.19    0.44 1.11       20
#>  sigma_resid               1.04    1.03 0.08 0.09    0.93    1.17 1.08       20
#>  sigma_phylo               0.15    0.13 0.12 0.12    0.01    0.35 1.21       20
#>  std_phylo_effects[1]     -0.01    0.09 1.17 1.32   -1.49    1.54 1.14       20
#>  std_phylo_effects[2]     -0.25   -0.25 1.04 0.85   -1.62    1.13 1.14       20
#>  std_phylo_effects[3]      0.08    0.06 0.78 0.69   -1.32    1.21 1.08       20
#>  std_phylo_effects[4]     -0.16   -0.42 1.16 1.29   -2.00    1.64 1.12       20
#>  std_phylo_effects[5]      0.00   -0.12 0.82 0.86   -1.18    1.31 1.01       20
#>  std_phylo_effects[6]      0.06    0.15 0.89 0.91   -1.30    1.35 1.26       20
#>  ess_tail
#>        20
#>        20
#>        20
#>        20
#>        20
#>        20
#>        20
#>        20
#>        20
#>        20
#> 
#>  # showing 10 of 405 rows (change via 'max_rows' argument or 'cmdstanr_max_rows' option)
#> 
#> $base_fit
#>                variable    mean  median   sd  mad      q5     q95 rhat ess_bulk
#>  lp__                   -148.26 -147.92 1.35 0.79 -150.93 -147.26 1.13       20
#>  centered_cov_intercept    0.09    0.07 0.11 0.09   -0.08    0.30 1.30       20
#>  sigma_resid               1.05    1.05 0.08 0.07    0.95    1.15 1.04       20
#>  intercept                 0.09    0.07 0.11 0.09   -0.08    0.30 1.30       20
#>  log_lik[1]               -0.97   -0.97 0.07 0.07   -1.06   -0.88 1.10       20
#>  log_lik[2]               -0.98   -1.00 0.07 0.08   -1.07   -0.89 1.11       20
#>  log_lik[3]               -1.25   -1.25 0.08 0.08   -1.38   -1.13 1.18       20
#>  log_lik[4]               -1.03   -1.03 0.08 0.06   -1.13   -0.96 1.08       20
#>  log_lik[5]               -1.43   -1.40 0.10 0.09   -1.66   -1.30 1.34       20
#>  log_lik[6]               -1.16   -1.14 0.08 0.05   -1.28   -1.07 1.17       20
#>  ess_tail
#>        20
#>        20
#>        20
#>        20
#>        20
#>        20
#>        20
#>        20
#>        20
#>        20
#> 
#>  # showing 10 of 204 rows (change via 'max_rows' argument or 'cmdstanr_max_rows' option)
#> 
#> $loo
#> $loo$pglmm_loo
#> 
#> Computed from 40 by 100 log-likelihood matrix.
#> 
#>          Estimate   SE
#> elpd_loo   -147.7  7.9
#> p_loo         3.3  0.6
#> looic       295.5 15.7
#> ------
#> MCSE of elpd_loo is NA.
#> MCSE and ESS estimates assume MCMC draws (r_eff in [0.6, 1.6]).
#> 
#> Pareto k diagnostic values:
#>                           Count Pct.    Min. ESS
#> (-Inf, 0.38]   (good)     67    67.0%   26      
#>    (0.38, 1]   (bad)      33    33.0%   <NA>    
#>     (1, Inf)   (very bad)  0     0.0%   <NA>    
#> See help('pareto-k-diagnostic') for details.
#> 
#> $loo$pglmm_ll_df
#> # A tibble: 40 × 100
#>    `log_lik[1]` `log_lik[2]` `log_lik[3]` `log_lik[4]` `log_lik[5]` `log_lik[6]`
#>           <dbl>        <dbl>        <dbl>        <dbl>        <dbl>        <dbl>
#>  1       -0.899       -0.900        -1.25       -0.949        -1.37        -1.08
#>  2       -0.912       -0.900        -1.32       -0.929        -1.30        -1.03
#>  3       -1.10        -1.09         -1.41       -1.10         -1.33        -1.16
#>  4       -1.11        -1.07         -1.60       -1.04         -1.18        -1.07
#>  5       -0.967       -0.945        -1.39       -0.956        -1.25        -1.03
#>  6       -0.888       -0.886        -1.26       -0.930        -1.35        -1.05
#>  7       -0.869       -0.888        -1.16       -0.963        -1.48        -1.12
#>  8       -0.997       -1.000        -1.27       -1.05         -1.38        -1.16
#>  9       -1.08        -1.15         -1.13       -1.22         -1.81        -1.42
#> 10       -1.04        -1.05         -1.27       -1.10         -1.45        -1.22
#> # ℹ 30 more rows
#> # ℹ 94 more variables: `log_lik[7]` <dbl>, `log_lik[8]` <dbl>,
#> #   `log_lik[9]` <dbl>, `log_lik[10]` <dbl>, `log_lik[11]` <dbl>,
#> #   `log_lik[12]` <dbl>, `log_lik[13]` <dbl>, `log_lik[14]` <dbl>,
#> #   `log_lik[15]` <dbl>, `log_lik[16]` <dbl>, `log_lik[17]` <dbl>,
#> #   `log_lik[18]` <dbl>, `log_lik[19]` <dbl>, `log_lik[20]` <dbl>,
#> #   `log_lik[21]` <dbl>, `log_lik[22]` <dbl>, `log_lik[23]` <dbl>, …
#> 
#> $loo$base_loo
#> 
#> Computed from 40 by 100 log-likelihood matrix.
#> 
#>          Estimate   SE
#> elpd_loo   -147.6  7.8
#> p_loo         2.1  0.4
#> looic       295.1 15.7
#> ------
#> MCSE of elpd_loo is NA.
#> MCSE and ESS estimates assume MCMC draws (r_eff in [0.2, 1.6]).
#> 
#> Pareto k diagnostic values:
#>                           Count Pct.    Min. ESS
#> (-Inf, 0.38]   (good)     38    38.0%   9       
#>    (0.38, 1]   (bad)      62    62.0%   <NA>    
#>     (1, Inf)   (very bad)  0     0.0%   <NA>    
#> See help('pareto-k-diagnostic') for details.
#> 
#> $loo$comparison
#>           elpd_diff se_diff
#> base_fit   0.0       0.0   
#> pglmm_fit -0.2       0.5   
#> 
#> 
```
