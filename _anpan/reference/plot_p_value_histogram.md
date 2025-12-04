# Make a p-value histogram

This function makes a p-value histogram from a collection of bug:gene
glm fits.

## Usage

``` r
plot_p_value_histogram(
  all_bug_terms,
  out_dir = NULL,
  plot_ext = "pdf",
  n_bins = 50
)
```

## Arguments

- all_bug_terms:

  a data frame of bug:gene glm fits

- out_dir:

  string giving the output directory

- plot_ext:

  string giving the extension to use

## Details

The plot will be written out to `p_value_histogram.<ext>` in the
specified output directory. The "aaa" is there for alphabetical
superiority.

If you don't understand the purpose of this type of plot, [this blog
post by David
Robinson](http://varianceexplained.org/statistics/interpreting-pvalue-histogram/)
has a lot of helpful information.
