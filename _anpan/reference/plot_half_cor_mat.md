# Plot a rotated half correlation matrix

Plot the lower triangle of a correlation matrix

## Usage

``` r
plot_half_cor_mat(cor_mat, border_width = 0)
```

## Arguments

- cor_mat:

  a correlation matrix (must have dimnames)

- border_width:

  width parameter of the border around each cell.

## Value

a ggplot of the lower triangle of the matrix.

## Details

If you see a thin, pixel-width grey border around each cell, try setting
border_width = 0.1 or so, depending on your output resolution.

No checks are made on the order of the columns. If you want the order to
line up with another plot, you'll need to check the input manually
beforehand.
