
<!-- README.md is generated from README.Rmd. Please edit that file -->

# aphylo2

<!-- badges: start -->
<!-- badges: end -->

The goal of aphylo2 is to …

## Installation

You can install the released version of aphylo2 from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("aphylo2")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("USCbiostats/aphylo2")
```

## Example

``` r
library(aphylo2)

# Preparing data
annotations <- list(
  c(0, 0, 1, 9, 9, 9),
  c(1, 1, 9, 0, 9, 9),
  c(1, 9, 0, 1, 9, 9)
)

geneid   <- c(0, 1, 2, 3, 4, 5)
parentid <- c(4, 4, 5, 5, 6, 6)

# Building the model
amodel <- new_model(annotations, geneid, parentid)

invisible({
  
  term_gains(amodel, 0:2)
  term_loss(amodel, 0:2)
  term_cogain(amodel, 0, 1)
  term_cogain(amodel, 0, 2)
  term_cogain(amodel, 1, 2)
  
  init(amodel)
  
})

# Testing
params <- c(
  # Gains
  .1, .1, .1,
  # Loss
  .1, .1, .1,
  # Co-gain
  .1, .1, .1,
  # Root probabilities
  .1, .1, .1
)

# Testing the likelihood
likelihood(amodel, params*0)
#> [1] 0.001953125

# Finding MLE
ans <- optim(params*0, function(f) {
  log(likelihood(amodel, f))
}, control = list(fnscale=-1, maxit = 2e3),
hessian = TRUE)

ans
#> $par
#>  [1]   29.13627   29.13600   29.13633   13.47175  259.89993   14.01999
#>  [7] -185.78661 -166.61386  150.07648  -95.66314  -60.97847  174.44247
#> 
#> $value
#> [1] -2.772593
#> 
#> $counts
#> function gradient 
#>     1365       NA 
#> 
#> $convergence
#> [1] 0
#> 
#> $message
#> NULL
#> 
#> $hessian
#>                [,1]          [,2]          [,3]          [,4]          [,5]
#>  [1,] -9.999998e-01  4.999999e-01  4.999999e-01  1.110223e-10  3.552714e-09
#>  [2,]  4.999999e-01 -4.999999e-01  0.000000e+00  0.000000e+00  0.000000e+00
#>  [3,]  4.999999e-01  0.000000e+00 -4.999999e-01  0.000000e+00  7.105427e-09
#>  [4,]  1.110223e-10  0.000000e+00  0.000000e+00 -2.820522e-06  0.000000e+00
#>  [5,]  3.552714e-09  0.000000e+00  7.105427e-09  0.000000e+00 -2.220446e-10
#>  [6,] -1.665335e-09 -2.664535e-09 -5.551115e-11  0.000000e+00  0.000000e+00
#>  [7,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
#>  [8,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
#>  [9,]  0.000000e+00  1.110223e-10  0.000000e+00  0.000000e+00  0.000000e+00
#> [10,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
#> [11,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
#> [12,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
#>                [,6] [,7] [,8]         [,9] [,10] [,11] [,12]
#>  [1,] -1.665335e-09    0    0 0.000000e+00     0     0     0
#>  [2,] -2.664535e-09    0    0 1.110223e-10     0     0     0
#>  [3,] -5.551115e-11    0    0 0.000000e+00     0     0     0
#>  [4,]  0.000000e+00    0    0 0.000000e+00     0     0     0
#>  [5,]  0.000000e+00    0    0 0.000000e+00     0     0     0
#>  [6,] -1.632916e-06    0    0 0.000000e+00     0     0     0
#>  [7,]  0.000000e+00    0    0 0.000000e+00     0     0     0
#>  [8,]  0.000000e+00    0    0 0.000000e+00     0     0     0
#>  [9,]  0.000000e+00    0    0 0.000000e+00     0     0     0
#> [10,]  0.000000e+00    0    0 0.000000e+00     0     0     0
#> [11,]  0.000000e+00    0    0 0.000000e+00     0     0     0
#> [12,]  0.000000e+00    0    0 0.000000e+00     0     0     0

# Root node probabilities
plogis(tail(ans$par, 3))
#> [1] 2.844634e-42 3.291434e-27 1.000000e+00

# Is it invertible?
diag(MASS::ginv(-ans$hessian))
#>  [1] -12999791.087 -12999789.216 -12999791.176    354544.332  -2496622.545
#>  [6]    612307.425         0.000         0.000      -275.916         0.000
#> [11]         0.000         0.000
# diag(solve(-ans$hessian, tol = 1e-100))
```

## Code of Conduct

Please note that the aphylo2 project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
