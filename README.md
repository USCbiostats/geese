
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
  c(1, 1, 0),
  c(9, 1, 0),
  c(0, 9, 1),
  c(1, 0, 9),
  c(9, 9, 9),
  c(9, 9, 9),
  c(9, 9, 9)
)

geneid   <- c(0, 1, 2, 3, 4, 5, 6)
parentid <- c(4, 4, 5, 5, 6, 6, -1)
duplication <- rep(TRUE, 7)

# Building the model
amodel <- new_model(annotations, geneid, parentid, duplication)

invisible({
  
  term_cogain(amodel, 0, 1)
  term_cogain(amodel, 0, 2)
  term_cogain(amodel, 1, 2)
  term_maxfuns(amodel, 2, 2)
  
  init(amodel)
  
})

# Testing
params <- c(
  # Cogain
  .1, .1, .1,
  # max funs
  .1, 
  # Root probabilities
  .1, .1, .1
)

# Testing the likelihood
likelihood(amodel, params*0)
#> [1] 9.106756e-312

# Finding MLE
ans <- optim(params*0, function(f) {
  log(likelihood(amodel, f))
}, control = list(fnscale=-1, maxit = 2e3),
hessian = TRUE)

ans
#> $par
#> [1]  40.58158 -44.56906 -21.95080  87.28769 -30.44540 -51.66042  21.95585
#> 
#> $value
#> [1] -714.459
#> 
#> $counts
#> function gradient 
#>      116       NA 
#> 
#> $convergence
#> [1] 0
#> 
#> $message
#> NULL
#> 
#> $hessian
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7]
#> [1,]    0    0    0    0    0    0    0
#> [2,]    0    0    0    0    0    0    0
#> [3,]    0    0    0    0    0    0    0
#> [4,]    0    0    0    0    0    0    0
#> [5,]    0    0    0    0    0    0    0
#> [6,]    0    0    0    0    0    0    0
#> [7,]    0    0    0    0    0    0    0

# Root node probabilities
plogis(tail(ans$par, 3))
#> [1] 5.994166e-14 3.665780e-23 1.000000e+00

# Is it invertible?
diag(MASS::ginv(-ans$hessian))
#> [1] 0 0 0 0 0 0 0
# diag(solve(-ans$hessian, tol = 1e-100))

# Simulating
ans <- sim_aphylo2(amodel, c(10, 10, -10, 10, -100, -100, -100)/5)
do.call(rbind, ans)
#>      [,1] [,2] [,3]
#> [1,]    0    1    1
#> [2,]    1    0    1
#> [3,]    1    1    0
#> [4,]    1    1    0
#> [5,]    1    1    0
#> [6,]    1    1    0
#> [7,]    0    0    0
```

## Code of Conduct

Please note that the aphylo2 project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
