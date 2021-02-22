
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

annotations <- list(
  c(0, 0, 1, 9, 9, 9),
  c(1, 1, 9, 0, 9, 9),
  c(1, 9, 0, 1, 9, 9)
)

geneid   <- c(0, 1, 2, 3, 4, 5)
parentid <- c(4, 4, 5, 5, 6, 6)

amodel <- new_model(annotations, geneid, parentid)

term_gains(amodel, 0:2)
#> [1] 0
term_loss(amodel, 0:2)
#> [1] 0
term_cogain(amodel, 0, 1)
#> [1] 0
term_cogain(amodel, 0, 2)
#> [1] 0
term_cogain(amodel, 1, 2)
#> [1] 0

init(amodel)
#> [1] 0

params <- c(
  # Main parameters
  .1, .1, .1, .1, .1, .1, .1, .1, .1,
  # Root probabilities
  .1, .1, .1
)

likelihood(amodel, params)
#> [1] 0.002040334
```

## Code of Conduct

Please note that the aphylo2 project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
