
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
n <- 100L
annotations <- replicate(n * 2 - 1, c(9, 9), simplify = FALSE)

# Random tree
set.seed(31)
tree <- aphylo::sim_tree(n)$edge - 1L

duplication <- rep(TRUE, n * 2 - 1)

# Reading the data in
amodel <- new_model(
  annotations = annotations,
  geneid = c(tree[, 2], n),
  parent = c(tree[, 1], -1),
  duplication = duplication
  )

# Preparing the model
invisible({
  
  term_gains(amodel, 0:1)
  term_loss(amodel, 0:1)
  term_maxfuns(amodel, 1, 1)
  
  init(amodel)
  
})

# Testing
params <- c(
  # Gains
  2, 1,
  # Loss
  -2, -1,
  # Max funs
  2, 
  # Root probabilities
  -10, -10
)

likelihood(amodel, params*0) # Equals 1 b/c all missings
#> [1] 1

# Simulating data
fake <- sim_aphylo2(p = amodel, par = params, seed = 1)

amodel <- new_model(
  annotations = fake,
  geneid = tree[, 2],
  parent = tree[, 1],
  duplication = duplication
  )

# Fitting the model

invisible({
  
  term_gains(amodel, 0:1)
  term_loss(amodel, 0:1)
  term_maxfuns(amodel, 1, 1)
  
  init(amodel)
  
})


# Finding MLE
ans <- aphylo2_mle(amodel, hessian = TRUE)
ans
#> $par
#> [1]  1.6611797  0.5034166 -0.4629501 -1.1796831  1.5409651 -0.2501636 -1.0789896
#> 
#> $value
#> [1] -117.3298
#> 
#> $counts
#> function gradient 
#>      284       NA 
#> 
#> $convergence
#> [1] 0
#> 
#> $message
#> NULL
#> 
#> $hessian
#>               [,1]         [,2] [,3]          [,4]          [,5]          [,6]
#> [1,] -3.2491893549  0.828509251    0  -2.176372295   1.403106516 -5.369882e-03
#> [2,]  0.8285092505 -4.856673442    0   4.991316501   1.944071759  7.987918e-03
#> [3,]  0.0000000000  0.000000000    0   0.000000000   0.000000000  0.000000e+00
#> [4,] -2.1763722948  4.991316501    0 -10.512375216  -7.317389530  1.401741e-03
#> [5,]  1.4031065163  1.944071759    0  -7.317389530 -16.171026651  3.249731e-03
#> [6,] -0.0053698823  0.007987918    0   0.001401741   0.003249731  4.711822e-04
#> [7,] -0.0009945342  0.002276106    0   0.001637474   0.001204327  6.339107e-05
#>               [,7]
#> [1,] -9.945342e-04
#> [2,]  2.276106e-03
#> [3,]  0.000000e+00
#> [4,]  1.637474e-03
#> [5,]  1.204327e-03
#> [6,]  6.339107e-05
#> [7,]  8.769163e-05

# Root node probabilities and odds ratios
plogis(tail(ans$par, 2))
#> [1] 0.4377832 0.2536973
exp(ans$par[1:5])
#> [1] 5.2655191 1.6543639 0.6294241 0.3073761 4.6690941

# Is it invertible?
diag(MASS::ginv(-ans$hessian))
#> [1]  5.033246e-01  4.258817e-01  6.320748e-26  3.972670e-01  1.418275e-01
#> [6] -2.184622e+03 -1.227846e+04
# diag(solve(-ans$hessian, tol = 1e-100))
```

``` r
set.seed(122)
ans_mcmc <- aphylo2_mcmc(
  amodel,
  initial = ans$par * 0,
  nsteps  = 10000,
  kernel  = fmcmc::kernel_ram(warmup = 2000)
  )
plot(window(ans_mcmc, start = 5000))
```

<img src="man/figures/README-mcmc-fit-1.png" width="100%" /><img src="man/figures/README-mcmc-fit-2.png" width="100%" />

``` r
summary(window(ans_mcmc, start = 5000))
#> 
#> Iterations = 5000:10000
#> Thinning interval = 1 
#> Number of chains = 1 
#> Sample size per chain = 5001 
#> 
#> 1. Empirical mean and standard deviation for each variable,
#>    plus standard error of the mean:
#> 
#>         Mean     SD Naive SE Time-series SE
#> par1  1.2080 0.8304 0.011743        0.17984
#> par2  0.6798 0.6630 0.009375        0.10187
#> par3  0.1542 1.6571 0.023433        0.45838
#> par4 -0.6455 0.5817 0.008225        0.11440
#> par5  1.2820 0.2937 0.004153        0.03244
#> par6  0.4064 1.3183 0.018641        0.25909
#> par7  0.8116 1.1591 0.016390        0.21856
#> 
#> 2. Quantiles for each variable:
#> 
#>           2.5%     25%      50%     75%  97.5%
#> par1 -0.007838  0.6083  1.12367  1.6672 3.3777
#> par2 -0.354193  0.2196  0.60095  1.0356 2.2473
#> par3 -3.437999 -0.6417  0.04904  1.1074 3.1886
#> par4 -1.855627 -1.0481 -0.62175 -0.2048 0.3188
#> par5  0.806523  1.0701  1.25124  1.4574 1.9751
#> par6 -2.645390 -0.3673  0.58166  1.2277 2.9185
#> par7 -2.040171  0.4205  0.86344  1.4287 3.0142
```

## Code of Conduct

Please note that the aphylo2 project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
