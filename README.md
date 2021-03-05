
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
  2, 1.5,
  # Loss
  -2, -1.5,
  # Max funs
  2, 
  # Root probabilities
  -10, -10
)
names(params) <- c("gain0", "gain1", "loss0", "loss1", "onefun", "root0", "root1")

likelihood(amodel, params*0) # Equals 1 b/c all missings
#> [1] 1

# Simulating data
fake <- sim_aphylo2(p = amodel, par = params, seed = 1110)
```

``` r
library(aphylo)
#> Loading required package: ape
ap <- new_aphylo(
  tree           = {set.seed(31);aphylo::sim_tree(n)},
  tip.annotation = as.data.frame(do.call(rbind, fake[1:n]))
)
plot(ap)
```

<img src="man/figures/README-viz-with-aphylo-1.png" width="100%" />

``` r
# Fitting the model
fake[101:199] <- replicate(99, c(9L,9L), simplify = FALSE)
amodel <- new_model(
  annotations = fake,
  geneid = tree[, 2],
  parent = tree[, 1],
  duplication = duplication
  )

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
#> [1]   5.885176   8.528545  -7.341997  -5.870773   7.153870 -13.242267   7.963540
#> 
#> $value
#> [1] -81.56194
#> 
#> $counts
#> function gradient 
#>      914       NA 
#> 
#> $convergence
#> [1] 0
#> 
#> $message
#> NULL
#> 
#> $hessian
#>               [,1]          [,2]          [,3]          [,4]          [,5]
#> [1,] -3.522599e+00  3.424469e-01  4.736160e+00 -3.355075e+00  4.561186e+00
#> [2,]  3.424469e-01 -8.947986e-01 -2.737647e-01 -5.364908e-02  2.249417e-01
#> [3,]  4.736160e+00 -2.737647e-01 -9.233159e+00  4.050577e+00 -9.644935e+00
#> [4,] -3.355075e+00 -5.364908e-02  4.050577e+00 -5.799351e+00  1.659960e+00
#> [5,]  4.561186e+00  2.249417e-01 -9.644935e+00  1.659960e+00 -1.277114e+01
#> [6,]  1.598721e-08 -9.947598e-08  6.750156e-08 -2.078338e-07 -5.329071e-08
#> [7,] -4.645173e-05 -5.136691e-05  1.401013e-05 -1.172520e-04 -5.419665e-06
#>               [,6]          [,7]
#> [1,]  1.598721e-08 -4.645173e-05
#> [2,] -9.947598e-08 -5.136691e-05
#> [3,]  6.750156e-08  1.401013e-05
#> [4,] -2.078338e-07 -1.172520e-04
#> [5,] -5.329071e-08 -5.419665e-06
#> [6,] -1.172396e-07  0.000000e+00
#> [7,]  0.000000e+00  5.663026e-05
# [1]  0.6881464  1.0305919 -0.9295311 -1.1126288  1.4606278 -0.1944617 -1.0264741
# Root node probabilities and odds ratios
plogis(tail(ans$par, 2))
#> [1] 1.774009e-06 9.996522e-01
exp(ans$par[1:5])
#> [1] 3.596660e+02 5.057081e+03 6.477554e-04 2.820692e-03 1.279046e+03

# Is it invertible?
diag(MASS::ginv(-ans$hessian))
#> [1]  7.275988e+03  7.276315e+03  7.275885e+03  7.275309e+03  7.275289e+03
#> [6]  1.223175e-04 -1.765650e+04
# diag(solve(-ans$hessian, tol = 1e-100))
```

``` r
set.seed(122)
ans_mcmc <- aphylo2_mcmc(
  amodel,
  nsteps  = 40000,
  kernel  = fmcmc::kernel_ram(warmup = 2000), 
  prior   = function(p) dlogis(p, scale = 2, log = TRUE)
  )
```

<img src="man/figures/README-mcmc-analysis-1.png" width="100%" />

    #> 
    #> Iterations = 15000:40000
    #> Thinning interval = 1 
    #> Number of chains = 1 
    #> Sample size per chain = 25001 
    #> 
    #> 1. Empirical mean and standard deviation for each variable,
    #>    plus standard error of the mean:
    #> 
    #>           Mean     SD Naive SE Time-series SE
    #> par1  3.249064 1.9591 0.012390        0.10719
    #> par2  2.081177 1.9647 0.012426        0.13206
    #> par3 -1.999795 1.0003 0.006326        0.04588
    #> par4 -2.401966 1.3356 0.008447        0.07514
    #> par5  2.562228 0.8483 0.005365        0.03982
    #> par6  0.007906 3.5372 0.022371        0.16822
    #> par7 -0.017456 3.5727 0.022595        0.15159
    #> 
    #> 2. Quantiles for each variable:
    #> 
    #>         2.5%     25%      50%    75%   97.5%
    #> par1  0.1628  1.9027  3.02451  4.313  8.0296
    #> par2 -0.8047  0.5023  1.78413  3.484  6.2260
    #> par3 -4.1927 -2.6005 -1.91258 -1.306 -0.2993
    #> par4 -5.1704 -3.3077 -2.34491 -1.400 -0.1107
    #> par5  1.2453  1.9455  2.44889  3.044  4.5671
    #> par6 -6.9166 -2.1409  0.07747  2.228  7.2453
    #> par7 -7.2802 -2.2749 -0.02967  2.200  7.2766

## Code of Conduct

Please note that the aphylo2 project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
