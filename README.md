
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
  nsteps  = 20000,
  kernel  = fmcmc::kernel_ram(warmup = 2000), 
  prior   = function(p) dlogis(p, scale = 2, log = TRUE)
  )
#> 
#> |0%                 |25%                |50%                |75%           100%|
#> --------------------------------------------------------------------------------
#> ////////////////////////////////////////////////////////////////////////////////
```

``` r
op <- par(mfrow = c(4, 2))
for (i in 1:ncol(ans_mcmc)) {
  tmpx <- window(ans_mcmc, start = 15000)[,i,drop=FALSE]
  plot(
    density(tmpx),
    main = names(params)[i]
    )
  abline(v = quantile(tmpx, .5), lty = 2, lwd = 2, col = "tomato")
  abline(v = params[i], lty=3, lwd=2, col = "steelblue")
}
par(op)
```

<img src="man/figures/README-mcmc-analysis-1.png" width="100%" />

``` r
summary(window(ans_mcmc, start = 15000))
#> 
#> Iterations = 15000:20000
#> Thinning interval = 1 
#> Number of chains = 1 
#> Sample size per chain = 5001 
#> 
#> 1. Empirical mean and standard deviation for each variable,
#>    plus standard error of the mean:
#> 
#>         Mean     SD Naive SE Time-series SE
#> par1  3.4737 2.3279  0.03292        0.37906
#> par2  2.8471 2.7343  0.03866        0.57655
#> par3 -2.0180 0.9663  0.01366        0.10986
#> par4 -2.2933 1.3051  0.01846        0.16502
#> par5  2.5153 0.8131  0.01150        0.09235
#> par6  0.5291 3.1101  0.04398        0.39694
#> par7 -0.5479 3.1468  0.04450        0.42594
#> 
#> 2. Quantiles for each variable:
#> 
#>         2.5%     25%     50%    75%   97.5%
#> par1  0.1964  2.1794  3.1351  4.276  9.9887
#> par2 -0.6170  0.7377  2.4391  4.068 10.0958
#> par3 -4.1523 -2.6471 -1.8872 -1.323 -0.4065
#> par4 -4.8967 -3.2830 -2.1809 -1.295 -0.1173
#> par5  1.2103  1.9248  2.4225  2.974  4.4613
#> par6 -5.6624 -1.5766  0.5858  2.324  7.4192
#> par7 -6.7773 -2.7457 -0.5531  1.516  5.7348
```

## Code of Conduct

Please note that the aphylo2 project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
