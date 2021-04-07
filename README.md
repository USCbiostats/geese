
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Integrative Methods of Analysis for Genetic
Epidemiology](https://raw.githubusercontent.com/USCbiostats/badges/master/tommy-image-badge.svg)](https://image.usc.edu)

# geese: *GE*ne-functional *E*volution using *S*uffici*E*ncy <img src="man/figures/logo.svg" align="right" width="180px"/>

<!-- badges: start -->
<!-- badges: end -->

The goal of geese is to …

## Installation

<!-- You can install the released version of geese from [CRAN](https://CRAN.R-project.org) with: -->
<!-- ``` r -->
<!-- install.packages("geese") -->
<!-- ``` -->

The development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("USCbiostats/geese")
```

# Examples

## Simulating annotations (two different sets)

``` r
library(geese)

# Preparing data
n <- 100L
annotations <- replicate(n * 2 - 1, c(9, 9), simplify = FALSE)

# Random tree
set.seed(31)
tree <- aphylo::sim_tree(n)$edge - 1L

duplication <- rep(TRUE, n * 2 - 1)

# Reading the data in
amodel <- new_geese(
  annotations = annotations,
  geneid = c(tree[, 2], n),
  parent = c(tree[, 1], -1),
  duplication = duplication
  )

# Preparing the model
term_gains(amodel, 0:1)
term_loss(amodel, 0:1)
term_maxfuns(amodel, 1, 1)
init_model(amodel)

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
#> [1] 0

# Simulating data
fake1 <- sim_geese(p = amodel, par = params, seed = 1110)
fake2 <- sim_geese(p = amodel, par = params)
```

We can now visualize either of the annotations using the
[aphylo](https://github.com/USCbiostats/aphylo) package.

``` r
library(aphylo)
#> Loading required package: ape
ap <- aphylo_from_data_frame(
  tree        = as.phylo(tree), 
  annotations = data.frame(
    id = c(tree[, 2], n),
    do.call(rbind, fake1)
    )
)
plot(ap)
```

<img src="man/figures/README-viz-with-aphylo-1.png" width="100%" />

## Model fitting MLE

``` r
# Creating the object
amodel <- new_geese(
  annotations = fake1,
  geneid      = c(tree[, 2], n),
  parent      = c(tree[, 1],-1),
  duplication = duplication
  )

# Adding the model terms
term_gains(amodel, 0:1)
term_loss(amodel, 0:1)
term_maxfuns(amodel, 1, 1)

# We need to initialize to do all the accountintg
init_model(amodel)

# Finding MLE
geese_mle(amodel, hessian = TRUE)
#> $par
#> [1]  1.3024291  1.4904183 -1.4668169 -0.9165704  1.8952822 -0.1918496 -0.5043625
#> 
#> $value
#> [1] -58.54824
#> 
#> $counts
#> function gradient 
#>      454       NA 
#> 
#> $convergence
#> [1] 0
#> 
#> $message
#> NULL
#> 
#> $hessian
#>               [,1]          [,2]          [,3]          [,4]          [,5]
#> [1,] -1.920339e+00  3.092661e-02  9.719647e-02 -5.826438e-01  1.230857e+00
#> [2,]  3.092661e-02 -2.155737e+00 -2.661554e-01  2.058254e+00  2.388066e+00
#> [3,]  9.719647e-02 -2.661554e-01 -1.562040e+00  1.126001e-01 -1.189657e+00
#> [4,] -5.826438e-01  2.058254e+00  1.126001e-01 -5.959471e+00 -3.427126e+00
#> [5,]  1.230857e+00  2.388066e+00 -1.189657e+00 -3.427126e+00 -2.441145e+01
#> [6,] -8.064660e-07  1.930012e-06 -1.187495e-06  3.366196e-06 -4.565237e-07
#> [7,]  6.354028e-06 -8.340448e-04  9.558576e-06 -6.562351e-04  1.888232e-04
#>               [,6]          [,7]
#> [1,] -8.064660e-07  6.354028e-06
#> [2,]  1.930012e-06 -8.340448e-04
#> [3,] -1.187495e-06  9.558576e-06
#> [4,]  3.366196e-06 -6.562351e-04
#> [5,] -4.565237e-07  1.888232e-04
#> [6,]  4.440892e-08  6.750156e-08
#> [7,]  6.750156e-08 -2.657963e-05
```

## Model fitting MCMC

``` r
set.seed(122)
ans_mcmc <- geese_mcmc(
  amodel,
  nsteps  = 40000,
  kernel  = fmcmc::kernel_ram(warmup = 2000), 
  prior   = function(p) dlogis(p, scale = 2, log = TRUE)
  )
```

We can take a look at the results like this:

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
    #>          Mean     SD Naive SE Time-series SE
    #> par1  1.26120 0.7629 0.004825       0.023220
    #> par2  1.93817 1.2494 0.007902       0.072333
    #> par3 -1.51192 0.8443 0.005340       0.026278
    #> par4 -0.86702 0.5525 0.003494       0.017003
    #> par5  1.99785 0.2394 0.001514       0.006938
    #> par6 -0.19538 3.6550 0.023116       0.175974
    #> par7 -0.02748 3.6706 0.023215       0.161282
    #> 
    #> 2. Quantiles for each variable:
    #> 
    #>         2.5%     25%      50%     75%   97.5%
    #> par1 -0.1486  0.7323  1.23488  1.7457 2.83010
    #> par2  0.1582  1.1853  1.75376  2.4721 4.84487
    #> par3 -3.2207 -2.0756 -1.47760 -0.9378 0.04801
    #> par4 -1.8700 -1.2413 -0.89249 -0.5323 0.30781
    #> par5  1.5569  1.8305  1.98664  2.1502 2.51727
    #> par6 -7.7258 -2.4187 -0.12163  2.0795 6.98715
    #> par7 -7.8283 -2.2335 -0.01834  2.2817 7.28420

## Using a flock

GEESE models can be grouped (pooled) into a flock.

``` r
flock <- new_flock()

# Adding first set of annotations
add_geese(
  flock,
  annotations = fake1,
  geneid      = c(tree[, 2], n),
  parent      = c(tree[, 1],-1),
  duplication = duplication  
)

# Now the second set
add_geese(
  flock,
  annotations = fake2,
  geneid      = c(tree[, 2], n),
  parent      = c(tree[, 1],-1),
  duplication = duplication  
)

# Adding the model terms
term_gains(flock, 0:1)
term_loss(flock, 0:1)
term_maxfuns(flock, 1, 1)

# We need to initialize to do all the accountintg
init_model(flock)
```

We can use the same program to fit the MCMC

``` r
set.seed(122)
ans_mcmc2 <- geese_mcmc(
  flock,
  nsteps  = 20000,
  kernel  = fmcmc::kernel_ram(warmup = 2000), 
  prior   = function(p) dlogis(p, scale = 2, log = TRUE)
  )
```

``` r
op <- par(
  mfrow = c(4, 2), #tcl=.5,
  las=1, mar = c(3,3,1,0),
  bty = "n", oma = rep(1,4)
  )
for (i in 1:ncol(ans_mcmc2)) {
  tmpx <- window(ans_mcmc2, start = 10000)[,i,drop=FALSE]
  
  coda::traceplot(
    tmpx, smooth = FALSE, ylim = c(-11,11), col = rgb(0, 128, 128, maxColorValue = 255), 
    main = names(params)[i]
    )
  abline(h = params[i], lty=3, lwd=4, col = "red")
}
par(op)
```

<img src="man/figures/README-viz-flock-1.png" width="100%" />

``` r
summary(window(ans_mcmc2, start = 10000))
#> 
#> Iterations = 10000:20000
#> Thinning interval = 1 
#> Number of chains = 1 
#> Sample size per chain = 10001 
#> 
#> 1. Empirical mean and standard deviation for each variable,
#>    plus standard error of the mean:
#> 
#>          Mean     SD Naive SE Time-series SE
#> par1  1.75029 0.5658 0.005658       0.031031
#> par2  1.15245 0.6418 0.006418       0.031898
#> par3 -1.27191 0.5489 0.005488       0.026757
#> par4 -1.37351 0.3635 0.003635       0.015714
#> par5  1.81130 0.1630 0.001630       0.007164
#> par6  0.20264 2.6196 0.026195       0.304291
#> par7 -0.09525 1.9988 0.019987       0.192587
#> 
#> 2. Quantiles for each variable:
#> 
#>          2.5%     25%     50%     75%   97.5%
#> par1  0.72459  1.3498  1.7308  2.1272  2.8592
#> par2  0.03824  0.7103  1.1054  1.5305  2.5470
#> par3 -2.37365 -1.6300 -1.2692 -0.9166 -0.1438
#> par4 -2.09160 -1.6242 -1.3736 -1.1207 -0.6764
#> par5  1.49584  1.6992  1.8025  1.9109  2.1439
#> par6 -4.78636 -1.4981  0.1150  1.8273  5.1856
#> par7 -4.12554 -1.3982 -0.1091  1.1896  4.0561
```

## Code of Conduct

Please note that the aphylo2 project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
