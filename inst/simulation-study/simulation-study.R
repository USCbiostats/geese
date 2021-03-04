#!/bin/sh
#SBATCH --account=pdthomas_136
#SBATCH --partition=thomas
#SBATCH --mail-user=vegayon@usc.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=aphylo2-sim

library(aphylo2)

n     <- 30
nsims <- 2e3
NJOBS <- 100

# Testing
params <- c(
  # Gains
  2, 1,
  # Loss
  -2, -1,
  # Maxfuns
  2,
  # Root probabilities
  -5, -5
)

# Preparing data
annotations <- replicate(n * 2 - 1, c(9, 9), simplify = FALSE)

# Random tree
set.seed(31)
tree <- aphylo::sim_tree(n)$edge - 1L

duplication <- rep(TRUE, n * 2 - 1)

# Reading the data in
amodel <- new_model(
  annotations = annotations,
  geneid = c(tree[, 2]),
  parent = c(tree[, 1]),
  duplication = duplication
)

# Preparing the model
invisible({

  term_gains(amodel, 0:1)
  term_loss(amodel, 0:1)
  term_maxfuns(amodel, 1, 1)

  init(amodel)

})



likelihood(amodel, params*0)

# Simulating
set_seed(amodel, 111)
ans <- replicate(nsims, {
  sim_aphylo2(amodel, params)
}, simplify = FALSE)

# Computing transition probabilities

# Checking distribution
idx <- tree[which(tree[,1] == n),2] + 1

last2 <- lapply(ans, `[`, idx)
last2 <- unlist(last2, recursive = FALSE)
last2 <- do.call(rbind, last2)
colMeans(last2) # It should be something like c(1, .5, .5)

# Finding MLE in each one of them
library(slurmR)
out <- slurmR::Slurm_lapply(ans, function(a) {
    # Building the model
    amodel <- new_model(
      a, geneid = tree[,2], parent = tree[,1], duplication)

    invisible({

      term_gains(amodel, 0:1)
      term_loss(amodel, 0:1)
      term_maxfuns(amodel, 1, 1)

      init(amodel)

    })

    # Fitting the model
    names(mu) <- c("gain0", "gain1", "loss0", "loss1", "onefun", "root0", "root1")
    ans_mcmc <- tryCatch(aphylo2_mcmc(
      amodel,
      initial = mu * 0,
      nsteps  = 20000,
      kernel  = fmcmc::kernel_ram(warmup = 2000),
      prior   = function(p) dlogis(p, scale = 2, log = TRUE)
    ), error = function(e) e)


    if (inherits(estimates, "error"))
      return(estimates)

    estimates

  },
  njobs = NJOBS,
  job_name = "aphylo2-lapply",
  tmp_path = "/scratch/vegayon/",
  sbatch_opt    = list(
    account     = "pdthomas_136",
    partition   = "thomas",
    "mail-user" = "vegayon@usc.edu",
    "mail-type" = "END,FAIL"
    ),
  mc.cores = 1L
)

# Saving the output
saveRDS(out, file = "simulation-study.rds")

