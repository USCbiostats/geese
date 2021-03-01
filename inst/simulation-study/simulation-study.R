library(aphylo2)

# Preparing data
geneid   <- c(0, 1, 2, 3, 4, 5, 6)
parentid <- c(4, 4, 5, 5, 6, 6, -1)
duplication <- rep(TRUE, 7)

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

# Simulating
ans <- replicate(100, {
  sim_aphylo2(amodel, c(10, 10, -10, 10, -100, -100, -100)/5)
}, simplify = FALSE)

# Checking distribution
last2 <- lapply(ans, "[", 5:6)
last2 <- unlist(last2, recursive = FALSE)
last2 <- do.call(rbind, last2)
colMeans(last2) # It should be something like c(1, .5, .5)

# Finding MLE in each one of them
iter <- 1
mles <- lapply(ans, function(a) {
  # Building the model
  amodel <- new_model(a, geneid, parentid, duplication)

  invisible({

    term_cogain(amodel, 0, 1)
    term_cogain(amodel, 0, 2)
    term_cogain(amodel, 1, 2)
    term_maxfuns(amodel, 2, 2)

    init(amodel)

  })

  # Fitting the model
  estimates <- tryCatch(optim(params*0, function(f) {
    log(likelihood(amodel, f))
  }, control = list(fnscale=-1, maxit = 2e3),
  hessian = TRUE), error = function(e) e)

  iter <<- iter + 1L
  if (!(iter %% 10))
    message(iter, " completed")

  if (inherits(estimates, "error"))
    return(estimates)

  with(estimates, list(par = par, hessian = hessian, counts = counts))
})

estimates <- sapply(mles, "[[", "par")
estimates <- do.call(rbind, estimates[sapply(estimates, length) > 0])

boxplot(estimates)
abline(h = 0, lwd=2, lty=2)
