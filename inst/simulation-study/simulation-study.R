library(aphylo2)

n <- 20

# Preparing data
annotations <- replicate(n * 2 - 1, c(9, 9, 9), simplify = FALSE)

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

  term_cogain(amodel, 0, 1)
  term_cogain(amodel, 0, 2)
  term_cogain(amodel, 1, 2)
  term_maxfuns(amodel, 2, 2)

  init(amodel)

})

# Testing
params <- c(
  # Cogain
  5, 5, -5,
  # max funs
  10,
  # Root probabilities
  -10, -10, -10
)

likelihood(amodel, params*0)

stop()

# Simulating
ans <- replicate(500, {
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

estimates <- lapply(mles, "[[", "par")
estimates <- do.call(rbind, estimates[sapply(estimates, length) > 0])

boxplot(estimates)
abline(h = 0, lwd=2, lty=2)
