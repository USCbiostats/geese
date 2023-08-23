library(aphylo)
library(geese)

# Preparing data
n <- 100L
annotations <- replicate(n * 2 - 1, c(9, 9), simplify = FALSE)

# Random tree
set.seed(31)
tree <- aphylo::sim_tree(n)$edge - 1L

# Tagging leaves
is_interior <- which(tree[, 2] %in% tree[, 1])

duplication <- sample.int(
  n = 2, size = n * 2 - 1, replace = TRUE, prob = c(.1, .9)
  ) == 1

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
term_maxfuns(amodel, 0, 1)
term_overall_changes(amodel, duplication = 0)
init_model(amodel)

# Testing
params <- c(
  # Gains
  2, 1.5,
  # Loss
  -2, -1.5,
  # Max funs
  2, 
  # Overall changes
  -5,
  # Root probabilities
  -10, -10
)

names(params) <- c("gain0", "gain1", "loss0", "loss1", "onefun", "root0", "root1")

# Simulating data
fake1 <- sim_geese(p = amodel, par = params, seed = 212)

fake1[is_interior] <- lapply(fake1[is_interior], \(x) c(NA, NA))

ap <- aphylo_from_data_frame(
  tree        = as.phylo(tree), 
  annotations = data.frame(
    id = c(tree[, 2], n),
    do.call(rbind, fake1)
    ),
  types = data.frame(
    id   = c(tree[, 2], n),
    dupl = ifelse(duplication, 0, 1)
    )
)
plot(ap)

ans_aphylo_mle <- aphylo_mle(ap ~ mu_d + mu_s + Pi)

amodel <- new_geese(
  annotations = fake1,
  geneid      = c(tree[, 2], n),
  parent      = c(tree[, 1],-1),
  duplication = duplication
  )

# Adding the model terms
term_overall_gains(amodel, duplication = 1)
term_overall_loss(amodel, duplication = 1)
term_overall_gains(amodel, duplication = 0)
term_overall_loss(amodel,  duplication = 0)
# We need to initialize to do all the accounting
init_model(amodel)

# Both methods should give the same loglikelihood ----------------------------------
params_aphylo <- c(.1, .01, .01, .001, .5)
params_geese  <- qlogis(c(params_aphylo, tail(params_aphylo, 1)))

ll_geese <- geese::likelihood(amodel, params_geese, as_log = TRUE)
ll_aphylo <- aphylo::LogLike(
  ap, psi = c(0,0),
  mu_d = params_aphylo[1:2], mu_s = params_aphylo[3:4],
  Pi = params_aphylo[5], eta = c(1,1))$ll

# Both methods should yield the same MLE -------------------------------------------

res <- geese_mle(amodel, lower = -10, upper = 10, method = "L-BFGS-B")
plogis(res$par)
coef(ans_aphylo_mle)

ans_geese <- predict_geese(amodel, res$par, only_annotated = TRUE) |>
  do.call(what = rbind)

# predict_geese_simulate(amodel, res$par, seed = 212, nsim = 1000) |>
#   do.call(what = rbind)

# The absolute difference should be below 0.01
differences <- abs(head(plogis(res$par), -2) - head(coef(ans_aphylo_mle), -1))/
    head(coef(ans_aphylo_mle), -1)

# prediction_score(ans_aphylo_mle)

# prediction_score(
#   x        = ans_geese,
#   expected = do.call(rbind, fake1)
# )


