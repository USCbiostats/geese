library(geese)

# Preparing data
n <- 1000L

annotations <- replicate(n * 2 - 1, c(9, 9), simplify = FALSE)

# Random tree
set.seed(31)
tree <- aphylo::sim_tree(n)$edge - 1L

# Sorting by the second column
tree <- tree[order(tree[, 2]), ]

duplication <- sample.int(
  n = 2, size = n * 2 - 1, replace = TRUE, prob = c(.4, .6)
  ) == 1

# Reading the data in
amodel <- new_geese(
  annotations = annotations,
  geneid = c(tree[, 2], n),
  parent = c(tree[, 1], -1),
  duplication = duplication
  )

# Preparing the model
term_gains(amodel, 0:1, duplication = 1)
term_loss(amodel, 0:1, duplication = 1)
term_gains(amodel, 0:1, duplication = 0)
term_loss(amodel, 0:1, duplication = 0)
term_maxfuns(amodel, 0, 1, duplication = 2)
init_model(amodel)

# Testing
params <- c(
  # Gains spe
  2, 1.5,
  # Loss
  -2, -1.5,
  # Gains spe
  -2, -1,
  # Loss spe
  -4, -4,
  # Max funs
  2, 
  # Root probabilities
  -10, -10
)
names(params) <- c(
  "gain0 dupl", "gain1 dupl",
  "loss0 dupl", "loss1 dupl",
  "gain0 spe", "gain1 spe",
  "loss0 spe", "loss1 spe",
  "onefun", 
  "root0", "root1"
  )


fake1 <- sim_geese(p = amodel, par = params, seed = 212)

# Adding some missings
set.seed(1231)
tomiss <- sample.int(
  n    = length(annotations),
  size = length(annotations) - 100, replace = FALSE
  )

fake1[tomiss] <- lapply(fake1[tomiss], function(x) {
  rep(NA, length(x))
  })

amodel <- new_geese(
  annotations = fake1,
  geneid      = c(tree[, 2], n),
  parent      = c(tree[, 1],-1),
  duplication = duplication
  )


# Adding the model terms
term_gains(amodel, 0:1, duplication = 1)
term_loss(amodel, 0:1, duplication = 1)
term_gains(amodel, 0:1, duplication = 0)
term_loss(amodel, 0:1, duplication = 0)
term_maxfuns(amodel, 0, 1, duplication = 2)
init_model(amodel)

# likelihood(amodel, params) 

predict_geese(amodel, params)
