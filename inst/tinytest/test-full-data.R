library(geese)


# Placeholder with simple test
params_jmor = c(-2.0, -4.0, -6.0)

par_gain = params_jmor[1]
par_loss = params_jmor[2]
par_root = 1.0/(1.0 + exp(-params_jmor[3]))

eta_0 = exp(2.0 * par_gain) + 2.0 * exp(par_gain) + 1.0
eta_1 = 1.0 + 2.0 * exp(par_loss) + exp(2.0 * par_loss)

ans0_jmor = 
    par_root * (
        exp(par_loss) / eta_1 * exp(par_gain) / eta_0 +
        exp(par_loss) / (eta_1^2.0)
    ) + 
    (1.0 - par_root) * (
        (exp(par_gain)/eta_0)^2.0 +
        exp(2.0 * par_gain) / eta_0 * exp(par_loss) / eta_1 
    )

# Creating the jmor tree
ann_jmor <- list(c(0), c(1), c(1), c(9), c(9))
geneid_jmor <- c(0, 1, 2, 3, 4)
parent_jmor <- c(3, 3, 4, 4, -1)
duplication_jmor <- c(TRUE, TRUE, TRUE, TRUE, TRUE)

jmor_model <- new_geese(
  annotations = ann_jmor,
  geneid      = geneid_jmor,
  parent      = parent_jmor,
  duplication = duplication_jmor
  )

term_gains(jmor_model, 0)
term_loss(jmor_model, 0)

init_model(jmor_model)

# Aphylo is equivalent in this case
ap <- aphylo_from_data_frame(
  tree        = as.phylo(cbind(parent_jmor[-5], geneid_jmor[-5])), 
  annotations = data.frame(
    id = geneid_jmor,
    do.call(rbind, ann_jmor)
    ), 
  types = data.frame(
    id   = c(geneid_jmor),
    dupl = ifelse(duplication_jmor, 0, 1)
  )
)

ans1_jmor <- likelihood(jmor_model, par = params_jmor, as_log = FALSE)

ans2_jmor <- aphylo::LogLike(
  ap, 
  psi = c(0,0),
  mu_d = plogis(params_jmor[1:2]),
  mu_s = c(0,0),
  Pi = plogis(params_jmor[3]),
  eta = c(1,1)
)$ll

ans2_jmor <- exp(ans2_jmor)


tinytest::expect_equivalent(
  ans0_jmor,
  ans1_jmor
)

tinytest::expect_equivalent(
  ans0_jmor,
  ans2_jmor
)

# Testing leave one out --------------------------------------------------------
loo_pred0 <- predict_geese(jmor_model, par = params_jmor, leave_one_out = TRUE)

preds_loo <- predict_geese(jmor_model, par = params_jmor, leave_one_out = FALSE)
for (i in 1:3) {

  tmp_ann <- ann_jmor
  tmp_ann[[i]] <- c(9L)

  tmp_model <- new_geese(
    annotations = tmp_ann,
    geneid      = geneid_jmor,
    parent      = parent_jmor,
    duplication = duplication_jmor
  )

  term_gains(tmp_model, 0)
  term_loss(tmp_model, 0)

  init_model(tmp_model)

  preds_loo[[i]][[1L]] <- predict_geese(
    tmp_model, par = params_jmor,
    leave_one_out = FALSE
    )[[i]]
    
}

expect_equivalent(
  loo_pred0,
  preds_loo
)

preds_sim <- predict_geese_simulate(
  jmor_model,
  par = params_jmor,
  nsim = 1000000
  )

expect_true(
  all(abs(unlist(loo_pred0) - unlist(preds_sim)) < 0.01) 
)
