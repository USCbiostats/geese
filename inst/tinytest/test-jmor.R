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
ans1_jmor <- likelihood(jmor_model, par = params_jmor, as_log = FALSE)

tinytest::expect_equivalent(
  ans0_jmor,
  ans1_jmor
)

predict_geese(jmor_model, par = params_jmor, leave_one_out = TRUE)
set.seed(112)
predict_geese_simulate(jmor_model, par = params_jmor, nsim = 10000)

# Fitting the model
ans <- geese_mle(jmor_model)
predict_geese(jmor_model, par = ans$par, leave_one_out = TRUE)
predict_geese_simulate(jmor_model, par = ans$par, nsim = 10000)

