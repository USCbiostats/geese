# Reading the MCMC estimates
x <- readRDS("inst/simulation-study/simulation-study.rds")

# DGP parameters
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

# Checking coverage
estimates <- lapply(x, "[[", "estimates")
coverage <- sapply(estimates, function(e) {
  (params >= e[1,]) & (params <= e[3,])
})

colMeans(t(coverage)) # Almost 99 pcent

# Checking bias
bias <- sapply(estimates, function(e) {
  e[2,]
})

rownames(bias) <- c(
  "Gain0", "Gain1",
  "Lose0", "Lose1",
  "OnlyOne", "Root0", "Root1"
)

# Pretty plot
boxplot(t(bias), border = "darkgray", col = "steelblue", outline = FALSE,
        ylim = c(-6, 3), ylab = "Parameter Value", xlab = "Parameter")
# grid()
abline(h = 0, lty = 2, lwd = 2, col = "gray")
points(
  x = 1:7, lwd = 2,
  y = params, cex = 2, pch=23, col = "black", bg = "tomato"
)
text(
  x = 1:7 + .3,
  y = params + .15,
  labels = sprintf("%.2f", params)
)

