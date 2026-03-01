simulate_logistic_data <- function(n = 200, beta0 = 1, beta1 = 0.3, seed = 123) {
  set.seed(seed)
  x <- rnorm(n, mean = 0, sd = 1)
  X <- cbind(1, x)  # intercept + x
  eta <- beta0 + beta1 * x
  p <- 1 / (1 + exp(-eta))
  y <- rbinom(n, size = 1, prob = p)
  list(
    n = n,
    beta_true = c(beta0, beta1),
    x = x,
    X = X,
    y = y
  )
}