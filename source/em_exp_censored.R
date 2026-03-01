# EM for censored exponential data
# Observed: y_i = min(t_i, c_i), delta_i = I(t_i <= c_i)
# Model: t_i ~ Exponential(lambda)  (rate parameter lambda > 0)

em_exp_censored <- function(y, delta, lambda_init = NULL, tol = 1e-10, max_iter = 10000) {
  y <- as.numeric(y)
  delta <- as.integer(delta)
  
  n <- length(y)
  d <- sum(delta)
  
  if (is.null(lambda_init)) {
    # a reasonable start
    lambda <- d / sum(y)
    if (!is.finite(lambda) || lambda <= 0) lambda <- 1 / mean(y)
  } else {
    lambda <- as.numeric(lambda_init)
  }
  
  # observed log-likelihood for censored exponential:
  # l(lambda) = sum delta_i log(lambda) - lambda * sum y_i
  loglik_obs <- function(lam) {
    if (lam <= 0) return(-Inf)
    d * log(lam) - lam * sum(y)
  }
  
  trace <- data.frame(iter = 0, lambda = lambda, loglik = loglik_obs(lambda))
  
  for (k in seq_len(max_iter)) {
    lambda_old <- lambda
    
    # E-step: compute expected complete sufficient statistic sum(t_i | y, delta, lambda_old)
    # If uncensored: t_i = y_i
    # If censored: t_i | (t_i > y_i) has E[t_i] = y_i + 1/lambda_old (memoryless property)
    S <- sum(y) + sum(1 - delta) * (1 / lambda_old)
    
    # M-step: for complete data exponential, MLE lambda = n / sum(t_i)
    lambda <- n / S
    
    ll <- loglik_obs(lambda)
    trace <- rbind(trace, data.frame(iter = k, lambda = lambda, loglik = ll))
    
    if (abs(lambda - lambda_old) < tol) break
  }
  
  list(
    lambda_hat = lambda,
    trace = trace,
    converged = (nrow(trace) - 1 < max_iter),
    iter = nrow(trace) - 1,
    d = d,
    n = n
  )
}

# Closed-form MLE for censored exponential (for comparison)
mle_censored_exp <- function(y, delta) {
  d <- sum(delta)
  d / sum(y)
}

# Wald CI using observed information: l''(lambda) = -d / lambda^2
# Var(lambda_hat) approx = lambda_hat^2 / d
# Use log-scale CI to keep positivity:
# log(lambda) +/- z / sqrt(d)
ci_lambda_logwald <- function(lambda_hat, d, level = 0.95) {
  z <- qnorm(1 - (1 - level)/2)
  se_log <- 1 / sqrt(d)
  lwr <- exp(log(lambda_hat) - z * se_log)
  upr <- exp(log(lambda_hat) + z * se_log)
  c(lwr = lwr, est = lambda_hat, upr = upr)
}