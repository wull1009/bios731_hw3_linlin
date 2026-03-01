logistic_pi <- function(eta) {
  1 / (1 + exp(-eta))
}

loglik_logistic <- function(beta, X, y) {
  eta <- as.vector(X %*% beta)
  # stable-ish form: sum(y*eta - log(1+exp(eta)))
  sum(y * eta - log1p(exp(eta)))
}

grad_logistic <- function(beta, X, y) {
  eta <- as.vector(X %*% beta)
  pi <- logistic_pi(eta)
  as.vector(t(X) %*% (y - pi))
}

hess_logistic <- function(beta, X, y) {
  eta <- as.vector(X %*% beta)
  pi <- logistic_pi(eta)
  w <- pi * (1 - pi)
  # Hessian of log-likelihood (concave): - X^T W X
  - crossprod(X * w, X)
}

vcov_from_hess <- function(beta_hat, X, y) {
  H <- hess_logistic(beta_hat, X, y)  # negative semidefinite
  # observed information is -H
  Iobs <- -H
  solve(Iobs)
}


wald_ci <- function(beta_hat, vcov_mat, level = 0.95) {
  z <- qnorm(1 - (1 - level)/2)
  se <- sqrt(diag(vcov_mat))
  cbind(
    estimate = beta_hat,
    se = se,
    lwr = beta_hat - z * se,
    upr = beta_hat + z * se
  )
}