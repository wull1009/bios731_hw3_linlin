newton_logistic <- function(X, y, beta_init = c(0, 0), tol = 1e-8, max_iter = 100) {
  beta <- as.numeric(beta_init)
  iter <- 0
   
  for (k in seq_len(max_iter)) {
    iter <- k
    g <- grad_logistic(beta, X, y)
    H <- hess_logistic(beta, X, y)
    
    step <- tryCatch(
      solve(H, g),  # H * step = g
      error = function(e) rep(NA_real_, length(beta))
    )
    
    if (anyNA(step)) break
    
    beta_new <- beta - step
    
    if (max(abs(beta_new - beta)) < tol) {
      beta <- beta_new
      break
    }
    beta <- beta_new
  }
  
  list(beta_hat = beta, iter = iter, converged = (iter < max_iter))
}