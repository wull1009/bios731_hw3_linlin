mm_logistic <- function(X, y, beta_init = c(0, 0), tol = 1e-8, max_iter = 200) {
  beta <- as.numeric(beta_init)
  n <- nrow(X)
  p <- ncol(X)  
  
  iter <- 0
  
  for (k in seq_len(max_iter)) {
    iter <- k
    beta_old <- beta
    
    eta_k <- as.vector(X %*% beta_old)
    a_i <- exp(eta_k) / (1 + exp(eta_k))  # = pi_k
    
    # update each coordinate j by solving the equation in (C)
    for (j in seq_len(p)) {
      xij <- X[, j]
      # constants for this coordinate update
      c_i <- a_i * xij * exp(-p * xij * beta_old[j])  # note: uses beta^{(k)} in the exp(-p Xij theta_j^{(k)})
      
      # We need solve: - sum_i c_i * exp(p xij * theta_j) + sum_i y_i xij = 0
      rhs <- sum(y * xij)
      
      f <- function(thetaj) {
        - sum(c_i * exp(p * xij * thetaj)) + rhs
      }
      
      # Try bracket for uniroot
      lo <- beta_old[j] - 10
      hi <- beta_old[j] + 10
      flo <- f(lo); fhi <- f(hi)
      
      if (is.finite(flo) && is.finite(fhi) && flo * fhi < 0) {
        beta[j] <- uniroot(f, lower = lo, upper = hi)$root
      } else {
        # fallback: 1D minimize squared equation
        obj <- function(thetaj) f(thetaj)^2
        beta[j] <- optim(par = beta_old[j], fn = obj, method = "Brent",
                         lower = beta_old[j] - 20, upper = beta_old[j] + 20)$par
      }
    } 
    
    if (max(abs(beta - beta_old)) < tol) break
  }
  
  list(beta_hat = beta, iter = iter, converged = (iter < max_iter))
}