run_all_methods <- function(X, y, beta_init = c(0, 0), tol = 1e-8, max_iter_newton = 100, max_iter_mm = 200) {
  out <- list()
  
  # ---------- Newton ----------
  t0 <- proc.time()
  fit_newton <- newton_logistic(X, y, beta_init = beta_init, tol = tol, max_iter = max_iter_newton)
  t1 <- proc.time()
  beta_n <- fit_newton$beta_hat
  vc_n <- vcov_from_hess(beta_n, X, y)
  ci_n <- wald_ci(beta_n, vc_n)
  out$Newton <- list(beta = beta_n, vcov = vc_n, ci = ci_n, time = (t1 - t0)[["elapsed"]], iter = fit_newton$iter)
  
  # ---------- MM ----------
  t0 <- proc.time()
  fit_mm <- mm_logistic(X, y, beta_init = beta_init, tol = tol, max_iter = max_iter_mm)
  t1 <- proc.time()
  beta_m <- fit_mm$beta_hat
  vc_m <- vcov_from_hess(beta_m, X, y)
  ci_m <- wald_ci(beta_m, vc_m)
  out$MM <- list(beta = beta_m, vcov = vc_m, ci = ci_m, time = (t1 - t0)[["elapsed"]], iter = fit_mm$iter)
  
  # ---------- glm ----------
  t0 <- proc.time()
  df <- data.frame(y = y, x = X[, 2])
  fit_glm <- glm(y ~ x, data = df, family = binomial())
  t1 <- proc.time()
  beta_g <- as.numeric(coef(fit_glm))
  vc_g <- as.matrix(vcov(fit_glm))
  ci_g <- wald_ci(beta_g, vc_g)
  # number of iterations used by IRLS:
  iter_g <- fit_glm$iter %||% NA_integer_
  out$glm <- list(beta = beta_g, vcov = vc_g, ci = ci_g, time = (t1 - t0)[["elapsed"]], iter = iter_g)
  
  # ---------- optim (BFGS) ----------
  negloglik <- function(beta) -loglik_logistic(beta, X, y)
  neggrad   <- function(beta) -grad_logistic(beta, X, y)
  
  t0 <- proc.time()
  fit_opt <- optim(par = beta_init, fn = negloglik, gr = neggrad,
                   method = "BFGS", hessian = TRUE, control = list(maxit = 1000, reltol = tol))
  t1 <- proc.time()
  beta_o <- fit_opt$par
  # optim returns Hessian of objective = Hessian of -loglik = observed information
  vc_o <- solve(fit_opt$hessian)
  ci_o <- wald_ci(beta_o, vc_o)
  out$optim_BFGS <- list(beta = beta_o, vcov = vc_o, ci = ci_o, time = (t1 - t0)[["elapsed"]], iter = fit_opt$counts[["function"]])
  
  # assemble a neat data.frame
  methods <- names(out)
  res <- do.call(rbind, lapply(methods, function(m) {
    ci <- out[[m]]$ci
    data.frame(
      method = m,
      beta0_hat = ci[1, "estimate"],
      beta0_lwr = ci[1, "lwr"],
      beta0_upr = ci[1, "upr"],
      beta1_hat = ci[2, "estimate"],
      beta1_lwr = ci[2, "lwr"],
      beta1_upr = ci[2, "upr"],
      time_sec = out[[m]]$time,
      iter = out[[m]]$iter
    )
  }))
  rownames(res) <- NULL
  
  list(detail = out, summary = res)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b