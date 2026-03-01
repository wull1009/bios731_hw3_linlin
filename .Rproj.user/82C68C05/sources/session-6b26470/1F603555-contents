run_problem4 <- function(save_prefix = "problem4", tol = 1e-10) {
  library(survival)
  library(ggplot2)
  library(dplyr)
  library(here)
  
  library(survival)
  
  if (exists("veteran", where = asNamespace("survival"), inherits = FALSE)) {
    veteran <- get("veteran", envir = asNamespace("survival"))
  } else if (exists("veteran", envir = .GlobalEnv)) {
    veteran <- get("veteran", envir = .GlobalEnv)
  } else {
    stop("Cannot find a `veteran` dataset/object. Your survival install doesn't include survival::veteran, and no object named `veteran` exists in your environment.")
  }
  
  if (!all(c("time", "status") %in% names(veteran))) {
    stop("`veteran` must contain columns named `time` and `status`.")
  }
  y <- veteran$time
  delta <- as.integer(veteran$status == 1)  
  
  if (sum(delta) == 0) stop("No events found: check status coding for veteran dataset.")
  
  # --- EM fit ---
  fit_em <- em_exp_censored(y, delta, tol = tol)
  lambda_em <- fit_em$lambda_hat
  ci_em <- ci_lambda_logwald(lambda_em, d = fit_em$d)
  
  # --- Closed-form MLE (should match EM limit) ---
  lambda_closed <- mle_censored_exp(y, delta)
  
  # --- AFT fit with exponential errors via survreg ---
  fit_aft <- survreg(Surv(veteran$time, delta) ~ 1, data = veteran, dist = "exponential")  
  mu_hat <- as.numeric(coef(fit_aft))
  se_mu  <- sqrt(vcov(fit_aft)[1, 1])
  lambda_aft <- exp(-mu_hat)
  
  z <- qnorm(0.975)
  mu_lwr <- mu_hat - z * se_mu
  mu_upr <- mu_hat + z * se_mu
  lambda_aft_ci <- c(lwr = exp(-mu_upr), est = lambda_aft, upr = exp(-mu_lwr))  
  # --- assemble summary table ---
  summary_tbl <- data.frame(
    method = c("EM", "Closed-form MLE", "AFT (survreg)"),
    lambda_hat = c(lambda_em, lambda_closed, lambda_aft),
    lambda_lwr = c(ci_em["lwr"], NA, lambda_aft_ci["lwr"]),
    lambda_upr = c(ci_em["upr"], NA, lambda_aft_ci["upr"]),
    iter = c(fit_em$iter, NA, NA)
  )
  
  # save CSV
  write.csv(summary_tbl, here("results", paste0(save_prefix, "_lambda_summary.csv")), row.names = FALSE)
  
  # --- EM trace plot ---
  p <- ggplot(fit_em$trace, aes(x = iter, y = lambda)) +
    geom_line() +
    geom_point(size = 1.5) +
    theme_bw() +
    labs(title = "EM trace for censored exponential (veteran data)",
         x = "Iteration", y = expression(lambda))
  
  ggsave(here("results", paste0(save_prefix, "_em_trace.png")), p, width = 7, height = 4, dpi = 300)
  
  list(
    y = y,
    delta = delta,
    fit_em = fit_em,
    fit_aft = fit_aft,
    summary = summary_tbl,
    trace_plot = p
  )
}