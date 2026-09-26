# Long-run diagnostic wrappers for the BCSTDTempVar rescue analysis.
# These functions retain the archived Toussaint/RPANDA likelihood and absolute
# rate transformation. Only optim() control and returned diagnostics are added.

fit_bd_diagnostics_maxit_v3 <- function(
    phylo, tot_time, f.lamb, f.mu, lamb_par, mu_par, f = 1,
    meth = "Nelder-Mead", cst.lamb = FALSE, cst.mu = FALSE,
    expo.lamb = FALSE, expo.mu = FALSE, fix.mu = FALSE,
    cond = "crown", control = list(maxit = 10000L)) {
  if (!inherits(phylo, "phylo")) stop("object phylo is not of class phylo")
  nobs <- Ntip(phylo)

  capture_optim <- function(init, objective) {
    warns <- character(0)
    fit <- withCallingHandlers(
      optim(init, objective, method = meth, control = control),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    list(fit = fit, warnings = unique(warns))
  }

  if (!fix.mu) {
    init <- c(lamb_par, mu_par)
    p <- length(init)
    n_lamb <- length(lamb_par)
    optimLH <- function(init) {
      lamb_par_now <- init[seq_len(n_lamb)]
      mu_par_now <- init[n_lamb + seq_along(mu_par)]
      f.lamb.par <- function(t) abs(f.lamb(t, lamb_par_now))
      f.mu.par <- function(t) abs(f.mu(t, mu_par_now))
      LH <- likelihood_bd(
        phylo, tot_time, f.lamb.par, f.mu.par, f,
        cst.lamb = cst.lamb, cst.mu = cst.mu,
        expo.lamb = expo.lamb, expo.mu = expo.mu, cond = cond
      )
      -LH
    }
    captured <- capture_optim(init, optimLH)
    temp <- captured$fit
    res <- list(
      model = "birth death", LH = -temp$value,
      aicc = 2 * temp$value + 2 * p + (2 * p * (p + 1)) / (nobs - p - 1),
      lamb_par = temp$par[seq_len(n_lamb)],
      mu_par = temp$par[n_lamb + seq_along(mu_par)]
    )
  } else {
    init <- c(lamb_par)
    p <- length(init)
    n_lamb <- length(lamb_par)
    optimLH <- function(init) {
      lamb_par_now <- init[seq_len(n_lamb)]
      f.lamb.par <- function(t) abs(f.lamb(t, lamb_par_now))
      f.mu.par <- function(t) abs(f.mu(t, mu_par))
      LH <- likelihood_bd(
        phylo, tot_time, f.lamb.par, f.mu.par, f,
        cst.lamb = cst.lamb, cst.mu = TRUE,
        expo.lamb = expo.lamb, cond = cond
      )
      -LH
    }
    captured <- capture_optim(init, optimLH)
    temp <- captured$fit
    res <- list(
      model = "birth death", LH = -temp$value,
      aicc = 2 * temp$value + 2 * p + (2 * p * (p + 1)) / (nobs - p - 1),
      lamb_par = temp$par[seq_len(n_lamb)]
    )
  }

  res$convergence <- temp$convergence
  res$message <- if (is.null(temp$message)) NA_character_ else temp$message
  res$counts <- temp$counts
  res$optim_value <- temp$value
  res$optim_par <- temp$par
  res$start <- init
  res$warnings <- captured$warnings
  res$method <- meth
  res$control <- control
  res
}

fit_env_bd_diagnostics_maxit_v3 <- function(
    phylo, env_data, tot_time, f.lamb, f.mu, lamb_par, mu_par,
    df = NULL, f = 1, meth = "Nelder-Mead", cst.lamb = FALSE,
    cst.mu = FALSE, expo.lamb = FALSE, expo.mu = FALSE,
    fix.mu = FALSE, cond = "crown", control = list(maxit = 10000L)) {
  if (is.null(df)) df <- smooth.spline(x = env_data[, 1], env_data[, 2])$df
  spline_result <- sm.spline(env_data[, 1], env_data[, 2], df = df)
  env_func <- function(t) predict(spline_result, t)

  lower_bound_control <- 0.10
  upper_bound_control <- 0.10
  lower_bound <- min(env_data[, 1])
  upper_bound <- max(env_data[, 1])
  time_tabulated <- seq(
    from = lower_bound * (1 - lower_bound_control),
    to = upper_bound * (1 + upper_bound_control), length.out = 1 + 1e6
  )
  env_tabulated <- env_func(time_tabulated)
  env_func_tab <- function(t) {
    b <- upper_bound * (1 + upper_bound_control)
    a <- lower_bound * (1 - lower_bound_control)
    n <- length(env_tabulated) - 1L
    index <- 1L + as.integer((t - a) * n / (b - a))
    env_tabulated[index]
  }
  f.lamb.env <- function(t, y) f.lamb(t, env_func_tab(t), y)
  f.mu.env <- function(t, y) f.mu(t, env_func_tab(t), y)
  res <- fit_bd_diagnostics_maxit_v3(
    phylo, tot_time, f.lamb.env, f.mu.env, lamb_par, mu_par, f,
    meth, cst.lamb, cst.mu, expo.lamb, expo.mu, fix.mu, cond, control
  )
  res$model <- "environmental birth death"
  res
}
