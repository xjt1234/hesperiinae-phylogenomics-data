# Diagnostic-only wrapper of Toussaint et al. (2025) fit_bd.R.
# The likelihood call, parameter transformations, model flags, optimizer and
# AICc expression are unchanged. Added return fields expose optim diagnostics.

fit_bd <- function(phylo, tot_time, f.lamb, f.mu, lamb_par, mu_par, f = 1,
                   meth = "Nelder-Mead", cst.lamb = FALSE, cst.mu = FALSE,
                   expo.lamb = FALSE, expo.mu = FALSE, fix.mu = FALSE,
                   cond = "crown") {
  if (!inherits(phylo, "phylo")) stop("object phylo is not of class phylo")
  nobs <- Ntip(phylo)

  capture_optim <- function(init, objective) {
    warns <- character(0)
    fit <- withCallingHandlers(
      optim(init, objective, method = meth),
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
      mu_par_now <- init[(n_lamb + 1L):length(init)]
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
      mu_par = temp$par[(n_lamb + 1L):length(init)]
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
  res
}
