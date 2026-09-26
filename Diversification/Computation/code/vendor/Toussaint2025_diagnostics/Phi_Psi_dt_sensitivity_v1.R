# Parameterized copies of the archived Toussaint/RPANDA .Phi and .Psi helpers.
#
# These constructors are for numerical-sensitivity diagnosis only.  For
# dt = 0.005 they reproduce the general-environment branches in the archived
# Phi.R and Psi.R.  The analytical constant/exponential branches are copied
# unchanged.  The formal model fit continues to use the archived helpers.

make_Phi_dt_sensitivity_v1 <- function(dt) {
  stopifnot(length(dt) == 1L, is.finite(dt), dt > 0)
  force(dt)
  function(t, f.lamb, f.mu, f, cst.lamb = FALSE, cst.mu = FALSE,
           expo.lamb = FALSE, expo.mu = FALSE) {
    if (cst.lamb && cst.mu) {
      lamb <- f.lamb(0)
      mu <- f.mu(0)
      r <- lamb - mu
      return(1 - r * exp(r * t) / (r / f + lamb * (exp(r * t) - 1)))
    }

    if (cst.lamb && expo.mu) {
      lamb0 <- f.lamb(0)
      mu0 <- f.mu(0)
      beta <- log(f.mu(1) / mu0)
      r.int <- function(x, y) lamb0 * (y - x) - mu0 / beta * (exp(beta * y) - exp(beta * x))
      g <- function(y) r.int(0, y)
      gvect <- function(y) mapply(g, y)
      r.int.0 <- function(y) exp(gvect(y)) * f.lamb(y)
      r.int.int <- function(x, y) .Integrate(r.int.0, x, y, stop.on.error = FALSE)
      return(1 - exp(r.int(0, t)) / (1 / f + r.int.int(0, t)))
    }

    if (expo.lamb && cst.mu) {
      lamb0 <- f.lamb(0)
      alpha <- log(f.lamb(1) / lamb0)
      mu0 <- f.mu(0)
      r.int <- function(x, y) lamb0 / alpha * (exp(alpha * y) - exp(alpha * x)) - mu0 * (y - x)
      g <- function(y) r.int(0, y)
      gvect <- function(y) mapply(g, y)
      r.int.0 <- function(y) exp(gvect(y)) * f.lamb(y)
      r.int.int <- function(x, y) .Integrate(r.int.0, x, y, stop.on.error = FALSE)
      return(1 - exp(r.int(0, t)) / (1 / f + r.int.int(0, t)))
    }

    if (expo.lamb && expo.mu) {
      lamb0 <- f.lamb(0)
      alpha <- log(f.lamb(1) / lamb0)
      mu0 <- f.mu(0)
      beta <- log(f.mu(1) / mu0)
      r.int <- function(x, y) {
        lamb0 / alpha * (exp(alpha * y) - exp(alpha * x)) -
          mu0 / beta * (exp(beta * y) - exp(beta * x))
      }
      g <- function(y) r.int(0, y)
      gvect <- function(y) mapply(g, y)
      r.int.0 <- function(y) exp(gvect(y)) * f.lamb(y)
      r.int.int <- function(x, y) .Integrate(r.int.0, x, y, stop.on.error = FALSE)
      return(1 - exp(r.int(0, t)) / (1 / f + r.int.int(0, t)))
    }

    ageMin <- 0
    ageMax <- t
    Nintervals <- 1L + as.integer((ageMax - ageMin) / dt)
    X <- seq(ageMin, ageMax, length.out = Nintervals + 1L)
    r <- function(z) f.lamb(z) - f.mu(z)
    r.int <- cumsum(r(X)) * (ageMax - ageMin) / Nintervals
    r.int.0 <- function(y) {
      index <- 1L + as.integer((y - ageMin) * Nintervals / (ageMax - ageMin))
      exp(r.int[index]) * f.lamb(y)
    }
    r.int.int.tab <- cumsum(r.int.0(X)) * (ageMax - ageMin) / Nintervals
    r.int.int <- function(x, y) {
      indy <- 1L + as.integer((y - ageMin) * Nintervals / (ageMax - ageMin))
      indx <- 1L + as.integer((x - ageMin) * Nintervals / (ageMax - ageMin))
      r.int.int.tab[indy] - r.int.int.tab[indx]
    }
    rit <- r.int[1L + Nintervals]
    ri0t <- r.int.int(0, t)
    1 - exp(rit) / (1 / f + ri0t)
  }
}

make_Psi_dt_sensitivity_v1 <- function(dt) {
  stopifnot(length(dt) == 1L, is.finite(dt), dt > 0)
  force(dt)
  function(s, t, f.lamb, f.mu, f, cst.lamb = FALSE, cst.mu = FALSE,
           expo.lamb = FALSE, expo.mu = FALSE) {
    if (cst.lamb && cst.mu) {
      lamb <- f.lamb(0)
      mu <- f.mu(0)
      r <- lamb - mu
      return(exp(r * (t - s)) *
        abs(1 + lamb * (exp(r * t) - exp(r * s)) /
          (r / f + lamb * (exp(r * s) - 1)))^(-2))
    }

    if (cst.lamb && expo.mu) {
      lamb0 <- f.lamb(0)
      mu0 <- f.mu(0)
      beta <- log(f.mu(1) / mu0)
      r.int <- function(x, y) lamb0 * (y - x) - mu0 / beta * (exp(beta * y) - exp(beta * x))
      g <- function(y) r.int(0, y)
      gvect <- function(y) mapply(g, y)
      r.int.0 <- function(y) exp(gvect(y)) * f.lamb(y)
      r.int.int <- function(x, y) .Integrate(r.int.0, x, y, stop.on.error = FALSE)
      return(exp(r.int(s, t)) * abs(1 + r.int.int(s, t) / (1 / f + r.int.int(0, s)))^(-2))
    }

    if (expo.lamb && cst.mu) {
      lamb0 <- f.lamb(0)
      alpha <- log(f.lamb(1) / lamb0)
      mu0 <- f.mu(0)
      r.int <- function(x, y) lamb0 / alpha * (exp(alpha * y) - exp(alpha * x)) - mu0 * (y - x)
      g <- function(y) r.int(0, y)
      gvect <- function(y) mapply(g, y)
      r.int.0 <- function(y) exp(gvect(y)) * f.lamb(y)
      r.int.int <- function(x, y) .Integrate(r.int.0, x, y, stop.on.error = FALSE)
      return(exp(r.int(s, t)) * abs(1 + r.int.int(s, t) / (1 / f + r.int.int(0, s)))^(-2))
    }

    if (expo.lamb && expo.mu) {
      lamb0 <- f.lamb(0)
      alpha <- log(f.lamb(1) / lamb0)
      mu0 <- f.mu(0)
      beta <- log(f.mu(1) / mu0)
      r.int <- function(x, y) {
        lamb0 / alpha * (exp(alpha * y) - exp(alpha * x)) -
          mu0 / beta * (exp(beta * y) - exp(beta * x))
      }
      g <- function(y) r.int(0, y)
      gvect <- function(y) mapply(g, y)
      r.int.0 <- function(y) exp(gvect(y)) * f.lamb(y)
      r.int.int <- function(x, y) .Integrate(r.int.0, x, y, stop.on.error = FALSE)
      return(exp(r.int(s, t)) * abs(1 + r.int.int(s, t) / (1 / f + r.int.int(0, s)))^(-2))
    }

    ageMin <- s
    ageMax <- t
    Nintervals <- 1L + as.integer((ageMax - ageMin) / dt)
    X <- seq(ageMin, ageMax, length.out = Nintervals + 1L)
    r <- function(z) f.lamb(z) - f.mu(z)
    r.int <- cumsum(r(X)) * (ageMax - ageMin) / Nintervals
    r.int.0 <- function(y) {
      index <- 1L + as.integer((y - ageMin) * Nintervals / (ageMax - ageMin))
      exp(r.int[index]) * f.lamb(y)
    }
    r.int.int.tab <- cumsum(r.int.0(X)) * (ageMax - ageMin) / Nintervals
    r.int.int <- function(x, y) {
      indy <- 1L + as.integer((y - ageMin) * Nintervals / (ageMax - ageMin))
      indx <- 1L + as.integer((x - ageMin) * Nintervals / (ageMax - ageMin))
      r.int.int.tab[indy] - r.int.int.tab[indx]
    }
    rst <- r.int[1L + Nintervals]
    rist <- r.int.int(s, t)
    ri0s <- r.int.int(0, s)
    exp(rst) * abs(1 + rist / (1 / f + ri0s))^(-2)
  }
}
