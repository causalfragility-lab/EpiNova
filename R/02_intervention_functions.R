#' @title Flexible Intervention (pi) Functions
#'
#' @description
#' eSIR restricted pi(t) to step functions or a single decaying
#' exponential.  EpiNova provides a richer library of transmission
#' modifiers and allows users to compose them or supply a completely
#' custom function.
#'
#' All builders return a plain R function \code{pi_fn(t)} in \eqn{[0,1]}.
#' @name EpiNova-interventions
#' @importFrom stats dnorm predict
#' @importFrom splines interpSpline
NULL

# ============================================================
# IDEA 4: Spline-based π(t)
#   Fits a monotone or unconstrained natural cubic spline through
#   user-supplied control knots.  Far more realistic than a step
#   function for policies that are phased in/out gradually.
# ============================================================

#' Build a natural cubic spline π(t)
#'
#' @param knot_times  Numeric vector of days (knot positions).
#' @param knot_values Numeric vector of pi values at each knot,
#'   all in \eqn{[0, 1]}.
#' @param extrapolate  How to handle t outside knot range:
#'   \code{"flat"} (default) or \code{"linear"}.
#'
#' @return A function \code{pi_fn(t)}.
#' @examples
#' pi_spline <- build_pi_spline(
#'   knot_times  = c(0, 10, 20, 40, 100),
#'   knot_values = c(1, 0.8, 0.5, 0.3, 0.3)
#' )
#' curve(pi_spline(x), 0, 120, ylab = expression(pi(t)))
#' @export
build_pi_spline <- function(knot_times, knot_values,
                            extrapolate = c("flat","linear")) {
  stopifnot(length(knot_times) == length(knot_values),
            all(knot_values >= 0 & knot_values <= 1))
  extrapolate <- match.arg(extrapolate)

  sp <- splines::interpSpline(knot_times, knot_values)

  function(t) {
    t_clamped <- if (extrapolate == "flat") {
      pmax(min(knot_times), pmin(max(knot_times), t))
    } else {
      t
    }
    val <- predict(sp, t_clamped)$y
    pmax(0, pmin(1, val))  # hard-clip to [0,1]
  }
}


# ============================================================
# IDEA 5: Gaussian Process (GP) prior on π(t) – for Bayesian fitting
#   Instead of fixing the shape, we let the data inform it via a GP.
#   This is used internally by the Stan fitting module.
#   Exposed here as a covariance builder for transparency.
# ============================================================

#' Build a squared-exponential covariance matrix for GP π(t)
#'
#' Used to construct the GP prior over transmission modifiers at
#' a discrete set of time points before passing to the Stan sampler.
#'
#' @param times  Numeric vector of time points.
#' @param l      Length-scale (controls smoothness).
#' @param sigma  Marginal standard deviation of the GP.
#'
#' @return A symmetric positive-definite matrix of dimension
#'   \code{length(times) x length(times)}.
#' @export
gp_cov_sqexp <- function(times, l = 14, sigma = 0.3) {
  D <- outer(times, times, FUN = function(a, b) (a - b)^2)
  sigma^2 * exp(-D / (2 * l^2)) +
    diag(1e-6, length(times))  # nugget for numerical stability
}


# ============================================================
# IDEA 6: Composite intervention – multiplicative combination
#   of multiple π functions (e.g., mask mandate × lockdown).
# ============================================================

#' Compose multiple π(t) functions multiplicatively
#'
#' Each component function represents an independent non-pharmaceutical
#' intervention (NPI).  Their effects are assumed multiplicative on
#' transmission.
#'
#' @param ...  Any number of functions \code{pi_fn(t)}, each in \eqn{[0,1]}.
#'
#' @return A composite function \code{pi_fn(t)} in \eqn{[0,1]}.
#' @examples
#' lockdown  <- build_pi_step(c(10, 60), c(1, 0.4, 0.7))
#' masks     <- build_pi_spline(c(0, 10, 20, 30), c(1, 0.92, 0.85, 0.75))
#' combined  <- compose_pi(lockdown, masks)
#' @export
compose_pi <- function(...) {
  fns <- list(...)
  stopifnot(all(sapply(fns, is.function)))
  function(t) {
    vals <- sapply(fns, function(f) f(t))
    prod(vals)
  }
}


#' Build a step-function π(t)  (reproduces eSIR Model 1 behaviour)
#'
#' @param change_times Numeric vector of change-point days.
#' @param pi_values    Numeric vector of length
#'   \code{length(change_times) + 1}.
#'
#' @return A function \code{pi_fn(t)}.
#' @export
build_pi_step <- function(change_times, pi_values) {
  stopifnot(length(pi_values) == length(change_times) + 1)
  function(t) {
    idx <- findInterval(t, change_times) + 1
    pi_values[idx]
  }
}


#' Build an exponential decay π(t) = exp(−λ t)
#'
#' @param lambda Decay rate (positive scalar).
#' @param t0     Start of decay (default 0).
#' @return A function \code{pi_fn(t)}.
#' @export
build_pi_exp <- function(lambda, t0 = 0) {
  function(t) pmax(0, exp(-lambda * pmax(0, t - t0)))
}


# ============================================================
# IDEA 7: φ(t) as a smooth Gaussian pulse instead of Dirac delta
#   eSIR's quarantine jumps were hard discrete steps.
#   A Gaussian pulse gives continuous differentiability for
#   gradient-based samplers (HMC in Stan).
# ============================================================

#' Build a smooth quarantine pulse φ(t)
#'
#' Approximates a Dirac delta at each change point with a narrow
#' Gaussian pulse.  The area under each pulse equals the quarantine
#' fraction \code{phi_values[i]}.
#'
#' @param change_times  Numeric vector of quarantine event days.
#' @param phi_values    Quarantine fractions in (0, 1).
#' @param bandwidth     Width (SD) of each Gaussian pulse (default 0.5 days).
#'
#' @return A function \code{phi_fn(t)}.
#' @export
build_phi_pulse <- function(change_times, phi_values,
                            bandwidth = 0.5) {
  stopifnot(length(change_times) == length(phi_values))
  function(t) {
    pulses <- mapply(function(mu, h) {
      h * dnorm(t, mean = mu, sd = bandwidth) /
        dnorm(0, mean = 0, sd = bandwidth)   # normalise peak to h
    }, change_times, phi_values)
    if (is.matrix(pulses)) rowSums(pulses) else sum(pulses)
  }
}
