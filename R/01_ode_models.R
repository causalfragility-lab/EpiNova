#' @title EpiNova: Flexible Epidemiological Compartmental Models
#' @description Core ODE system for multiple compartmental model types
#' @name EpiNova-ode
#' @importFrom deSolve ode
#' @importFrom dplyr mutate select
#' @importFrom stats predict setNames
NULL

# ============================================================
# IDEA 1: Unified compartmental model dispatcher
#   eSIR only had SIR. EpiNova supports a hierarchy of models
#   chosen by a single `model` argument.
# ============================================================

#' Solve a compartmental ODE model
#'
#' Dispatches to the appropriate ODE system based on the model type.
#' Supports SIR, SEIR, SEIRD, SVEIR, SVEIRD, and age-stratified
#' variants.  All models share a common calling convention so that
#' downstream functions (fitting, forecasting, plotting) are
#' model-agnostic.
#'
#' @param params  Named list of epidemiological parameters.
#' @param init    Named numeric vector of initial compartment
#'   proportions (must sum to 1).
#' @param times   Numeric vector of time points (days).
#' @param model   Character string: one of \code{"SIR"}, \code{"SEIR"},
#'   \code{"SEIRD"}, \code{"SVEIR"}, \code{"SVEIRD"},
#'   \code{"age_SEIR"}.
#' @param pi_fn   Optional function \code{pi_fn(t)} returning a
#'   transmission modifier in \eqn{[0,1]} at time \code{t}.  If \code{NULL},
#'   the unmodified transmission rate is used.
#' @param phi_fn  Optional function \code{phi_fn(t)} returning the
#'   instantaneous quarantine removal rate at time \code{t} (Dirac
#'   delta approximated by a narrow Gaussian pulse). If \code{NULL},
#'   no quarantine pulse is applied.
#'
#' @return A data frame with one row per time point and one column per
#'   compartment, plus a column \code{time}.
#'
#' @examples
#' p <- list(beta = 0.3, gamma = 0.1, sigma = 0.2,
#'           delta = 0.005, vax_rate = 0.002, wane = 0.01, ve = 0.85)
#' y0 <- c(S = 0.990, E = 0.005, I = 0.004, R = 0.001, D = 0, V = 0)
#' t  <- seq(0, 200, by = 1)
#' result <- solve_model(p, y0, t, model = "SVEIRD")
#'
#' @export
solve_model <- function(params, init, times,
                        model  = c("SIR","SEIR","SEIRD",
                                   "SVEIR","SVEIRD","age_SEIR"),
                        pi_fn  = NULL,
                        phi_fn = NULL) {

  model <- match.arg(model)

  # Default identity modifier if none supplied
  if (is.null(pi_fn))  pi_fn  <- function(t) 1
  if (is.null(phi_fn)) phi_fn <- function(t) 0

  ode_fn <- switch(model,
    SIR     = .ode_SIR,
    SEIR    = .ode_SEIR,
    SEIRD   = .ode_SEIRD,
    SVEIR   = .ode_SVEIR,
    SVEIRD  = .ode_SVEIRD,
    age_SEIR = .ode_age_SEIR
  )

  out <- deSolve::ode(
    y     = init,
    times = times,
    func  = ode_fn,
    parms = c(params, list(pi_fn = pi_fn, phi_fn = phi_fn)),
    method = "lsoda"
  )

  as.data.frame(out)
}


# ============================================================
# IDEA 2: SVEIRD – adds Vaccination (V) and explicit Death (D)
#   compartments.  eSIR lumped deaths into "Removed".
#   Waning vaccine immunity (V -> S) and direct vaccination
#   of susceptibles (S -> V) are both modelled.
# ============================================================

.ode_SVEIRD <- function(t, state, parms) {
  with(as.list(c(state, parms)), {

    pi_t  <- pi_fn(t)   # transmission modifier in [0,1]
    phi_t <- phi_fn(t)  # quarantine pulse

    force_of_infection <- pi_t * beta * I

    dS <- -force_of_infection * S - vax_rate * S + wane * V - phi_t * S
    dV <-  vax_rate * S - wane * V - force_of_infection * V * (1 - ve)
    dE <-  force_of_infection * (S + V * (1 - ve)) - sigma * E
    dI <-  sigma * E - gamma * I - delta * I
    dR <-  gamma * I
    dD <-  delta * I

    list(c(dS, dV, dE, dI, dR, dD))
  })
}

# ve  = vaccine efficacy against infection
# delta = disease-induced mortality rate
# wane  = rate of waning vaccine-induced immunity


.ode_SIR <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    pi_t  <- pi_fn(t)
    phi_t <- phi_fn(t)
    lambda <- pi_t * beta * I
    dS <- -lambda * S - phi_t * S
    dI <-  lambda * S - gamma * I
    dR <-  gamma * I + phi_t * S
    list(c(dS, dI, dR))
  })
}

.ode_SEIR <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    pi_t  <- pi_fn(t)
    phi_t <- phi_fn(t)
    lambda <- pi_t * beta * I
    dS <- -lambda * S - phi_t * S
    dE <-  lambda * S - sigma * E
    dI <-  sigma * E - gamma * I
    dR <-  gamma * I + phi_t * S
    list(c(dS, dE, dI, dR))
  })
}

.ode_SEIRD <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    pi_t  <- pi_fn(t)
    phi_t <- phi_fn(t)
    lambda <- pi_t * beta * I
    dS <- -lambda * S - phi_t * S
    dE <-  lambda * S - sigma * E
    dI <-  sigma * E - (gamma + delta) * I
    dR <-  gamma * I + phi_t * S
    dD <-  delta * I
    list(c(dS, dE, dI, dR, dD))
  })
}

.ode_SVEIR <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    pi_t  <- pi_fn(t)
    phi_t <- phi_fn(t)
    lambda <- pi_t * beta * I
    dS <- -lambda * S - vax_rate * S + wane * V - phi_t * S
    dV <-  vax_rate * S - wane * V - lambda * V * (1 - ve)
    dE <-  lambda * (S + V * (1 - ve)) - sigma * E
    dI <-  sigma * E - gamma * I
    dR <-  gamma * I + phi_t * S
    list(c(dS, dV, dE, dI, dR))
  })
}


# ============================================================
# IDEA 3: Age-stratified SEIR
#   n_age age groups, contact matrix C, group-specific params
# ============================================================

#' @keywords internal
.ode_age_SEIR <- function(t, state, parms) {
  with(parms, {
    n <- n_age
    S <- state[1:n]
    E <- state[(n+1):(2*n)]
    I <- state[(2*n+1):(3*n)]
    R <- state[(3*n+1):(4*n)]

    pi_t <- pi_fn(t)

    # Force of infection: C is the n x n contact matrix
    lambda <- pi_t * beta_vec * as.vector(C %*% I)

    dS <- -lambda * S
    dE <-  lambda * S - sigma_vec * E
    dI <-  sigma_vec * E - gamma_vec * I
    dR <-  gamma_vec * I

    list(c(dS, dE, dI, dR))
  })
}
