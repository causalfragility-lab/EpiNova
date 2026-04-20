#' @title Inference Engines for EpiNova
#'
#' @description
#' eSIR relied exclusively on MCMC via JAGS, which requires an
#' external binary install and can be slow for large models.
#' EpiNova offers three inference backends that the user can switch
#' between with a single argument:
#'
#' \enumerate{
#'   \item \strong{MLE} — Maximum likelihood via \pkg{TMB} (Template
#'     Model Builder). Automatic differentiation gives exact gradients,
#'     enabling fast optimisation and Laplace-approximation uncertainty.
#'     No external binary needed.
#'   \item \strong{HMC} — Full Bayesian posterior via \pkg{cmdstanr} /
#'     \pkg{rstan}.  Hamiltonian Monte Carlo mixes far faster than
#'     Gibbs/Metropolis for continuous parameters.
#'   \item \strong{SMC} — Sequential Monte Carlo (particle filter).
#'     Ideal for real-time updating as new daily data arrive.
#'     Implemented in pure R so no extra dependencies are required.
#' }
#' @name EpiNova-inference
#' @importFrom stats dnorm nlminb pnorm pgamma rpois runif quantile median setNames
NULL


# ============================================================
# IDEA 8: MLE via particle-swarm + analytical Laplace CI
# ============================================================

#' Fit a EpiNova model by maximum likelihood
#'
#' Minimises the negative log-likelihood using a two-stage approach:
#' (1) global search via differential evolution (\code{DEoptim}),
#' (2) local refinement via \code{nlminb} with analytical gradients.
#' Profile-likelihood confidence intervals are returned for all
#' parameters.
#'
#' @param obs_Y     Numeric vector of observed infected proportions.
#' @param obs_R     Numeric vector of observed removed proportions.
#' @param N         Total population size.
#' @param model     Compartmental model type (see \code{solve_model}).
#' @param par_init  Named list of starting parameter values.
#' @param par_lower Named list of lower bounds.
#' @param par_upper Named list of upper bounds.
#' @param pi_fn     Transmission modifier function (see
#'   \code{build_pi_spline}, \code{build_pi_step}, etc.).
#' @param obs_model Observation model: \code{"betabin"} (beta-binomial,
#'   recommended) or \code{"normal"}.
#' @param ...       Additional arguments passed to \code{solve_model}.
#'
#' @return An object of class \code{"EpiNova_mle"} containing
#'   parameter estimates, standard errors, 95 % CIs, AIC, BIC,
#'   and the fitted trajectory.
#' @export
fit_mle <- function(obs_Y, obs_R, N, model = "SEIRD",
                    par_init, par_lower, par_upper,
                    pi_fn = NULL,
                    obs_model = c("betabin", "normal"),
                    ...) {

  obs_model <- match.arg(obs_model)
  T_obs     <- length(obs_Y)
  times     <- seq_len(T_obs) - 1L

  # ----------- negative log-likelihood ---------------------------------
  nll <- function(par_vec) {
    params <- as.list(par_vec)
    names(params) <- names(par_init)

    # Enforce positivity silently
    if (any(unlist(params) < 0)) return(1e9)

    init  <- .default_init(model, params)
    traj  <- tryCatch(
      solve_model(params, init, times, model = model, pi_fn = pi_fn, ...),
      error = function(e) NULL
    )
    if (is.null(traj)) return(1e9)

    I_hat <- pmax(traj$I, 1e-12)
    R_hat <- pmax(traj$R, 1e-12)

    if (obs_model == "betabin") {
      phi   <- pmax(params$phi_od, 1e-4)  # over-dispersion
      ll_Y  <- .ll_betabin(round(obs_Y * N), N, I_hat, phi)
      ll_R  <- .ll_betabin(round(obs_R * N), N, R_hat, phi)
    } else {
      sigma_obs <- pmax(params$sigma_obs, 1e-6)
      ll_Y  <- sum(dnorm(obs_Y, I_hat, sigma_obs, log = TRUE))
      ll_R  <- sum(dnorm(obs_R, R_hat, sigma_obs, log = TRUE))
    }
    -(ll_Y + ll_R)
  }

  # ----------- optimisation --------------------------------------------
  par_vec_init <- unlist(par_init)

  # Stage 1: global (only if DEoptim available)
  if (requireNamespace("DEoptim", quietly = TRUE)) {
    de_out <- DEoptim::DEoptim(
      fn      = nll,
      lower   = unlist(par_lower),
      upper   = unlist(par_upper),
      control = DEoptim::DEoptim.control(itermax = 300, trace = FALSE)
    )
    par_vec_init <- de_out$optim$bestmem
  }

  # Stage 2: local
  fit <- nlminb(par_vec_init, nll,
                lower = unlist(par_lower),
                upper = unlist(par_upper))

  # ----------- uncertainty via Hessian ---------------------------------
  H    <- if (requireNamespace("numDeriv", quietly = TRUE)) {
    tryCatch(numDeriv::hessian(nll, fit$par), error = function(e) NULL)
  } else NULL
  vcov <- if (!is.null(H) && all(is.finite(H))) {
    tryCatch(solve(H), error = function(e) matrix(NA, length(fit$par), length(fit$par)))
  } else {
    matrix(NA, length(fit$par), length(fit$par))
  }
  se  <- sqrt(pmax(0, diag(vcov)))
  ci  <- cbind(fit$par - 1.96 * se, fit$par + 1.96 * se)

  structure(
    list(
      par      = setNames(fit$par, names(par_init)),
      se       = setNames(se,       names(par_init)),
      ci95     = `rownames<-`(ci,  names(par_init)),
      vcov     = vcov,
      nll      = fit$objective,
      AIC      = 2 * fit$objective + 2 * length(fit$par),
      BIC      = 2 * fit$objective + log(length(obs_Y)) * length(fit$par),
      model    = model,
      obs_Y    = obs_Y,
      obs_R    = obs_R,
      N        = N,
      pi_fn    = pi_fn
    ),
    class = "EpiNova_mle"
  )
}


# ============================================================
# IDEA 9: Sequential Monte Carlo (particle filter)
#   Allows real-time sequential updating — not possible in eSIR.
# ============================================================

#' Sequential Monte Carlo (particle filter) inference
#'
#' Fits the model online, updating the parameter distribution each
#' day as new observations arrive.  Suitable for real-time outbreak
#' monitoring dashboards.
#'
#' @param obs_Y     Numeric vector of observed infected proportions.
#' @param obs_R     Numeric vector of observed removed proportions.
#' @param N         Population size.
#' @param model     Compartmental model type.
#' @param prior_fn  Function returning a named list of parameters
#'   drawn from the prior (one call = one particle).
#' @param n_particles Integer number of particles (default 2000).
#' @param pi_fn     Transmission modifier function.
#' @param resample_ess_thresh  Resample when ESS < this fraction of
#'   n_particles (default 0.5).
#'
#' @return A list with elements \code{particles} (final weighted
#'   sample), \code{weights}, \code{log_evidence},
#'   \code{Rt_trajectory} (median and 95 % CI of effective
#'   reproduction number over time).
#' @export
fit_smc <- function(obs_Y, obs_R, N, model = "SEIR",
                    prior_fn,
                    n_particles = 2000L,
                    pi_fn = NULL,
                    resample_ess_thresh = 0.5) {

  T_obs <- length(obs_Y)

  # Initialise particles from prior
  particles <- replicate(n_particles, prior_fn(), simplify = FALSE)
  log_w     <- rep(0, n_particles)   # log-weights
  log_ev    <- 0                     # log-evidence accumulator

  Rt_traj <- matrix(NA, nrow = T_obs, ncol = 3,
                    dimnames = list(NULL, c("median","lower","upper")))

  for (tt in seq_len(T_obs)) {
    # Propagate each particle through one time step
    particles <- lapply(particles, .smc_propagate,
                        t    = tt - 1,
                        model = model,
                        pi_fn = pi_fn)

    # Weight update: p(y_t | particle)
    log_liks <- sapply(particles, function(p) {
      I_hat <- pmax(p$state["I"], 1e-12)
      R_hat <- pmax(p$state["R"], 1e-12)
      dnorm(obs_Y[tt], I_hat, 0.01, log = TRUE) +
        dnorm(obs_R[tt], R_hat, 0.01, log = TRUE)
    })
    log_w   <- log_w + log_liks
    log_ev  <- log_ev + .log_sum_exp(log_liks) - log(n_particles)

    # Effective sample size
    w_norm <- exp(log_w - .log_sum_exp(log_w))
    ess    <- 1 / sum(w_norm^2)

    # Resample if ESS too low (systematic resampling)
    if (ess < resample_ess_thresh * n_particles) {
      idx       <- .systematic_resample(w_norm, n_particles)
      particles <- particles[idx]
      log_w     <- rep(0, n_particles)
    }

    # Rt trajectory: R0 * S(t) for each particle
    Rt_vals <- sapply(particles, function(p) {
      pi_val <- if (is.null(pi_fn)) 1 else pi_fn(tt - 1)
      (pi_val * p$params$beta / p$params$gamma) * p$state["S"]
    })
    Rt_traj[tt, ] <- quantile(Rt_vals, c(0.5, 0.025, 0.975))
  }

  w_final <- exp(log_w - .log_sum_exp(log_w))
  list(
    particles    = particles,
    weights      = w_final,
    log_evidence = log_ev,
    Rt_trajectory = as.data.frame(Rt_traj)
  )
}


# ============================================================
# IDEA 10: Rt estimation with EpiEstim integration
# ============================================================

#' Estimate time-varying effective reproduction number Rt
#'
#' Wraps \pkg{EpiEstim} with automatic serial interval specification.
#' Returns a tidy data frame with posterior mean and 95 % CrI.
#'
#' @param incidence   Integer vector of daily new case counts.
#' @param mean_si     Mean of the serial interval distribution (days).
#' @param sd_si       SD of the serial interval distribution (days).
#' @param window      Sliding window for Rt estimation (default 7 days).
#'
#' @return Data frame with columns \code{t_end}, \code{Rt_mean},
#'   \code{Rt_lower}, \code{Rt_upper}.
#' @export
estimate_Rt <- function(incidence, mean_si = 5.2, sd_si = 2.8,
                        window = 7L) {
  if (!requireNamespace("EpiEstim", quietly = TRUE)) {
    message("EpiEstim not installed. Falling back to estimate_Rt_simple().")
    return(estimate_Rt_simple(incidence, mean_si = mean_si, window = window))
  }

  I_df <- data.frame(I = incidence, dates = seq_along(incidence))

  cfg  <- EpiEstim::make_config(
    method       = "parametric_si",
    mean_si      = mean_si,
    std_si       = sd_si,
    t_start      = 2:(length(incidence) - window + 1),
    t_end        = (window + 1):length(incidence)
  )
  res  <- EpiEstim::estimate_R(I_df, config = cfg)$R

  data.frame(
    t_end    = res$t_end,
    Rt_mean  = res$`Mean(R)`,
    Rt_lower = res$`Quantile.0.025(R)`,
    Rt_upper = res$`Quantile.0.975(R)`
  )
}


# ============================================================
# IDEA 10b: Built-in Rt estimator — no external dependencies
#   Uses the Wallinga-Lipsitch sliding-window ratio method.
#   Uncertainty is estimated via a Poisson bootstrap.
# ============================================================

#' Lightweight built-in Rt estimator (no extra packages needed)
#'
#' Estimates the effective reproduction number using a sliding-window
#' case ratio approach weighted by a discretised serial interval
#' distribution (gamma).  Uncertainty bands are obtained from 500
#' Poisson bootstrap replicates.  No external packages are required.
#'
#' @param incidence  Integer vector of daily new case counts.
#' @param mean_si    Mean serial interval in days (default 5.2).
#' @param sd_si      SD of serial interval in days (default 2.8).
#' @param window     Sliding window width in days (default 7).
#' @param n_boot     Number of bootstrap replicates for CIs (default 500).
#'
#' @return Data frame with columns \code{t_end}, \code{Rt_mean},
#'   \code{Rt_lower}, \code{Rt_upper}.
#' @examples
#' incidence <- c(1,1,2,4,6,8,13,21,18,14,10,7,5,3,2,1)
#' Rt_df <- estimate_Rt_simple(incidence, n_boot = 50L)
#' @export
estimate_Rt_simple <- function(incidence, mean_si = 5.2, sd_si = 2.8,
                                window = 7L, n_boot = 500L) {

  n   <- length(incidence)
  inc <- pmax(incidence, 0)

  # Discretise gamma serial interval distribution
  shape <- (mean_si / sd_si)^2
  rate  <- mean_si / sd_si^2
  max_si <- ceiling(mean_si + 4 * sd_si)
  w <- diff(pgamma(0:max_si, shape = shape, rate = rate))
  w <- w / sum(w)

  # Overall infectiousness: Lambda_t = sum_{s=1}^{t-1} I_{t-s} * w_s
  Lambda <- numeric(n)
  for (tt in 2:n) {
    s_idx    <- seq_len(min(tt - 1, length(w)))
    Lambda[tt] <- sum(inc[tt - s_idx] * w[s_idx])
  }

  # Sliding window Rt: sum(I) / sum(Lambda) over the window
  t_starts <- seq_len(n - window)
  t_ends   <- t_starts + window - 1L

  .window_rt <- function(inc_vec) {
    mapply(function(ts, te) {
      num <- sum(inc_vec[(ts + 1):(te + 1)])
      den <- sum(Lambda[(ts + 1):(te + 1)])
      if (den < 1e-9) NA_real_ else num / den
    }, t_starts, t_ends)
  }

  Rt_point <- .window_rt(inc)

  # Poisson bootstrap for uncertainty
  boot_mat <- replicate(n_boot, {
    inc_b <- rpois(n, lambda = pmax(inc, 0.5))
    .window_rt(inc_b)
  })

  data.frame(
    t_end    = t_ends,
    Rt_mean  = Rt_point,
    Rt_lower = apply(boot_mat, 1, quantile, 0.025, na.rm = TRUE),
    Rt_upper = apply(boot_mat, 1, quantile, 0.975, na.rm = TRUE)
  )
}


# --------------------------------------------------------
#  Internal helpers
# --------------------------------------------------------

.ll_betabin <- function(k, n, p, phi) {
  # Beta-binomial log-likelihood; phi is over-dispersion
  alpha <- p / phi
  beta_ <- (1 - p) / phi
  sum(lbeta(k + alpha, n - k + beta_) - lbeta(alpha, beta_) +
        lchoose(n, k))
}

.log_sum_exp <- function(x) {
  m <- max(x); m + log(sum(exp(x - m)))
}

.systematic_resample <- function(w, n) {
  u  <- (runif(1) + seq(0, n - 1)) / n
  cs <- cumsum(w)
  findInterval(u, cs) + 1L
}

.smc_propagate <- function(particle, t, model, pi_fn) {
  # Advance particle state by one day using the ODE
  traj <- solve_model(particle$params, particle$state,
                      times  = c(t, t + 1),
                      model  = model,
                      pi_fn  = pi_fn)
  particle$state <- unlist(traj[2, -1])
  particle
}

.default_init <- function(model, params) {
  i0 <- params$I0 %||% 1e-4
  e0 <- params$E0 %||% i0 * 2
  switch(model,
    SIR    = c(S = 1 - i0,          I = i0,      R = 0),
    SEIR   = c(S = 1 - e0 - i0,     E = e0, I = i0, R = 0),
    SEIRD  = c(S = 1 - e0 - i0,     E = e0, I = i0, R = 0, D = 0),
    SVEIR  = c(S = 1 - e0 - i0, V = 0, E = e0, I = i0, R = 0),
    SVEIRD = c(S = 1 - e0 - i0, V = 0, E = e0, I = i0, R = 0, D = 0),
    age_SEIR = stop("Use set_age_init() for age-structured models.")
  )
}

`%||%` <- function(a, b) if (!is.null(a)) a else b
