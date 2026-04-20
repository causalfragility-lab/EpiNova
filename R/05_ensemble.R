#' @title Ensemble Forecasting and Uncertainty Quantification
#'
#' @description
#' eSIR produced credible intervals only for its fixed model class.
#' EpiNova provides:
#' \itemize{
#'   \item Model-averaged forecasts (Bayesian model averaging over
#'     different compartmental structures).
#'   \item Scenario-based projection (what-if analysis).
#'   \item Proper scoring rules to evaluate forecast calibration.
#' }
#' @name EpiNova-ensemble
#' @importFrom stats median quantile setNames rnorm
NULL

utils::globalVariables(c('median', 'quantile'))

# ============================================================
# IDEA 13: Bayesian Model Averaging across compartmental models
# ============================================================

#' Ensemble forecast via Bayesian Model Averaging (BMA)
#'
#' Fits multiple compartmental models and combines their forecasts
#' weighted by their marginal likelihoods (approximated by AIC weights).
#'
#' @param obs_Y      Numeric vector of observed infected proportions.
#' @param obs_R      Numeric vector of observed removed proportions.
#' @param N          Population size.
#' @param models     Character vector of model types to average over.
#' @param par_inits  Named list of initial parameter lists (one per model).
#' @param par_bounds List with elements \code{lower} and \code{upper},
#'   each a named list of bounds applying to all models.
#' @param pi_fn      Common transmission modifier (or \code{NULL}).
#' @param T_forecast Integer number of days to forecast ahead.
#' @param n_sim      Number of simulation draws per model.
#'
#' @return A list with \code{ensemble_forecast} (data frame),
#'   \code{model_weights}, and \code{individual_fits}.
#' @export
ensemble_forecast <- function(obs_Y, obs_R, N,
                               models = c("SEIR","SEIRD","SVEIRD"),
                               par_inits, par_bounds,
                               pi_fn = NULL,
                               T_forecast = 60L,
                               n_sim = 500L) {

  fits <- vector("list", length(models))
  aics <- numeric(length(models))

  for (i in seq_along(models)) {
    message("Fitting model: ", models[i])
    fits[[i]] <- tryCatch(
      fit_mle(obs_Y, obs_R, N,
              model     = models[i],
              par_init  = par_inits[[i]],
              par_lower = par_bounds$lower,
              par_upper = par_bounds$upper,
              pi_fn     = pi_fn),
      error = function(e) { message("  Failed: ", e$message); NULL }
    )
    aics[i] <- if (!is.null(fits[[i]])) fits[[i]]$AIC else Inf
  }

  # AIC weights
  delta_aic <- aics - min(aics)
  w_raw     <- exp(-0.5 * delta_aic)
  weights   <- w_raw / sum(w_raw)

  # Simulate from each model
  T_total <- length(obs_Y) + T_forecast
  times   <- seq(0, T_total - 1)

  all_sims <- lapply(seq_along(models), function(i) {
    if (is.null(fits[[i]])) return(NULL)
    fit <- fits[[i]]
    .parametric_bootstrap(fit, times, n_sim, pi_fn)
  })

  # Weighted ensemble quantiles
  ensemble <- .bma_quantiles(all_sims, weights, T_total, n_sim)

  list(
    ensemble_forecast = ensemble,
    model_weights     = setNames(weights, models),
    individual_fits   = setNames(fits, models)
  )
}


# ============================================================
# IDEA 14: Scenario projections ("what-if" analysis)
#   User defines alternative intervention scenarios and compares
#   projected outcomes — not available in eSIR.
# ============================================================

#' Project scenarios under alternative intervention strategies
#'
#' @param fit        A fitted model object (\code{EpiNova_mle}).
#' @param scenarios  Named list of π(t) functions, one per scenario.
#' @param T_forecast Integer forecast horizon in days.
#' @param n_sim      Simulation draws for uncertainty.
#'
#' @return A tidy data frame with columns \code{scenario},
#'   \code{time}, \code{I_median}, \code{I_lower}, \code{I_upper},
#'   \code{peak_day}, \code{peak_I}, \code{total_infected}.
#' @export
project_scenarios <- function(fit, scenarios, T_forecast = 120L,
                               n_sim = 500L) {

  T_obs   <- length(fit$obs_Y)
  times   <- seq(0, T_obs + T_forecast - 1)

  results <- lapply(names(scenarios), function(sc_name) {
    pi_fn_sc <- scenarios[[sc_name]]
    sims     <- .parametric_bootstrap(fit, times, n_sim, pi_fn_sc)

    I_mat <- do.call(rbind, lapply(sims, `[[`, "I"))

    data.frame(
      scenario    = sc_name,
      time        = times,
      I_median    = apply(I_mat, 2, median),
      I_lower     = apply(I_mat, 2, quantile, 0.025),
      I_upper     = apply(I_mat, 2, quantile, 0.975),
      peak_day    = times[which.max(apply(I_mat, 2, median))],
      total_inf   = apply(I_mat, 2, function(x) sum(x) * fit$N)
    )
  })

  do.call(rbind, results)
}


# ============================================================
# IDEA 15: Proper scoring rules for forecast evaluation
# ============================================================

#' Evaluate forecast calibration with proper scoring rules
#'
#' Computes the Continuous Ranked Probability Score (CRPS) and
#' interval coverage for each time point.
#'
#' @param forecast_df Data frame with columns \code{I_median},
#'   \code{I_lower}, \code{I_upper} (from \code{project_scenarios}
#'   or \code{ensemble_forecast}).
#' @param actual_Y    Numeric vector of actual observed values.
#'
#' @return Data frame with \code{CRPS}, \code{coverage_50},
#'   \code{coverage_95}, and \code{MAE}.
#' @export
score_forecast <- function(forecast_df, actual_Y) {
  n   <- length(actual_Y)
  med <- forecast_df$I_median[seq_len(n)]
  lo  <- forecast_df$I_lower[seq_len(n)]
  hi  <- forecast_df$I_upper[seq_len(n)]

  # Gaussian approximation for CRPS
  sigma_approx <- (hi - lo) / (2 * 1.96)
  crps <- sigma_approx * (
    (actual_Y - med) / sigma_approx *
      (2 * pnorm((actual_Y - med) / sigma_approx) - 1) +
      2 * dnorm((actual_Y - med) / sigma_approx) - 1 / sqrt(pi)
  )

  data.frame(
    CRPS        = mean(crps, na.rm = TRUE),
    coverage_95 = mean(actual_Y >= lo & actual_Y <= hi, na.rm = TRUE),
    MAE         = mean(abs(actual_Y - med), na.rm = TRUE)
  )
}


# --------------------------------------------------------
# Internal helpers
# --------------------------------------------------------

.parametric_bootstrap <- function(fit, times, n_sim, pi_fn) {
  # Draw parameters from multivariate normal around MLE
  if (all(is.finite(fit$vcov))) {
    par_draws <- if (requireNamespace("MASS", quietly = TRUE)) {
      MASS::mvrnorm(n_sim, fit$par, fit$vcov)
    } else {
      matrix(fit$par, nrow = n_sim, ncol = length(fit$par), byrow = TRUE) + 
        matrix(rnorm(n_sim * length(fit$par), 0, 1e-4), n_sim)
    }
  } else {
    par_draws <- matrix(fit$par, nrow = n_sim, ncol = length(fit$par),
                        byrow = TRUE)
  }
  par_draws <- pmax(par_draws, 1e-8)

  lapply(seq_len(n_sim), function(j) {
    params <- as.list(par_draws[j, ])
    names(params) <- names(fit$par)
    init   <- .default_init(fit$model, params)
    tryCatch(
      solve_model(params, init, times,
                  model = fit$model, pi_fn = pi_fn),
      error = function(e) NULL
    )
  })
}

.bma_quantiles <- function(all_sims, weights, T_total, n_sim) {
  # Weighted quantile across ensemble
  I_list <- lapply(seq_along(all_sims), function(i) {
    sims <- all_sims[[i]]
    if (is.null(sims)) return(NULL)
    w <- weights[i]
    valid <- Filter(Negate(is.null), sims)
    mat   <- do.call(rbind, lapply(valid, function(s) s$I))
    list(mat = mat, w = w)
  })
  I_list <- Filter(Negate(is.null), I_list)

  # Stack weighted samples
  I_all <- do.call(rbind, lapply(I_list, function(x) x$mat * x$w))

  data.frame(
    time     = seq(0, T_total - 1),
    I_median = apply(I_all, 2, median),
    I_lower  = apply(I_all, 2, quantile, 0.025),
    I_upper  = apply(I_all, 2, quantile, 0.975)
  )
}
