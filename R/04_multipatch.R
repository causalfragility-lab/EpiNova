#' @title Multi-patch Spatial SIR Models
#'
#' @description
#' eSIR was a single-population model.  EpiNova adds a network of
#' patches (cities, regions, countries) connected by a mobility /
#' commuting matrix.  Each patch has its own compartmental dynamics;
#' patches exchange infectious individuals according to a movement
#' kernel.
#'
#' This is particularly useful for modelling cross-regional seeding,
#' airport-based importations, and the effect of travel restrictions.
#' @name EpiNova-multipatch
NULL

# ============================================================
# IDEA 11: Metapopulation / multi-patch SEIR
# ============================================================

#' Build a mobility-coupled multi-patch SEIR ODE system
#'
#' @param n_patches    Integer number of patches.
#' @param M            n_patches x n_patches mobility matrix.
#'   \code{M[i,j]} is the daily fraction of population in patch i
#'   that travels to patch j (rows sum to <= 1).
#' @param beta_vec     Length-n_patches vector of transmission rates.
#' @param gamma_vec    Length-n_patches vector of recovery rates.
#' @param sigma_vec    Length-n_patches vector of incubation rates.
#' @param pi_fn_list   List of n_patches pi(t) functions (one per patch).
#'   Use \code{NULL} for no intervention in that patch.
#' @param N_vec        Population sizes for each patch.
#'
#' @return A function suitable for \code{deSolve::ode}.
#' @export
build_multipatch_SEIR <- function(n_patches, M,
                                   beta_vec, gamma_vec, sigma_vec,
                                   pi_fn_list = NULL,
                                   N_vec = rep(1, n_patches)) {

  if (is.null(pi_fn_list))
    pi_fn_list <- replicate(n_patches, function(t) 1, simplify = FALSE)

  # State vector layout: [S1..Sn, E1..En, I1..In, R1..Rn]
  function(t, state, parms) {
    S <- state[1:n_patches]
    E <- state[(n_patches + 1):(2 * n_patches)]
    I <- state[(2 * n_patches + 1):(3 * n_patches)]
    R <- state[(3 * n_patches + 1):(4 * n_patches)]

    pi_t <- sapply(pi_fn_list, function(f) f(t))

    # Local transmission
    lambda <- pi_t * beta_vec * I

    # Mobility: effective infectious pressure from visitors
    # Lagrangian mobility (residents travel, carry disease home)
    I_eff <- as.vector(t(M) %*% I + diag(1 - rowSums(M)) %*% I)
    lambda_mob <- pi_t * beta_vec * I_eff

    dS <- -lambda_mob * S
    dE <-  lambda_mob * S - sigma_vec * E
    dI <-  sigma_vec * E  - gamma_vec * I
    dR <-  gamma_vec * I

    list(c(dS, dE, dI, dR))
  }
}


#' Solve a multi-patch SEIR model
#'
#' @param ode_fn   ODE function from \code{build_multipatch_SEIR}.
#' @param init_mat n_patches x 4 matrix of initial conditions
#'   (columns: S, E, I, R).
#' @param times    Time points.
#' @param n_patches Number of patches.
#'
#' @return A data frame with columns \code{time} and
#'   \code{S_i, E_i, I_i, R_i} for each patch i.
#' @export
solve_multipatch <- function(ode_fn, init_mat, times, n_patches) {
  y0  <- as.vector(init_mat)  # column-major: S1..Sn, E1..En, ...
  out <- deSolve::ode(y = y0, times = times, func = ode_fn,
                      parms = list(), method = "lsoda")
  df  <- as.data.frame(out)
  compartments <- c("S","E","I","R")
  for (ci in seq_along(compartments)) {
    for (pi in seq_len(n_patches)) {
      col_idx <- 1 + (ci - 1) * n_patches + pi  # +1 for time column
      df[[paste0(compartments[ci], "_", pi)]] <- df[[col_idx]]
    }
  }
  df[, c("time", paste0(rep(compartments, each = n_patches),
                         "_", seq_len(n_patches)))]
}


# ============================================================
# IDEA 12: Gravity model for mobility matrix
#   When empirical commuting data are unavailable, a gravity model
#   lets users generate a plausible M from population sizes and
#   geographic distances.
# ============================================================

#' Build a gravity-model mobility matrix
#'
#' Uses the functional form
#' \deqn{M_{ij} = \kappa \frac{N_i^\alpha N_j^\beta}{d_{ij}^\gamma}}
#' normalised so that each row sums to at most \code{max_travel}.
#'
#' @param N_vec      Length-n population sizes.
#' @param dist_mat   n x n symmetric matrix of pairwise distances.
#' @param kappa      Proportionality constant.
#' @param alpha,beta Population exponents (default 1).
#' @param grav_gamma Distance decay exponent (default 2).
#' @param max_travel Maximum daily travel fraction per patch (0–1).
#'
#' @return A row-normalised n x n mobility matrix.
#' @export
gravity_mobility <- function(N_vec, dist_mat,
                              kappa = 1e-5,
                              alpha = 1, beta = 1,
                              grav_gamma = 2,
                              max_travel = 0.1) {
  n <- length(N_vec)
  M <- matrix(0, n, n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i != j && dist_mat[i, j] > 0) {
        M[i, j] <- kappa *
          N_vec[i]^alpha * N_vec[j]^beta /
          dist_mat[i, j]^grav_gamma
      }
    }
  }
  # Normalise rows so they sum to at most max_travel
  row_sums <- rowSums(M)
  scale    <- pmin(1, max_travel / pmax(row_sums, 1e-12))
  M * scale
}
