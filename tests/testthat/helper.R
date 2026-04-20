# Helper: small SEIRD trajectory used across tests
.make_seird_traj <- function(times = 0:60) {
  params <- list(beta = 0.3, gamma = 0.1, sigma = 0.2, delta = 0.005)
  init   <- c(S = 0.9989, E = 0.0006, I = 0.0005, R = 0, D = 0)
  solve_model(params, init, times, model = "SEIRD")
}
