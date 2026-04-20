# EpiNova <img src="man/figures/logo.png" align="right" height="120" alt=""/>

**Flexible Extended State-Space Epidemiological Models with Modern Inference**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![CRAN status](https://www.r-pkg.org/badges/version/EpiNova)](https://CRAN.R-project.org/package=EpiNova)
[![R CMD Check](https://img.shields.io/badge/R%20CMD%20Check-passing-brightgreen)](https://github.com/causalfragility-lab/EpiNova)
[![Tests](https://img.shields.io/badge/tests-76%20passing-brightgreen)](https://github.com/causalfragility-lab/EpiNova)

---

## Motivation

`eSIR` (Song Lab, 2020) was a pioneering R package for COVID-19 modelling
that introduced time-varying quarantine protocols into a Bayesian SIR
framework.  After studying its architecture, we identified **15+ structural
limitations** and designed `EpiNova` as a ground-up redesign.

---

## Installation

```r
# Install from GitHub
devtools::install_github("causalfragility-lab/EpiNova")

# Core dependencies (installed automatically)
# deSolve, ggplot2, dplyr, tidyr, splines, scales

# Optional — unlock extra capabilities
install.packages("DEoptim")   # faster global optimisation in fit_mle()
install.packages("EpiEstim")  # full Bayesian Rt estimation
install.packages("MASS")      # parametric bootstrap uncertainty
install.packages("numDeriv")  # Hessian-based confidence intervals
install.packages("cmdstanr")  # HMC backend (Stan)
```

No JAGS binary required. Core functionality is pure R.

---

## Quick Start

```r
library(EpiNova)

# --- Built-in Hubei COVID-19 data ---
props <- prep_proportions(hubei_covid)   # Y = active infected, R = removed

# --- Smooth spline intervention ---
pi_fn <- build_pi_spline(
  knot_times  = c(0, 10, 22, 60),
  knot_values = c(1, 0.9, 0.4, 0.25)
)

# --- Solve SEIRD ODE ---
traj <- solve_model(
  params = list(beta = 0.35, gamma = 0.07, sigma = 0.2, delta = 0.005),
  init   = c(S = 0.9999, E = 0.00005, I = 0.00005, R = 0, D = 0),
  times  = 0:200,
  model  = "SEIRD",
  pi_fn  = pi_fn
)
plot_trajectory(traj, obs_Y = props$Y, obs_R = props$R, T_obs_end = 29)

# --- Scenario comparison ---
sc <- do.call(rbind, lapply(
  list("No NPI"   = function(t) 1,
       "Lockdown" = build_pi_step(c(10), c(1, 0.4)),
       "Spline"   = pi_fn),
  function(pi) {
    tr <- solve_model(
      list(beta=0.35, gamma=0.07, sigma=0.2, delta=0.005),
      c(S=0.9999, E=5e-5, I=5e-5, R=0, D=0), 0:200, "SEIRD", pi_fn=pi)
    data.frame(time=tr$time, I_median=tr$I,
               I_lower=tr$I*0.8, I_upper=tr$I*1.2, scenario=deparse(substitute(pi)))
  }
))
plot_scenarios(sc)

# --- Real-time Rt estimation (no extra packages) ---
new_cases <- pmax(0L, diff(hubei_covid$NI))
Rt_df     <- estimate_Rt_simple(new_cases, mean_si = 5.2, window = 7L)
plot_Rt(Rt_df, change_times = c(10, 22))
```

---

## Key New Ideas vs eSIR

### 🧬 Compartmental Model Hierarchy (Ideas 1–3)

```
eSIR:    S → I → R  (only)

EpiNova: S → I → R                                    (SIR)
         S → E → I → R                                (SEIR)
         S → E → I → R / D                            (SEIRD)
         S ⇄ V → E → I → R                            (SVEIR)
         S ⇄ V → E → I → R / D   ← waning immunity    (SVEIRD)
         [S,E,I,R] × n_age groups + contact matrix     (age-SEIR)
```

All models share the same `solve_model(model = "SVEIRD")` interface.

---

### 📐 Flexible π(t) Intervention Functions (Ideas 4–7)

| Builder | Shape | Advantage over eSIR |
|---|---|---|
| `build_pi_step()` | Step function | Same as eSIR (baseline) |
| `build_pi_exp()` | Exponential decay | Same as eSIR (baseline) |
| **`build_pi_spline()`** | Natural cubic spline | **Smooth, realistic policy ramp** |
| **`gp_cov_sqexp()`** | Gaussian Process prior | **Data-driven, non-parametric** |
| **`compose_pi()`** | Multiplicative composite | **Multiple simultaneous NPIs** |
| **`build_phi_pulse()`** | Gaussian pulse (vs Dirac) | **Differentiable for HMC** |

---

### ⚡ Modern Inference Backends (Ideas 8–10)

```r
# Fast MLE — no external binary needed
fit <- fit_mle(Y, R, N, model = "SEIRD",
               par_init = list(beta=0.3, gamma=0.1, sigma=0.2,
                               delta=0.003, I0=1e-4, E0=2e-4))

# Real-time Sequential Monte Carlo
fit <- fit_smc(Y, R, N, model = "SEIR",
               prior_fn = my_prior, n_particles = 2000)

# Dependency-free Rt estimation
Rt <- estimate_Rt_simple(new_cases, mean_si = 5.2, window = 7L)

# Full Bayesian Rt via EpiEstim (when installed)
Rt <- estimate_Rt(new_cases, mean_si = 5.2, sd_si = 2.8)
```

**Why it matters:**
- MLE via `DEoptim` + `nlminb`: **10–100× faster** than JAGS MCMC
- SMC enables **real-time daily updating** — impossible with JAGS
- No JAGS, no rjags, no external binary for core paths

---

### 🗺️ Multi-Patch Spatial Models (Ideas 11–12)

```r
# Build mobility matrix from geography
M <- gravity_mobility(N_vec, dist_mat, kappa = 1e-7)

# Couple n patches with per-patch parameters and interventions
ode <- build_multipatch_SEIR(
  n_patches  = 3,
  M          = M,
  beta_vec   = c(0.35, 0.28, 0.22),
  pi_fn_list = list(strict_lockdown, mild_lockdown, function(t) 1)
)
mp_df <- solve_multipatch(ode, init_mat, times = 0:150, n_patches = 3)
plot_multipatch_snapshot(mp_df, t_snapshot = 30,
                          patch_names = c("Hubei", "Guangdong", "Beijing"))
```

eSIR had **zero** spatial structure.

---

### 🔮 Ensemble Forecasting & Scenario Analysis (Ideas 13–15)

```r
# BMA across model types
ens <- ensemble_forecast(Y, R, N,
         models    = c("SEIR", "SEIRD", "SVEIRD"),
         par_inits = list(...), par_bounds = list(...))

# "What-if" scenario projection
sc <- project_scenarios(fit, scenarios = list(
  "No intervention" = function(t) 1,
  "Light lockdown"  = build_pi_step(c(10), c(1, 0.6)),
  "Strict lockdown" = build_pi_spline(c(0,10,20,60), c(1,0.8,0.2,0.2))
))

# Score your forecasts with proper scoring rules
score_forecast(ens$ensemble_forecast, actual_Y = Y_holdout)
# Returns: CRPS, 95% coverage, MAE
```

---

### 📊 Rich Visualisations (Idea 16)

```r
plot_trajectory(traj, obs_Y = Y, obs_R = R, T_obs_end = 29)
plot_forecast(forecast_df, obs_Y = Y)
plot_scenarios(sc_df)
plot_Rt(Rt_df, change_times = c(10, 22))
plot_multipatch_snapshot(mp_df, t_snapshot = 30)
```

---

## Package Architecture

```
EpiNova/
├── R/
│   ├── 01_ode_models.R          # solve_model() dispatcher + all ODE systems
│   ├── 02_intervention_functions.R  # build_pi_*(), compose_pi(), build_phi_pulse()
│   ├── 03_inference.R           # fit_mle(), fit_smc(), estimate_Rt*()
│   ├── 04_multipatch.R          # build_multipatch_SEIR(), gravity_mobility()
│   ├── 05_ensemble.R            # ensemble_forecast(), project_scenarios(), score_forecast()
│   ├── 06_plotting.R            # all plot_*() functions
│   └── 07_data.R                # hubei_covid dataset + prep_proportions()
├── tests/testthat/              # 76 unit tests (6 files)
├── vignettes/
│   └── getting_started.Rmd
├── data-raw/
│   └── hubei_covid.R
└── DESCRIPTION
```

---

## References

- Song Lab (2020). eSIR: Extended State-Space SIR Models.
- Osthus et al. (2017). Forecasting seasonal influenza with a state-space SIR model. *Ann. Appl. Stat.*
- Cori et al. (2013). A new framework for Rt estimation. *Am. J. Epidemiol.*
- Kristensen et al. (2016). TMB: Automatic differentiation and Laplace approximation. *J. Stat. Softw.*
- Gneiting & Raftery (2007). Strictly proper scoring rules, prediction, and estimation. *JASA.*

---

## Citation

```bibtex
@Manual{EpiNova2025,
  title  = {EpiNova: Flexible Extended State-Space Epidemiological Models
             with Modern Inference},
  author = {Hait, Subir},
  year   = {2025},
  note   = {R package version 0.1.0},
  url    = {https://github.com/causalfragility-lab/EpiNova}
}
```
