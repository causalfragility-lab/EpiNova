# EpiNova <img src="man/figures/logo.png" align="right" height="120"/>

**Flexible Extended State-Space Epidemiological Models with Modern Inference**

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

---

## Motivation

`eSIR` (Song Lab, 2020) was a pioneering R package for COVID-19 modelling
that introduced time-varying quarantine protocols into a Bayesian SIR
framework.  After studying its architecture, we identified **15+ structural
limitations** and designed `EpiNova` as a ground-up redesign.

---

## Key New Ideas vs eSIR

### 🧬 Compartmental Model Hierarchy (Idea 1–3)
```
eSIR:    S → I → R  (only)

EpiNova: S → I → R                           (SIR)
         S → E → I → R                       (SEIR)
         S → E → I → R / D                   (SEIRD)
         S ⇄ V → E → I → R                   (SVEIR)
         S ⇄ V → E → I → R / D               (SVEIRD) ← waning immunity
         [S,E,I,R] × n_age groups + contact matrix  (age-SEIR)
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
# Option A: Fast MLE (no external binary needed)
fit <- fit_mle(Y, R, N, model = "SEIRD",
               par_init = list(beta=0.3, gamma=0.1, ...))

# Option B: Full Bayesian via Stan/HMC
fit <- fit_hmc(Y, R, N, model = "SEIRD",   # uses cmdstanr
               prior_list = list(beta = c(0.3, 0.1), ...))

# Option C: Real-time Sequential Monte Carlo
fit <- fit_smc(Y, R, N, model = "SEIR",
               prior_fn = my_prior,
               n_particles = 2000)

# Option D: Rt estimation (EpiEstim integration)
Rt  <- estimate_Rt(new_cases, mean_si = 5.2, sd_si = 2.8)
```

**Why it matters:**
- MLE via `nlminb` + `DEoptim` global search: **10-100× faster** than JAGS MCMC
- HMC mixes dramatically faster than Gibbs for continuous parameters
- SMC enables **real-time daily updating** without re-running the full chain
- No JAGS binary install required for MLE and SMC paths

---

### 🗺️ Multi-Patch Spatial Models (Ideas 11–12)

```r
# Build mobility matrix from geography (gravity model)
M <- gravity_mobility(N_vec, dist_mat, kappa = 1e-7)

# Couple n patches with different beta, intervention
ode <- build_multipatch_SEIR(
  n_patches  = 5,
  M          = M,
  beta_vec   = c(0.35, 0.28, 0.22, 0.20, 0.18),
  pi_fn_list = list(lockdown_A, lockdown_B, no_npi, no_npi, no_npi)
)
```
eSIR had **zero** spatial structure.

---

### 🔮 Ensemble Forecasting & Scenario Analysis (Ideas 13–15)

```r
# BMA across model types
ens <- ensemble_forecast(Y, R, N,
         models    = c("SEIR","SEIRD","SVEIRD"),
         par_inits = list(...), par_bounds = list(...))

# "What-if" scenarios
sc <- project_scenarios(fit, scenarios = list(
  "No intervention"    = function(t) 1,
  "Light lockdown"     = build_pi_step(c(10), c(1, 0.6)),
  "Strict lockdown"    = build_pi_spline(c(0,10,20), c(1,0.8,0.2))
))

# Score your forecasts
score_forecast(ens$ensemble_forecast, actual_Y = Y_holdout)
# Returns: CRPS, 95% coverage, MAE
```

---

### 📊 Rich Visualisations (Idea 16)

```r
plot_trajectory(traj, obs_Y = Y, obs_R = R)   # all compartments
plot_forecast(forecast_df, obs_Y = Y)          # ribbon plot
plot_scenarios(sc_df)                          # scenario comparison
plot_Rt(Rt_df, change_times = c(10, 22))       # Rt with interventions
plot_multipatch_snapshot(mp_df, t = 30)        # spatial bar chart
```

---

## Installation

```r
# Install dependencies
install.packages(c("deSolve","ggplot2","dplyr","tidyr","splines"))

# Optional but recommended
install.packages("DEoptim")        # global optimisation
install.packages("EpiEstim")       # Rt estimation
install.packages("cmdstanr")       # HMC backend

# Install EpiNova
devtools::install_github("yourname/EpiNova")
```

No JAGS, no rjags, no external binary required for core functionality.

---

## Quick Start

```r
library(EpiNova)

# 1. Define a smooth intervention
pi_fn <- build_pi_spline(
  knot_times  = c(0, 10, 22, 60),
  knot_values = c(1, 0.9, 0.4, 0.25)
)

# 2. Solve SEIRD ODE
traj <- solve_model(
  params = list(beta=0.35, gamma=0.07, sigma=0.2, delta=0.005),
  init   = c(S=0.9999, E=0.00005, I=0.00005, R=0, D=0),
  times  = 0:200,
  model  = "SEIRD",
  pi_fn  = pi_fn
)

# 3. Fit to data
fit <- fit_mle(obs_Y, obs_R, N = 58.5e6,
               model    = "SEIRD",
               par_init = list(beta=0.3, gamma=0.08, sigma=0.15,
                               delta=0.003, I0=1e-4, E0=2e-4))

# 4. Compare intervention scenarios
sc <- project_scenarios(fit, scenarios = list(
  "Baseline"    = pi_fn,
  "Relaxed"     = build_pi_step(c(10), c(1, 0.6)),
  "Counterfactual (no NPI)" = function(t) 1
))

plot_scenarios(sc, obs_Y = obs_Y)
```

---

## Architecture

```
EpiNova/
├── R/
│   ├── 01_ode_models.R          # solve_model() dispatcher + ODE systems
│   ├── 02_intervention_functions.R  # build_pi_*(), compose_pi(), build_phi_pulse()
│   ├── 03_inference.R           # fit_mle(), fit_smc(), estimate_Rt()
│   ├── 04_multipatch.R          # build_multipatch_SEIR(), gravity_mobility()
│   ├── 05_ensemble.R            # ensemble_forecast(), project_scenarios(), score_forecast()
│   └── 06_plotting.R            # plot_trajectory(), plot_forecast(), plot_scenarios(), plot_Rt()
├── inst/stan/                   # Stan model files for HMC backend
├── vignettes/
│   └── getting_started.Rmd
└── DESCRIPTION
```

---

## References

- Song Lab (2020). eSIR: Extended State-Space SIR Models.
- Osthus et al. (2017). Forecasting seasonal influenza with a state-space SIR model. *Ann. Appl. Stat.*
- Cori et al. (2013). A new framework for Rt estimation. *Am. J. Epidemiol.*
- Kristensen et al. (2016). TMB: Automatic differentiation and Laplace approximation. *J. Stat. Softw.*
- Gneiting & Raftery (2007). Strictly proper scoring rules, prediction, and estimation. *JASA.*
