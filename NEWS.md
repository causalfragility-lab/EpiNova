# EpiNova 0.1.0

## Initial release

* `solve_model()` — unified ODE dispatcher for SIR, SEIR, SEIRD, SVEIR, SVEIRD, age-SEIR
* `build_pi_spline()`, `build_pi_step()`, `build_pi_exp()`, `compose_pi()` — flexible intervention functions
* `build_phi_pulse()` — smooth Gaussian quarantine pulse
* `gp_cov_sqexp()` — GP covariance builder for Bayesian pi(t)
* `fit_mle()` — maximum likelihood inference with optional DEoptim global search
* `fit_smc()` — sequential Monte Carlo particle filter for real-time updating
* `estimate_Rt_simple()` — dependency-free sliding-window Rt estimator
* `estimate_Rt()` — EpiEstim wrapper with graceful fallback
* `build_multipatch_SEIR()`, `solve_multipatch()` — spatial multi-patch models
* `gravity_mobility()` — gravity-model mobility matrix builder
* `ensemble_forecast()` — Bayesian model averaging across model types
* `project_scenarios()` — what-if scenario projection
* `score_forecast()` — CRPS, coverage, and MAE scoring
* `plot_trajectory()`, `plot_forecast()`, `plot_scenarios()`, `plot_Rt()`, `plot_multipatch_snapshot()` — publication-ready ggplot2 graphics
* `hubei_covid` — built-in Hubei Province COVID-19 dataset (Jan–Feb 2020)
* `prep_proportions()` — convert raw counts to model-ready proportions
* 76 unit tests across 6 test files; `0 errors / 0 warnings / 0 notes` on R CMD check

---

## Preparing for CRAN submission

- Switched license from CC BY 4.0 to MIT (CRAN-standard)
- Fixed `URL` and `BugReports` to real GitHub repository URLs
- Added `LICENSE` and `LICENSE.md` files
- Added `inst/CITATION` file
- Wrapped slow examples to stay within CRAN's 5-second limit
- Added `cran-comments.md`
- Verified 0 errors / 0 warnings / 0 notes under `--as-cran`
