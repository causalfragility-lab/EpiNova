#' @title Visualisation Layer for EpiNova
#'
#' @description
#' Publication-ready ggplot2 graphics. All plot functions return
#' ggplot2 objects that can be further customised.
#' @name EpiNova-plotting
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon geom_point
#'   geom_vline geom_hline geom_col scale_color_manual scale_fill_manual
#'   scale_fill_gradient labs theme_minimal theme element_text facet_wrap
#' @importFrom tidyr pivot_longer
#' @importFrom scales hue_pal
#' @importFrom stats median quantile setNames
NULL

utils::globalVariables(c(
  "time", "proportion", "compartment",
  "I_median", "I_lower", "I_upper",
  "Rt_mean", "Rt_lower", "Rt_upper", "t_end",
  "scenario", "patch", "I"
))

# palette consistent across all plots
.EpiNova_palette <- c(
  S = "#4C72B0", E = "#DD8452", I = "#C44E52",
  R = "#55A868", D = "#8C8C8C", V = "#9467BD",
  Rt = "#E377C2"
)


#' Plot model trajectory with observed data
#'
#' @param traj_df   Data frame from \code{solve_model}.
#' @param obs_Y     Observed infected proportions (optional).
#' @param obs_R     Observed removed proportions (optional).
#' @param T_obs_end Last observed day (draws a vertical cutoff line).
#' @param title     Plot title.
#'
#' @return A ggplot2 object.
#' @export
plot_trajectory <- function(traj_df, obs_Y = NULL, obs_R = NULL,
                             T_obs_end = NULL,
                             title = "Epidemic Trajectory") {

  compartments <- setdiff(names(traj_df), "time")
  long_df <- tidyr::pivot_longer(traj_df,
                                  cols      = tidyr::all_of(compartments),
                                  names_to  = "compartment",
                                  values_to = "proportion")

  p <- ggplot2::ggplot(long_df,
         ggplot2::aes(x = time, y = proportion, color = compartment)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::scale_color_manual(values = .EpiNova_palette,
                                na.value = "grey50") +
    ggplot2::labs(title = title, x = "Days",
                  y = "Population proportion", color = "Compartment") +
    ggplot2::theme_minimal(base_size = 13)

  if (!is.null(obs_Y)) {
    obs_df_I <- data.frame(time = seq_along(obs_Y) - 1,
                            proportion = obs_Y, compartment = "I_obs")
    p <- p + ggplot2::geom_point(data = obs_df_I, shape = 16,
                                  size = 2,
                                  color = .EpiNova_palette["I"])
  }
  if (!is.null(obs_R)) {
    obs_df_R <- data.frame(time = seq_along(obs_R) - 1,
                            proportion = obs_R, compartment = "R_obs")
    p <- p + ggplot2::geom_point(data = obs_df_R, shape = 17,
                                  size = 2,
                                  color = .EpiNova_palette["R"])
  }
  if (!is.null(T_obs_end))
    p <- p + ggplot2::geom_vline(xintercept = T_obs_end,
                                  linetype = "dashed",
                                  color = "steelblue", linewidth = 0.7)
  p
}


#' Plot forecast with uncertainty ribbon
#'
#' @param forecast_df Data frame with \code{time}, \code{I_median},
#'   \code{I_lower}, \code{I_upper}.
#' @param obs_Y       Observed infected proportions.
#' @param T_obs_end   Last observed day.
#' @param title       Plot title.
#'
#' @return A ggplot2 object.
#' @export
plot_forecast <- function(forecast_df, obs_Y = NULL,
                           T_obs_end = NULL,
                           title = "Infection Forecast") {
  p <- ggplot2::ggplot(forecast_df, ggplot2::aes(x = time)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = I_lower, ymax = I_upper),
                         fill = .EpiNova_palette["I"], alpha = 0.25) +
    ggplot2::geom_line(ggplot2::aes(y = I_median),
                       color = .EpiNova_palette["I"], linewidth = 1) +
    ggplot2::labs(title = title, x = "Days",
                  y = "Infected proportion") +
    ggplot2::theme_minimal(base_size = 13)

  if (!is.null(obs_Y))
    p <- p + ggplot2::geom_point(
      data = data.frame(time = seq_along(obs_Y) - 1, I_median = obs_Y),
      ggplot2::aes(y = I_median), size = 2,
      color = .EpiNova_palette["I"])

  if (!is.null(T_obs_end))
    p <- p + ggplot2::geom_vline(xintercept = T_obs_end,
                                  linetype = "dashed", linewidth = 0.7,
                                  color = "steelblue")
  p
}


#' Plot scenario comparison
#'
#' @param scenario_df   Data frame from \code{project_scenarios} or a
#'   manually constructed data frame with columns \code{time},
#'   \code{I_median}, \code{I_lower}, \code{I_upper}, \code{scenario}.
#' @param obs_Y         Optional observed data to overlay.
#'
#' @return A ggplot2 object.
#' @export
plot_scenarios <- function(scenario_df, obs_Y = NULL) {
  n_sc   <- length(unique(scenario_df$scenario))
  colors <- scales::hue_pal()(n_sc)
  names(colors) <- unique(scenario_df$scenario)

  p <- ggplot2::ggplot(scenario_df,
         ggplot2::aes(x = time, color = scenario, fill = scenario)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = I_lower, ymax = I_upper),
                         alpha = 0.18) +
    ggplot2::geom_line(ggplot2::aes(y = I_median), linewidth = 0.95) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_fill_manual(values  = colors) +
    ggplot2::labs(title = "Scenario Comparison",
                  x = "Days", y = "Infected proportion",
                  color = "Scenario", fill = "Scenario") +
    ggplot2::theme_minimal(base_size = 13)

  if (!is.null(obs_Y))
    p <- p + ggplot2::geom_point(
      data = data.frame(time = seq_along(obs_Y) - 1,
                         I_median = obs_Y, scenario = "Observed"),
      ggplot2::aes(y = I_median), size = 2, color = "black")
  p
}


#' Plot effective reproduction number Rt over time
#'
#' @param Rt_df   Data frame with \code{t_end}, \code{Rt_mean},
#'   \code{Rt_lower}, \code{Rt_upper}.
#' @param change_times Optional numeric vector of intervention change
#'   points to mark as vertical lines.
#'
#' @return A ggplot2 object.
#' @export
plot_Rt <- function(Rt_df, change_times = NULL) {
  p <- ggplot2::ggplot(Rt_df, ggplot2::aes(x = t_end)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = Rt_lower, ymax = Rt_upper),
                         fill = .EpiNova_palette["Rt"], alpha = 0.3) +
    ggplot2::geom_line(ggplot2::aes(y = Rt_mean),
                       color = .EpiNova_palette["Rt"], linewidth = 1) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed",
                        color = "darkred", linewidth = 0.7) +
    ggplot2::labs(
      title = expression("Effective Reproduction Number " * R[t]),
      x = "Day", y = expression(R[t])) +
    ggplot2::theme_minimal(base_size = 13)

  if (!is.null(change_times))
    p <- p + ggplot2::geom_vline(xintercept = change_times,
                                  linetype = "dotted", color = "grey40")
  p
}


#' Multi-patch bar chart of infected proportion by patch
#'
#' @param multipatch_df   Data frame from \code{solve_multipatch}.
#' @param t_snapshot      Day to visualise.
#' @param patch_names     Character vector of patch names.
#'
#' @return A ggplot2 bar chart.
#' @export
plot_multipatch_snapshot <- function(multipatch_df, t_snapshot,
                                      patch_names = NULL) {
  row    <- multipatch_df[multipatch_df$time == t_snapshot, ]
  I_cols <- grep("^I_", names(row), value = TRUE)
  df <- data.frame(
    patch = if (!is.null(patch_names)) patch_names
            else gsub("I_", "Patch ", I_cols),
    I     = as.numeric(row[, I_cols])
  )
  df$patch <- factor(df$patch,
                      levels = df$patch[order(df$I, decreasing = TRUE)])

  ggplot2::ggplot(df, ggplot2::aes(x = patch, y = I, fill = I)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_gradient(low = "#FFF3B0", high = "#C44E52") +
    ggplot2::labs(
      title = paste0("Infected proportion by patch (day ", t_snapshot, ")"),
      x = "Patch", y = "Infected proportion", fill = "I") +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35,
                                                        hjust = 1))
}
