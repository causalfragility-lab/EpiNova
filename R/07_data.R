#' Hubei Province COVID-19 data (Jan 13 - Feb 11, 2020)
#'
#' Cumulative daily confirmed and removed (recovered + deceased) COVID-19
#' case counts for Hubei Province, China, covering 30 days from
#' 13 January 2020 to 11 February 2020. Same dataset used in the
#' original eSIR paper.
#'
#' @format A named list with six elements:
#' \describe{
#'   \item{NI}{Integer vector (length 30). Cumulative confirmed cases.}
#'   \item{RI}{Integer vector (length 30). Cumulative removed cases
#'     (recovered + deceased).}
#'   \item{N}{Numeric scalar. Hubei population size (58,500,000).}
#'   \item{begin_date}{Character. First observation date ("2020-01-13").}
#'   \item{end_date}{Character. Last observation date ("2020-02-11").}
#'   \item{description}{Character. Plain-text description of the dataset.}
#' }
#' @source DXY.cn daily COVID-19 situation reports.
#' @examples
#' Y <- hubei_covid$NI / hubei_covid$N - hubei_covid$RI / hubei_covid$N
#' R <- hubei_covid$RI / hubei_covid$N
#' plot(Y, type = "l", ylab = "Infected proportion", xlab = "Day")
#' @export
hubei_covid <- list(
  NI = c(41L, 41L, 41L, 45L, 62L, 131L, 200L, 270L, 375L, 444L,
         549L, 729L, 1052L, 1423L, 2714L, 3554L, 4903L, 5806L,
         7153L, 9074L, 11177L, 13522L, 16678L, 19665L, 22112L,
         24953L, 27100L, 29631L, 31728L, 33366L),
  RI = c(1L, 1L, 7L, 10L, 14L, 20L, 25L, 31L, 34L, 45L, 55L,
         71L, 94L, 121L, 152L, 213L, 252L, 345L, 417L, 561L,
         650L, 811L, 1017L, 1261L, 1485L, 1917L, 2260L, 2725L,
         3284L, 3754L),
  N           = 58.5e6,
  begin_date  = "2020-01-13",
  end_date    = "2020-02-11",
  description = paste(
    "Cumulative COVID-19 cases in Hubei Province, China",
    "(Jan 13 - Feb 11 2020). Source: DXY.cn."
  )
)


#' Prepare population proportions from a hubei_covid-style list
#'
#' Converts raw cumulative counts into the proportion vectors \code{Y}
#' (active infected) and \code{R} (removed) expected by
#' \code{solve_model} and the fitting functions.
#'
#' @param dat A list with elements \code{NI}, \code{RI}, and \code{N}.
#'
#' @return A list with numeric vectors \code{Y} and \code{R}.
#' @examples
#' props <- prep_proportions(hubei_covid)
#' str(props)
#' @export
prep_proportions <- function(dat) {
  stopifnot(all(c("NI", "RI", "N") %in% names(dat)))
  list(
    Y = dat$NI / dat$N - dat$RI / dat$N,
    R = dat$RI / dat$N
  )
}
