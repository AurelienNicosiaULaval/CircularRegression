#' Sandhopper Escape Trials
#'
#' Repeated escape orientations from \link{Sandhopper}, with one row per
#' animal and trial. The original angles are converted from degrees to radians.
#' Wind speed retains its original units and missing values are retained.
#'
#' @format A data frame with 360 rows and six variables.
#' \describe{
#'   \item{Anim}{Animal identifier.}
#'   \item{trial}{Trial number within animal, from 1 to 5.}
#'   \item{escape}{Escape orientation relative to north, in radians.}
#'   \item{sun}{Sun azimuth at release, in radians.}
#'   \item{wind}{Wind direction, in radians.}
#'   \item{wind_speed}{Wind speed, in the original recorded units.}
#' }
#' @source \link{Sandhopper}. The script \code{data-raw/sandhopper_trials.R}
#'   in the package repository reshapes the data and converts the angles.
#' @seealso \code{\link{angular_re}}, \code{\link{Sandhopper}}
#' @keywords datasets
"sandhopper_trials"
