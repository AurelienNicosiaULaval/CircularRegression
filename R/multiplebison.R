#' Two plains bison tracked by GPS (Prince Albert National Park, 2013)
#'
#' Hourly GPS-derived movement metrics for two plains bison (IDs `1044-a` and
#' `1045-a`) monitored in Prince Albert National Park (Saskatchewan, Canada)
#' during July-October 2013. From the hourly locations, the dataset includes
#' per-animal direction of movement, step length (distance between consecutive
#' fixes), latitude/longitude, the inter-animal distance at each hour, and
#' turning-angles (difference between consecutive movement directions).
#'
#' The file results from merging individual hourly GPS tracks (one per animal),
#' retaining rows where the time gap between successive fixes is exactly one hour,
#' and reshaping to a "wide" format with one row per timestamp holding both
#' animals' variables. Turning-angles are computed from consecutive directions.
#'
#' @format A data frame with \eqn{n = 1475} rows (hours kept after filtering) and the
#' following columns (angles in radians unless stated otherwise):
#' \describe{
#'   \item{datetime_round_next}{POSIXct. Timestamp of the retained hourly record.}
#'   \item{heure}{integer. Hour of day (0-23).}
#'   \item{bearing_ang}{numeric. Bearing from bison 1044-a to 1045-a.}
#'   \item{diff}{numeric. Difference between the animals' movement directions (1044-a minus 1045-a).}
#'   \item{direction_1044-a}{numeric. Movement direction of bison 1044-a.}
#'   \item{direction_1045-a}{numeric. Movement direction of bison 1045-a.}
#'   \item{distance_1044-a}{numeric. Step length (m) of bison 1044-a between consecutive fixes.}
#'   \item{distance_1045-a}{numeric. Step length (m) of bison 1045-a between consecutive fixes.}
#'   \item{distance}{numeric. Inter-animal distance (m) at the timestamp.}
#'   \item{lat_1044-a}{numeric. Latitude of bison 1044-a (decimal degrees).}
#'   \item{lat_1045-a}{numeric. Latitude of bison 1045-a (decimal degrees).}
#'   \item{long_1044-a}{numeric. Longitude of bison 1044-a (decimal degrees).}
#'   \item{long_1045-a}{numeric. Longitude of bison 1045-a (decimal degrees).}
#'   \item{turn_ang_1044}{numeric. Turning angle of bison 1044-a between consecutive directions.}
#'   \item{turn_ang_1045}{numeric. Turning angle of bison 1045-a between consecutive directions.}
#' }
#'
#' @details
#' - Angles are measured around the circle; if stored in degrees in your copy,
#'   convert to radians with \code{circular::rad(x, units = "degrees")}, or
#'   \code{circular::deg2rad()} before modeling.
#' - Step lengths and distances were computed from consecutive GPS fixes; rows
#'   are kept only when the time difference equals one hour.
#' - The turning-angle is the signed circular difference between consecutive
#'   movement directions (current minus previous, wrapped to \eqn{(-\pi,\pi]}).
#'
#' @references
#' Mardia, K. V., & Jupp, P. E. (2000). \emph{Directional Statistics}. Wiley.
#'
#' @examples
#' ## Basic glimpse
#' head(multiplebison)
#'
#'
#' ## Example: inter-animal proximity rate (< 500 m)
#' # mean(multiplebison$Distance < 500, na.rm = TRUE)
#'
"multiplebison"
