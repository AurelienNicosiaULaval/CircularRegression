#' Bison movement data
#'
#' A subset of the bison movement dataset used for illustrating angular
#' regression models of the form \eqn{y = x_1 + x_2:z_2 + \cdots}, where `x.*`
#' denote angular covariates (target directions) and `z.*` the corresponding
#' positive covariates (weights) in the product \code{x:z}. The data originate
#' from GPS-tracked bison trajectories and environmental variables describing
#' landscape features such as meadows and canopy gaps.
#'
#' @format A data frame with \code{5,696} rows and 7 columns:
#' \describe{
#'   \item{y.dir}{Direction of movement at time \eqn{t}.}
#'   \item{y.prec}{Direction of movement at time \eqn{t-1}.}
#'   \item{y.prec2}{Direction of movement at time \eqn{t-2}.}
#'   \item{x.meadow}{Angle toward the nearest meadow.}
#'   \item{z.meadow}{Distance toward the nearest meadow associated with \code{x.meadow}.}
#'   \item{x.gap}{Angle toward the nearest canopy gap.}
#'   \item{z.gap}{Distance toward the nearest canopy gap associated with \code{x.gap}.}
#' }
#'
#' @details
#' This reduced dataset is designed to work directly with the functions of the
#' \emph{CircularRegression} package implementing the general angular regression
#' model of Rivest et al. (2016). Each \code{x:z} pair represents a weighted
#' target direction contributing to the mean resultant vector
#' \eqn{\sum_j \beta_j z_j (\cos x_j, \sin x_j)}, where \code{z.*} modulates the
#' attraction toward the angular covariate \code{x.*}. Variables
#' \code{y.prec} and \code{y.prec2} can be included to model directional
#' persistence between successive movement steps.
#'
#' Ecologically, these variables describe the orientation of bison movement
#' relative to attractive features in the landscape: open meadows and canopy
#' gaps. The dataset serves as a clean example for model testing, illustration,
#' and teaching applications involving directional movement analysis.
#'
#' @references
#' Rivest, L.-P., Duchesne, T., Nicosia, A., and Fortin, D. (2016).
#' A general angular regression model for the analysis of data on animal movement in ecology.
#' \emph{Journal of the Royal Statistical Society: Series C (Applied Statistics)}, 65(3), 445–463.
#'
#' @seealso \code{\link{angular}}, \code{\link{consensus}}, \code{\link{pick_reference_angle}}
#'
#' @examples
#' # View structure
#' str(bison)
#'
#'
"bison"
