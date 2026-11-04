#' Create autoregressive dataset
#'
#' This function augments a data frame with lagged versions of an angular response
#' variable so that autoregressive circular regression models can be fitted.
#'
#' @param data A data frame containing at least the response variable.
#' @param response A character string specifying the name of the angular response variable.
#' @param p Integer giving the order of the autoregressive structure.
#' @param ID Optional. A character string specifying the name of a grouping variable. If provided, lags are computed within each group.
#'
#' @return The input \code{data} with \code{p} additional columns named \code{response_tmk} for \eqn{k=1,\ldots,p} containing the lagged values of \code{response}.
#' @details Lags follow the existing row order within each group. Sort the data
#'   by time first if necessary. Groups shorter than a requested lag receive
#'   missing values, and missing group identifiers are rejected. This utility
#'   does not model dependence or provide serially robust inference.
#' @examples
#' data(bison)
#' head(autoregressivedata(bison[1:6, ], "y.dir", p = 2))
#'
#' @export
autoregressivedata <- function(data, response, p, ID = NULL) {
  if (!is.data.frame(data))
    stop("'data' must be a data frame.")
  if (!is.character(response) || length(response) != 1L || is.na(response) || !(response %in% names(data)))
    stop("The response variable specified in 'response' is not present in 'data'.")
  if (!is.numeric(p) || length(p) != 1 || !is.finite(p) || p < 1 ||
      p != floor(p) || p >= .Machine$integer.max)
    stop("'p' must be a positive integer.")
  if (!is.null(ID) && (!is.character(ID) || length(ID) != 1L || is.na(ID) || !(ID %in% names(data))))
    stop("The grouping variable specified in 'ID' is not present in 'data'.")

  if (!is.null(ID) && anyNA(data[[ID]])) stop("Grouping identifiers must not be missing.")
  p <- as.integer(p)
  groups <- if (is.null(ID)) list(seq_len(nrow(data))) else split(seq_len(nrow(data)), data[[ID]])
  for (k in seq_len(p)) {
    colname <- paste0(response, "_tm", k)
    data[[colname]] <- data[[response]][rep(NA_integer_, nrow(data))]
    for (idx in groups) {
      if (length(idx) > k) {
        data[[colname]][idx[(k + 1L):length(idx)]] <- data[[response]][idx[seq_len(length(idx) - k)]]
      }
    }
  }
  data
}
