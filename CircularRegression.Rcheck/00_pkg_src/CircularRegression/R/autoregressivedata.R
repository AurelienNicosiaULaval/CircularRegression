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
#'
#' @export
autoregressivedata <- function(data, response, p, ID = NULL) {
  if (!is.data.frame(data))
    stop("'data' must be a data frame.")
  if (!(response %in% names(data)))
    stop("The response variable specified in 'response' is not present in 'data'.")
  if (!is.numeric(p) || length(p) != 1 || p < 1)
    stop("'p' must be a positive integer.")
  if (!is.null(ID) && !(ID %in% names(data)))
    stop("The grouping variable specified in 'ID' is not present in 'data'.")

  p <- as.integer(p)
  for (k in seq_len(p)) {
    colname <- paste0(response, "_tm", k)
    data[[colname]] <- NA
    if (is.null(ID)) {
      data[[colname]] <- c(rep(NA, k), head(data[[response]], -k))
    } else {
      for (id in unique(data[[ID]])) {
        idx <- which(data[[ID]] == id)
        x <- data[[response]][idx]
        data[[colname]][idx] <- c(rep(NA, k), head(x, -k))
      }
    }
  }
  data
}
