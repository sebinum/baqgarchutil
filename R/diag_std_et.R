################################################################################
#' Standardized Univariate Series \eqn{e_{t}}{e_t}
#'
#' Transform a \eqn{k}-dimensional series \eqn{a_{t}}{a_t} to
#' the standardized univariate series \eqn{e_{t}}{e_t}.
#'
#' @param x A \code{matrix} / \code{data.frame} / \code{numeric vector} of
#'   (multivariat) financial time series. Each column contains a series, each
#'   row an observation of the series.
#' @return The standardized univariate series \eqn{e_{t}}{e_t} of \code{x} as
#'   a \code{matrix}.
#'
#' @section Details: {
#'
#' #' A \eqn{k}-dimensional series \eqn{a_{t}}{a_t} can be transformed to a
#' standardized univariate series \eqn{e_{t}}{e_t}:
#'
#' \deqn{e_{t} = a_{t}^{\prime} \sum^{-1} a_{t} - k}{e_t = a'_t * \sum^-1 *
#' a_t - k}
#'
#' where \eqn{\sum} denotes the unconditional covariance matrix of the
#' \eqn{k}-dimensional series \eqn{a_{t}}{a_t}.
#'
#' }
#'
#' @references {
#'
#'   Dufour J. M. & Roy, R. (1985). The \eqn{t} copula and related copulas.
#'   Working Paper. Department of Mathematics, Federal Institute of Technology.
#'
#'   Dufour J. M. & Roy, R. (1986). Generalized portmanteau statistics and tests
#'   of randomness. Communications in Statistics-Theory and Methods, 15:
#'   2953-2972.
#'
#'   Tsay, R. S. (2014). Multivariate Time Series Analysis with R
#'   and Financial Applications. John Wiley. Hoboken, NJ.
#'
#'   Tsay, R. S. (2015). MTS: All-Purpose Toolkit for Analyzing
#'   Multivariate Time Series (MTS) and Estimating Multivariate Volatility
#'   Models. R package version 0.33.
#'
#' }
#'
#' @seealso \code{\link{diag_dufour_roy}} & \code{\link{diag_ljung_box}} for Test
#'   Statistics for the transformed standardized series \code{diag_std_et}
#'
#' @examples
#' # create heteroscedastic data
#' dat <- mgarchBEKK::simulateBEKK(3, 500)
#' eps <- data.frame(eps1 = dat$eps[[1]], eps2 = dat$eps[[2]],
#'                   eps3 = dat$eps[[3]])
#'
#' # transform to standardized univariate series e_t
#' diag_std_et(eps)
#'
#' @export
diag_std_et <- function(x) {

  if(!is.data.frame(x) & !is.matrix(x) & !is.numeric(x)) {
    stop("x needs to be an object of class data.frame, matrix or numeric")
  }

  # save names
  if(is.null(names(x))) {
    x_name <- "et"} else {
    x_name <- paste0("et_of_", paste0(names(x), collapse = "&"))
  }

  # force as matrix
  x <- as.matrix(x)

  at <- scale(x, scale = FALSE)
  et <- as.matrix(apply((at %*% solve(stats::cov(x)) * at), 1, sum) - ncol(x))

  dimnames(et) <- list(NULL, x_name)

  structure(et, class = "std_et")
}


#' Reports whether x is a std_et object
#' @param x An object to test
#' @keywords internal
#' @export
is.std_et <- function(x) inherits(x, "std_et")
