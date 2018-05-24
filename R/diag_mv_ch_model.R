################################################################################
#' Multivariate Conditional Heteroscedasticity Model Checking
#'
#' Apply four portmanteau test statistics to check the validity of a fitted
#' multivariate volatility model.
#'
#' @param x A \code{list} of two \code{data.frame}s / \code{matrices}
#'  containting: \cr
#'   Element 1 - eps: A T-by-k matrix of residuals for a k-dimensional asset
#'   return series (return series minus conditional mean / residuals of the
#'   mean-equation). \cr
#'   Element 2 - cnd_h: The fitted volatility matrices (estimated
#'   conditional covariance matrix). The dimension is a T-by-k^2 matrix. \cr
#'   Handles \code{baq_nif} class objects on its own.
#' @param lags The number of lags used in the tests. Defaults to 10.
#' @param baq_err_cor Logical switch if potential negative variances should be
#'   handled for \code{baq_nif} object. This is achieved by removing all
#'   observations up to the last negative variance in descending order from eps
#'   and ccovm.
#' @return Various test statistics and their p-values.
#'
#' @details {
#'
#' For the four test statistics employed check the details section of
#' \code{link{mv_ch_tests}}. For the transformation of \code{eps} & \code{cnd_h}
#' check the details section of \code{\link{diag_std_et_cnd}}. For an
#' comprehensive explanation see the references, especially Tsay (2014).
#'
#'}
#'
#' @references {
#'
#'   Ljung G. & Box G. E. P. (1978). On a measure of lack of fit in time series
#'   models. Biometrika 66: 67-72.
#'
#'   Dufour, J. M. & Roy R. (1985). The \eqn{t} copula and related copulas.
#'   Working Paper. Department of Mathematics, Federal Institute of Technology.
#'
#'   Dufour, J. M. & Roy R. (1986). Generalized portmanteau statistics and tests
#'   of randomness. Communications in Statistics-Theory and Methods, 15:
#'   2953-2972.
#'
#'   Ling, S. & Li, W. K. (1997). Diagnostic checking of nonlinear multivariate
#'   time series with multivariate ARCH errors. Journal of Time Series
#'   Analysis, 18: 447–464.
#'
#'   Tse, Y. K. (2002). Residual-based diagnostics for conditional
#'   heteroscedasticity models. Econometric Journal, 5: 358–373.
#'
#'   Tsay, R. S. (2014). Multivariate Time Series Analysis with R and Financial
#'   Applications. John Wiley. Hoboken, NJ.
#'
#'   Tsay, R. S. (2015). MTS: All-Purpose Toolkit for Analyzing Multivariate
#'   Time Series (MTS) and Estimating Multivariate Volatility Models.
#'   R package version 0.33.
#'
#' }
#'
#' @seealso \code{\link{diag_std_et_cnd}}, \code{\link{diag_dufour_roy}},
#'   \code{\link{diag_ljung_box}}
#'
#' @examples
#' # create heteroscedastic data
#' dat <- mgarchBEKK::simulateBEKK(2, 150)
#' eps <- data.frame(eps1 = dat$eps[[1]], eps2 =dat$eps[[2]])
#'
#' # fit a GARCH model
#' gjr <- mGJR(eps[, 1], eps[, 2])
#'
#' # conditional covariance matrices
#' cnd_h <- matrix(unlist(gjr$H.estimated), ncol = 4, byrow = TRUE)
#'
#' #diagnostics on estimated conditional covariance matrices
#' diag_mv_ch_model(x = list(as.matrix(eps), cnd_h))
#'
#' #alternative for a fitted baq_nif object
#' baq <- baq_nifunction(gjr)
#' diag_mv_ch_model(baq)
#'
#' @export
diag_mv_ch_model <- function (x, lags = c(8, 10, 12), baq_err_cor = TRUE)  {

  if(inherits(x, "baq_nif")) {
    #get conditional covariance matrix/df
    baq_h <- as.matrix(x$baq_h)

    # save active setting for warnings
    par_warn <- getOption("warn")

    # deactivate warnings for if there are now negative values
    options(warn = -1)

    # get max row for negative values, which are to be filtered out
    err_row <- max(max(which(baq_h[, 1] < 0), max(which(baq_h[, 4] < 0))))

    # reset warnings to previous value
    options(warn = par_warn)

    # remove all observations from eps / baq_h up to the last negative
    # variance found (err_row)
    if (all(err_row != -Inf & baq_err_cor)) {
      eps <- x$eps[-(1:err_row), ]
      cnd_h <- baq_h[-(1:err_row), ]
      bln_err_cor <- TRUE
    } else {
      eps <- as.matrix(x$eps)
      cnd_h <- as.matrix(x$baq_h)
      bln_err_cor <- FALSE
    }

  } else {
      eps <- as.matrix(x[[1]])
      cnd_h <- as.matrix(x[[2]])
      bln_err_cor <- FALSE
  }

  # transform a multivariate series a_t and its conditional covariance matrices
  # to a  standardized univariate (scalar) series e_t and a standardized
  # multivariate series
  et <- diag_std_et_cnd(eps = eps, cnd_h = cnd_h)

  # extract series
  et_uv <- et[[1]]
  et_mv <- et[[2]]


  # Ljung box on et univariate
  test_lb_uv <- diag_ljung_box(et_uv, lags = lags)
  # Rank-Based Roy & Dufour on et uv
  test_rd_uv <- diag_dufour_roy(et, lags = lags)
  # multivariate Ljung Box on et multivariate
  test_lb_mv <- diag_ljung_box(et_mv, lags = lags, squared = TRUE)

  # mv robust Ljung Box on et mv / trimmed by upper 5 percent of et uv
  x_r <- et_mv[(et_uv[, 1] <= stats::quantile(et_uv[, 1], 0.95)), ]
  test_lbr_mv <- diag_ljung_box(x_r, lags = lags, squared = TRUE)

  test_names <- c(
    "Q(m) - uv Ljung-Box Test on std_et_cnd_uv",
    "QR(m) - uv Rank-Based Test (Dufour & Roy) on rank of std_et_cnd_uv",
    "QK(m) - mv Ljung-Box Test on std_et_cnd_mv",
    "QKr(m) - mv robust Ljung-Box Test on trimmed std_et_cnd_mv"
  )

  structure(list(test_lb_uv, test_rd_uv, test_lb_mv, test_lbr_mv),
    "names" = test_names)
}
