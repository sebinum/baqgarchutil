################################################################################
#' Standardized Series \eqn{e_{t}}{e_t} for Diagnostics
#'
#' Transform a (multivariate) \eqn{k}-dimensional series \eqn{a_{t}}{a_t} and
#' its corresponding conditional variance to a standardized univariate series
#' \eqn{e_{t}}{e_t} and a standardized multivariate series
#' \eqn{e_{t}^{k}}{e_t^k}.
#'
#' @param eps A \code{matrix} / \code{data.frame} / \code{numeric vector} of
#'   (multivariat) financial time series. Each column contains a series, each
#'   row an observation of the series.
#' @param cnd_h A \code{matrix} / \code{data.frame} with the conditional
#'   variance of eps. The dimensions of \code{cnd_h} need to fullfill the
#'   following conditions: \code{nrow(cnd_h) == nrow(eps)} and
#'   \code{ncol(cnd_h) == ncol(eps)^2}.
#' @return A list containing the standardized univariate conditional series
#'   \eqn{e_{t}}{e_t} and the (marginally) standardized multivariate
#'   conditional series \eqn{\hat{e}_{t}}{e_t}.
#'
#' @section Details: {
#'
#' Ling & Li (1997) proposed model diagnostics based on a standardized scalar
#' series \eqn{\hat{e}_{t}}{ê_t}.
#'
#' The residuals of a \eqn{k}-dimensional financial time series
#' \eqn{\hat{a}_{t} = z_{t} - \hat{\mu}}{a_t = z_t - \mu} (where
#' \eqn{z_{t}}{z_t} stands for the return series and \eqn{\hat{\mu}}{\mu} for
#' the return series conditional mean) and it's time conditional covariance
#' matrices can be transformed to a quadratic standardized residual series
#' \eqn{\hat{e}_{t}}{ê_t}:
#'
#' \deqn{\hat{e}_{t} = \hat{a}_{t}^{\prime} \widehat{\textstyle{\sum}}_{t}^{-1}
#' \hat{a}_{t} - k}{ê_t = â'_t * \sum_t^-1 * â_t}
#'
#' where \eqn{\hat{\sum}_{t}}{\sum_t} denotes the estimated conditional
#' covariance matrices of the \eqn{k}-dimensional series \eqn{a_{t}}{a_t}. The
#' lag-\eqn{\ell}{l} autocorrelation of \eqn{\hat{e}_{t}}{ê_t} is denoted by
#'
#' \deqn{\hat{\rho}_{\ell} = \frac{\sum_{t=\ell+1}^{T} (\hat{e}_{t} - k)
#' (\hat{e}_{t-\ell} - k)}{\sum_{t=1}^{T} (\hat{e}_{t} - k)^{2}}}{\rho_l =
#' (\sum_{t=l+1}^T * (ê_t - k) * (ê_{t-l} - k)) / (\sum_{t=1}^T (ê_t - k)²}
#'
#' where \eqn{E(\hat{a}_{t}^{\prime} \hat{\sum}_{t}^{-1} \hat{a}_{t})
#' = k}{E(â'_t \sum_t^{-1} â_t) = k} (for more details see the references). The
#' series returned is \eqn{\hat{e}_{t} - k}{ê_t - k} and thus can be used to
#' compute the autocorrelation.
#'
#' An approach focusing on the squared elements of a (marginally) multivariate
#' standardized series was proposed by Tse (2002), where the \eqn{i}th
#' standardized residual is denoted by
#'
#' \deqn{\hat{\eta}_{it} = \frac{\hat{a}_{it}}{\sqrt{\hat{\sigma}_{ii,t}}},
#' i = 1, ..., k}{\eta_{it} = â_{it} / \sqrt{\sigma_{ii,t}},
#' i = 1, ..., k}
#'
#' \eqn{\hat{\sigma}_{ii,t}}{\sigma_{ii,t}} stands for the \eqn{(i,i)}th element
#' of the time depenent conditional covariance matrices
#' \eqn{\hat{sum}_{t}}{sum_t} and  \eqn{\hat{a}_{t} = z_{t} - \hat{\mu}}{a_t =
#'  z_t - \mu} (where again \eqn{z_{t}}{z_t} stands for the return series and
#'  \eqn{\hat{\mu}}{\mu} for the return series conditional mean).
#'
#' }
#'
#' @references {
#'
#'   Ling, S. & Li, W. K. (1997). Diagnostic checking of nonlinear multivariate
#'   time series with multivariate ARCH errors. Journal of Time Series
#'   Analysis, 18: 447–464.
#'
#'   Tse, Y. K. (2002). Residual-based diagnostics for conditional
#'   heteroscedasticity models. Econometric Journal, 5: 358–373.
#'
#'   Tsay, R. S. (2014). Multivariate Time Series Analysis with R
#'   and Financial Applications. John Wiley. Hoboken, NJ.
#'
#'   Tsay, R. S (2015). MTS: All-Purpose Toolkit for Analyzing
#'   Multivariate Time Series (MTS) and Estimating Multivariate Volatility
#'   Models. R package version 0.33.
#'
#' }
#'
#' @seealso \code{\link{diag_dufour_roy}}, \code{\link{diag_ljung_box}},
#'   \code{\link{diag_mv_ch_model}} for Test Statistics which can be used on
#'   the series \code{diag_std_et_cnd}
#'
#' @examples
#'
#' # create data
#' eps <- mgarchBEKK::simulateBEKK(2, 150)
#'
#' # fit the model
#' gjr <- mgarchBEKK::mGJR(eps$eps[[1]], eps$eps[[2]])
#'
#' # apply the news impact function to the model
#' nif <- baq_nifunction(gjr)
#'
#' # get the standardized series
#' et_cnd <- diag_std_et_cnd(eps = nif$eps, cnd_h = nif$baq_h)
#'
#' @export
diag_std_et_cnd <- function(eps, cnd_h) {

  if(!is.data.frame(eps) & !is.matrix(eps) & !is.numeric(eps)) {
    stop("eps needs to be an object of class data.frame, matrix or numeric")
  }

  if(!is.data.frame(cnd_h) & !is.matrix(cnd_h) & !is.numeric(cnd_h)) {
    stop("cnd_h needs to be an object of class data.frame, matrix or numeric")
  }
  # force to matrix
  eps <- as.matrix(eps)
  cnd_h <- as.matrix(cnd_h)

  # dimensions
  n <- nrow(eps)
  k <- ncol(eps)

  if(nrow(cnd_h) != n | ncol(cnd_h) != k^2) {
    stop("\ncondition nrow(eps) == nrow(cnd_h) or\n",
         "condition ncol(eps)^2 == ncol(cnd_h) not fullfilled")
  }

  # save names
  eps_names <- dimnames(eps)[[2]]
  cndh_names <- dimnames(cnd_h)[[2]]

  et_univariate <- function(i) {
    hinv <- solve(matrix(cnd_h[i, ], ncol = k))
    meps <- matrix(eps[i, ])
    t(meps) %*% hinv %*% meps - k
  }

  et_uv <- as.matrix(apply(matrix(1:n), 1, et_univariate))

  et_multivariate <- function(i) {
    heig <- eigen(matrix(cnd_h[i, ], ncol = k))
    meps <- matrix(eps[i, ], nrow = 1)
    deig <- diag(1/sqrt(heig$values))
    meps %*% (heig$vectors %*% deig %*% t(heig$vectors))
  }

  et_mv <- t(apply(matrix(1:n), 1, et_multivariate))

  structure(list(et_uv, et_mv), class = "std_et_cnd")
}

#' Reports whether x is a std_et_cnd object
#' @param x An object to test
#' @keywords internal
#' @export
is.std_et_cnd <- function(x) inherits(x, "std_et_cnd")
