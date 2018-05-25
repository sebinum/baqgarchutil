################################################################################
#' Rank-Based Test for Serial Correlation
#'
#' The Rank-Based Test by Dufour & Roy (1985, 1986) for serial correlation.
#'
#' @param x An object of class \code{std_et} or \code{std_et_cnd}. See
#' \code{\link{diag_std_et}} and \code{\link{diag_std_et_cnd}}.
#' @param lags The number of lags of cross-correlation matrices used in the
#'   tests. Can take multiple values. Defaults to \code{lags = c(8, 10, 12)}.
#' @return The results of the Rank-Based Test by Dufour & Roy as a
#'   \code{data.frame}.
#'
#' @section Details: {
#'
#' Extreme observations in return series (heavy tails) can have pronounced
#' effects on the results of \eqn{Q^{*}(m)}{Q*(m)}. One approach to circumvent
#' the heavy tails problem is the Rank-Based test on the rank series of
#' \eqn{e_{t}}{e_t} by Dufour & Roy (1985, 1986). With \eqn{R_{t}}{R_t} being
#' the rank of \eqn{e_{t}}{e_t}, the lag-\eqn{\ell}{l} rank autocorrelation of
#' \eqn{e_{t}}{e_t} can be defined as
#'
#' \deqn{\tilde{\rho_{\ell}} = \frac{\sum_{t=\ell+1}^{T}(R_{t} - \bar{R})
#' (R_{t-\ell} - \bar{R})}{\sum_{t=1}^{T} (R_{t} - \bar{R})^{2}},
#' \ell = 1, 2, ..., }{\rho_l = (\sum_{t=l+1}^T * (R_t - R) * (R_{t-l} - R)) /
#' (\sum_{t=1}^T * (R_t - R)²) for l = 1, 2, ..., }
#'
#' where
#'
#' \deqn{\bar{R} = \sum_{i=1}^{T}R_{T}/T = (T + 1) / 2,}{R = \sum_{i=1}^T *
#' R_T/ T = (T + 1) / 2,}
#' \deqn{\sum_{t=1}^{T}(R_{t} - \bar{R}^{2}) = T(T^{2} - 1) / 12.}{\sum_{t=1}^T
#' * (R_t - R)² = T * (T² - 1) / 12.}
#'
#' It can be shown that
#'
#' \deqn{E(\tilde{\rho}_{\ell}) = -(T - \ell) / [T(T - 1)]}{E(\rho_l) = -(T - l)
#' / [T(T - 1)]}
#' \deqn{
#' Var(\tilde{\rho}_{\ell}) = \frac{5T^{4} - (5\ell + 9)T^{3} + 9(\ell -
#' 2)T^{2} + 2\ell(5\ell + 8)T + 16\ell^{2}}
#' {5(T -1)^{2}T^{2}(T + 1)}}{
#' Var(\rho_l) = (5T^4 - (5l + 9)T³ + 9 * (l - 2) * T² + 2l * (5l + 8) * T +
#' 16l²) / (5(T - 1)² * T² * (T + 1))}
#'
#' The Test Statistic
#'
#' \deqn{Q_{R}(m) = \displaystyle\sum_{i=1}^{m}\frac{[\tilde{\rho}_{i} -
#' E(\tilde{\rho}_{i})]^{2}}{Var(\tilde{\rho}_{i})}}{Q_R(m) = \sum_{i=1}^m *
#' [([\rho_{i} - E(\rho_i)]²] / Var(\rho_i)]}
#'
#' is distributed as \eqn{\chi^{2}_{m}}{\chi²_m} asymptotically if
#' \eqn{e_{t}}{e_t} is serially uncorrelated.
#'
#' }
#'
#' @references {
#'
#'   Dufour, J. M. & Roy R. (1985). The \eqn{t} copula and related copulas.
#'   Working Paper. Department of Mathematics, Federal Institute of Technology.
#'
#'   Dufour, J. M. & Roy R. (1986). Generalized portmanteau statistics and tests
#'   of randomness. Communications in Statistics-Theory and Methods, 15:
#'   2953-2972.
#'
#'   Li, W. K. (2004). Diagnostic Checks in Time Series. Chapman & Hall / CRC.
#'   Boca Raton, FL.
#'
#'   Tsay, R. S. (2014). Multivariate Time Series Analysis with R
#'   and Financial Applications. John Wiley. Hoboken, NJ.
#'
#'   Tsay, R. S. (2015). MTS: All-Purpose Toolkit for Analyzing Multivariate
#'   Time Series (MTS) and Estimating Multivariate Volatility Models.
#'   R package version 0.33.
#'
#' }
#'
#' @seealso \code{\link{diag_std_et}} & \code{\link{diag_std_et_cnd}} for the
#'   transformation to the input-series for \code{x}
#'
#' @examples
#' # create heteroscedastic data
#' dat <- mgarchBEKK::simulateBEKK(3, 500)
#' eps <- data.frame(eps1 = dat$eps[[1]], eps2 = dat$eps[[2]],
#'                   eps3 = dat$eps[[3]])
#'
#' # transform to standardized residuals
#' et <- diag_std_et(eps)
#'
#' # perform rank based test (dufour & roy)
#' diag_dufour_roy(et)
#'
#' @export
diag_dufour_roy <- function(x, lags = c(8, 10, 12)) {

  if(!is.std_et(x) & !is.std_et_cnd(x)) {
    stop("x needs to be an object of class std_et or std_et_cnd")
  } else if (is.std_et(x)) {
    # return of diag_std_et is a matrix
    et <- x
  } else {
    # return of diag_std_et_cnd contains et_uv and et_mv, here et_uv is needed
    et <- x[[1]]
  }

  # n observations
  n = nrow(et)
  # maximum lag specified
  m = max(lags)

  mean_etr <- function(n, l) {-(n - l) / (n*(n - 1))}
  var_etr <- function(n, l) {
    (5 * n^4 - (5 * l + 9) * n^3 + 9 * (l - 2) * n^2 + 2 * l *
       (5 * l + 8) * n + 16 * l^2) / (5 * (n - 1)^2 * n^2 * (n + 1))
  }

  etr_m <- mean_etr(n, 1:m)
  etr_var <- var_etr(n, 1:m)
  etr_acf <- stats::acf(rank(et), m, plot = FALSE)$acf[-1]

  q_stat <- cumsum((etr_acf - etr_m)^2/etr_var)

  # degrees of freedom
  df <- lags
  # remove negative degrees of freedom
  df[which(df < 0)] <- NA

  # p-value for the given lags m
  p_val <- 1 - stats::pchisq(q_stat[lags], df)

  data.frame(test = rep("Dufour-Roy", length(m)),
             lag = lags,
             statistic = q_stat[lags],
             dof = df,
             pvalue = p_val)

}
