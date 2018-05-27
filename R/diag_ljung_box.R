################################################################################
#' Ljung-Box Test for Serial Correlation
#'
#' The Ljung-Box Test for serial correlation.
#'
#' @param x A \code{matrix} / \code{data.frame} / \code{numeric vector} of
#'   (multivariat) financial time series. Each column contains a series, each
#'   row an observation of the series.
#' @param lags The number of lags of cross-correlation matrices used in the
#'   tests. Can take multiple values. Defaults to \code{lags = c(8, 10, 12)}.
#' @param order If the test is performed on residuals of a fitted time series
#'   model such as (V)ARMA or GARCH \code{order} generally equals the count of
#'   coefficients from the fitted model. Defaults to \code{order = 0}.
#' @param squared Logical switch, if the squared series should be tested.
#'   Defaults to \code{squared = FALSE}.
#' @return The Ljung-Box Test statistic, it's \eqn{p}-value and further
#'   parameters as a \code{data.frame}.
#'
#' @section Details: {
#'
#' The univariate Ljung-Box (1978) Test is denoted by
#'
#' \deqn{Q^{*}(m) = T(T + 2) \sum^{m}_{i=1} \hat{\rho}^{2}_{i} / (T - i)}{
#' Q*(m) = T * (T + 2) * \sum^m_i=1 * \rho²_i / (T - i)}
#'
#' where \eqn{T} stands for the sample size and \eqn{\hat{\rho}_{i}}{\rho_i}
#' for the lag-\eqn{i} sample autocorrelation of \eqn{e_{t}}{e_t}.
#'
#' The Hypothesis \eqn{H_{0} : \rho_{1} = ... = \rho_{m} = 0}{H0 :
#' \rho_1 = ... = \rho_m = 0} is tested against \eqn{H_{1} : \rho_{i} \neq
#' 0}{H1 : \rho_i != 0} for \eqn{i = (1 \leq i \leq m)}{i = (1 \le i \le m)}.
#' Under the null hypothesis of no conditional heteroscedasticity in
#' \eqn{a_{t}}{a_t}, the test statistic \eqn{Q^{*}(m)}{Q*(m)} is
#' asymptotically distributed as \eqn{\chi^{2}_{m}}{\chi²_m}.
#'
#' The generalization of the Ljung-Box Test to the multivariate case,
#' see e.g. Hosking (1980, 1981) is denoted by
#'
#' \deqn{Q^{*}_{k}(m) = T^{2} \displaystyle\sum_{i=1}^{m} b_{i}^{\prime}
#' (\hat{\rho}^{-1}_{0} \otimes \hat{\rho}^{-1}_{0}) b_{i}}{Q*_k(m) = T² *
#' \sum_{i=1}^m * b'_i(\rho_0^-1 \otimes \rho_0^-1) * b_i}
#'
#' where \eqn{T} stands for the sample size, \eqn{k} for the dimension of
#' \eqn{a_{t}}{a_t} and \eqn{b_{i} = vec(\hat{\rho}^{\prime}_{i})}{b_i =
#' vec(\rho'_i)} with \eqn{\hat{\rho}_{j}}{\rho_j} being the
#' lag-\eqn{\hat{\rho}_{j}}{\rho_j} cross-correlation matrix of
#' \eqn{a^{2}_{t}}{a²_t}. As in the univariate case, under the null hypothesis
#' of no conditional heteroscedasticity in \eqn{a_{t}}{a_t},
#' \eqn{Q_{k}^{*}(m)}{Q*_k(m)} is asymptoically distributed as
#' \eqn{\chi^{2}_{k^{2}m}}{\chi²_{k²m}}.
#'
#' }
#'
#' @references {
#'
#'   Ljung, G. & Box, G. E. P  (1978). On a measure of lack of fit in time
#'   series models. Biometrika 66: 67-72.
#'
#'   Hosking, J. R. M. (1980). The multivariate portmanteau statistic. Journal
#'   of the American Statistical Association, 75: 602–607.
#'
#'   Hosking, J. R. M. (1981). Lagrange-multiplier tests of multivariate time
#'   series model. Journal of the Royal Statistical Society,
#'   Series B, 43: 219–230.
#'
#'   Li, W. K. (2004). Diagnostic Checks in Time Series. Chapman & Hall / CRC.
#'   Boca Raton, FL.
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
#' @seealso \code{\link{diag_std_et}} for the transformation of a multivariate
#'   financial time series to a standardized scalar series which can be tested
#'   for conditional heteroscedasticity (ARCH effect) with the Ljung-Box Test,
#'   \code{\link{mv_ch_tests}} for different varieties of ARCH tests,
#'   \code{\link{diag_std_et_cnd}} for the transformation of a multivariate
#'   financial time series to a standardized scalar series and a multivariate
#'   (marginally) standardized series based on fitted conditional covariace
#'   matrices.
#'
#' @examples
#' # create heteroscedastic data
#' dat <- mgarchBEKK::simulateBEKK(3, 500)
#' eps <- data.frame(eps1 = dat$eps[[1]], eps2 = dat$eps[[2]],
#'                   eps3 = dat$eps[[3]])
#'
#' # perform multivariate arch test
#' diag_ljung_box(eps)
#'
#' @export
diag_ljung_box <- function(x, lags = c(8, 10, 12) , order = 0,
                           squared = FALSE) {

  if(!is.matrix(x) & !is.data.frame(x) & !is.numeric(x)) {
    stop("x needs to be an object of class data.frame, matrix or numeric")
  }
  # force to matrix so dims can be calculated
  x <- as.matrix(x)
  # n observations
  n = nrow(x)
  # k series
  k = ncol(x)
  # maximum lag specified
  m = max(lags)

  if(m > (n - 1)) {
    stop("maximum lag specified is bigger than (n-observations - 1): the
         maximum lag allowed is (n-observations - 1)")
  } #should

  if(squared) {x <- x^2}

  # autocorrelation of the squared k series
  ac <- stats::acf(x, m, plot = FALSE)$acf

  # kronecker product of the inverse auto-cross correlation matrix at
  # lag 0 (p_0^-1 otimes p_0^-1 )
  kron <- kronecker(solve(ac[1, , ]), solve(ac[1, , ]))

  # 'vectorized' autocorrelation for each lag (column = lag)
  # drop = F for when m = 1
  b <- t(matrix(ac, ncol = k^2)[-(1), , drop = FALSE])

  lbsum <- function(lag) {
    bi <- b[, lag, drop = FALSE]
    1 / (n - lag) * crossprod(bi, crossprod(kron, bi))}

  eq <- apply(matrix(1:m), 1, lbsum)

  if(k == 1) {
    q_stat <- n * (n + 2) * cumsum(eq)
  } else {
    q_stat <- n^2 * cumsum(eq)
  }

  # degrees of freedom
  df <- k^2 * (lags - order)
  # remove negative degrees of freedom
  df[which(df < 0)] <- NA
  # p-value for the given lags m
  p_val <- 1 - stats::pchisq(q_stat[lags], df)

  # return results as data.frame
  data.frame(test = rep("Ljung-Box", length(m)),
             lag = lags,
             statistic = q_stat[lags],
             dof = df,
             pvalue = p_val)
}
