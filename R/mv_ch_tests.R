################################################################################
#' Multivariate Conditional Heteroscedasticity (ARCH) Tests
#'
#' Performs tests to check whether conditional heteroscedasticity in a
#' multivariate  time series vector is statistically significant. This is a
#' wrapper function for \code{diag_ljung_box} and
#' \code{diag_dufour_roy}.
#'
#' @param x A \code{matrix} / \code{data.frame} / \code{numeric vector} of
#'   (multivariat) financial time series. Each column contains a series, each
#'   row an observation of the series.
#' @param lags The number of lags of cross-correlation matrices used in the
#'   tests. Can take multiple values. Defaults to \code{lags = c(8, 10, 12)}.
#' @return Four different test statistics and their p-values to determine
#'   multivariate ARCH-effect as a \code{data.frame}. For more information
#'   see the Details.
#'
#' @section Details: {
#'
#' The four test statistics are different approaches to detect conditional
#' heteroscedasticity (ARCH-effect) in multivariate time-series as employed
#' by \strong{Ruey S. Tsay (2014)} in \emph{Multivariate Time Series Analysis
#' with R and Financial Applications}.
#'
#' The \eqn{k}-dimensional series \eqn{a_{t}}{a_t} can be transformed to a
#' standardized univariate series \eqn{e_{t}}{e_t}:
#'
#' \deqn{e_{t} = a_{t}^{\prime} \sum^{-1} a_{t} - k}{e_t = a'_t * \sum^-1 *
#' a_t - k}
#'
#' where \eqn{\sum} denotes the unconditional covariance matrix of the
#' \eqn{k}-dimensional series \eqn{a_{t}}{a_t}.
#'
#' \strong{1. \eqn{Q^{*}(m)}{Q*(m)}: univariate Ljung-Box Test on the
#' standardized series \eqn{e_{t}}{e_t}}
#'
#'   The univariate series \eqn{e_{t}}{e_t} is the basis for the univariate
#'   Ljung-Box Test
#'
#'   \deqn{Q^{*}(m) = T(T + 2) \sum^{m}_{i=1} \hat{\rho}^{2}_{i} / (T - i)}{
#'   Q*(m) = T * (T + 2) * \sum^m_i=1 * \rho²_i / (T - i)}
#'
#'   where \eqn{T} stands for the sample size and \eqn{\hat{\rho}_{i}}{\rho_i}
#'   for the lag-\eqn{i} sample autocorrelation of \eqn{e_{t}}{e_t}.
#'
#'   The Hypothesis \eqn{H_{0} : \rho_{1} = ... = \rho_{m} = 0}{H0 :
#'   \rho_1 = ... = \rho_m = 0} is tested against \eqn{H_{1} : \rho_{i} \neq
#'   0}{H1 : \rho_i != 0} for \eqn{i = (1 \leq i \leq m)}{i = (1 \le i \le m)}.
#'   Under the null hypothesis of no conditional heteroscedasticity in
#'   \eqn{a_{t}}{a_t}, the test statistic \eqn{Q^{*}(m)}{Q*(m)} is
#'   asymptotically distributed as \eqn{\chi^{2}_{m}}{\chi²_m}.
#'
#' \strong{2. \eqn{Q_{R}(m)}{Q_R(m)}: Rank-Based Test on the the ranked
#' standardized series \eqn{e_{t}}{e_t}}
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
#' \eqn{e_{t}}{e_t} shows no signs of serial dependency.
#'
#'
#' \strong{3. \eqn{Q_{k}^{*}(m)}{Q*_k(m)}: multivariate Ljung-Box Test on the
#' \eqn{k}-dimensional series \eqn{a_{t}}{a_t}}
#'
#' The multivariate Ljung-Box Statistic is:
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
#' \strong{4. \eqn{Q_{k}^{r}(m)}{Q_k^r(m)}: robust multivariate Ljung-Box Test on
#' the \eqn{k}-dimensional series \eqn{a_{t}}{a_t}}
#'
#' \eqn{Q^{*}_{k}(m)}{Q*_k(m)} may not perform very well when \eqn{a_{t}}{a_t}
#' has heavy tails. To make the test results more robust, the heavy tails from
#' \eqn{a_{t}}{a_t} are trimmed. This is achieved by removing 5\% of the
#' observations of \eqn{a_{t}}{a_t} corresponding to the upper 5\% quantile from
#' the univariate standardized series \eqn{e_{t}}{e_t} (see the details to
#' \eqn{Q^{*}(m)}{Q*(m)}). \eqn{Q_{k}^{*}(m)}{Q*_k(m)} performed on the
#' trimmed upper 5\% quantile series \eqn{a_{t}}{a_t} is denoted by
#' \eqn{Q_{k}^{r}(m)}{Q_k^r(m)}.
#'
#' }
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
#' @seealso \code{\link{diag_std_et}} for the transformation of a multivariate
#'   financial time series to a standardized scalar series which can be tested
#'   for conditional heteroscedasticity, \code{\link{diag_ljung_box}} for the
#'   Ljung-Box Test statistic, \code{\link{diag_dufour_roy}} for the Rank-Based
#'   Test for serial correlation
#'
#' @examples
#' # create heteroscedastic data
#' dat <- mgarchBEKK::simulateBEKK(3, 150)
#' eps <- data.frame(eps1 = dat$eps[[1]], eps2 = dat$eps[[2]],
#'                   eps3 = dat$eps[[3]])
#'
#' # perform multivariate arch test
#' mv_ch_tests(eps)
#'
#' @export
mv_ch_tests <- function(x, lags = c(8, 10, 12)) {

  if(!is.matrix(x) & !is.data.frame(x) & !is.numeric(x)) {
    stop("x needs to be an object of class data.frame, matrix or numeric")
  }
  # force to matrix so dims can be calculated
  x <- as.matrix(x)

  # transform a multivariate series a_t to a standardized univariate series e_t
  et <- diag_std_et(x)

  # Ljung box on et
  test_lb_uv <- diag_ljung_box(et, lags = lags)
  # Rank-Based Roy & Dufour on et
  test_rd_uv <- diag_dufour_roy(et, lags = lags)
  # multivariate Ljung Box on at
  test_lb_mv <- diag_ljung_box(x, lags = lags, squared = TRUE)

  # mv robust Ljung Box on at / trimmed by upper 5 percent of et
  x_r <- x[(et[, 1] <= stats::quantile(et[, 1], 0.95)), ]
  test_lbr_mv <- diag_ljung_box(x_r, lags = lags, squared = TRUE)

  test_names <- c("Q(m) - uv Ljung-Box Test on standardized series e_t",
                  "QR(m) - uv Rank-Based Test (Dufour & Roy) on rank of e_t",
                  "QK(m) - mv Ljung-Box Test on a_t",
                  "QKr(m) - robust mv Ljung-Box Test on trimmed a_t")

  structure(list(test_lb_uv, test_rd_uv, test_lb_mv, test_lbr_mv),
            "names" = test_names)
}
