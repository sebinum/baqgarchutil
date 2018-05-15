################################################################################
#' Multivariate ARCH test
#'
#' Performs tests to check the conditional heteroscedasticity in a vector
#' time series.
#'
#' @param zt A \code{matrix} of multivariat financial time series. Each
#'   column contains a series, each row an observation of the series.
#' @param m The number of lags of cross-correlation matrices used in the
#'   tests.
#' @return Four different test statistics and their p-values to determine
#'   multivariate ARCH-effect as a \code{data.frame}. For more information
#'   see the Details.
#'
#' @section Details: {
#' The four test statistics are different approaches to detect conditional
#' heteroscedasticity (ARCH-effect) in multivariate time-series as employed
#' by \strong{Ruey S. Tsay (2014)} in \emph{Multivariate Time Series Analysis
#' with R and Financial Applications}. The code is a copy of the function
#' \code{\link[MTS:MarchTest]{MTS::MarchTest}} but returns the results as a
#' \code{data.frame} instead of (only) printing the results to the console.
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
#' \deqn{\sum_{t=1}^{T}(R_{t} - \bar{R}^{2} = T(T^{2} - 1) / 12.}{\sum_{t=1}^T
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
#' \deqn{Q^{*}_{k}(m) = T^{2} \displaystyle\sum_{i=1}^{m} (\hat{\rho}^{-1}_{0}
#'       \otimes \hat{\rho}^{-1}_{0}) b_{i}}{Q*_k(m) = T² * \sum_{i=1}^m *
#'       (\rho_0^-1 \otimes \rho_0^-1) * b_i}
#' where \eqn{T} stands for the sample size, \eqn{k} for the dimension of
#' \eqn{a_{t}}{a_t} and \eqn{b_{i} = vec(\hat{\rho}^{\prime}_{i})}{b_i =
#' vec(\rho'_i)} with \eqn{\hat{\rho}_{j}}{\rho_j} being the
#' lag-\eqn{\hat{\rho}_{j}}{\rho_j} cross-correlation matrix of
#' \eqn{a^{2}_{t}}{a²_t}. As in the univariate case, under the null hypothesis
#' of no conditional heteroscedasticity in \eqn{a_{t}}{a_t},
#' \eqn{Q_{k}^{*}(m)}{Q*_k(m)} is asymptoically distributed as
#' \eqn{\chi^{2}_{k^{2}m}}{\chi²_{k²m}}.
#'
#' \strong{4. \eqn{Q_{k}^{r}(m)}{Q_k^r(m)}}
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
#'   J. M. Dufour & R. Roy (1985). The \eqn{t} copula and related copulas.
#'   Working Paper. Department of Mathematics, Federal Institute of Technology.
#'
#'   J. M. Dufour & R. Roy (1986). Generalized portmanteau statistics and tests
#'   of randomness. Communications in Statistics-Theory and Methods, 15:
#'   2953-2972.
#'
#'   W. K. Li (2004). Diagnostic Checks in Time Series. Chapman & Hall / CRC.
#'   Boca Raton, FL.
#'
#'   R. S. Tsay (2014, Chapter 7). Multivariate Time Series Analysis with R
#'   and Financial Applications. John Wiley. Hoboken, NJ.
#'
#'   R. S. Tsay (2015, Page 26-28). MTS: All-Purpose Toolkit for Analyzing
#'   Multivariate Time Series (MTS) and Estimating Multivariate Volatility
#'   Models. R package version 0.33.
#'
#' }
#'
#' @examples
#' # create heteroscedastic data
#' dat <- mgarchBEKK::simulateBEKK(2, 500)
#' eps <- data.frame(eps1 = dat$eps[[1]], eps2 = dat$eps[[2]])
#'
#' # perform multivariate arch test
#' march_test(eps)
#'
#' @export
march_test = function (zt, m = 10){
  # the MarchTest function from Tsays MTS package
  # modified to output the results as a DF
  if (!is.matrix(zt)) {zt = as.matrix(zt)}

  # dimension of zt (k = series, nT = observations)
  nT = dim(zt)[1]
  k = dim(zt)[2]
  # uncondtional stats::covariance-matrix of zt
  C0 = stats::cov(zt)
  # center at mean = 0
  zt1 = scale(zt, center = TRUE, scale = FALSE)
  # inverse of covariance matrix
  C0iv = solve(C0)
  # matrix-multiplication: centered zt with inverse covariance matrix
  wk = zt1 %*% C0iv
  # element wise multiplication of wk and zt1
  wk = wk * zt1
  rt2 = apply(wk, 1, sum) - k
  m1 = acf(rt2, lag.max = m, plot = F)
  acf = m1$acf[2:(m + 1)]
  c1 = c(1:m)
  deno = rep(nT, m) - c1
  Q = sum(acf^2/deno) * nT * (nT + 2)
  pv1 = 1 - stats::pchisq(Q, m)

  rk = rank(rt2)
  m2 = acf(rk, lag.max = m, plot = F)
  acf = m2$acf[2:(m + 1)]
  mu = -(rep(nT, m) - c(1:m))/(nT * (nT - 1))
  v1 = rep(5 * nT^4, m) - (5 * c(1:m) + 9) * nT^3 + 9 *
    (c(1:m) - 2) * nT^2 + 2 * c(1:m) * (5 * c(1:m) +
        8) * nT + 16 * c(1:m)^2
  v1 = v1/(5 * (nT - 1)^2 * nT^2 * (nT + 1))
  QR = sum((acf - mu)^2/v1)
  pv2 = 1 - stats::pchisq(QR, m)

  x = zt^2
  g0 = stats::var(x)
  ginv = solve(g0)
  qm = 0
  df = 0
  for (i in 1:m) {
    x1 = x[(i + 1):nT, ]
    x2 = x[1:(nT - i), ]
    g = stats::cov(x1, x2)
    g = g * (nT - i - 1)/(nT - 1)
    h = t(g) %*% ginv %*% g %*% ginv
    qm = qm + nT * nT * sum(diag(h))/(nT - i)
    df = df + k * k
  }
  Qkm = qm
  pv3 = 1 - stats::pchisq(qm, df)

  cut1 = stats::quantile(rt2, 0.95)
  idx = c(1:nT)[rt2 <= cut1]
  x = zt[idx, ]^2
  eT = length(idx)
  g0 = stats::var(x)
  ginv = solve(g0)
  qm = 0
  df = 0
  for (i in 1:m) {
    x1 = x[(i + 1):eT, ]
    x2 = x[1:(eT - i), ]
    g = stats::cov(x1, x2)
    g = g * (eT - i - 1)/(eT - 1)
    h = t(g) %*% ginv %*% g %*% ginv
    qm = qm + eT * eT * sum(diag(h))/(eT - i)
    df = df + k * k
  }
  Qkr = qm
  pv4 = 1 - stats::pchisq(qm, df)

  df_results = data.frame(Test = c("Q*(m)",
                                   "Q_R(m)",
                                   "Q*_k(m)",
                                   "Q^r_k(m)"),
                          m = rep(m, 4),
                          Test_statistic = c(Q, QR, Qkm, Qkr),
                          pvalue = c(pv1, pv2, pv3, pv4))
  df_results
}
