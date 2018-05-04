################################################################################
#' Multivariate ARCH test
#'
#' Performs tests to check the conditional heteroscedasticity in a vector
#' time series.
#'
#' @param zt  nT-by-k data matrix of a k-dimensional financial time series,
#'   each column contains a series.
#' @param m The number of lags of cross-correlation matrices used in the
#'   tests.
#' @return Various test statistics and their p-values.
#'
#' @references {
#'   Tsay (2014, Chapter 7). Multivariate Time Series Analysis with R and
#'   Financial Applications. John Wiley. Hoboken, NJ.
#'
#'   Ruey S. Tsay (2015, Page 26-28). MTS: All-Purpose Toolkit for Analyzing
#'   Multivariate Time Series (MTS) and Estimating Multivariate Volatility
#'   Models. R package version 0.33.
#'
#' }
#'
#' @examples
#' # create heteroscedastic data
#' dat <- mgarchBEKK::simulateBEKK(2, 500)
#' eps <- data.frame(eps1 = dat$eps[[1]], eps2 =dat$eps[[2]])
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

  df_results = data.frame(Test = c("Q(m) of squared series(LM test)",
                                   "Q_R(m) Rank-based test",
                                   "Q_k(m) of squared series",
                                   "Q_k^r Robust Test(5% upper tail trimming)"),
                          m = rep(m, 4),
                          Test_statistic = c(Q, QR, Qkm, Qkr),
                          pvalue = c(pv1, pv2, pv3, pv4))
  df_results
}
