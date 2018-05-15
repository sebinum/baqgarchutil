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
#'   Element 2 - est_ccovm: The fitted volatility matrices (estimated
#'   conditional covariance matrix). The dimension is a T-by-k^2 matrix. \cr
#'   Handles \code{baq_nif} class objects on its own.
#' @param m The number of lags used in the tests. Defaults to 10.
#' @param baq_err_cor Logical switch if potential negative variances should be
#'   handled for \code{baq_nif} object. This is achieved by removing all
#'   observations up to the last negative variance in descending order from eps
#'   and ccovm.
#' @param quiet Logical switch if information about applying the mch_diag
#'   function should be shown in the console.
#' @return Various test statistics and their p-values.
#'
#' @details This function is a copy of \code{\link[MTS]{MCHdiag}}. It deviates
#'   from the originial in two ways:
#'
#'  1) it displays the result as a \code{data.frame} and
#'
#'  2) thus makes it possible to save the results to a variable.
#'
#'  The test-statistics are from Tsays(2014) Multivariate Time Series as
#'  follows:
#'
#'  Q(m):
#'
#' @references {
#'   Ruey S. Tsay (2014, Chapter 7). Multivariate Time Series Analysis with
#'   R and Financial Applications. John Wiley. Hoboken, NJ.
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
#' # fit a GARCH model
#' gjr <- mGJR(eps[, 1], eps[, 2])
#'
#' # conditional covariance matrices
#' ccovm <- as.data.frame(matrix(unlist(gjr$H.estimated),
#'                               ncol = 4, byrow = TRUE))
#'
#' #diagnostics on estimated conditional covariance matrices
#' mch_diag(x = list(eps, ccovm))
#'
#' #alternative for a fitted baq_nif object
#' baq <- baq_nifunction(gjr)
#' mch_diag(baq)
#'
#' @export
mch_diag <- function (x, m = 10, baq_err_cor = TRUE, quiet = FALSE)  {

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
      est_ccovm <- baq_h[-(1:err_row), ]
      bln_err_cor <- TRUE
    } else {
    eps <- x$eps
    est_ccovm <- x$baq_h
    bln_err_cor <- FALSE
    }
  } else {
    eps <- x[[1]]
    est_ccovm <- x[[2]]
    bln_err_cor <- FALSE
  }

  if (!is.matrix(eps)) {
    eps <- as.matrix(eps)
  }
  if (!is.matrix(est_ccovm)) {
    est_ccovm <- as.matrix(est_ccovm)
  }

  nT <- dim(eps)[1]
  k <- dim(eps)[2]
  nT1 <- dim(est_ccovm)[1]
  k1 <- dim(est_ccovm)[2]
  if ((nT != nT1) || (k1 != k^2)) {
    stop("Inconsistency in dimensions. Following conditions need to be met:\n",
         "dim(eps)[1] == dim(est_ccovm)[1] & dim(eps)[2]^2 == dim(est_ccovm)[2]")
  }

  et <- NULL
  etat <- NULL
  for (i in 1:nT) {
    Vt <- matrix(est_ccovm[i, ], k, k)
    Vtinv <- solve(Vt)
    x <- matrix(eps[i, ], 1, k)
    tmp <- x %*% Vtinv %*% t(x) - k
    et <- c(et, tmp)
    m1 <- eigen(Vt)
    P <- m1$vectors
    lam <- m1$values
    d1 <- diag(1/sqrt(lam))
    Vthalf <- P %*% d1 %*% t(P)
    wk <- x %*% Vthalf
    etat <- rbind(etat, wk)
  }
  m1 <- acf(et, m, plot = FALSE)
  acf <- m1$acf[2:(m + 1)]
  tmp <- acf^2/c(rep(nT, m) - c(1:m))
  Q1 <- sum(tmp) * nT * (nT + 2)
  pv1 <- 1 - stats::pchisq(Q1, m)
  lag <- m
  mu <- -(rep(nT, lag) - c(1:lag))/(nT * (nT - 1))
  v1 <- rep(5 * nT^4, lag) - (5 * c(1:lag) + 9) * nT^3 +
    9 * (c(1:lag) - 2) * nT^2 + 2 * c(1:lag) * (5 * c(1:lag) +
        8) * nT + 16 * c(1:lag)^2
  v1 <- v1/(5 * (nT - 1)^2 * nT^2 * (nT + 1))
  ret <- rank(et)
  m2 <- acf(ret, m, plot = FALSE)
  acf <- m2$acf[2:(m + 1)]
  Qr <- sum((acf - mu)^2/v1)
  pv2 <- 1 - stats::pchisq(Qr, m)
  x <- etat^2
  g0 <- stats::var(x)
  ginv <- solve(g0)
  qm <- 0
  for (i in 1:lag) {
    x1 <- x[(i + 1):nT, ]
    x2 <- x[1:(nT - i), ]
    g <- stats::cov(x1, x2)
    g <- g * (nT - i - 1)/(nT - 1)
    h <- t(g) %*% ginv %*% g %*% ginv
    qm <- qm + nT * nT * sum(diag(h))/(nT - i)
  }
  QKm <- qm
  pv3 <- 1 - stats::pchisq(QKm, k^2 * m)

  q95 <- stats::quantile(et, 0.95)
  idx <- c(1:nT)[et <= q95]
  x <- etat[idx, ]^2
  eT <- length(idx)
  g0 <- stats::var(x)
  ginv <- solve(g0)
  qm <- 0
  for (i in 1:lag) {
    x1 <- x[(i + 1):eT, ]
    x2 <- x[1:(eT - i), ]
    g <- stats::cov(x1, x2)
    g <- g * (eT - i - 1)/(eT - 1)
    h <- t(g) %*% ginv %*% g %*% ginv
    qm <- qm + eT * eT * sum(diag(h))/(eT - i)
  }
  Qrm <- qm
  pv4 <- 1 - stats::pchisq(Qrm, k^2 * m)

  df_results <- data.frame(
    Test_name = c("Q(m) of et", "Rank-based test", "Qk(m) of et", "Robust Qk(m)"),
    Test_statistic = c(Q1, Qr, QKm, Qrm),
    m = rep(m, 4),
    pvalue = c(pv1, pv2, pv3, pv4)

    )
  if(bln_err_cor) {
    df_results <- cbind(df_results , exempt_rows = rep(paste0("1:",err_row), 4))
  }

  df_results
}
