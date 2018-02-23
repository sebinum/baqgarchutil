################################################################################
#' Apply a News Impact Function
#'
#' Apply a News Impact Function to a mGJR class object.
#'
#' @param mgjrclass A mGJR class object.
#' @param epsnames  Custom names for the variables eps1 and eps2. Character
#'   vector of length two. Defaults to c("eps1", "eps2").
#' @param er_grid   The minimum and maximum value for the hypothetical returns
#'   in the News Impact Function. Numerical vector of length one.
#' @param er_grid_by Steps by which to construct the sequence of returns.
#' Numerical vector of length one.
#' @return The applied News Impact Function packaged as a \code{baq_nif} class
#'   object. The News Impact Functions' conditional variance/correlation in
#'   long- (\code{data.frame}) and and wide-format (\code{list} of
#'   \code{matrices}).
#'  The values are defined as:
#' \describe{
#'   \item{eps}{\code{data.frame} with eps1 and eps2 inherited from the
#'   \code{mGJR} class object.}
#'   \item{est.params}{\code{list} of \code{matrices} with the baqGARCH
#'     parameters needed for the News Impact Function.}
#'   \item{baqH_long}{\code{data.frame} containing the News Impact on the
#'   conditional variance/correlation of eps1/eps2 in long-format.}
#'   \item{baqH_wide}{A \code{list} of \code{matrices} with the News Impact on
#'    the conditional variance of eps1/eps2 and conditional correlation in
#'    wide-format.}
#' }
#'
#' @references {
#'   H. Schmidbauer & A. Roesch. Volatility Spillovers Between Crude Oil
#'   Prices. International Conference on Policy Modeling. EcoMod,
#'   Berlin, 2008.
#'
#'   H. Schmidbauer & A. Roesch. Volatility Spillovers Between Crude Oil
#'   Prices and Us Dollar To Euro Exchange Rates. 4th IAEE Asian Conference,
#'   Beijing, 2014.
#' }
#'
#' @examples {
#' # create data
#' eps <- mgarchBEKK::simulateBEKK(2, 100)
#'
#' # fit the model
#' gjr <- mgarchBEKK::mGJR(eps$eps[[1]], eps$eps[[2]])
#'
#' # apply the news impact function to the model
#' nif <- baq_nifunction(gjr)
#' }
#' @export
baq_nifunction <- function(mgjrclass, epsnames = c("eps1", "eps2"),
                            er_grid = 10, er_grid_by = 0.2) {
  # check if input has appropriate class
  if(!inherits(mgjrclass, "mGJR")) {
    stop("mgjrclass needs to be a mGJR class object.")
  }

  # check if input order is 1, 1, 1
  if(!all(mgjrclass$order == c(1, 1, 1))) {
    stop("This function is only implemented for class mGJR objects with
      order = c(1,1,1).")
  }

  # check if type/length of epsnames is correct format
  # if not: override with default
  if(!all(is.character(epsnames) && length(epsnames) == 2)) {
    cat("epsnames needs to be of type character with length = 2.\n")
    cat("typeof(epsnames): ", typeof(epsnames), "\n")
    cat("length(epsnames): ", length(epsnames), "\n")
    epsnames <- c("eps1", "eps2")
    cat("epsnames set to default ", epsnames, "\n")
  }

  # check if er_grid and er_grid_by are correct format
  if(!all(is.numeric(er_grid) | is.numeric(er_grid_by) |
      er_grid<0 | er_grid_by<0)) {
    stop("er_grid and er_grid_by need to be positive values of type numeric.")
  }

  eps <- structure(
    data.frame(eps1 = mgjrclass$eps1, eps2 = mgjrclass$eps2),
    names = epsnames
    )

  p <- structure(
    c(mgjrclass$est.params, list(stats::cov(eps))),
    names = c("C", "A", "B", "G", "w", "uccov")
  )

  # prepare coefficients for nif
  pC <- t(p[["C"]]) %*% p[["C"]]                   # C, constant
  pA <- p[["A"]]                                   # A, ARCH-coefficient
  pB <- t(p[["B"]]) %*% p[["uccov"]] %*% p[["B"]]  # B, GARCH-term
  tG <- p[["G"]]                                   # G, Gamma-coefficient
  tw <- p[["w"]]                                   # w, angle for S-function

  H_nif = function(xy, x2) {
    # xy/x2 stand for the modelled returns of series 1/2
    x <- c(xy, x2)
    S <-  1-0.5*((cos(pi/4+tw)*x[1]+sin(pi/4+tw)*x[2])/sqrt(x[1]^2+x[2]^2)+1)

    pC + t(pA) %*% x %*% t(x) %*% pA + pB + S * (t(tG) %*% x %*% t(x) %*% tG)
  }

  # create xy grid
  xy <- seq(-er_grid, er_grid, er_grid_by)
  xg <- expand.grid(x=xy, y=xy)

  # Alternative: mapply
  # df_H <- as.data.frame(t(mapply(H_nif, xg$x, xg$y))
  df_H <- as.data.frame(t(apply(xg, 1, function(x) H_nif(x[1], x[2]))))

  # bind names and results, remove redundant covariance
  df_H <- cbind(xg, df_H[, -(3)])

  # function to calculate the conditional correlation
  ccor <- function(vec) {
    vec[4]/sqrt(vec[3]*vec[5])
  }

  # replace conditional covariance with conditional correlation
  df_H[, 4] <- ccor(df_H)

  # replace NaNs (at xy/x2 = 0) with the minimum cond. variance for nicer plots
  df_H[(length(df_H[, 1])-1)/2+1, ] <- c(
    0,                         # column: x
    0,                         # column: y
    min(df_H[, 3], na.rm = T), # column: conditional variance series 1
    0,                         # column: conditional correlation
    min(df_H[, 5], na.rm = T)  # column: conditional variance series 2
    )

  # generic names for the conditional variance/correlation
  names(df_H) <- c("x", "y", "z_var1", "z_cor", "z_var2")

  # function to convert df_H from long to wide for r-base plots
  long_to_wide <- function(df, rcnames) {
    output <- as.matrix(unname(stats::reshape(df, idvar = "x",
      timevar = "y", direction = "wide")[, -(1)]))
    colnames(output) <- round(rcnames, 2)
    rownames(output) <- round(rcnames, 2)
    return(output)
  }

  z <- list(
    x_y = xy,
    z_var1 = long_to_wide(df = df_H[, 1:3], rcnames = xy),
    z_var2 = long_to_wide(df = df_H[, c(1:2, 5)], rcnames = xy),
    z_cor = long_to_wide(df = df_H[, c(1:2, 4)], rcnames = xy)
  )

  cat("News Impact Function successfully applied.\n")

  output <- list(
    eps = eps,
    est.params = p,
    baqH_long = df_H,
    baqH_wide = z
  )
  class(output) <- "baq_nif"

  cat("Class attributes are accesible via the following names:\n")
  cat(names(output), "\n")

  return(output)
}

#' Reports whether x is a baq_nif object
#' @param x An object to test
#' @keywords internal
#' @export
is.baq_nif <- function(x) inherits(x, "baq_nif")

#' @export
print.baq_nif <- function(x, ...) {
  n1 <- names(x$eps[1])
  n2 <- names(x$eps[2])
  ni <- structure(
    x$baqH_long[, c(3, 5, 4)],
    names = c(paste0("cvar ", n1), paste0("cvar ", n2), "ccor")
  )

  cat("News Impact Function based on a fitted baqGARCH model.\n\n")
  cat("------------------------------------------------------\n\n")
  cat("Series 1 (eps1): ", n1, "\n")
  cat("Series 2 (eps2): ", n2, "\n")
  cat("------------------------------------------------------\n\n")
  cat("Input parameters for the News Impact Function:\n\n")
  print(x$est.params)
  cat("------------------------------------------------------\n\n")
  cat("Summary statistics:\n\n")
  print(summary(ni))
  cat("------------------------------------------------------\n\n")
  cat("Class attributes are accesible via the following names:\n")
  cat(names(x), "\n")
}
