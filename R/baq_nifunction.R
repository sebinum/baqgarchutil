################################################################################
#' Apply a news impact Function
#'
#' Apply a news impact Function to a mGJR class object.
#'
#' @param x A mGJR class object.
#' @param epsnames  Custom names for the variables eps1 and eps2. Character
#'   vector of length two. Defaults to c("eps1", "eps2").
#' @param er_grid   The minimum and maximum value for the hypothetical returns
#'   in the news impact Function. Numerical vector of length one.
#' @param er_grid_by Steps by which to construct the sequence of returns.
#' Numerical vector of length one.
#' @param ni_long Logical switch if the conditional variance after news impact
#'   should be returned in long format (as \code{data.frame}).
#' @param ni_wide Logical switch if the conditional variance after news impact
#'   should be returned in wide format (as \code{matrix}).
#' @param quiet Logical switch if information about applying baq_nifunction
#'   should be shown in the console.
#'
#' @return The applied news impact Function packaged as a \code{baq_nif} class
#'   object. The news impact Functions' conditional variance/correlation
#'   (optional) as long- (\code{data.frame}) and and wide-format (\code{list}
#'   of \code{matrices}).
#'  The values are defined as:
#' \describe{
#'   \item{series_names}{character \code{vector} with the names for eps1 and
#'     eps2}
#'   \item{eps}{\code{data.frame} with eps1 and eps2 inherited from the
#'     \code{mGJR} class object.}
#'   \item{baq_h}{\code{data.frame} with the estimated conditional covariance
#'     matrices inherited from the \code{mGJR} class object in a more accessible
#'     format.}
#'   \item{coef}{\code{data.frame} with the baqGARCH coefficients
#'     needed for the news impact Function.}
#'   \item{coef_se}{\code{data.frame} with the baqGARCH coefficients
#'     Standard Errors.}
#'   \item{coef_tval}{\code{data.frame} with the baqGARCH coefficients
#'     T-values.}
#'   \item{baq_ni_ccovm_long}{\code{data.frame} containing the news impact on
#'     the conditional variance/correlation of eps1/eps2 in long-format.}
#'   \item{baq_ni_ccovm_wide}{A \code{list} of \code{matrices} with the news
#'     impact on the conditional variance of eps1/eps2 and conditional
#'     correlation in wide-format.}
#' }
#' @section Details:{
#'
#' The news impact on the conditional volatility of a fitted baqGARCH model can
#' be analysed by letting the conditional volatility matrices \eqn{H = (h_{ij})}
#' depend on \eqn{x = (x_{1}, x_{2})}{x = (x_1, x_2)}:
#'
#' \deqn{x \mapsto H(x) = C^{\prime}C + A^{\prime}xx^{\prime}A +
#' B^{\prime}\sum B + S_{w}(x) * \Gamma^{\prime}xx^{\prime}\Gamma}{x -> H(x) =
#' C'C + A'xx'A + B'\sumB + S_w(x) * \Gamma'xx'\Gamma}
#'
#' where \eqn{\sum} is the unconditional covariance matrices of the bivariate
#' time series and \eqn{x} the vector of potential innovations (i.e. returns) in
#' the bivariate series affecting its' conditional volatility.
#'
#' The contour lines (conditional variance after news impact) are based on the
#' functions:
#'
#' \deqn{x \mapsto h_{11}(x), x \mapsto h_{22}(x), x \mapsto h_{12}(x)/
#' \sqrt{h_{11}(x)h_{22}(x)}}{x -> h_{11}(x), x -> h_{22}(x), x -> h_{12}
#' (x)/\sqrt(h_{11}(x) * h_{22}(x)),}
#'
#' where the function \eqn{x_{11}} stands for the news impact on the next day's
#' conditional variance of returns on series 1, \eqn{x_{22}} stands for the news
#' impact on the next day's conditional variance of returns on series 2 and
#' \eqn{h_{12}(x)/\sqrt(h_{11}(x) * h_{22}(x))} for the conditional correlation
#' of returns (series 1 & 2).
#' }
#'
#'
#' @references {
#'   H. Schmidbauer & A. Roesch (2008). Volatility Spillovers Between Crude Oil
#'   Prices. International Conference on Policy Modeling. EcoMod,
#'   Berlin.
#'
#'   H. Schmidbauer & A. Roesch (2014). Volatility Spillovers Between Crude Oil
#'   Prices and Us Dollar To Euro Exchange Rates. 4th IAEE Asian Conference,
#'   Beijing.
#' }
#'
#'
#' @examples
#' # create data
#' eps <- mgarchBEKK::simulateBEKK(2, 100)
#'
#' # fit the model
#' gjr <- mgarchBEKK::mGJR(eps$eps[[1]], eps$eps[[2]])
#'
#' # apply the news impact function to the model
#' nif <- baq_nifunction(gjr)
#'
#' @export
baq_nifunction <- function(x, epsnames = c("series1", "series2"),
                            er_grid = 10, er_grid_by = 0.2,
                            ni_long = TRUE, ni_wide = TRUE, quiet = FALSE) {
  # check if input has appropriate class
  if(!inherits(x, "mGJR")) {
    stop("x needs to be an object of class mGJR.")
  }

  # check if input order is 1, 1, 1
  if(!all(x$order == c(1, 1, 1))) {
    stop("This function is only implemented for objects of class mGJR  with
      order = c(1,1,1).")
  }

  # check if type/length of epsnames is correct format
  # if not: override with default
  if(!all(is.character(epsnames) && length(epsnames) == 2)) {
    cat("epsnames needs to be of type character with length = 2.\n")
    cat("typeof(epsnames): ", typeof(epsnames), "\n")
    cat("length(epsnames): ", length(epsnames), "\n")
    epsnames <- c("series1", "series2")
    cat("epsnames set to default ", epsnames, "\n")
  }

  # check if er_grid and er_grid_by are correct format
  if(!all(is.numeric(er_grid) | is.numeric(er_grid_by) |
      er_grid<0 | er_grid_by<0)) {
    stop("er_grid and er_grid_by need to be positive values of type numeric.")
  }

  # replace space with underscore for df names
  m_names <- gsub(" ", "_", epsnames)

  eps <- structure(
    data.frame(eps1 = x$eps1, eps2 = x$eps2),
    names = paste0(m_names, "_eps")
    )

  res <- structure(
    data.frame(res1 = x$resid1, res2 = x$resid2),
    names = paste0(m_names, "_res")
  )

  par_names <- c("C", "A", "B", "G", "w")

  p <- structure(
    c(x$est.params, list(stats::cov(eps))),
    names = c(par_names, "uccov")
  )

  # prepare coefficients for nif
  pC <- t(p[["C"]]) %*% p[["C"]]                   # C, constant
  pA <- p[["A"]]                                   # A, ARCH-coefficient
  pB <- t(p[["B"]]) %*% p[["uccov"]] %*% p[["B"]]  # B, GARCH-term
  tG <- p[["G"]]                                   # G, Asym.-coefficient
  tw <- p[["w"]]                                   # w, angle for S-function

  H_nif = function(ret_series1, ret_series2) {
    # input is returns / news impact of series 1 & 2
    x <- c(ret_series1, ret_series2)
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

  # function to calculate the conditional correlation for each row of df_H
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
  # function to convert a list of matrices to a data.frame
  par_list_to_df <- function(x, colname_ext = "",
    colname_ext_sep = " ") {
    # replace parameter w with replicated matrix so each list element has the same dims
    x[[5]] <-  matrix(rep(x[[5]], times = 4), nrow = 2, ncol = 2)


    # reformat list of matrices to data frame
    df <- as.data.frame(matrix(unlist(x), nrow = 4, byrow = FALSE))

    #
    len_df <- length(df)

    # rownames are the matrix location of each list element
    df_rn <- c("[1,1]", "[2,1]", "[1,2]", "[2,2]")
    # colnames
    df_cn <- c("C", "A", "B", "G", "w")

    # add name for unconditional covariance, if there are 6 columns
    if (len_df == 6) {
      df_cn <- c(df_cn, "uc_cov")
    }
    # add an extension to the column name, if default = "" was overwritten
    if (colname_ext != "") {
      df_cn <- paste0(df_cn, colname_ext_sep, colname_ext)
    }

    structure(
      df,
      row.names = df_rn,
      names = df_cn
    )
  }
  baq_h <- as.data.frame(matrix(unlist(x$H.estimated), ncol = 4, byrow = T))
  names(baq_h) <- c(paste0("cvar_", m_names[1]),
                    rep("ccovar", 2), paste0("cvar_", m_names[1]))
  coef <- par_list_to_df(p)
  coef_se <- par_list_to_df(x$asy.se.coef,
                                    colname_ext = "SE")
  coef_tval <- coef[, -(6)]/coef_se
  names(coef_tval) <- paste0(par_names, " ", "T-value")

  output <- list(
    series_names = epsnames,
    eps = eps,
    baq_h = baq_h,
    res = res,
    coef = coef,
    coef_se = coef_se,
    coef_tval = coef_tval
  )

  # add conditional variance after news impact in wide / long format
  # if logical switch is set to TRUE (ni_long can use up lots of space)
  if (ni_long) output[["baq_ni_ccovm_long"]] = df_H
  if (ni_wide) output[["baq_ni_ccovm_wide"]] = z

  # set class to baq_nif
  class(output) <- "baq_nif"
  if (!quiet) {
    cat("News impact function successfully applied.\n")
    cat("Class attributes are accesible via the following names:\n")
    cat(names(output), "\n")
  }

  output
}

#' Reports whether x is a baq_nif object
#' @param x An object to test
#' @keywords internal
#' @export
is.baq_nif <- function(x) inherits(x, "baq_nif")

#' Reports whether x is a baq_nif object
#' @param x An object to \code{\link[base]{print}}
#' @keywords internal
#' @export
print.baq_nif <- function(x, ...) {
  n1 <- names(x$eps[1])
  n2 <- names(x$eps[2])
  ni <- structure(
    x$baq_ni_ccovm_long[, c(3, 5, 4)],
    names = c(paste0("cvar ", n1), paste0("cvar ", n2), "ccor")
  )

  cat("news impact Function based on a fitted baqGARCH model.\n")
  cat("------------------------------------------------------\n")
  cat("Series 1 (eps1): ", n1, "\n")
  cat("Series 2 (eps2): ", n2, "\n")
  cat("------------------------------------------------------\n\n")
  cat("GARCH coefficients for the news impact Function:\n")
  print(x$coef)
  cat("\nGARCH coefficients Standard Errors:\n")
  print(x$coef_se)
  cat("\nGARCH coefficients absolute T-values:\n")
  print(abs(x$coef_tval))
  cat("------------------------------------------------------\n\n")
  cat("Class attributes are accesible via the following names:\n")
  cat(names(x), "\n")
}
