################################################################################
#' Plot News Impact
#'
#' Plot the News Impact on the conditional variance / correlation of a baq_nif
#' class object.
#'
#' @param baq_nif A baq_nif class object.
#' @param ni_plots A character vector specifying which conditional correlation
#'   / variance to plot. Can contain any combination of "var1", "var2" and
#'   "cor". Defaults to c("var1", "var2", "cor").
#' @param plotdir A character vector of length one to determine if the plots
#'   should be created in a row ("horizontal") or a column ("vertical").
#'   Defaults to "horizontal".
#' @return A plot of the News Impact on the conditional variance / correlation.
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
#' # plot the results
#' baq_niplot(nif)
#'
#' # plot only the the first series
#' baq_niplot(nif, ni_plots = "var1")
#'
#' \dontrun{
#'
#' # save plot as pdf
#' pdf("baqplot_example.pdf", width = 15, height = 7.2)
#' baq_niplot(nif)
#' dev.off()
#' }
#' @export
baq_niplot <- function(baq_nif, ni_plots = c("var1", "var2", "cor"),
  plotdir = "horizontal") {

  if(!inherits(baq_nif, "baq_nif")) {
    stop("baq_nif needs to be a baq_nif class object.")
  }

  plotpar <- c("var1", "var2", "cor")
  pcount <- length(ni_plots)
  pdir <- list(
    "horizontal" = c(1, pcount),
    "vertical" = c(pcount, 1)
    )

  if(!all(ni_plots %in% plotpar)) {
    stop("Only 'var1', 'var2' and 'cor' are valid entries for the parameter
      ni_plots.")
  }

  if(!all(plotdir %in% names(pdir) && length(plotdir) == 1)) {
    cat("Only 'horizontal' OR 'vertical' are valid entries for the parameter
      plotdir.\n")
    cat("Input ignored and set to 'horizontal'.\n")
    plotdir <- "horizontal"
  }

  oldoptions <- options(stringsAsFactors = FALSE)
  on.exit(options(oldoptions), add = TRUE)

  epsnames <- names(baq_nif$eps)

  # the color palette from the tim.colors function from the 'fields' package
  tim.colors <- c(
      "#00008F", "#00009F", "#0000AF", "#0000BF", "#0000CF",
      "#0000DF", "#0000EF", "#0000FF", "#0010FF", "#0020FF",
      "#0030FF", "#0040FF", "#0050FF", "#0060FF", "#0070FF",
      "#0080FF", "#008FFF", "#009FFF", "#00AFFF", "#00BFFF",
      "#00CFFF", "#00DFFF", "#00EFFF", "#00FFFF", "#10FFEF",
      "#20FFDF", "#30FFCF", "#40FFBF", "#50FFAF", "#60FF9F",
      "#70FF8F", "#80FF80", "#8FFF70", "#9FFF60", "#AFFF50",
      "#BFFF40", "#CFFF30", "#DFFF20", "#EFFF10", "#FFFF00",
      "#FFEF00", "#FFDF00", "#FFCF00", "#FFBF00", "#FFAF00",
      "#FF9F00", "#FF8F00", "#FF8000", "#FF7000", "#FF6000",
      "#FF5000", "#FF4000", "#FF3000", "#FF2000", "#FF1000",
      "#FF0000", "#EF0000", "#DF0000", "#CF0000", "#BF0000",
      "#AF0000", "#9F0000", "#8F0000", "#800000")

  graphics::par(pty = 's',
    mfrow = pdir[[plotdir]]
  )

  xy <- baq_nif$baqH_wide[["x_y"]]
  xy_length <- length(xy)

  # lookup table for the titles of the News Impact plots
  pmain <- c(
    "var1" = paste0("News Impact on the variance of ", epsnames[1]),
    "var2" = paste0("News Impact on the variance of ", epsnames[2]),
    "cor" = paste("News Impact on the correlation of", epsnames[1],
      epsnames[2], sep = " ")
    )

  # lookup table for the colors used in the News Impact plots
  pcol <- list(
    "var1" = rev(grDevices::heat.colors(xy_length)),
    "var2" = rev(grDevices::heat.colors(xy_length)),
    "cor" = pcol <- tim.colors
  )

  # create the j specified plots in ni_plots
  for (j in ni_plots) {
    z <- paste0("z_", j)
    graphics::image(x = xy, y = xy, z = baq_nif$baqH_wide[[z]], col = pcol[[j]],
      main = pmain[[j]],
      xlab = epsnames[1], ylab = epsnames[2])
    graphics::abline(h = 0, v = 0, lty = 1, col = "black")
    graphics::contour(x = xy, y = xy, z = baq_nif$baqH_wide[[z]], add = TRUE)
    graphics::box()
  }
}
