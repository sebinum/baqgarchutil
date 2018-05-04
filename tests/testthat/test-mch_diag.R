context("mch_diag")

# generate heteroscedastic data
seps <- mgarchBEKK::simulateBEKK(2, 150)
ueps <- data.frame(eps1 = seps$eps[[1]], eps2 = seps$eps[[2]])

# create mGJR class object
show_failure(invisible(capture.output(
  my_mgjr <- mgarchBEKK::mGJR(ueps$eps1, ueps$eps2)
)))
#conditional covariance matrices
ccovm <- as.data.frame(matrix(unlist(my_mgjr$H.estimated),
                       ncol = 4, byrow = T))

# object to perform diagnostics on
x <- list(ueps, ccovm)

# apply news impact function
my_baq <- baq_nifunction(my_mgjr, quiet = T)

x2 <- list(my_baq$eps, my_baq$baq_h)


################################################################################
test_that("mch_diag handles baq_nif class the same way it does normal input", {
  expect_identical(mch_diag(x), mch_diag(x2))
})


# guarantee variances are positive
ccovm[1, ] <- abs(ccovm[1, ])
ccovm[4, ] <- abs(ccovm[4, ])
