context("baq_nifunction")

# generate heteroscedastic data
seps <- mgarchBEKK::simulateBEKK(2, 150)

# create mGJR class object
invisible(capture.output(
  my_mgjr <- mgarchBEKK::mGJR(eps1 = seps$eps[[1]], eps2 = seps$eps[[2]])
  ))

# test 1
test_that("baq_nifunction only handles mGJR class objects", {
  expect_error(baq_nifunction(seps))  # data.frame
  expect_error(baq_nifunction(1))     # numeric
  expect_error(baq_nifunction("ten")) # character
  expect_equal(class(baq_nifunction(my_mgjr, quiet = T)), "baq_nif")
})

test_that("baq_nifunctions quiet parameter works as intended", {
  expect_silent(baq_nifunction(my_mgjr, quiet = T))
  expect_output(baq_nifunction(my_mgjr))
})
