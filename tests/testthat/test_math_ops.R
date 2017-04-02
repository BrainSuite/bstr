context("Test math ops")

test_that("test math ops", {

  expect_equal(log10_transform(1), -log10(1))
  expect_equal(log10_transform(2), -log10(2))
  expect_equal(log10_transform(-5.6), log10(5.6))
  expect_equal(log10_transform(0), -log10(.Machine$double.eps))

})

