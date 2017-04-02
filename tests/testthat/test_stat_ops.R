context("Test statistical vectorized operations")

test_that("unpaired ttest_vec is same as R t.test", {

  # Unpaired t-test
  X1 <- rnorm(n=30, m=1, sd=1)
  dim(X1) <- c(30, 1)
  X2 <- rnorm(n=40, m=2, sd=5)
  dim(X2) <- c(40, 1)
  res_ttest_r <- t.test(X1, X2)
  res_ttest_vec <- ttest_vec(X1, X2)

  expect_equal(res_ttest_vec$tvalues, res_ttest_r$statistic[["t"]])
  expect_equal(res_ttest_vec$pvalues, res_ttest_r$p.value)
})

test_that("paired ttest_vec is same as R t.test", {

  # Paired t-test
  X1 <- rnorm(n=30, m=1, sd=1)
  dim(X1) <- c(30, 1)
  X2 <- rnorm(n=30, m=2, sd=5)
  dim(X2) <- c(30, 1)
  res_ttest_r <- t.test(X1, X2, paired = TRUE)
  res_ttest_vec <- ttest_vec(X1, X2, paired = TRUE)

  expect_equal(res_ttest_vec$tvalues, res_ttest_r$statistic[["t"]])
  expect_equal(res_ttest_vec$pvalues, res_ttest_r$p.value)
})

test_that("corr_vec is same as R cor.test", {

  X1 <- rnorm(10,sd=2)
  X2 <- 2*X1 + rnorm(10,sd=0.5)
  dim(X1) <- c(10, 1)

  expect_equal(corr_vec(X1, X2)$pvalues, cor.test(X1, X2)$p.value)
  expect_equal(corr_vec(X1, X2)$tvalues, cor.test(X1, X2)$statistic[["t"]])
  expect_equal(corr_vec(X1, X2)$corr_coeff, cor.test(X1, X2)$estimate[["cor"]])

})


