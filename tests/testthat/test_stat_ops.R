# BrainSuite Statistics Toolbox in R (bstr)
# Copyright (C) 2017 The Regents of the University of California
# Creator: Shantanu H. Joshi, Department of Neurology, Ahmanson Lovelace Brain Mapping Center, UCLA
#
# This program is free software; you can redistribute it and/or modify it under the terms
# of the GNU General Public License as published by the Free Software Foundation; version 2.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License version 2 for more details.
#
# You should have received a copy of the GNU General Public License along with this program;
# if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

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

test_that("lm_vec is same as R lm", {

  data <- read.csv(system.file("extdata/testdata", "design.csv", package="bstr"))
  main_effect <- "Age"
  covariates <- "Sex + Height"

  # Fit model using lm
  lm_full <- lm(formula(sprintf('V1000 ~ %s', paste(main_effect, '+', covariates))), data = data)

  # Fit model using lm_vec
  bst_data <- new("BstrData", getwd(), system.file("extdata/testdata", "design.csv", package = 'bstr'), exclude_col="")

  bst_data@data_array <- as.matrix(data[, "V1000"])
  bst_model <- lm_vec(main_effect = main_effect, covariates = covariates, bst_data)

  # Test for equality of tvalue, beta coefficient, and pvalue
  expect_equal(bst_model@tvalues, summary(lm_full)$coefficients[main_effect, "t value"])
  expect_equal(bst_model@beta_coeff[[main_effect, 1]], summary(lm_full)$coefficients[main_effect, "Estimate"])
  expect_equal(bst_model@pvalues, summary(lm_full)$coefficients[main_effect, "Pr(>|t|)"])
})

test_that("bstr_anova is same as R model comparison using anova", {

  data <- read.csv(system.file("extdata/testdata", "design.csv", package="bstr"))
  main_effect <- "Age"
  covariates <- "Sex + Height"

  # Fit full model using R first
  lm_full <- lm(formula(sprintf('V1000 ~ %s', paste(main_effect, '+', covariates))), data = data)
  # Fit null model using R first
  lm_null <- lm(formula(sprintf('V1000 ~ %s', covariates)), data = data)
  # Compare full and null in R
  R_model_cmp <- anova(lm_full, lm_null)

  bst_data <- new("BstrData", getwd(), system.file("extdata/testdata", "design.csv", package = 'bstr'), exclude_col="")
  bst_data@data_array <- as.matrix(data[, "V1000"])
  # Fit full model using bstr first
  bst_lm_full <- lm_vec(main_effect = main_effect, covariates = covariates, bstr_data = bst_data)
  # Fit null model using bstr first
  bst_lm_null <- lm_vec(main_effect = "", covariates = covariates, bstr_data = bst_data)
  # Compare full and null in bstr
  bst_model <- anova_vec(bst_lm_full, bst_lm_null, bst_data)

  # Test for equality of pvalue, Fstat, and RSS
  expect_equal(bst_model@pvalues, R_model_cmp$`Pr(>F)`[2])
  expect_equal(bst_model@Fstat, R_model_cmp$F[2])
  expect_equal(bst_lm_full@rss, R_model_cmp$RSS[1])
  expect_equal(bst_lm_null@rss, R_model_cmp$RSS[2])
})


