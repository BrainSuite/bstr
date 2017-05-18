#' Perform analysis of variance (ANOVA) for brain imaging data.
#'
#' This function accepts a \code{main_effect} and a set of covariates (using the R formula notation) and uses an
#' F-test to compare the full model including the \code{main_effect + covariates} with the reduced (null) model
#' that only includes the \code{covariates}.
#'
#' Slightly different from the standard R anova function, \code{bss_anova} currently does not directly accept the
#' results from \code{lm_vec}. This could be accomodated in the future versions.

#' @param main_effect Character string containing an independent variable whose effect you want to measure.
#' It could be disease status, age, gender etc. This should strictly be a single variable. This can be
#' either a categorical or a continuous variable.
#' @param covariates Character string containing a set of other predictors (variables) in the model. If more than
#' one covariates are included, they should be separated by a \code{+} operator similar to an R formula.
#' @param  bss_data Object of type \code{\link{BssData}}
#' @param  mult_comp method for multiple comparisons correction. The default method is "fdr". See \code{\link{bss_p_adjust}} for valid values.
#' @seealso \code{\link{lm_vec}} for linear regression, \code{\link{bss_ttest}} for independent sample and paired t-tests.
#'
#' @export
bss_anova <- function(main_effect="", covariates="", bss_data, mult_comp="fdr", niter=5000) {

  if (class(bss_data) == "BssROIData") {
    return(bss_roi_anova(main_effect = main_effect, covariates = covariates, bss_data = bss_data))
  }

  message('Running the statistical model. This may take a while...', appendLF = FALSE)
  bss_lm_full <- lm_vec(main_effect = main_effect, covariates = covariates, bss_data = bss_data)
  bss_lm_null <- lm_vec(main_effect = "", covariates = covariates, bss_data = bss_data)
  bss_model <- anova_vec(bss_lm_full, bss_lm_null, bss_data)

  switch(mult_comp,
    perm={
      cl <- parallel::makeCluster(parallel::detectCores())
      registerDoParallel(cl)
      options(warn=-1)
      pvalue_and_nulldist <- maxTperm(main_effect = main_effect, covariates = covariates, bss_data = bss_data, niter)
      bss_model@pvalues <- pvalue_and_nulldist[[1]]
      bss_model@pvalues[is.nan(bss_model@pvalues)] <- 1
      bss_model@pvalues <- bss_model@pvalues*bss_model@tvalues_sign
      bss_model@tvalues[abs(bss_model@pvalues) >= 0.05] <- 0
      nulldist <- pvalue_and_nulldist[[2]]
      bss_model@pvalues_adjusted <- perm_p_adjust(main_effect = main_effect, covariates = covariates, bss_data, nulldist)
      options(warn=0)
      stopCluster(cl)
    },
    fdr={
      bss_model@pvalues[is.nan(bss_model@pvalues)] <- 1
      bss_model@pvalues <- bss_model@pvalues*bss_model@tvalues_sign
      bss_model@tvalues[abs(bss_model@pvalues) >= 0.05] <- 0
      bss_model@pvalues_adjusted <- bss_p_adjust(bss_model@pvalues, mult_comp)
    }
  )

  message('Done.')
  return(bss_model)
}

#' A vectorized version of analysis of variance (ANOVA).
#'
#' This function compares results of model fitting after \code{\link{lm_vec}}.
#' It accepts a full model and and a reduced model and compares them using an F-test.
#' For most scenarios, the user does not need to call this function directly. This function
#' will be called internally from \code{\link{bss_anova}}
#'
#' @param bss_lm_full An object of type \code{\link{BssModel}} returned from \code{\link{lm_vec}}.
#' This is a full model including both the main effect and covariates.
#' @param bss_lm_null An object of type \code{\link{BssModel}} returned from \code{\link{lm_vec}}.
#' This is a null model including only the covariates.
#' @param  bss_data Object of type \code{\link{BssData}}
#'
#' @seealso \code{\link{bss_anova}} for most commonly used function for ANOVA, \code{\link{lm_vec}} for vectorized linear regression, \code{\link{ttest_vec}} for
#' vectorized independent sample and paired t-tests.
#'
#' @export
anova_vec <- function(bss_lm_full, bss_lm_null, bss_data) {

  N <- nrow(bss_data@data_array)
  Fstat <- (bss_lm_null@rss - bss_lm_full@rss)/bss_lm_full@rss * (N - bss_lm_full@Npfull - 1)/(bss_lm_full@Npfull - bss_lm_null@Npnull)  # F statistic

  model_unique_idx <- which(bss_lm_null@unique %in% bss_lm_null@fullvars) + 1    # Add 1, because the first column in the design matrix is the intercept

  se_full_unique <- sqrt(diag(solve(t(bss_lm_full@X_design_full) %*% bss_lm_full@X_design_full)))[model_unique_idx] *
    sqrt(bss_lm_full@rss / (N - bss_lm_full@Npfull - 1))

  tvalues <- bss_lm_full@beta_coeff[model_unique_idx, ]/(se_full_unique + .Machine$double.eps)
  pvalues <- 1 - pf(Fstat, bss_lm_full@Npfull - bss_lm_null@Npnull, N - bss_lm_full@Npfull - 1)
  tvalues_sign <- (bss_lm_full@beta_coeff[model_unique_idx, ] + .Machine$double.eps)/(abs(bss_lm_full@beta_coeff[model_unique_idx, ]) + .Machine$double.eps)

  bss_model <- new("BssModel", model_type="bss_anova", main_effect = bss_lm_full@main_effect, covariates = bss_lm_full@covariates,
                   demographics = bss_data@demographics, mspec_file="")
  bss_model@pvalues <- pvalues
  bss_model@tvalues <- tvalues
  bss_model@tvalues_sign <- tvalues_sign
  bss_model@se <- se_full_unique
  bss_model@Fstat <- Fstat
  return(bss_model)
}

# linear regression
bss_lm <- function(main_effect="", covariates="", bss_data) {

  # if (class(bss_data) == "BssROIData") {
  #   return(bss_roi_lm(main_effect = main_effect, covariates = covariates, bss_data = bss_data))
  # }
  #
  # # Check the model type and call the appropriate method
  # bss_model <- new("BssModel", model_type="bss_lm", main_effect = main_effect, covariates = covariates,
  #                  demographics = bss_data@demographics, mspec_file="")
  # message('Running the statistical model. This may take a while...', appendLF = FALSE)
  # bss_lm_full <- lm_vec(main_effect = main_effect, covariates = covariates, bss_data = bss_data)
  # browser()
  # bss_lm_null <- lm_vec(main_effect = "", covariates = covariates, bss_data = bss_data)
  #
  # Xtemp <- solve(t(bss_model@X_design_full) %*% bss_model@X_design_full) %*% t(bss_model@X_design_full) # Pre Hat matrix
  # beta_full <- Xtemp %*% bss_data@data_array  # beta coefficients
  # y_full <- bss_model@X_design_full %*% beta_full  # Predicted response
  # RSS_full <- colSums((bss_data@data_array - y_full)^2)
  #
  # Xtemp <- solve(t(bss_model@X_design_null) %*% bss_model@X_design_null) %*% t(bss_model@X_design_null) # Pre Hat matrix
  # beta_null <- Xtemp %*% bss_data@data_array  # beta coefficients
  # y_null <- bss_model@X_design_null %*% beta_null  # Predicted response
  # RSS_null <- colSums((bss_data@data_array - y_null)^2)
  #
  # N <- nrow(bss_data@data_array)
  # Fstat <- (RSS_null - RSS_full)/RSS_full * (N - bss_model@Npfull - 1)/(bss_model@Npfull - bss_model@Npnull)  # F statistic
  # model_unique_idx <- which(bss_model@unique %in% bss_model@fullvars) + 1    # Add 1, because the first column in the design matrix is the intercept
  #
  # se_full_unique <- sqrt(diag(solve(t(bss_model@X_design_full) %*% bss_model@X_design_full)))[model_unique_idx] *
  #   sqrt(RSS_full / (N - bss_model@Npfull - 1))
  #
  # tvalue_sign <- (beta_full[model_unique_idx, ] + .Machine$double.eps)/(abs(beta_full[model_unique_idx, ]) + .Machine$double.eps)
  #
  # pvalues <- 1 - pf(Fstat, bss_model@Npfull - bss_model@Npnull, N - bss_model@Npfull - 1)
  #
  # pvalues[is.nan(pvalues)] <- 1
  #
  # pvalues <- pvalues*tvalue_sign
  # tvalues <- beta_full[model_unique_idx, ]/(se_full_unique + .Machine$double.eps)
  # bss_model@pvalues <- pvalues
  # bss_model@tvalues <- tvalues
  # bss_model@tvalues[abs(pvalues) >= 0.05] <- 0
  # bss_model@pvalues_adjusted <- p.adjust(bss_model@pvalues, 'BH')
  # message('Done.')
  # return(bss_model)
}

#' Vectorized linear regression for brain imaging phenotypes.
#'
#' This function accepts a \code{main_effect} and a set of covariates (using the R formula notation) and performs
#' a linear regression including \code{main_effect + covariates}.
#'
#' Slightly different from the standard R \code{lm} function, \code{lm_vec} currently does not directly accept an R formula.
#' This could be accomodated in the future versions.
#' Also currently, this function returns the p-values and the t-statistics for the \code{main_effect}
#' only. Returning the statistics for all variables could be accomodated in the future versions.
#' @param main_effect Character string containing an independent variable whose effect you want to measure.
#' It could be disease status, age, gender etc. This should strictly be a single variable. This can be
#' either a categorical or a continuous variable.
#' @param covariates Character string containing a set of other predictors (variables) in the model. If more than
#' one covariates are included, they should be separated by a \code{+} operator similar to an R formula.
#' @param  bss_data Object of type \code{\link{BssData}}
#'
#' @export
lm_vec <- function(main_effect = "", covariates = "", bss_data) {

  bss_model <- new("BssModel", model_type="bss_lm", main_effect = main_effect, covariates = covariates,
                   demographics = bss_data@demographics, mspec_file="")

  # lm_formula <- formula(sprintf('~ %s', paste(main_effect, '+', covariates)))
  lm_formula <- bss_model@lm_formula

  # Fit model
  N <- dim(bss_data@data_array)[1]
  Np <- length(unlist(strsplit(as.character(lm_formula)[2], '\\+')))

  X <- model.matrix(lm_formula, data = bss_data@demographics)
  X_hat <- solve(t(X) %*% X) %*% t(X) # pre hat matrix
  beta_coeff <- X_hat %*% bss_data@data_array  # beta coefficients
  Y <- X %*% beta_coeff  # predicted response
  rss <- colSums((bss_data@data_array - Y)^2) # residual sum of squares

  if (main_effect == "")
    main_effect = "(Intercept)" # If main_effect is empty, return the parameters of the Intercept

  se <- sqrt(diag(solve(t(X) %*% X)))[[main_effect]] * sqrt(rss / (N-Np-1)) # standard error
  tvalues <- as.numeric(beta_coeff[main_effect, ]/(se + .Machine$double.eps)) # tvalue
  pvalues <- 2*pt(abs(tvalues), N-Np-1, lower.tail = FALSE) # pvalue
  residuals <- bss_data@data_array - Y
  bss_model@pvalues <- pvalues
  bss_model@tvalues <- tvalues
  bss_model@beta_coeff <- beta_coeff
  bss_model@rss <- rss
  bss_model@residuals <- residuals
  return(bss_model)
}

bss_roi_anova <- function(main_effect="", covariates="", bss_data=bss_data) {
  # Check the model type and call the appropriate method
  bss_model <- new("BssModel", model_type="bss_lm", main_effect = main_effect, covariates = covariates,
                   demographics = bss_data@demographics, mspec_file="")
  message('Running the statistical model. This may take a while...', appendLF = FALSE)

  bss_data@demographics[paste('ROI_', as.character(bss_data@roiid), sep = '' )]

  cmd1 <- sprintf("lm_full <- lm(%s, data = bss_data@demographics)",
                  paste('ROI_', as.character(bss_data@roiid), ' ~ ', bss_model@fullmodel, sep = ''))
  cmd2 <- sprintf("lm_null <- lm(%s, data = bss_data@demographics)",
                  paste('ROI_', as.character(bss_data@roiid), ' ~ ', bss_model@nullmodel, sep = ''))
  cmd3 <- "pander::pander(anova(lm_full, lm_null))"

  stats_commands <- c(cmd1, cmd2, cmd3)

  for (cmd in stats_commands) {
    eval(parse(text = cmd))
  }
  bss_model@stats_commands <- stats_commands

  return(bss_model)

}

bss_roi_lm <- function(main_effect="", covariates="", bss_data=bss_data) {
  # Check the model type and call the appropriate method
  bss_model <- new("BssModel", model_type="bss_lm", main_effect = main_effect, covariates = covariates,
                   demographics = bss_data@demographics, mspec_file="")
  message('Running the statistical model. This may take a while...', appendLF = FALSE)

  bss_data@demographics[paste('ROI_', as.character(bss_data@roiid), sep = '' )]

  cmd1 <- sprintf("lm_full <- lm(%s, data = bss_data@demographics)",
                  paste('ROI_', as.character(bss_data@roiid), ' ~ ', bss_model@fullmodel, sep = ''))
  cmd2 <- sprintf("lm_null <- lm(%s, data = bss_data@demographics)",
                  paste('ROI_', as.character(bss_data@roiid), ' ~ ', bss_model@nullmodel, sep = ''))
  cmd3 <- "pander::pander(anova(lm_full, lm_null))"

  stats_commands <- c(cmd1, cmd2, cmd3)

  for (cmd in stats_commands) {
    eval(parse(text = cmd))
  }
  bss_model@stats_commands <- stats_commands

  return(bss_model)

}

#' Test for Correlation between a variable \code{corr_var} and a brain imaging phenotype.
#'
#' Test for correlation between a brain imaging phenotype (cortical thickness, determinant
#' of the jacobian matrix) and \code{corr_var} using the Pearson's product moment correlation
#' coefficient. The brain imaging phenotype is automatically selected from the type of \code{bss_data}.
#' @param corr_var Character variable name. This should be present in the demographics csv file associated
#' with \code{bss_data}.
#' @param  bss_data Object of type \code{\link{BssData}}.
#' @param  mult_comp method for multiple comparisons correction. The default method is "fdr". See \code{\link{bss_p_adjust}} for valid values.
#' @details
#' \code{bss_data} can be of the type "cbm", "tbm", or "roi".
#'
#' @export
bss_corr <- function(corr_var, bss_data, mult_comp="fdr") {

  message('Running correlations...', appendLF = FALSE)
  bss_model <- new("BssModel", model_type="bss_corr", corr_var = corr_var,
                   demographics = bss_data@demographics, mspec_file="")

  corr_result <- corr_vec(bss_data@data_array, bss_data@demographics[[corr_var]])
  corr_coeff <- corr_result$corr_coeff
  bss_model@tvalues <- corr_result$tvalues
  bss_model@pvalues <- corr_result$pvalues

  bss_model@pvalues <- sign(corr_coeff)*bss_model@pvalues
  bss_model@pvalues[is.na(bss_model@pvalues)] <- 1 # Set the p-values with the NA correlations to 1
  bss_model@corr_values <- corr_coeff
  bss_model@pvalues_adjusted <- bss_p_adjust(bss_model@pvalues, mult_comp)
  # bss_model@pvalues_adjusted <- p.adjust(bss_model@pvalues, 'BH')
  bss_model@corr_values[abs(bss_model@pvalues) >= 0.05] <- 0
  message('Done.')
  return(bss_model)
}

#' Vectorized correlation between a variable and brain imaging data.
#'
#' For most scenarios, the user does not need to call this function directly.
#' Instead call \code{\link{bss_corr}} which calls this function internally.
#' @param X matrix of dimensions (\eqn{N x T}), where \eqn{N} = number of subjects and \eqn{T} = number of vertices/voxels.
#' @param Y vector of length \eqn{N}.
#' @export
corr_vec <- function(X, Y) {

  N <- length(Y)
  X_dev <- sweep(X, 2, colMeans(X))
  Y_dev <- Y - mean(Y)
  corr_coeff <- as.numeric((Y_dev %*% X_dev)/sqrt(colSums(X_dev^2)*sum(Y_dev^2)))
  tvalues <- corr_coeff * sqrt((N-2)/(1-corr_coeff^2 + .Machine$double.eps))
  pvalues <- 2*pt(abs(tvalues), N-2, lower.tail = FALSE)
  return(list("tvalues"=tvalues, "pvalues"=pvalues, "corr_coeff"=corr_coeff))
}

#' T-test for for brain imaging phenotypes
#'
#' Perform independent sample and paired sample t-tests for differences between means of brain
#' imaging phenotypes for a categorical variable.
#' @param group_var Categorical variable name. This should be present in the demographics csv file associated
#' with \code{bss_data}.
#' @param  bss_data Object of type \code{\link{BssData}}.
#' @param  paired logical; is TRUE if \code{group_var} contains matching (dependent) samples. The default value is \code{FALSE}.
#' @param  mult_comp method for multiple comparisons correction. The default method is "fdr". See \code{\link{bss_p_adjust}} for valid values.
#' @details
#' The degrees of freedom are calculated using the Welch–Satterthwaite approximation by default.
#' \code{bss_data} can be of the type "cbm", "tbm", or "roi".
#'
#' @export
bss_ttest <- function(group_var, bss_data, paired = FALSE, mult_comp="fdr") {

  if (paired == FALSE)
    bss_model <- new("BssModel", model_type="unpairedttest", group_var = group_var,
                     demographics = bss_data@demographics, mspec_file="")
  else
    bss_model <- new("BssModel", model_type="pairedttest", group_var = group_var,
                     demographics = bss_data@demographics, mspec_file="")

  group1 <- levels(bss_data@demographics[[group_var]])[1]
  group2 <- levels(bss_data@demographics[[group_var]])[2]
  idx_group1 <- which(bss_data@demographics[[group_var]] == group1)
  idx_group2 <- which(bss_data@demographics[[group_var]] == group2)

  message(sprintf("The 2 groups are %s and %s", group1, group2), appendLF = TRUE)
  message('Running t-tests...', appendLF = FALSE)

  test_result <- ttest_vec(bss_data@data_array[idx_group1,], bss_data@data_array[idx_group2,], paired)
  pvalues <- test_result$pvalues
  tvalues <- test_result$tvalues

  pvalues[is.na(pvalues)] <- 1
  pvalues <- pvalues*sign(tvalues)
  bss_model@pvalues <- pvalues
  bss_model@tvalues <- tvalues
  bss_model@tvalues[abs(pvalues) >= 0.05] <- 0
  bss_model@pvalues_adjusted <- bss_p_adjust(bss_model@pvalues, mult_comp)
  # bss_model@pvalues_adjusted <- p.adjust(abs(bss_model@pvalues), 'BH')

  message('Done.')
  return(bss_model)
}

#' Vectorized t-test for for brain imaging phenotypes
#'
#' Perform independent sample and paired sample t-tests between two numerfor differences between means of brain
#' imaging phenotypes for a categorical variable.
#' For most scenarios, the user does not need to call this function directly. This function
#' will be called internally from \code{\link{bss_ttest}}
#' @param X1 matrix of dimensions (\eqn{N1 x T}), where \eqn{N1} = number of subjects and \eqn{T} = number of vertices/voxels.
#' @param  X2 matrix of dimensions (\eqn{N2 x T}), where \eqn{N2} = number of subjects and \eqn{T} = number of vertices/voxels.
#' @param  paired logical; is TRUE if \code{group_var} contains matching (dependent) samples. The default value is \code{FALSE}.
#' @details
#' For an independent samples t-test \eqn{N1} not equal to \eqn{N2}.
#' For a dependent (paired) samples t-test, \eqn{N1 = N2}.
#' The degrees of freedom are calculated using the Welch–Satterthwaite approximation by default.
#' \code{bss_data} can be of the type "cbm", "tbm", or "roi".
#'
#' @export
ttest_vec <- function(X1, X2, paired=FALSE) {

  n1 <- dim(X1)[1]
  n2 <- dim(X2)[1]

  if (paired == TRUE) {
    D = X1 - X2
    D_mean <- colMeans(D)
    D_dev <- sweep(D, 2, D_mean)
    s1 <- sqrt(colSums(D_dev^2)/(n1-1))
    tvalues <- D_mean/(s1/sqrt(n1))
    pvalues <- 2*pt(abs(tvalues),n1-1, lower.tail = FALSE)
  }
  else {
    X1_mean <- colMeans(X1)
    X2_mean <- colMeans(X2)

    X1_dev <- sweep(X1, 2, X1_mean)
    X2_dev <- sweep(X2, 2, X2_mean)

    SS1 <- colSums(X1_dev^2) # sum of squares of group 1
    SS2 <- colSums(X2_dev^2) # sum of squares of group 2

    s1_sq <- SS1/(n1-1)
    s2_sq <- SS2/(n2-1)

    se_diff <- sqrt(s1_sq/n1 + s2_sq/n2)
    tvalues <- (X1_mean - X2_mean + .Machine$double.eps)/(se_diff + .Machine$double.eps)
    # Calculate the degrees of freedom using the Welch–Satterthwaite approximation
    deg <- (s1_sq/n1 + s2_sq/n2)^2/(s1_sq^2/(n1^2*(n1-1)) + s2_sq^2/(n2^2*(n2-1)))
    pvalues <- 2*pt(abs(tvalues), deg, lower.tail = FALSE)
  }
  return(list("tvalues"=tvalues, "pvalues"=pvalues))
}

#' Adjust p-values for multiple comparisons testing
#'
#' Perform multiple comparisons correction for mass univariate tests.
#' @param pvalues numeric vector of p values.
#' @param method character string specifying the method for correction. The default method is 'fdr'. Supported methods are
#' all types given in \code{\link{p.adjust.methods}}
#'
#' @export
bss_p_adjust <- function(pvalues, method='fdr') {

  valid_methods <- c(p.adjust.methods, "perm")
  if ( !(method  %in% valid_methods) ) {
    warning(sprintf("%s is not a valid multiple comparisons method. Using no correction.", method), call. = FALSE)
    return(pvalues)
  }

  if (method %in% p.adjust.methods && method != "perm") {
    return(p.adjust(abs(pvalues), method))
  }

  # if (method == "perm") {
  #   tvalues_null <- bss_perm()
  # }


}



# bss_maxTperm <- function(main_effect="", covariates="", bss_data, num_of_perm) {
# 
#   if (class(bss_data) == "BssROIData") {
#     return(bss_roi_anova(main_effect = main_effect, covariates = covariates, bss_data = bss_data))
#   }
# 
#   message('Running the statistical model and permutation tests. This may take a while...', appendLF = FALSE)
#   bss_lm_full <- lm_vec(main_effect = main_effect, covariates = covariates, bss_data = bss_data)
#   bss_lm_null <- lm_vec(main_effect = "", covariates = covariates, bss_data = bss_data)
#   bss_model <- anova_vec(bss_lm_full, bss_lm_null, bss_data)
# 
#   ## TODO: maybe give user option to set number of cores?
#   cl <- parallel::makeCluster(parallel::detectCores())
#   registerDoParallel(cl)
# 
#   null_distribution_t <- maxTperm(main_effect = main_effect, covariates = covariates, bss_data = bss_data, num_of_perm)
#   pvalues_adjusted <- perm_p_adjust(main_effect = main_effect, covariates = covariates, bss_data, null_distribution_t)
# 
#   #TODO: add stopCluster(cl)
#   bss_model@pvalues[is.nan(bss_model@pvalues)] <- 1
#   bss_model@pvalues <- bss_model@pvalues*bss_model@tvalues_sign
#   bss_model@tvalues[abs(bss_model@pvalues) >= 0.05] <- 0
#   bss_model@pvalues_adjusted <- pvalues_adjusted
#   # bss_model@pvalues_adjusted <- p.adjust(bss_model@pvalues, 'BH')
#   message('Done.')
#   return(bss_model)
# }


maxTperm <- function(main_effect = "", covariates = "", bss_data, num_of_perm){
  
  N <- dim(bss_data@data_array)[1]
  
  bss_lm_full <- lm_vec(main_effect = main_effect, covariates = covariates, bss_data = bss_data)
  bss_lm_null <- lm_vec(main_effect = "", covariates = covariates, bss_data = bss_data)
  
  # (1) compute full model t-statistic
  T0 <- bss_lm_full@tvalues
  maxT0 <- T0[ which.max( abs(T0) ) ]
  
  # (2) compute estimated gamma_hat and estimated residuals from reduced model
  gamma_hat <- bss_lm_null@beta_coeff
  residuals_null <- bss_lm_null@residuals
  
  # (3) compute a set of permuted data Y
  bss_data_cp <- bss_data
  
  tvalues_all <- foreach(j=1:(num_of_perm-1), .export=c('T0', 'maxT0', 'N', 'residuals_null', 'bss_lm_null', 'gamma_hat','bss_data_cp',
                                                     'main_effect', 'covariates', 'lm_vec'), .packages=c('Matrix', 'bit')) %dopar% {
                                                       set.seed(j)
                                                       pmatrix <- as(sample(N), "pMatrix")
                                                       Y_j <- (pmatrix %*% residuals_null) + (bss_lm_null@X_design_null %*% gamma_hat)
                                                       
                                                       # (4) regress permuted data Y_j against the full model
                                                       # bss_data_cp <- bss_data
                                                       bss_data_cp@data_array <- Y_j
                                                       bss_lm_full_perm <- lm_vec(main_effect = main_effect, covariates = covariates, bss_data = bss_data_cp )
                                                       
                                                       ## binarize vector after comparing permuted T and observed T
                                                       idx <- which( abs(T0) <= abs(bss_lm_full_perm@tvalues) )
                                                       
                                                       t_bin <- as.bit(rep(FALSE, dim(bss_data_cp@data_array)[2]))
                                                       t_bin[idx] <- TRUE
                                                       
                                                       tvalue_j <- bss_lm_full_perm@tvalues[ which.max( abs(bss_lm_full_perm@tvalues) ) ]
                                                       
                                                       return(list(t_bin, tvalue_j))
                                                       
                                                     }
  tvalues_all[[num_of_perm]] <- list(as.bit(rep(TRUE, dim(bss_data@data_array)[2])), maxT0)
  
  pvalues <- foreach(n=1:dim(bss_data@data_array)[2], .export=c('tvalues_all', 'num_of_perm')) %dopar% {
    count <- sum(sapply( seq(1, num_of_perm), function(x) (tvalues_all[[x]][[1]][n])) == TRUE)
    pvalue <- as.double(count/num_of_perm)
    pvalue_widx <- array(c(pvalue, n), dim=c(1,2))
    # pvalue_widx <- list(pvalue, n)
  }
  
  p_sort <- do.call('rbind', pvalues)
  pvalues_sort <- p_sort[order(p_sort[,2]), ]
  
  tvalues_null <- sapply( seq(1, num_of_perm), function(x) tvalues_all[[x]][[2]])
  
  return(list(pvalues_sort[,1], tvalues_null))
  
}


perm_p_adjust <- function(main_effect = "", covariates = "", bss_data, tvalues_null){
  bss_lm_full <- lm_vec(main_effect = main_effect, covariates = covariates, bss_data = bss_data)

  tvalues <- bss_lm_full@tvalues

  pvalues <- foreach(i=1:length(tvalues), .export=c()) %dopar% {
    p <- sum(abs(tvalues_null) >= abs(tvalues[i])) / length(tvalues_null)
    pvalue_widx <- array(c(p,i), dim = c(1,2))
  }

  p_sort <- do.call('rbind', pvalues)
  pvalues_sort <- p_sort[order(p_sort[,2]), ]
  return(pvalues_sort[,1])
}


