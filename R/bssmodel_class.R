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


#' S4 class for representing the statistical model
#' @slot mspec_file modelspec file
#' @slot main_effect character string containing an independent variable whose effect you want to measure.
#' It could be disease status, age, gender etc. This should strictly be a single variable. This can be
#' either a categorical or a continuous variable.
#' @slot covariates character string containing a set of other predictors (variables) in the model. If more than
#' one covariates are included, they should be separated by a `+` operator similar to an R formula.
#' @slot corr_var character variable name. This should be present in the demographics csv file associated
#' with `bstr_data`.
#' @slot corr_values numeric vector to store correlation coefficients
#' @slot corr_values_masked_adjusted numeric vector storing the masked correlation coefficients corresponding to the adjusted p-values
#' @slot group_var Categorical variable name. This should be present in the demographics csv file associated
#' with `bstr_data`.
#' @slot model_type character string denoting the type of model. Should be one of `"bstr_anova"`,
#' `"bstr_corr"`, `"bstr_corr"`, `"pairedttest"`, `"unpairedttest"` or `"bstr_lmer"`
#' @slot fullmodel character string like an R formula denoting the full model including both
#' the main effect and covariates.
#' @slot nullmodel character string like an R formula denoting the null model including covariates
#' only
#' @slot fullvars character list of individual variables in the full model
#' @slot nullvars character list of individual variables in the null model
#' @slot lm_formula R [formula()]
#' @slot X_design_full a design matrix of the type [model.matrix()] for the full model
#' @slot X_design_null a design matrix of the type [model.matrix()] for the null model
#' @slot Npfull number of variables in the full model
#' @slot Npnull number of variables in the null model
#' @slot unique unique variable obtained as a `setdiff(fullvars, nullvars)`
#' @slot pvalues numeric vector storing the p-values
#' @slot tvalues numeric vector storing the t-statistics
#' @slot tvalues_sign numeric vector storing the sign of the t-statistics
#' @slot Fstat numeric vector storing the F-statistics
#' @slot beta_coeff numeric vector storing the beta coefficients
#' @slot rss numeric vector storing the residual sum of squares
#' @slot residuals numeric vector storing residuals # TODO: Check if this variable can be eliminated
#' @slot se numeric vector storing the standard error
#' @slot mult_comp character type of multiple comparison adjustment. Takes values of "fdr"
#' @slot pvalues_adjusted numeric vector storing the adjusted p-values
#' @slot tvalues_adjusted numeric vector storing the t-values corresponding to the adjusted p-values
#' @slot stats_commands list of R commands (primarily for ROI analysis)
#' @slot load_data_command character string for the command used to load the data
#'
#' @export
BstrModel <- setClass(
  "BstrModel",
  slots = list(
    mspec_file = "character",
    main_effect = "character",
    covariates = "character",
    corr_var = "character",
    corr_values = "numeric",
    corr_values_masked_adjusted = "numeric",
    group_var = "character",
    mult_comp = "character",
    model_type = "character",
    fullmodel = "character",
    fullvars = "character",
    nullvars = "character",
    nullmodel = "character",
    lm_formula = "formula",
    X_design_full = "matrix",
    X_design_null = "matrix",
    Npfull = "integer",
    Npnull = "integer",
    unique = "character",
    pvalues = "numeric",
    tvalues = "numeric",
    tvalues_sign = "numeric",
    Fstat = "numeric",
    beta_coeff = "matrix",
    rss = "numeric",
    residuals = "matrix",
    se = "numeric",
    pvalues_adjusted = "numeric",
    tvalues_adjusted = "numeric",
    stats_commands = "vector",
    load_data_command = "character"
  )
)
#' Returns an error if necessary elements are missing in the model
#' @param main_effect string designating which column of demographics is the main effect
#' @param covariates string designating which column of demographics is the covariates
#' @param corr_var string designating which column of demographics is the correlation variable
#' @param group_var string designating which column of demographics is the group variable
#' @param model_type string designating the type of model
#' @param demographics data frame of the demographics
#'

parse_lm <- function(main_effect="", covariates="", corr_var="", group_var = "", model_type="", demographics) {
  main_effect_present <- !(main_effect == "")
  covariates_present <- !(covariates == "")
  if (!main_effect_present && !covariates_present)
    stop('Either the main effect or the covariates should be specified.', call. = FALSE)

  if (main_effect_present) {
    if (grepl(main_effect, covariates)) {
      stop(sprintf("Main effect *%s* also occurs in the list of covariates *%s*.\nMain effect and covariates should be disjoint. %s.\n",
                     main_effect, covariates), call. = FALSE)
    }
  }

  if (main_effect_present) {
    if (!is.null(main_effect)) {
      if (grepl("\\+", main_effect)) {
        stop(sprintf("Main effect should only contain a single variable. You specified it as %s.\n", main_effect))
      }
    }
  }
  return(list("main_effect_present"=main_effect_present, "covariates_present"=covariates_present,
              "corr_var_present"=FALSE, "group_var_present"=FALSE))

}

#' Returns a list denoting which elements are present in the model
#' @param main_effect string designating which column of demographics is the main effect
#' @param covariates string designating which column of demographics is the covariates
#' @param corr_var string designating which column of demographics is the correlation variable
#' @param group_var string designating which column of demographics is the group variable
#' @param model_type string designating the type of model
#' @param demographics data frame of the demographics
#'

parse_model <- function(main_effect="", covariates="", corr_var="", group_var = "", model_type="", demographics) {

  if (! model_type %in% model_type_list) {
    stop(sprintf('model_type should be one of the following: %s',
                 paste(unlist(model_type_list, use.names = FALSE), collapse = ', ')), call. = FALSE)
  }
  main_effect_present <- !(main_effect == "")
  covariates_present <- !(covariates == "")
  corr_var_present <- !(corr_var == "")
  group_var_present <- !(group_var == "")

  if (main_effect_present & covariates_present & corr_var_present & group_var_present)
    stop('Only the main effect and covariates or corr_var or group_var should be specified separately.', call. = FALSE)

  if ( (main_effect_present & !covariates_present) | (covariates_present & !main_effect_present) )
    stop('main_effect and covariates should be specified together.', call. = FALSE)

  if (!main_effect_present & !covariates_present & !corr_var_present & !group_var_present)
    stop('Either the main effect and covariates or corr_var or group_var should be specified.', call. = FALSE)

  if (main_effect_present & covariates_present) {
    if (!is.null(main_effect)) {
      if (grepl("\\+", main_effect)) {
        stop(sprintf("Main effect should only contain a single variable. You specified it as %s.\n", main_effect))
      }
    }

    if (!is.null(main_effect) && !is.null(covariates)) {
      if (grepl(main_effect, covariates)) {
        stop(sprintf("Main effect *%s* also occurs in the list of covariates *%s*.\nMain effect and covariates should be disjoint. %s.\n",
                     main_effect, covariates), call. = FALSE)
      }
    }

    if (!main_effect %in% colnames(demographics)) {
      stop(sprintf("Main effect *%s* doesn't occur in the demographics csv file.\n", main_effect), call. = FALSE)
    }
    return(list("main_effect_present"=TRUE, "covariates_present"=TRUE, "corr_var_present"=FALSE))
  }

  if (!corr_var == "") {
    if(!corr_var %in% colnames(demographics)) {
      stop(sprintf("corr_var *%s* doesn't occur in the demographics csv file.\n", corr_var), call. = FALSE)
    }
    return(list("main_effect_present"=FALSE, "covariates_present"=FALSE, "corr_var_present"=TRUE))
  }

  if (!group_var == "") {
    if(!group_var %in% colnames(demographics)) {
      stop(sprintf("group_var *%s* doesn't occur in the demographics csv file.\n", group_var), call. = FALSE)
    }
    demographics[[group_var]] <- as.factor(demographics[[group_var]])

    # If model_type is pairedttest group_var should have is a factor having exactly 2 levels
    if( model_type == "pairedttest") {
      # Check if group_var is a factor having exactly 2 levels
      if( nlevels(demographics[[group_var]]) != 2)
        stop("group_var should be a factor having exactly 2 levels.\n", call. = FALSE)
      group1 <- levels(demographics[[group_var]])[1]
      group2 <- levels(demographics[[group_var]])[2]
      group1_elems <- demographics[[group_var]][demographics[[group_var]] == group1]
      group2_elems <- demographics[[group_var]][demographics[[group_var]] == group2]
      if (! length(group1_elems) == length(group2_elems))
        stop(sprintf("For a paired design, there should be equal number of subjects for the two levels: %s. \nPlease check for missing or duplicate data.\n",
                     paste(group1, group2, sep=', ')), call. = FALSE)
    }

    return(list("main_effect_present"=FALSE, "covariates_present"=FALSE, "corr_var_present"=FALSE, "group_var_present"=TRUE))
  }

  # TODO: Validate covariates
}

# TODO: Call read_modelspec from within initialize
setMethod("initialize", valueClass = "BstrModel", signature = "BstrModel",
          function(.Object, model_type, main_effect="", covariates="", corr_var="", group_var="", mult_comp="", demographics, mspec_file) {

          if (model_type == "bstr_lm" || model_type == "bstr_anova")
            parse_model_result <- parse_lm(main_effect, covariates, corr_var, group_var, model_type, demographics)
          else
            parse_model_result <- parse_model(main_effect, covariates, corr_var, group_var, model_type, demographics)

          .Object@main_effect <- main_effect
          .Object@covariates <- covariates
          .Object@corr_var <- corr_var
          .Object@group_var <- group_var
          .Object@model_type <- model_type
          .Object@mult_comp <- mult_comp

          if (parse_model_result$main_effect_present || parse_model_result$covariates_present) {
            .Object <- initialize_lm(.Object, main_effect, covariates, demographics)
            return (.Object)
          }

          .Object@mspec_file <- mspec_file
          return(.Object)
})


setGeneric("initialize_lm", valueClass = "BstrModel", function(.Object, main_effect, covariates, demographics) {
  standardGeneric("initialize_lm")
})


setMethod("initialize_lm", signature("BstrModel", "character", "character", "data.frame"), function(.Object, main_effect, covariates, demographics) {
  .Object@fullmodel <- paste(main_effect, '+', covariates)
  .Object@nullmodel <- paste(covariates)
  # Design matrix for full model
  .Object@X_design_full <- model.matrix(formula(sprintf('~ %s', .Object@fullmodel)), data = demographics)
  # Design matrix for null model
  .Object@X_design_null <- model.matrix(formula(sprintf('~ %s', .Object@nullmodel)), data = demographics)
  .Object@Npfull <- length(unlist(strsplit(.Object@fullmodel, '\\+')))
  .Object@Npnull <- length(unlist(strsplit(.Object@nullmodel, '\\+')))

  .Object@fullvars <- unlist(lapply(unlist(strsplit(.Object@fullmodel, '\\+')), function (x) {gsub("\\s+", '', x)}))
  .Object@nullvars <- unlist(lapply(unlist(strsplit(.Object@nullmodel, '\\+')), function (x) {gsub("\\s+", '', x)}))
  .Object@unique <- setdiff(.Object@fullvars, .Object@nullvars)
  if (main_effect == "")
    .Object@lm_formula <- formula(sprintf('~ %s', covariates))
  else
    .Object@lm_formula <- formula(sprintf('~ %s', paste(main_effect, '+', covariates)))

  return(.Object)
})

model_type_list <- list(
  bstr_anova = 'bstr_anova',
  bstr_lm = 'bstr_lm',
  bstr_lmer = 'bstr_lmer',
  bstr_corr = 'bstr_corr',
  pairedttest = 'pairedttest',
  unpairedttest = 'unpairedttest'
)
