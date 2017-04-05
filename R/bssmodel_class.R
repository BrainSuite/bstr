#' Defines an S4 class for the statistical model
#' @export

BssModel <- setClass(
  "BssModel",
  slots = list(
    mspec_file = "character",
    main_effect = "character",
    covariates = "character",
    corr_var = "character",
    corr_values = "numeric",
    group_var = "character",
    model_type = "character",
    fullmodel = "character",
    fullvars = "character",
    nullvars = "character",
    nullmodel = "character",
    X_design_full = "matrix",
    X_design_null = "matrix",
    Npfull = "integer",
    Npnull = "integer",
    unique = "character",
    pvalues = "numeric",
    tvalues = "numeric",
    beta_coeff = "matrix",
    pvalues_adjusted = "numeric",
    stats_commands = "vector"
  )
)

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
    # Check if group_var is a factor having exactly 2 levels
    if( nlevels(demographics[[group_var]]) != 2)
      stop("group_var should be a factor having exactly 2 levels.\n", call. = FALSE)

    # If model_type is pairedttest group_var should have is a factor having exactly 2 levels
    if( model_type == "pairedttest") {
      group1 <- levels(bss_data@demographics[[group_var]])[1]
      group2 <- levels(bss_data@demographics[[group_var]])[2]
      group1_elems <- bss_data@demographics[[group_var]][bss_data@demographics[[group_var]] == group1]
      group2_elems <- bss_data@demographics[[group_var]][bss_data@demographics[[group_var]] == group2]
      if (! length(group1_elems) == length(group2_elems))
        stop(sprintf("For a paired design, there should be equal number of subjects for the two levels: %s. \nPlease check for missing or duplicate data.\n",
                     paste(group1, group2, sep=', ')), call. = FALSE)
    }

    return(list("main_effect_present"=FALSE, "covariates_present"=FALSE, "corr_var_present"=FALSE, "group_var_present"=TRUE))
  }

  # TODO: Validate covariates
}

# TODO: Call read_modelspec from within initialize
setMethod("initialize", valueClass = "BssModel", signature = "BssModel",
          function(.Object, model_type, main_effect="", covariates="", corr_var="", group_var="", demographics, mspec_file) {

          parse_model_result <- parse_model(main_effect, covariates, corr_var, group_var, model_type, demographics)
          .Object@main_effect <- main_effect
          .Object@covariates <- covariates
          .Object@corr_var <- corr_var
          .Object@group_var <- group_var
          .Object@model_type <- model_type

          if (parse_model_result$main_effect_present & parse_model_result$covariates_present) {
            .Object <- initialize_lm(.Object, main_effect, covariates, demographics)
            return (.Object)
          }

          .Object@mspec_file <- mspec_file
          return(.Object)
})


setGeneric("initialize_lm", valueClass = "BssModel", function(.Object, main_effect, covariates, demographics) {
  standardGeneric("initialize_lm")
})


setMethod("initialize_lm", signature("BssModel", "character", "character", "data.frame"), function(.Object, main_effect, covariates, demographics) {
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
  return(.Object)
})

setGeneric("run", valueClass = "BssModel", function(bss_model, bss_data) {
  standardGeneric("run")
})

#' @export
setMethod("run", signature = "BssModel", function(bss_model, bss_data) {
  return(dispatch(bss_model, bss_data))
})

setMethod("run", signature = c("BssModel", "BssROIData"), function(bss_model, bss_data) {
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
})

setGeneric("dispatch", valueClass = "BssModel", function(bss_model, bss_data) {
  standardGeneric("dispatch")
})

setMethod("dispatch", signature = "BssModel", function(bss_model, bss_data) {
  # Check the model type and call the appropriate method
  switch(bss_model@model_type,
         bss_lm = { bss_model <- bss_lm(bss_model@main_effect, bss_model@covariates, bss_data) }
         # bss_corr = { bss_model <- bss_corr(bss_model, bss_data) }
  )
  return(bss_model)
})

model_type_list <- list(
  bss_lm = 'bss_lm',
  bss_corr = 'bss_corr',
  roi = 'pairedttest',
  dbm = 'unpairedttest'
)
