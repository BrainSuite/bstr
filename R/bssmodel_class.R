#' Defines an S4 class for the statistical model
#' @export

BssModel <- setClass(
  "BssModel",
  slots = list(
    mspec_file = "character",
    main_effect = "character",
    covariates = "character",
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
    pvalues_adjusted = "numeric",
    stats_commands = "vector"
  )
)

parse_model <- function(main_effect, covariates, demographics) {

  if (grepl("\\+", main_effect)) {
    stop(sprintf("Main effect should only contain a single variable. You specified it as %s.\n", main_effect))
  }

  if (grepl(main_effect, covariates)) {
    stop(sprintf("Main effect *%s* also occurs in the list of covariates *%s*.
                 Main effect and covariates should be disjoint. %s.\n", main_effect, covariates))
  }

  if (!main_effect %in% colnames(demographics)) {
    stop(sprintf("Main effect *%s* doesn't occur in the demographics csv file.\n", main_effect))
  }

  # TODO: Validate covariates
  }

# TODO: Call read_modelspec from within initialize
setMethod("initialize", valueClass = "BssModel", signature = "BssModel",
          function(.Object, model_type, main_effect, covariates, demographics, mspec_file) {
            parse_model(main_effect, covariates, demographics)
            .Object@main_effect <- main_effect
            .Object@covariates <- covariates
            .Object@model_type <- model_type
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
            .Object@mspec_file <- mspec_file

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
         bss_lm = { bss_model <- bss_lm(bss_model, bss_data) },
         bss_corr = { bss_model <- bss_corr(bss_model, bss_data) }
  )
  return(bss_model)
})


setGeneric("bss_lm", valueClass = "BssModel", function(bss_model, bss_data) {
  standardGeneric("bss_lm")
})

setMethod("bss_lm", signature = "BssModel", function(bss_model, bss_data) {
  # Check the model type and call the appropriate method
  message('Running the statistical model. This may take a while...', appendLF = FALSE)
  Xtemp <- solve(t(bss_model@X_design_full) %*% bss_model@X_design_full) %*% t(bss_model@X_design_full) # Pre Hat matrix
  beta_full <- Xtemp %*% bss_data@data_array  # beta coefficients
  y_full <- bss_model@X_design_full %*% beta_full  # Predicted response
  RSS_full <- colSums((bss_data@data_array - y_full)^2)

  Xtemp <- solve(t(bss_model@X_design_null) %*% bss_model@X_design_null) %*% t(bss_model@X_design_null) # Pre Hat matrix
  beta_null <- Xtemp %*% bss_data@data_array  # beta coefficients
  y_null <- bss_model@X_design_null %*% beta_null  # Predicted response
  RSS_null <- colSums((bss_data@data_array - y_null)^2)

  N <- nrow(bss_data@data_array)
  Fstat <- (RSS_null - RSS_full)/RSS_full * (N - bss_model@Npfull - 1)/(bss_model@Npfull - bss_model@Npnull)  # F statistic
  model_unique_idx <- which(bss_model@unique %in% bss_model@fullvars) + 1    # Add 1, because the first column in the design matrix is the intercept

  se_full_unique <- sqrt(diag(solve(t(bss_model@X_design_full) %*% bss_model@X_design_full)))

  se_full_unique <- sqrt(diag(solve(t(bss_model@X_design_full) %*% bss_model@X_design_full)))[model_unique_idx] *
    sqrt(RSS_full / (N - bss_model@Npfull - 1))

  tvalue_sign <- (beta_full[model_unique_idx, ] + .Machine$double.eps)/(abs(beta_full[model_unique_idx, ]) + .Machine$double.eps)

  pvalues <- 1 - pf(Fstat, bss_model@Npfull - bss_model@Npnull, N - bss_model@Npfull - 1)

  pvalues[is.nan(pvalues)] <- 1

  pvalues <- pvalues*tvalue_sign
  tvalues <- beta_full[model_unique_idx, ]/(se_full_unique + .Machine$double.eps)
  bss_model@pvalues <- pvalues
  bss_model@tvalues <- tvalues
  bss_model@pvalues_adjusted <- p.adjust(bss_model@pvalues, 'BH')
  message('Done.')
  return(bss_model)
})

setGeneric("bss_corr", valueClass = "BssModel", function(bss_model, bss_data) {
  standardGeneric("bss_corr")
})

setMethod("bss_corr", signature = "BssModel", function(bss_model, bss_data) {
  return(bss_model)

})
