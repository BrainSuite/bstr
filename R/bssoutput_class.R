#' Bss output class
#' Defines an S4 class for statistical outputs
#' @export

check_files <- function(object){
  if (!dir.exists(object@subjdir)) {
    stop(sprintf("Subjects directory %s does not exist.\n", object@subjdir), call. = FALSE)
  }

  if (!file.exists(object@csv)) {
    stop(sprintf("Demographics csv file %s does not exist.\n", object@csv), call. = FALSE)
  }
}

BssOutput <- setClass(
  "BssOutput",
  slots = list(
    outdir = "character"
  ),
  validity = check_files
)

setMethod("initialize", valueClass = "BssOutput", signature = "BssOutput", function(.Object, outdir) {
  if (!dir.exists(outdir)) {
    dir.create(outdir)
    .Object@outdir <- outdir
  }
  else {
    .Object@outdir <- outdir
    warning(sprintf("The output directory %s already exists. Will overwrite it's contents.", outdir),
            immediate. = TRUE, call.=FALSE)
  }
  return(.Object)
})

#' @export
setGeneric("save_out", valueClass = "BssOutput", function(bss_out, bss_data, bss_model) {
  standardGeneric("save_out")
})

#' @export
setMethod("save_out", valueClass = "BssOutput", signature = "BssOutput", function(bss_out, bss_data, bss_model) {
  print("in save")
})

BssCBMOutput <- setClass(
  "BssCBMOutput",
  contains = "BssOutput"
)

BssTBMOutput <- setClass(
  "BssTBMOutput",
  contains = "BssOutput"
)

BssROIOutput <- setClass(
  "BssROIOutput",
  contains = "BssOutput"
)

setMethod("save_out", valueClass = "BssCBMOutput", signature = "BssCBMOutput", function(bss_out, bss_data, bss_model) {

  log_pvalues <- log10_transform(bss_model@pvalues)
  s1 <- bss_data@atlas_surface

  s1$attributes <- log_pvalues
  bss_cmap <- new("BssColormap", "logpvalue", "Spectral", log_pvalues)
  s1$vColor <- get_colors(bss_cmap)
  s1$vColor <- matrix(s1$vColor, nrow=3, ncol=s1$hdr$nVertices, byrow = TRUE)
  writedfs(file.path(bss_out@outdir, 'log_pvalues.dfs'), s1)

  log_pvalues <- log10_transform(p.adjust(bss_model@pvalues, method = 'BH'))
  s1$attributes <- log_pvalues
  bss_cmap <- new("BssColormap", "logpvalue", "Spectral", log_pvalues)
  s1$vColor <- get_colors(bss_cmap)
  s1$vColor <- matrix(s1$vColor, nrow=3, ncol=s1$hdr$nVertices, byrow = TRUE)
  writedfs(file.path(bss_out@outdir, 'log_pvalues_adjusted.dfs'), s1)


  bss_model@tvalues[abs(log_pvalues) <= -1*log10(0.05)] <- 0
  bss_cmap <- new("BssColormap", "tvalue", "Spectral", bss_model@tvalues)
  s1$attributes <- bss_model@tvalues
  s1$vColor <- get_colors(bss_cmap)
  s1$vColor <- matrix(s1$vColor, nrow=3, ncol=s1$hdr$nVertices, byrow = TRUE)
  writedfs(file.path(bss_out@outdir, 'tvalues.dfs'), s1)

  return(bss_out)
  }
)

setMethod("save_out", valueClass = "BssTBMOutput", signature = "BssTBMOutput", function(bss_out, bss_data, bss_model) {

  log_pvalues <- rep(1, length(bss_data@atlas_image))
  log_pvalues[bss_data@mask_idx] <- log10_transform(bss_model@pvalues)
  dim(log_pvalues) <- dim(bss_data@atlas_image)
  RNifti::writeNifti(log_pvalues, file.path(bss_out@outdir, 'log_pvalues.nii.gz'), template = bss_data@atlas_image)

  return(bss_out)
  }
)

setMethod("save_out", valueClass = "BssROIOutput", signature = "BssROIOutput", function(bss_out, bss_data, bss_model) {

  outdir = bss_out@outdir


  nb_header <- "### BrainSuite ROI statistical analysis report"
  nb_libraries <-"```{r librar_cmds, echo=FALSE}
  library('bssr')
  ```"

  nb_load_data <- sprintf("```{r echo=FALSE, message=FALSE, warning=FALSE, load_data}\n")
  # nb_load_data <- paste(nb_load_data, bss_data@load_data_command, sep = "")
  nb_load_data <- paste(nb_load_data, "\nDT::datatable(bss_data@demographics)\n", sep = "")
  nb_load_data <- paste(nb_load_data, "\n```\n", sep = "")


  nb_commands <- "```{r warning=FALSE, run_command}\n"

  for (i in bss_model@stats_commands) {
    nb_commands <- paste(nb_commands, i, "\n", sep = "")
  }
  nb_commands <- paste(nb_commands, "```\n", sep = "")
  nb_commands <- paste(nb_commands, sprintf("\n#### Main effect of %d %s on %s controlling for %s
                                            ", bss_data@roiid, bss_data@roimeas, bss_model@main_effect,
                                            bss_model@covariates ), sep = "")

  rmdfileconn<-file(file.path(outdir, "report.Rmd"))
  writeLines(c(nb_header, nb_libraries, nb_load_data, nb_commands), rmdfileconn)
  close(rmdfileconn)

  # Render the markdown
  rmarkdown::render(file.path(outdir, "report.Rmd"), output_file=file.path(outdir, "report.html"), quiet = TRUE)

  return(bss_out)
  }
)



